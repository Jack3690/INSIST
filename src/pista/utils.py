"""This module contains additional functions for the package"""
from pathlib import Path
from astropy.modeling import models
from astropy.coordinates import angular_separation
from astropy.table import Table
import astropy.units as u
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

data_path = Path(__file__).parent.joinpath()

sb.set_style('white')
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['figure.figsize'] = (10, 10)


def bandpass(wav, flux, inputs, plot=True, fig=None, ax=None):
    """
    Function to convolve response functions

    Parameters
    ----------
    wav : numpy.ndarray
            wavelenth in angstrom
    flux: numpy.ndarray
            flux normalized to [0,1]
    plot: bool,
            If true shows plots with input and convolved response functions

    fig : matplotlib.pyplot.figure
            User defined figure
    ax  : matplotlib.pyplot.axes
            User defined axes

    Returns
    -------

    fig, ax, data, params

    data : tuple,
            (wavelenth array, flux_array, convolved flux array)

    params: tuple,
            (effective wavelength, integrated flux, Effective Width)

    """
    lambda_ = wav
    flux_AB = flux

    if plot:
        if fig is None or ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        ax.plot(lambda_, flux_AB/flux_AB.max(), label=r'$F(\lambda)$',
                alpha=0.7)

    R_eff = 1

    for i in inputs:
        file_name = i.split(',')[0]
        n = float(i.split(',')[1])
        f_max = float(i.split(',')[2])

        filt_dat = np.loadtxt(file_name)
        wav = filt_dat[:, 0]
        flux = filt_dat[:, 1]

        if np.amax(flux) > 1:
            flux /= f_max

        indices = np.where((wav > lambda_[0]) & (wav < lambda_[-1]))
        wav_new = wav[indices]
        flux_new = flux[indices]

        wav_new = np.concatenate([[lambda_[0]], [wav_new[0] - 1], wav_new,
                                  [wav_new[-1] + 1], [lambda_[-1]]])

        flux_new = np.concatenate([[0], [0], flux_new, [0], [0]])

        flux_out = np.interp(lambda_, wav_new, flux_new)

        R_eff *= flux_out**n

        if plot:
            ax.plot(lambda_, flux_out/flux_out.max(),
                    label=f"{file_name.split('/')[-1][:-4]}x{n}", alpha=0.7)

    # Wavelength units
    conv_flux = R_eff*flux_AB
    int_flux = np.trapz(lambda_*conv_flux, lambda_)/np.trapz(lambda_*R_eff,
                                                             lambda_)
    W_eff = np.trapz(R_eff, lambda_)/R_eff.max()
    lambda_phot = np.trapz(lambda_**2*conv_flux,
                           lambda_)/np.trapz(lambda_*conv_flux, lambda_)

    c1 = lambda_ >= lambda_phot-W_eff/2
    c2 = lambda_ <= lambda_phot+W_eff/2

    R_sq = np.where((c1 & c2), 1, 0)
    flux_ratio = np.trapz(R_eff, lambda_)/np.trapz(R_sq, lambda_)

    if plot:
        ax.plot(lambda_, R_sq,
                label="Square Filter", alpha=0.7)

    # Frequency units
    R_eff_Jy = R_eff*lambda_**2*3.34e4
    flux_AB = flux_AB*lambda_**2*3.34e4
    nu = 3e18/lambda_

    conv_flux_Jy = R_eff_Jy*flux_AB
    int_flux_Jy = np.trapz(nu*conv_flux_Jy, nu)/np.trapz(nu*R_eff_Jy, nu)

    # Comparing to a square filter with same width

    data = lambda_, conv_flux, R_eff
    params = lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio

    if plot:
        ax.plot(lambda_, conv_flux/conv_flux.max(), label='Convloved Flux',
                linewidth=5)
        y = np.linspace(0, 1)
        x = y*0 + lambda_phot
        label = r'$\lambda_{phot} = $' + f'{round(lambda_phot,3)}' + r' $\AA$'
        ax.plot(x, y, '--', color='black', label=label)

        ax.set_xlabel(r'$\AA$')
        ax.set_ylabel(r'Normalized Flux')
        fig.suptitle('Bandpass', fontsize=20, y=0.95)
        ax.legend()

    return fig, ax, data, params


def generate_psf(npix, params, function='Gaussian'):

    """
    Function for generating user defined PSF

    npix : int,
           number of pixels along one axis for pixel array

    sigma: float,
           standard deviation of the PSF in pixels

    function: str,
               type of PSF function

    Returns
    -------

    numpy.ndarray

    """
    x = np.linspace(0, npix - 1, npix)
    y = x
    yy, xx = np.meshgrid(x, y)
    if function == 'Gaussian':
        sigma_x = params['sigma_x']
        sigma_y = params['sigma_y']
        psf = models.Gaussian2D(1, npix//2, npix//2,
                                sigma_x, sigma_y)(xx, yy)
        psf /= psf.sum()

    elif function == 'Moffat':
        gamma = params['gamma']
        alpha = params['alpha']
        psf = models.Moffat2D(1, npix//2, npix//2,
                              gamma, alpha)(xx, yy)
        psf /= psf.sum()

    np.save('user_defined_psf.npy', psf)
    return psf


def Xmatch(df1, df2, r=1):
    """Function for crossmatching two
    catalogs using RAs and Decs"""
    if isinstance(df1, Table):
        df1 = df1.to_pandas()

    if isinstance(df2, Table):
        df2 = df2.to_pandas()

    # define the catalogs
    matched = []
    for i, row in df1.iterrows():
        dist = angular_separation(row['ra'], row['dec'],
                                  df2['ra'].values, df2['dec'].values)
        dist = dist*u.radian
        dist = dist.to(u.arcsec).value

        if dist.min() <= r:
            index = np.where(dist == dist.min())[0][0]
            matched.append([i, index, dist.min()])
        else:
            matched.append([i, None, None])

    return matched
