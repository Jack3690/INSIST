"""This module contains additional functions for the package"""
from pathlib import Path
from astropy.modeling import models
from astropy.coordinates import angular_separation
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
from astropy.coordinates import SkyCoord, Distance, ICRS
from scipy.spatial import cKDTree
import astropy.units as u
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

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
    int_flux = np.trapezoid(lambda_*conv_flux, lambda_)/np.trapezoid(lambda_*R_eff,
                                                             lambda_)
    W_eff = np.trapezoid(R_eff, lambda_)/R_eff.max()
    lambda_phot = np.trapezoid(lambda_**2*conv_flux,
                           lambda_)/np.trapezoid(lambda_*conv_flux, lambda_)

    c1 = lambda_ >= lambda_phot-W_eff/2
    c2 = lambda_ <= lambda_phot+W_eff/2

    R_sq = np.where((c1 & c2), 1, 0)
    flux_ratio = np.trapezoid(R_eff, lambda_)/np.trapezoid(R_sq, lambda_)

    if plot:
        ax.plot(lambda_, R_sq,
                label="Square Filter", alpha=0.7)

    # Frequency units
    R_eff_Jy = R_eff*lambda_**2*3.34e4
    flux_AB = flux_AB*lambda_**2*3.34e4
    nu = 3e18/lambda_

    conv_flux_Jy = R_eff_Jy*flux_AB
    int_flux_Jy = np.trapezoid(nu*conv_flux_Jy, nu)/np.trapezoid(nu*R_eff_Jy, nu)

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


def redshift_corr(df):
    """
    Function for redshift correction of input data

    Parameters
    ----------

    df: astropy.table.Table or pandas.DataFrame
        table with columns 'wav', 'flux', 'z1' and 'z2'

        wav: numpy.ndarray
            wavelength in Angstrom
        flux: numpy.ndarray
                flux in ergs/s/cm2/A

        z1: float,
            Original redshift to star
        z2: float,
            New redshift to star

    Returns
    -------

    astropy.table.Table or pandas.DataFrame
    """
    if 'z1' in df.keys() and 'z2' in df.keys():
        red_corr = df['z2'].value.reshape(-1, 1)-df['z1'].value.reshape(-1, 1)
        df['wav'] = df['wav'].value*(1 + red_corr)
        d1 = cosmo.luminosity_distance(df['z1']).value
        d2 = cosmo.luminosity_distance(df['z2']).value
        flux_corr = (d1/d2)**2

    elif 'd1' in df.keys() and 'd2' in df.keys():
        flux_corr = (df['d1'].value/df['d2'].value)**2

    df['flux'] = df['flux'].value*flux_corr.reshape(-1, 1)

    return df


def spectra_to_mags_df(df, inputs):
    """
    Function to convert spectra to magnitude using telescope response functions
    Using Astropy Table

    Parameters
    ----------

    df: astropy.table.Table
        table with columns 'wav', 'flux', 'd1' and 'd2' or 'z1' or 'z2'

        wav: numpy.ndarray
            wavelength in Angstrom
        flux: numpy.ndarray
                flux in ergs/s/cm2/A

        d1: float,
            Original distance to star.
        d2: float
            New distance to star

        z1: float,
            Original redshift to star
        z2: float,
            New redshift to star

    inputs: list,
            list of path to response functions

    Returns
    -------

    astropy.table.Table

    """
    mags = []
    c1 = ('z1' in df.keys() and 'z2' in df.keys())
    c2 = ('d1' in df.keys() and 'd2' in df.keys())

    if c1 or c2:
        df = redshift_corr(df)
    else:
        KeyError("z1 and z2 or d1 and d2")

    fluxes = df['flux'].value
    wavs = df['wav'].value

    for wav, flux in zip(wavs, fluxes):
        out = bandpass(wav, flux, inputs=inputs,
                       plot=False)
        params = out[3]
        int_flux_Jy = params[2]
        ABmag = -2.5*np.log10(int_flux_Jy/3631)

        if ABmag == np.nan:
            ABmag = 100
        mags.append(ABmag)

        df['mag'] = mags

    return df


def spectra_to_mags(wav, flux, inputs, z1=None, z2=None, d1=None, d2=None):
    """
    Function to convert spectra to magnitude using telescope response functions
    Using Astropy Table

    Parameters
    ----------

    wav: list of numpy.ndarray
        wavelength in Angstrom
    flux: list of numpy.ndarray
            flux in ergs/s/cm2/A

    inputs: list,
        list of path to response functions

    d1: float,
        Original distance to star.
    d2: float
        New distance to star

    z1: float,
        Original redshift to star
    z2: float,
        New redshift to star

    Returns
    -------
    ABmag: numpy.ndarray
            magnitude of stars from spectra
    """
    c1 = z1 is not None and z2 is not None
    c2 = d1 is not None and d2 is not None

    if c1:

        wav = wav*(1 + z1 - z2)
        d1 = cosmo.luminosity_distance(z1).value
        d2 = cosmo.luminosity_distance(z2).value
        flux_corr = (d1/d2)**2

    elif c2:
        flux_corr = (d1/d2)**2

    flux *= flux_corr
    out = bandpass(wav, flux, inputs=inputs,
                   plot=False)

    params = out[3]
    int_flux_Jy = params[2]
    ABmag = -2.5*np.log10(int_flux_Jy/3631)

    if ABmag == np.nan:
        ABmag = 100

    return ABmag


def calc_mos_size(L, B, PA, lw):
    """
    Function to calculate minimum size of a Multi-object Spectrometer (MOS)
    detector based on it's orientation and shape

    Parameters
    ----------

    L: float,
        Length of MOS Field of View (FoV).

    B: float,
        Breadth of MOS FoV.

    PA: float,
            Orientation of MOS FoV with respect to sky plane. (degrees)

    lw: int,
            Number of pixels in dispersion direction.

    Returns
    -------

    x_size, y_size: tuple
                        (int,int)
    """
    PA *= np.pi/180
    x_size = lw + B*np.sin(PA) + L*np.cos(PA)
    y_size = B*np.cos(PA) + L*np.sin(PA)
    x_size = int(np.round(x_size, 0))
    y_size = int(np.round(y_size, 0))

    return x_size, y_size


def count_sources_within_radius(catalog, radius):
    """
        Function for counting the number of source around each stat
        within a given radius using KDTree

        Parameters
        ----------

        catalog: (np.ndarray, np.ndarray),
                source coordinates

        radius: float,
                radius of circle to count sources within

        Returns
        -------

        counts: list,
                number density of sources around each star

    """
    tree = cKDTree(catalog)
    counts = []

    for star in catalog:
        indices = tree.query_ball_point(star, radius)
        count = len(indices) - 1  # Exclude the star itself from the count
        counts.append(count)

    return counts


def select_mos(df, cen_x, cen_y, radius=10, min_sep=10, ny=1000):
    """
    Function to select sources from a given database such that overlap of
    spectra is minimum
    Parameters
    ----------

        df1: astropy.table.Table,
        table with columns 'ra' and 'dec'

        ra: float,
            Right Ascension in degrees

        dec: float,
            Declnation in degrees

        cen_x, float,
                Reference x position in pixel coordinates

        cen_y: float,
                Reference x position in pixel coordinates

        radius: float,
                Radius of circle to count stars within.

        min_sep: float,
                minimum separation between stars in the spatial axis,


        ny: int,
            length of spatial axis in pixels


    Returns
    -------

    mos_df: astropy.table.Table,
                table containing list of selected sources
    res_df: astropy.table.Table,
                table containing list of remaining sources


    """
    catalog = df[['x', 'y']].values
    result = count_sources_within_radius(catalog, radius)
    df['n_density'] = result

    top = np.flip(np.arange(cen_y, 0, -min_sep))
    bottom = np.arange(cen_y + min_sep, ny, min_sep)

    y_pos = np.concatenate([top, bottom])
    mos_df = []
    df['x_sep'] = abs(df['x'] - cen_x)
    for i in y_pos:
        p = df[(df['y'] > i-2) & (df['y'] < i+2)]
        p = p[p['n_density'] == 0]

        if len(p) > 0:
            p = p[p['x_sep'] == p['x_sep'].min()]
            mos_df.append(p.values[0])

    mos_df = pd.DataFrame(mos_df, columns=df.keys())

    res_df = mos_df[0:0]
    for id in df['objid']:
        if id not in mos_df['objid'].values:
            t = df[df['objid'] == id]
            res_df = pd.concat([res_df, t])

    return mos_df, res_df, cen_x, cen_y


def distance_transform(ras, decs, cen_ra, cen_dec, d1, d2):
    """
        Function for distance transformation of sources
        scattered around a point.

        Parameters
        ----------

        ras: np.ndarray,
            Right Ascension of sources.

        decs: np.ndarray,
            Declination of sources.

        cen_ra: float,
                Reference Right Ascension

        cen_dec: float,
                Reference Declination

        d1: astropy.units.ly, astropy.units.pc, astropy.units.Mpc

        d2: astropy.units.ly, astropy.units.pc, astropy.units.Mpc

        Returns
        -------

    new_star_coords.ra.value, new_star_coords.ra.value: np.ndarray, np.ndarray

    """

    # Define the coordinates of the galaxy in the original catalog
    galaxy_ra = cen_ra*u.deg
    galaxy_dec = cen_dec*u.deg

    # Define the coordinates of the stars in the original catalog
    star_ra = ras*u.deg
    star_dec = decs*u.deg

    # Create a SkyCoord object for the galaxy in the original catalog
    galaxy_coord = SkyCoord(ra=galaxy_ra, dec=galaxy_dec,
                            distance=Distance(d1), frame=ICRS())

    # Create a SkyCoord object for the stars in the original catalog
    star_coords = SkyCoord(ra=star_ra, dec=star_dec, frame=ICRS())

    # Calculate the factor by which to scale the coordinates
    scale_factor = d1/d2

    # Calculate the separation between the stars and the galaxy
    separation = star_coords.separation(galaxy_coord)

    # Calculate the new position angle using the original coordinates
    position_angle = star_coords.position_angle(galaxy_coord)

    # Scale the separation by the scale factor
    scaled_separation = np.arctan(np.tan(separation)*scale_factor)

    # Calculate the new star coordinates using the scaled separation
    #  and position angle
    new_star_coords = galaxy_coord.directional_offset_by(position_angle,
                                                         scaled_separation)

    return new_star_coords.ra.value, new_star_coords.dec.value
