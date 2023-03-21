"""This module contains classes that can be utilized for visualizing and
    analyzing the simulated images and spectra"""
from photutils import aperture as aper
from matplotlib import colors as col
import matplotlib.pyplot as plt

from photutils.detection import DAOStarFinder
from photutils.psf import DAOPhotPSFPhotometry, FittableImageModel

from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy import units as u

import numpy as np


class Analyzer():
    """Analyzer class is used for visualizing the simulated images and for
    peforming photometry"""

    def __init__(self):
        """
        A class to visualize and analyze the simulated image
        Parameters
        ----------
        Imager.init()
        Returns
        -------
        None.

        """

    def __call__(self, df=None, wcs=None, data=None,
                 photometry=None, detect_sources=False, fwhm=3, sigma=3,
                 zero_point=None):
        """
        Performs sim simulation and sim Photometry

        Imager.call()
        do_photometry : Bool, Default : True
                        Do Aperture Photometry
        """
        self.photometry_type = photometry
        if photometry == 'Aper':
            self.aper_photometry(data, wcs, df, fwhm, sigma, detect_sources,
                                 zero_point)
        elif photometry == 'PSF':
            self.psf_photometry(data, wcs, df, fwhm, sigma, zero_point)

    def aper_photometry(self, data, wcs, input_df, fwhm, sigma, detect,
                        zero_point):
        """This function performs PSF photometry using Photutils DAOPhot
        Method"""
        # Calculate zero point
        if zero_point is None:
            zero_p_flux = self.zero_flux + self.sky_bag_flux
            zero_p_flux += np.mean(self.DR)*self.exp_time
            zero_p_flux += self.det_params['NF'] + self.det_params['bias']

            zero_p_flux *= self.gain*0.9499142715255932
        else:
            zero_p_flux = zero_point

        # if detect flag is set to True, detect sources in the image
        if detect:

            # calculate the mean, median and standard deviation of the data
            mean, median, std = sigma_clipped_stats(data, sigma=sigma)
            # create DAOStarFinder object to detect sources
            daofind = DAOStarFinder(fwhm=fwhm, threshold=median + sigma*std)
            # detect sources in the image
            sources = daofind(data)
            # get the source positions
            positions = np.transpose(
                (sources['xcentroid'], sources['ycentroid']))
        else:

            # create SkyCoord object from ra and dec values in the dataframe
            coords = SkyCoord(input_df['ra'], input_df['dec'], unit=u.deg)
            # convert the sky coordinates to pixel coordinates
            pix = wcs.world_to_array_index(coords)
            positions = [(i, j) for i, j in zip(pix[1], pix[0])]

        # create circular aperture object
        self.aps = aper.CircularAperture(positions, r=2*fwhm)
        # count number of pixels within the aperture
        ap_pix = np.count_nonzero(self.aps.to_mask()[0])
        # create circular annulus object
        self.bags = aper.CircularAnnulus(positions, r_in=3*fwhm,
                                         r_out=5*fwhm)
        # count number of pixels within the annulus
        bag_pix = np.count_nonzero(self.bags.to_mask()[0])

        # perform aperture photometry on the data
        phot_table = aper.aperture_photometry(data, [self.aps, self.bags])

        # calculate sky flux
        phot_table['sky_flux'] = phot_table['aperture_sum_1']*(ap_pix/bag_pix)
        # calculate source flux
        phot_table['flux'] = phot_table['aperture_sum_0'].value - \
            phot_table['sky_flux'].value
        # calculate error on the source flux
        phot_table['flux_err'] = np.sqrt(
            phot_table['flux'].value + phot_table['sky_flux'].value)

        # calculate signal to noise ratio
        phot_table['SNR'] = phot_table['flux'].value / \
            phot_table['flux_err'].value

        if not detect and zero_point is None and len(phot_table) > 2:
            phot_table['mag_in'] = input_df['mag'].values
            phot_table['ra'] = input_df['ra'].values
            phot_table['dec'] = input_df['dec'].values
            if len(phot_table) > 3:
                zero_p_flux = 0

                for i in range(3):
                    flux = pow(10, -0.4*phot_table['mag_in'].value[i])
                    zero_p_flux += phot_table['flux'].value[i]/flux
                zero_p_flux /= 3

        phot_table['mag_out'] = -2.5*np.log10(phot_table['flux']/zero_p_flux)
        phot_table['mag_err'] = 1.082/phot_table['SNR']
        self.header['ZP'] = zero_p_flux

        self.header['EXPTIME'] = self.exp_time
        self.header['BUNIT'] = 'DN'

        coords = np.array(wcs.pixel_to_world_values(positions))
        phot_table['ra'] = coords[:, 0]
        phot_table['dec'] = coords[:, 1]
        self.phot_table = phot_table

    def psf_photometry(self, data, wcs, input_df, fwhm, sigma, ZP):

        mean, median, std = sigma_clipped_stats(data, sigma=3)

        psf_model = FittableImageModel(self.image_g_sub)
        self.psf_model = psf_model

        photometry = DAOPhotPSFPhotometry(crit_separation=2,
                                          threshold=mean + sigma*std,
                                          fwhm=fwhm,
                                          aperture_radius=3,
                                          psf_model=psf_model,
                                          fitter=LevMarLSQFitter(),
                                          fitshape=(11, 11),
                                          niters=1)

        result_tab = photometry(image=data)
        positions = np.array([result_tab['x_fit'], result_tab['y_fit']]).T
        coords = np.array(wcs.pixel_to_world_values(positions))

        result_tab['ra'] = coords[:, 0]
        result_tab['dec'] = coords[:, 1]
        result_tab['SNR'] = result_tab['flux_fit']/result_tab['flux_unc']

        if ZP is None:
            zero_p_flux = self.zero_flux + self.sky_bag_flux
            zero_p_flux += np.mean(self.DR)*self.exp_time + \
                self.det_params['NF'] + self.det_params['bias']
            zero_p_flux *= self.gain
        else:
            zero_p_flux = ZP

        result_tab['mag_out'] = -2.5 * \
            np.log10(result_tab['flux_fit']/zero_p_flux)
        result_tab['mag_err'] = 1.082/result_tab['SNR']

        self.header['ZP'] = zero_p_flux

        self.phot_table = result_tab

    def show_field(self, figsize=(10, 10), marker='.', cmap='jet'):
        """
        Function for creating a scatter plot of sources within the FoV

        Parameters
        ----------
        figsize : tuple,
                  Figure size

        Returns
        -------
        fig, ax
        """
        # Cropping data frame to show source within n_x x n_y

        x_min_cut = (self.df['x'] > self.n_pix_sub - 1)
        x_max_cut = (self.df['x'] < self.n_y_sim - self.n_pix_sub - 1)

        img_df = self.df[x_min_cut & x_max_cut]

        y_min_cut = (self.df['y'] > self.n_pix_sub - 1)
        y_max_cut = (self.df['y'] < self.n_x_sim - self.n_pix_sub - 1)

        img_df = img_df[y_min_cut & y_max_cut]

        fov_x = (self.n_x*self.pixel_scale)/3600
        fov_y = (self.n_y*self.pixel_scale)/3600

        fov_x = np.round(fov_x, 4)
        fov_y = np.round(fov_y, 4)

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        x = img_df['ra']
        y = img_df['dec']
        c = img_df['mag']

        img = ax.scatter(x, y, c=c, marker=marker, cmap=cmap)
        cb = plt.colorbar(img)
        cb.set_label('mag (ABmag)')
        ax.set_title(f"""Requested Center : {self.name} | {len(img_df)} sources
        Fov(RA) : {fov_x} (deg) | Fov(Dec) : {fov_y} (deg)""")
        ax.invert_xaxis()

        ax.set_xlabel('RA (Degrees)')
        ax.set_ylabel('Dec (Degrees)')

        return fig, ax

    def show_image(self, source='Digital', fig=None, ax=None, cmap='jet',
                   figsize=(15, 10), download=False, show_wcs=True,
                   overlay_apertures=False):
        """
        Function for plotting the simulated field image

        Source: str,
                Choose from
                            'Digital' : Final digial image
                            'Charge'  : electrons, Light(Source + sky)
                                        + Dark Current + Noises
                            'Source'  : Source + Sky + Noises
                            'Sky'     : Sky + shot_noise
                            'DC'      : Dark Current + DNFP
                            'QE'      : Quantum efficiency fluctuation
                                          across detector
                            'Bias'    : Charge offset
                            'PRNU'    : Photon Response Non-Uniformity
                            'DNFP'    : Dark Noise Fixed Pattern
                            'QN'      : Quantization Noise


        fig : matplotlib.pyplot.figure
              User defined figure
        ax  : matplotlib.pyplot.axes
              User defined axes
        cmap : str,
              matplotlib.pyplot colormap
        figsize : tuple
        download : bool
        show_wcs : bool
                  If true adds WCS projection to the image
        Returns
        -------
        Image

        fig, ax
        """
        if np.all(self.image) is not None:
            if fig is None or ax is None:
                fig = plt.figure(figsize=figsize)
                if show_wcs:
                    ax = fig.add_subplot(projection=self.wcs)
                else:
                    ax = fig.add_subplot()

            norm = None
            if source == 'Digital':
                data = self.digital
                norm = col.LogNorm()
            elif source == 'Charge':
                data = self.charge
                norm = col.LogNorm()
            elif source == 'Source':
                data = self.light_array
                norm = col.LogNorm()
            elif source == 'Sky':
                data = self.sky_photons
            elif source == 'DC':
                data = self.DC_array
            elif source == 'Bias':
                data = self.bias_array
            elif source == 'PRNU':
                data = self.PRNU_array
            elif source == 'DNFP':
                norm = col.LogNorm()
                data = self.DNFP_array
            elif source == 'QN':
                data = self.QN_array
            else:
                print("Invalid Input")
                return None

            if data.min() < 0:
                print('Negative values in image. Increase Bias')
                data += data.min()

            img = ax.imshow(data, cmap=cmap, norm=norm)
            ax.grid(False)
            cb = plt.colorbar(img, ax=ax)
            cb.set_label('DN')
            ax.set_title(f'{source} \nRequested center : {self.name}')
            ax.grid(False)

            if overlay_apertures and self.photometry_type == "Aper":
                for aperture in self.aps:
                    if aperture is not None:
                        aperture.plot(ax=ax, color='red', lw=1.5)
                for aperture in self.bags:
                    if aperture is not None:
                        aperture.plot(ax=ax, color='yellow', lw=1.5)

            if download:
                fig.savefig(f"{source}.png", format='png')
            return fig, ax

        else:
            print("Run Simulation")

    def show_hist(self, source='Digital', bins=None,
                  fig=None, ax=None, figsize=(15, 8)):
        """
        Function for plotting histogram of various stages of simulation

        Parameters
        ----------

        Source: str,
                Choose from
                          'Digital' : Final digial image
                          'Charge'  : electrons, Light(Source + sky)
                                      + Dark Current + Noises
                          'Source'  : Source + Sky + Noises
                          'Sky'     : Sky + shot_noise
                          'DC'      : Dark Current + DNFP
                          'QE'      : Quantum efficiency fluctuation
                                       across detector
                          'Bias'    : Charge offset
                          'PRNU'    : Photon Response Non-Uniformity
                          'DNFP'    : Dark Noise Fixed Pattern
                          'QN'      : Quantization Noise

        bins : numpy.array,
              bins for making histogram
        fig : matplotlib.pyplot.figure
              User defined figure
        ax  : matplotlib.pyplot.axes
              User defined axes
        figsize : tuple
        """

        if np.all(self.image) is not None:
            if fig is None or ax is None:
                fig, ax = plt.subplots(1, 1, figsize=figsize)

            if source == 'Digital':
                data = self.digital.ravel()
            elif source == 'Charge':
                data = self.charge.ravel()
            elif source == 'Source':
                data = self.light_array
            elif source == 'Sky':
                data = self.sky_photons.ravel()
            elif source == 'DC':
                data = self.DC_array.ravel()
            elif source == 'Bias':
                data = self.bias_array.ravel()
            elif source == 'PRNU':
                data = self.PRNU_array.ravel()
            elif source == 'DNFP':
                data = self.DNFP_array.ravel()
            elif source == 'QN':
                data = self.QN_array.ravel()

            if bins is None:
                bins = np.linspace(data.min(), data.max(), 20)
            ax.hist(data, bins=bins)
            ax.set_title(f'{source} histogram')
            ax.set_ylabel('Count')
            ax.set_yscale('log')
            return fig, ax
        else:
            print("Run Simulation")

    def get_image(self, source='Digital'):
        """
        Function of retrieving image array at different stages of simulation.

        Parameters
        ----------

        Source: str,
            Choose from
                          'Digital' : Final digial image
                          'Charge'  : electrons, Light(Source + sky)
                                      + Dark Current + Noises
                          'Source'  : Source + Sky + Noises
                          'Sky'     : Sky + shot_noise
                          'DC'      : Dark Current + DNFP
                          'QE'      : Quantum efficiency fluctuation
                                      across detector
                          'Bias'    : Charge offset
                          'PRNU'    : Photon Response Non-Uniformity
                          'DNFP'    : Dark Noise Fixed Pattern
                          'QN'      : Quantization Noise
        """

        if np.all(self.image) is not None:
            if source == 'Digital':
                data = self.digital
            elif source == 'Charge':
                data = self.charge
            elif source == 'Sky':
                data = self.sky_photoelec
            elif source == 'DC':
                data = self.DC_array
            elif source == 'QE':
                data = self.qe_array
            elif source == 'Bias':
                data = (self.bias_array + self.DC_array)
            elif source == 'PRNU':
                data = self.PRNU_array
            elif source == 'DNFP':
                data = self.DNFP_array
            elif source == 'QN':
                data = self.QN_array
            else:
                data = 0
            return data
        else:
            print("Run Simulation")

    def writeto(self, name, source='Digital', user_source=None):
        """
        Function for downloading a fits file of simulated field image

        Parameters
        ----------
        name : str
              filename, Example : simulation.fits

        Source: str,
            Choose from
                          'Digital' : Final digial image
                          'Charge'  : electrons, Light(Source + sky)
                                          + Dark Current
                                          + Noises
                          'Source'  : Source + Sky + Noises
                          'Sky'     : Sky + shot_noise
                          'DC'      : Dark Current + DNFP
                          'Bias'    : Charge offset
                          'PRNU'    : Photon Response Non-Uniformity
                          'DNFP'    : Dark Noise Fixed Pattern
                          'QN'      : Quantization Noise

        user_source : numpy.ndarray
                      2D numpy array user wants to save as FITS
        """

        if np.all(self.image) is not None:
            if user_source is not None:
                if isinstance(user_source) == np.ndarray:
                    data = user_source
                else:
                    raise Exception("Input is not an array")

            elif source == 'Digital':
                data = self.digital
            elif source == 'Charge':
                data = self.charge
            elif source == 'Source':
                data = self.light_array
            elif source == 'Sky':
                data = self.sky_photoelec
            elif source == 'DC':
                data = self.DC_array
            elif source == 'Bias':
                data = self.bias_array
            elif source == 'PRNU':
                data = self.PRNU_array
            elif source == 'DNFP':
                data = self.DNFP_array
            elif source == 'QN':
                data = self.QN_array

            else:
                print(f"{source} is not a valid source")

            hdu = fits.PrimaryHDU(data, header=self.header)
            hdu.wcs = self.wcs
            hdul = fits.HDUList([hdu])
            hdul.writeto(f'{name}', overwrite=True)
        else:
            print("Run Simulation")
