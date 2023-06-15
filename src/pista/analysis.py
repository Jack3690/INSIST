"""This module contains classes that can be utilized for visualizing and
    analyzing the simulated images and spectra"""
from photutils import aperture as aper
from matplotlib import colors as col
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from photutils.detection import DAOStarFinder
from photutils.psf import DAOPhotPSFPhotometry, FittableImageModel

from astropy.io import fits
from astropy.io.fits import CompImageHDU
from astropy.modeling.fitting import LevMarLSQFitter

from astropy.stats import sigma_clipped_stats
from .utils import Xmatch

import numpy as np


class Analyzer(object):
    def __init__(self):
        """
        A class to visualize and analyze the simulated image

        Parameters
        ----------


        Returns
        -------
        None.

        """
    def __call__(self, df=None, wcs=None, data=None,
                 photometry=None, detect_sources=False, fwhm=3, sigma=13,
                 ZP=None):
        """
        Performs sim simulation and sim Photometry

        Imager.call()

        do_photometry : Bool, Default : True
                        Do Aperture Photometry
        """
        self.photometry_type = photometry
        if photometry == 'Aper':
            self.aper_photometry(data, wcs, df, fwhm, sigma, detect_sources,
                                 ZP)
        elif photometry == 'PSF':
            self.psf_photometry(data, wcs, df, fwhm, sigma, ZP)

    def aper_photometry(self, data, wcs, df, fwhm, sigma, detect, ZP):
        """
        Function to perform Aperture photometry

        Parameters
        ----------

        data: np.ndarray,
            image to perform photometry on

        wcs: astropy.wcs.WCS
            WCS object of the image

        df: pandas.DataFrame,
            Source catalog of source in the image from simulation for
            reference

        fwhm : float, pixels
                    During aperture photometry,
                    fwhm corresponds to FWHM circular aperture for
                    aperture photometry
                    During PSF photometry,
                    fwhm corresponds FWHM kernel to use for PSF photometry
        sigma: float,
                The numbers of standard deviations above which source has to be
                detected
        detect: bool,
                If true, DARStarFinder is used to detect sources for aperture
                photometry

                if false, input catalog is used for getting positions
                of sources for aperture photometry
        ZP    : float,
                zero point of the telescope.

        Returns
        -------

        phot_table: astropy.table.Table
                    table containing photometry of the souces

                    Columns
                    'x-centeroid'
                    'y-centeroid'
                    'sky'
                    'flux'
                    'mag_in'
                    'mag_out'
                    'mag_err'
                    'SNR'

        """

        # if detect flag is set to True, detect sources in the image
        if detect:

            # Calculate the mean, median and standard deviation of the data
            mean, median, std = sigma_clipped_stats(data, sigma=sigma)
            # Create DAOStarFinder object to detect sources
            daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma*std)
            # Detect sources in the image
            sources = daofind(data-median)
            # Get the source positions
            positions = np.transpose((sources['xcentroid'],
                                      sources['ycentroid']))

        else:

            # create SkyCoord object from ra and dec values in the dataframe
            coords = np.array([df['ra'], df['dec']])
            # convert the sky coordinates to pixel coordinates
            pix = wcs.world_to_pixel_values(coords.T)
            positions = np.array(pix)

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
        phot_table['flux_err'] = np.sqrt(phot_table['flux'].value +
                                         phot_table['sky_flux'].value)

        # calculate signal to noise ratio
        phot_table['SNR'] = phot_table['flux']/phot_table['flux_err']

        if not detect:

            phot_table['ra'] = df['ra'].values
            phot_table['dec'] = df['dec'].values
            phot_table['mag_in'] = df['mag'].values

        else:
            coords = np.array(wcs.pixel_to_world_values(positions))

            phot_table['ra'] = coords[:, 0]
            phot_table['dec'] = coords[:, 1]

            matched = Xmatch(phot_table, df, 3*fwhm)
            mag_in = []
            sep = []
            for i, j, k in matched:
                if j is not None:
                    mag_in.append(df['mag'].values[j])
                    sep.append(k)
                else:
                    mag_in.append(np.nan)
                    sep.append(np.nan)

            phot_table['mag_in'] = mag_in
            phot_table['sep'] = sep

        phot_table['mag_out'] = -2.5*np.log10(phot_table['flux']) + ZP
        phot_table['mag_err'] = 1.082/phot_table['SNR']

        coords = np.array(wcs.pixel_to_world_values(positions))
        phot_table['ra'] = coords[:, 0]
        phot_table['dec'] = coords[:, 1]
        self.phot_table = phot_table

        def psf_photometry(self, data, wcs, df, fwhm, sigma, ZP):
            """
            Function to perform PSF photometry

            Parameters
            ----------

            data: np.ndarray,
                image to perform photometry on

            wcs: astropy.wcs.WCS
                WCS object of the image

            df: pandas.DataFrame,
                Source catalog of source in the image from simulation for
                reference

            fwhm : float, pixels
                        During aperture photometry,
                        fwhm corresponds to FWHM circular aperture for
                        aperture photometry
                        During PSF photometry,
                        fwhm corresponds FWHM kernel to use for PSF photometry
            sigma: float,
                    The numbers of standard deviations above which source
                      has to be detected
            detect: bool,
                    If true, DARStarFinder is used to detect sources
                    for aperture photometry

                    if false, input catalog is used for getting positions
                    of sources for aperture photometry
            ZP    : float,
                    zero point of the telescope.

            Returns
            -------

            phot_table: astropy.table.Table
                        table containing photometry of the souces

                        Columns
                        'x-centeroid'
                        'y-centeroid'
                        'sky'
                        'flux'
                        'mag_in'
                        'mag_out'
                        'mag_err'
                        'SNR'
            """
            mean, median, std = sigma_clipped_stats(data, sigma=3)

            psf_model = FittableImageModel(self.psf)
            self.psf_model = psf_model

            photometry = DAOPhotPSFPhotometry(crit_separation=2,
                                              threshold=mean + sigma*std,
                                              fwhm=fwhm,
                                              aperture_radius=3,
                                              psf_model=psf_model,
                                              fitter=LevMarLSQFitter(),
                                              fitshape=(11, 11),
                                              niters=3)

            phot_table = photometry(image=data)
            positions = np.array([phot_table['x_fit'], phot_table['y_fit']]).T
            coords = np.array(wcs.pixel_to_world_values(positions))

            phot_table['ra'] = coords[:, 0]
            phot_table['dec'] = coords[:, 1]
            phot_table['SNR'] = phot_table['flux_fit']/phot_table['flux_unc']

            phot_table['mag_out'] = -2.5*np.log10(phot_table['flux_fit'])
            phot_table['mag_out'] += ZP
            phot_table['mag_err'] = 1.082/phot_table['SNR']

            matched = Xmatch(phot_table, df, 3*fwhm)
            mag_in = []
            sep = []
            for i, j, k in matched:
                if j is not None:
                    mag_in.append(df['mag'].values[j])
                    sep.append(k)
                else:
                    mag_in.append(np.nan)
                    sep.append(np.nan)

            phot_table['mag_in'] = mag_in
            phot_table['sep'] = sep

            self.phot_table = phot_table

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
        scale = self.n_x/self.n_y
        figsize = figsize[0]*scale, figsize[1]

        # Cropping Dataframe based on FoV
        left = (self.n_x_sim - self.n_x)//2
        right = left + self.n_x

        df = self.sim_df

        x_min_cut = (df['x'] > left)
        x_max_cut = (df['x'] < right)

        df = df[x_min_cut & x_max_cut]

        bottom = (self.n_y_sim - self.n_y)//2
        top = bottom + self.n_y

        y_min_cut = (df['y'] > bottom)
        y_max_cut = (df['y'] < top)

        fov_x = (self.n_x*self.pixel_scale)/3600
        fov_y = (self.n_y*self.pixel_scale)/3600

        fov_x = np.round(fov_x, 4)
        fov_y = np.round(fov_y, 4)

        df = df[y_min_cut & y_max_cut]

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        if hasattr(self, 'spec_bins'):
            x = df['x']
            y = df['y']
            c = None
            cmap = None
            color = None

            if 'mag' in df.keys():
                c = df['mag']
                cmap = 'jet'
                color = None

            img = ax.scatter(x, y, c=c, color=color, cmap=cmap,
                             marker=marker)
            cb = plt.colorbar(img, fraction=1/scale, pad=0.01)
            cb.set_label('mag (ABmag)')

            ax.set_xlabel('x (pix)')
            ax.set_ylabel('y (pix)')

            if hasattr(self, 'L'):
                B = self.B
                L = self.L
                lw = len(self.spec_bins)
                PA = self.PA

                delta = np.arctan(B/L)
                d = np.sqrt(L**2 + B**2)/2
                om = PA - delta

                x_corr = self.n_x_sim//2 - d*np.cos(om) - lw/2
                y_corr = self.n_y_sim//2 - d*np.sin(om) - B*np.cos(PA)

                start = B*np.sin(PA) + lw/2 + x_corr
                end = start + L*np.cos(PA)
                x = np.linspace(start, end, 100)
                y = np.tan(PA)*(x - B*np.sin(PA) - lw/2 - x_corr) + y_corr

                ax.plot(x, y, color='red')

                start = lw/2 + x_corr
                end = start + L*np.cos(PA)
                x = np.linspace(start, end, 100)
                y = np.tan(PA)*(x - lw/2 - x_corr) + B*np.cos(PA) + y_corr

                ax.plot(x, y, color='red')

                start = y_corr
                end = start + B*np.cos(PA)
                y = np.linspace(start, end, 100)
                x = -np.tan(PA)*(y - B*np.cos(PA) - y_corr) + lw/2 + x_corr

                ax.plot(x, y, color='red')

                start = y_corr + L*np.sin(PA)
                end = start + B*np.cos(PA)
                y = np.linspace(start, end, 100)
                x = -np.tan(PA)*(y - B*np.cos(PA) - L*np.sin(PA) - y_corr)\
                    + L*np.cos(PA) + lw/2 + x_corr

                ax.plot(x, y, color='red')
                ax.set_xlim(0, self.n_x_sim)
                ax.set_ylim(0, self.n_y_sim)

        else:

            x = df['ra']
            y = df['dec']
            c = df['mag']

            img = ax.scatter(x, y, c=c, marker=marker, cmap=cmap)
            cb = plt.colorbar(img)
            cb.set_label('mag (ABmag)')

            ax.invert_xaxis()
            ax.set_xlabel('RA (Degrees)')
            ax.set_ylabel('Dec (Degrees)')
            ax.invert_xaxis()
            ax.set_xlim(self.ra+fov_x/2, self.ra-fov_x/2)
            ax.set_ylim(self.dec-fov_y/2, self.dec+fov_y/2)

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

            ax.tick_params(which='both', width=2, direction="in",
                           top=True, right=True,
                           bottom=True, left=True)

        title1 = f'Requested Center : {self.name} {len(df)} sources'
        title2 = f'Fov(RA) : {fov_x} (deg) | Fov(Dec) : {fov_y} (deg)'
        ax.set_title(f"""{title1}
                            {title2}""")

        return fig, ax

    def show_image(self, source='Digital', fig=None, ax=None, cmap='jet',
                   figsize=(15, 10), download=False, show_wcs=True,
                   overlay_apertures=False):
        """
        Function for plotting the simulated field image

        Source: str,
                Choose from
                            'Digital' : Final digial image
                            'Charge'  : electrons, Light(Source + sky) +
                                        Dark Current + Noises
                            'Source'  : Source + Sky + Noises
                            'Sky'     : Sky + shot_noise
                            'DC'      : Dark Current + DNFP
                            'QE'      : Quantum efficiency fluctuation
                                            across
                                        detector
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
        if hasattr(self, 'digital'):
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

            if hasattr(self, 'L'):
                B = self.B
                L = self.L
                lw = len(self.spec_bins)
                PA = self.PA

                delta = np.arctan(B/L)
                d = np.sqrt(L**2 + B**2)/2
                om = PA - delta

                x_corr = self.n_x//2 - d*np.cos(om) - lw/2
                y_corr = self.n_y//2 - d*np.sin(om) - B*np.cos(PA)

                start = B*np.sin(PA) + lw/2 + x_corr
                end = start + L*np.cos(PA)
                x = np.linspace(start, end, 100)
                y = np.tan(PA)*(x - B*np.sin(PA) - lw/2 - x_corr) + y_corr

                ax.plot(x, y, color='red')

                start = lw/2 + x_corr
                end = start + L*np.cos(PA)
                x = np.linspace(start, end, 100)
                y = np.tan(PA)*(x - lw/2 - x_corr) + B*np.cos(PA) + y_corr

                ax.plot(x, y, color='red')

                start = y_corr
                end = start + B*np.cos(PA)
                y = np.linspace(start, end, 100)
                x = -np.tan(PA)*(y - B*np.cos(PA) - y_corr) + lw/2 + x_corr

                ax.plot(x, y, color='red')

                start = y_corr + L*np.sin(PA)
                end = start + B*np.cos(PA)
                y = np.linspace(start, end, 100)
                x = -np.tan(PA)*(y - B*np.cos(PA) - L*np.sin(PA) - y_corr)\
                    + L*np.cos(PA) + lw/2 + x_corr

                ax.plot(x, y, color='red')

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
                        'QE'      : Quantum efficiency fluctuation across
                                    detector
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

        if hasattr(self, 'digital'):
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

    def getImage(self, source='Digital'):
        """
        Function of retrieving image array at different
        stages of simulation.

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
                        'QE'      : Quantum efficiency fluctuation across
                                    detector
                        'Bias'    : Charge offset
                        'PRNU'    : Photon Response Non-Uniformity
                        'DNFP'    : Dark Noise Fixed Pattern
                        'QN'      : Quantization Noise
        """

        if hasattr(self, 'digital'):
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

    def writeto(self, name, source='Digital', user_source=None,
                with_dark_flat=False):

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
                                    + Dark Current + Noises
                        'Source'  : Source + Sky + Noises
                        'Sky'     : Sky + shot_noise
                        'DC'      : Dark Current + DNFP
                        'Bias'    : Charge offset
                        'PRNU'    : Photon Response Non-Uniformity
                        'DNFP'    : Dark Noise Fixed Pattern
                        'QN'      : Quantization Noise

        user_source : numpy.ndarray
                    2D numpy array user wants to save as FITS
        with_dark_flat : bool
                        True : Output fits will provide dark
                        and flat frames
        """
        if hasattr(self, 'digital'):
            c1 = user_source is not None
            c2 = isinstance(user_source, np.ndarray)
            if c1 and c2:
                data = user_source
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
            hdus = [hdu]

            if with_dark_flat:
                dark = self.dark_frame/self.exp_time
                header = self.header
                header['Frame'] = 'Dark'

                hdu = fits.ImageHDU(dark, header=header)
                hdus.append(hdu)

                flat = self.flat_frame
                header['Frame'] = 'Flat'

                hdu = fits.ImageHDU(flat, header=header)
                hdus.append(hdu)

            hdul = fits.HDUList(hdus)
            hdul.writeto(f'{name}', overwrite=True, checksum=True)
        else:
            print("Run Simulation")

    def writecomp(self, name):

        hdu = CompImageHDU(data=self.digital, header=self.header)
        hdus = [fits.PrimaryHDU(), hdu]
        hdul = fits.HDUList(hdus)
        hdul.writeto(name, overwrite=True, checksum=True)
