"""This modules contains classes for simulating Imaging, Mosaicing, and
    Spectroscopic modes"""
from pathlib import Path
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord

import numpy as np

import cv2
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs

from .utils import bandpass
from .analysis import Analyzer

data_path = Path(__file__).parent.joinpath() / 'data/'


class Imager(Analyzer):
    """Imager class uses dataframe containing position and magntidue
    information to simulate image based on user defined telescope
    and detector characteristics
    """
    def __init__(self, df, coords=None, tel_params=None, n_x=1000,
                 n_y=1000, exp_time=100, plot=False, user_profiles=None,
                 **kwargs):
        """
        Parameters
        ----------
        df     : pd.DataFrame,
                 Pandas dataframe with source catalog
        coords : (float, float),
                 (RA, Dec) in degrees
        tel_params : dict,
                     {'aperture'       : float,  cm
                      'pixel_scale'    : float,  arcsecs/pixels
                      'sim_file'       : fits,npy
                      'response_funcs' : list, [filename.dat, n] where n is
                                              number of times to multiply
                                              filter
                                              profile
                      'coeffs'         : float, filter coefficients if not
                                              response_funcs
                      }
        n_x      : int,
                   number of pixels along RA direction
        n_y      : int,
                   number of pixels along Dec direction
        exp_time : float
                   Exposure time in seconds
        """
        super().__init__()

        # Flags
        self.shot_noise = True
        self.QE = True
        self.sky = True
        self.PRNU = True
        self.DC = True
        self.DCNU = True
        self.DNFP = True
        self.QN = True

        # TBD
        self.cosmic_rays = False

        # Telescope and Detector Parameters
        psf_file = f'{data_path}/PSF/INSIST/off_axis_hcipy.npy'
        sky_resp = f'{data_path}/Sky_mag.dat'

        self.tel_params = {'aperture': 100,  # cm
                           'pixel_scale': 0.1,
                           'psf_file': psf_file,
                           'response_funcs': [],
                           'sky_resp': sky_resp,
                           'coeffs': 1,
                           'theta': 0
                           }

        self.user_profiles = {
                              'sky': None,
                              'PRNU': None,
                              'QE': None,
                              'T': None,
                              'DC': None,
                              'DNFP': None,
                              'Bias': None,
                             }
        if user_profiles is not None:
            self.user_profiles.update(user_profiles)
        if tel_params is not None:
            self.tel_params.update(tel_params)
        self.tel_params.update(kwargs)

        self.det_params = {
                          'shot_noise': 'Gaussian',
                          'M_sky':  27,
                          'qe_response':  [],         # Wavelength dependence
                          'qe_mean':  0.95,        # Effective QE
                          'bias':  35,         # electrons
                          'G1':  1,
                          'bit_res':  14,
                          'RN':  5,          # elec/pix
                          'PRNU_frac':  0.25/100,   # PRNU sigma
                          'T':  223,        # K
                          'DFM':  1.424e-2,   # 14.24 pA
                          'pixel_area':  1e-6,       # cm2
                          'DCNU':  0.1/100,    # percentage
                          'DNFP':  0.,          # electrons
                          'NF':  0.,          # electrons
                          'FWC':  1.4e5,      # electrons
                          'C_ray_r':  2/50        # hits/second
                          }

        self.det_params.update(kwargs)

        self.df = df.copy()
        self.n_x = n_x
        self.n_y = n_y
        self.pixel_scale = self.tel_params['pixel_scale']
        self.theta = self.tel_params['theta']*np.pi/180

        self.response_funcs = self.tel_params['response_funcs']
        self.coeffs = self.tel_params['coeffs']

        # DN/electrons
        self.gain = pow(2, self.det_params['bit_res'])/self.det_params['FWC']
        self.gain *= self.det_params['G1']

        self.tel_area = np.pi*(self.tel_params['aperture']/2)**2

        self.exp_time = exp_time  # seconds
        self.sim_run = False

        self.psf_file = self.tel_params['psf_file']

        self.check_df()

        if coords is None:
            self.ra = np.median(self.df['ra'])
            self.dec = np.median(self.df['dec'])
        else:
            self.ra = coords[0]
            self.dec = coords[1]

        ra_n = np.round(self.ra,  3)
        dec_n = np.round(self.dec, 3)
        self.name = f" RA : {ra_n} degrees, Dec : {dec_n} degrees"

        self.generate_sim_field(plot)

        # Init Cosmic Rays
        # area = ((n_x*self.pixel_scale)/3600)*((n_y*self.pixel_scale)/3600)

        # eff_area = (0.17/area)*self.exp_time

        # self.n_cosmic_ray_hits = int(eff_area*self.det_params['C_ray_r'])

    def generate_sim_field(self, plot):
        """This function creates array with FoV a bit wider
        than user defined size for flux conservation"""
        if self.df is not None:
            self.calc_zp(plot=plot)
            self.init_psf_patch()

            # Cropping df to sim_field
            x_left = self.n_pix_sub - 1
            x_right = self.n_x_sim - self.n_pix_sub - 1
            y_left = self.n_pix_sub - 1
            y_right = self.n_y_sim - self.n_pix_sub - 1

            self.sim_df = self.init_df(df=self.df,
                                       n_x=self.n_x_sim, n_y=self.n_y_sim,
                                       x_left=x_left, x_right=x_right,
                                       y_left=y_left, y_right=y_right)
            if len(self.sim_df) < 1:
                print("Not Enough sources inside FoV. Increase n_x\
                                and n_y")
        else:
            print("df cannot be None")

    def check_df(self):
        # Input Dataframe
        if 'mag' not in self.df.keys():
            raise Exception("'mag' column not found input dataframe")

        if 'ra' not in self.df or 'dec' not in self.df.keys():
            if 'x' in self.df.keys() and 'y' in self.df.keys():
                print("Converting xy to ra-dec")
                self.df.x = self.df.x + 2
                self.df.y = self.df.y + 2
                self.df = self.xy_to_radec(self.df, self.n_x, self.n_y,
                                           self.pixel_scale)
            else:
                raise Exception("'ra','dec','x',or 'y', \
                 columns not found in input dataframe ")

    def calc_zp(self, plot=False):
        if len(self.response_funcs) > 0:
            wav = np.linspace(1000, 10000, 10000)
            flux = 3631/(3.34e4*wav**2)   # AB flux

            fig, ax, _, params = bandpass(wav, flux, self.response_funcs,
                                          plot=plot)

            lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

            self.lambda_phot = lambda_phot
            self.int_flux = int_flux
            self.W_eff = W_eff
            self.int_flux_Jy = int_flux_Jy
            self.flux_ratio = flux_ratio

            filt_dat = np.loadtxt(self.tel_params['sky_resp'])
            
            wav = filt_dat[:, 0]
            flux = filt_dat[:, 1]

            _, _, _, params = bandpass(wav, flux, self.response_funcs,
                                       plot=False)

            int_flux = params[1]
            self.det_params['M_sky'] = int_flux

        else:

            print("Response functions not provided. Using default values")
            self.int_flux_Jy = 3631
            self.W_eff = 1000
            self.lambda_phot = 2250
            self.flux_ratio = 1

        self.photons = 1.51e3*self.int_flux_Jy*(self.W_eff/self.lambda_phot)
        self.photons *= self.flux_ratio

        self.zero_flux = self.exp_time*self.tel_area*self.photons
        self.zero_flux *= self.coeffs
        self.M_sky_p = self.det_params['M_sky'] \
            - 2.5*np.log10(self.pixel_scale**2)

    def init_psf_patch(self, return_psf=False):
        """Creates PSF array from NPY or fits files"""
        ext = self.psf_file.split('.')[-1]

        if ext == 'npy':
            image = np.load(self.psf_file)
        elif ext == 'fits':
            image = fits.open(self.psf_file)[0].data

        image /= image.sum()  # Flux normalized to 1
        self.image_g_sub = image
        self.sky_bag_flux = self.zero_flux*pow(10, -0.4*self.M_sky_p)

        self.n_pix_sub = self.image_g_sub.shape[0]

        # Defining shape of simulation field
        self.n_x_sim = self.n_x + 2*(self.n_pix_sub-1)
        self.n_y_sim = self.n_y + 2*(self.n_pix_sub-1)

        if return_psf:
            return image*self.zero_flux

    def init_df(self, df, n_x, n_y, x_left, x_right, y_left, y_right):
        """Bounds sources to boundary defined by x and y limits"""
        wcs = self.create_wcs(n_x, n_y, self.ra, self.dec,
                              self.pixel_scale, self.theta)

        c = SkyCoord(df['ra'], df['dec'], unit=u.deg)
        pix = wcs.world_to_array_index(c)
        df['x'] = pix[1]
        df['y'] = pix[0]

        # Cropping Dataframe based on FoV
        x_min_cut = (df['x'] > x_left)
        x_max_cut = (df['x'] < x_right)

        df = df[x_min_cut & x_max_cut]

        y_min_cut = (df['y'] > y_left)
        y_max_cut = (df['y'] < y_right)

        df = df[y_min_cut & y_max_cut]

        return df

    def init_image_array(self, return_img=False):
        """
        Creates a base image array for adding photons

        Parameters
        ----------
        return_img : bool, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        numpy.ndarray
            if return_img is true return base image array

        """
        self.image = np.zeros((self.n_y_sim, self.n_x_sim))
        self.wcs = self.create_wcs(self.n_x_sim, self.n_y_sim,
                                   self.ra, self.dec, self.pixel_scale,
                                   self.theta)
        if return_img:
            return self.image, self.wcs

    def xy_to_radec(self, df, n_x, n_y, pixel_scale):

        w = WCS(naxis=2)
        w.wcs.crpix = [n_x//2, n_y//2]
        w.wcs.cdelt = np.array([-pixel_scale/3600, pixel_scale/3600])
        w.wcs.crval = [10, 10]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        pos = np.array([df.x, df.y])
        coords = np.array(w.pixel_to_world_values(pos.T))
        df['ra'] = coords[:, 0]
        df['dec'] = coords[:, 1]

        return df

    def create_wcs(self, n_x, n_y, ra, dec, pixel_scale, theta=0):
        """
        Parameters
        ----------
        n_x : int
            number of pixels in RA direction
        n_y : int
            number of pixels in Dec direction
        ra : float (degrees)
            right ascension of center of image.
        dec : float (degrees)
            declination of center of image.
        pixel_scale : floats
            arcsecs/pixel.

        Returns
        -------
        w : wcs object

        """
        w = WCS(naxis=2)
        w.wcs.crpix = [(n_x)//2, (n_y)//2]
        w.wcs.cdelt = np.array([-pixel_scale/3600, self.pixel_scale/3600])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.pc = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta),  np.cos(theta)]])
        return w

    def compute_shot_noise(self, array, type_='Poisson'):
        """
        Parameters
        ----------
        array : numpy.ndarray
            input array
        type_ : str, optional
             The default is 'Poisson'.
        Returns
        -------
        shot_noise : numpy.ndarray
             Return array with shot noise
        """
        if type(array) == np.float64 or type(array) == float:
            n_x = self.n_y
            n_y = self.n_x
        else:
            n_x = array.shape[0]
            n_y = array.shape[1]
        if type_ == 'Gaussian':
            shot_noise = np.random.normal(loc=array, scale=np.sqrt(array),
                                          size=(n_x, n_y))
        elif type_ == 'Poisson':
            shot_noise = np.random.poisson(lam=array,
                                           size=(n_x, n_y)
                                           ).astype(np.float64)
        else:
            print('Invalid type')
        return shot_noise

    def compute_coeff_arrays(self):
        """

        Computed coefficients based on input parameters
        Returns
        -------
        None.
        """
        n_x = self.n_y
        n_y = self.n_x

        if self.sky:
            if self.user_profiles['sky'] is not None:
                if self.user_profiles['sky'].shape == (n_x, n_y):
                    self.sky_photons = self.user_profiles['sky']
                else:
                    raise Exception(f"""User defined sky array shape: \
                    {self.user_profiles['sky'].shape} \
                    is not same as detector shape {(n_x,n_y)}""")
            else:
                self.sky_photons = self.compute_shot_noise(self.sky_bag_flux)
        else:
            self.sky_photons = 0

        if self.QE:
            if len(self.det_params['qe_response']) != 0:
                wav = np.linspace(1000, 10000, 10000)
                flux = 3631/(3.34e4*wav**2)   # AB flux

                _, _, _, params = bandpass(wav, flux,
                                           self.det_params['qe_response'],
                                           plot=False)

                _, _, _, _, flux_ratio = params
                self.det_params['qe_mean'] = flux_ratio
        else:
            self.det_params['qe_mean'] = 1

        if self.user_profiles['Bias'] is not None:
            if self.user_profiles['Bias'].shape == (n_x, n_y):
                self.bias_array = self.user_profiles['Bias']
            else:
                raise Exception(f"""User defined Bias array shape: \
                {self.user_profiles['Bias'].shape} \
                is not same as detector shape {(n_x,n_y)}""")
        else:
            self.bias_array = np.random.normal(loc=self.det_params['bias'],
                                               scale=self.det_params['RN'],
                                               size=(n_x, n_y))
        if self.PRNU:
            if self.user_profiles['PRNU'] is not None:
                if self.user_profiles['PRNU'].shape == (n_x, n_y):
                    self.PRNU_array = self.user_profiles['PRNU']
                else:
                    raise Exception(f"""User defined PRNU array shape: \
                    {self.user_profiles['PRNU'].shape} \
                    is not same as detector shape {(n_x,n_y)}""")
            else:
                scale = self.det_params['PRNU_frac']
                self.PRNU_array = np.random.normal(loc=0,
                                                   scale=scale,
                                                   size=(n_x, n_y))
        else:
            self.PRNU_array = 0

        if self.DC:
            if self.user_profiles['DC'] is not None:
                if self.user_profiles['DC'].shape == (n_x, n_y):
                    self.DR = self.user_profiles['DC']
                else:
                    raise Exception(f"""User defined DC array shape:
                    {self.user_profiles['DC'].shape}
                    is not same as detector shape {(n_x,n_y)}""")
            else:
                if self.user_profiles['T'] is not None:
                    if self.user_profiles['T'].shape == (n_x, n_y):
                        area = self.det_params['pixel_area']
                        self.DR = self.dark_current(self.user_profiles['T'],
                                                    self.det_params['DFM'],
                                                    area)
                    else:
                        raise Exception(f"""User defined DC array shape:
                        {self.user_profiles['DC'].shape}
                        is not same as detector shape {(n_x,n_y)}""")
                else:
                    self.DR = self.dark_current(self.det_params['T'],
                                                self.det_params['DFM'],
                                                self.det_params['pixel_area'])
            # Dark Current Non-uniformity
            if self.DCNU:
                sigma = self.exp_time*self.DR*self.det_params['DCNU']
                self.DCNU_array = np.random.lognormal(mean=0,
                                                      sigma=sigma,
                                                      size=(n_x, n_y))
                self.DR *= self.DCNU_array
            self.DC_array = self.compute_shot_noise(self.DR*self.exp_time)
            # Dark Current Fixed Pattern
            if self.DNFP:
                if self.user_profiles['DNFP'] is not None:
                    if self.user_profiles['DNFP'].shape == (n_x, n_y):
                        self.DNFP_array = self.user_profiles['DNFP']
                    else:
                        raise Exception(f"""User defined DNFP array shape:
                                            {self.user_profiles['DNFP'].shape}
                                is not same as detector shape {(n_x,n_y)}""")
                else:
                    arr = self.compute_shot_noise(self.det_params['DNFP'])
                    self.DNFP_array = arr
                self.DC_array += self.DNFP_array
        else:
            self.DC_array = 0
        # Quantization Noise
        if self.QN:
            # electrons
            A = self.det_params['FWC']
            B = pow(2, self.det_params['bit_res'])*np.sqrt(12)
            self.QN_value = (A/B)
            self.QN_array = self.QN_value*np.random.randint(-1, 2,
                                                            size=(n_x, n_y))
        else:
            self.QN_array = 0

    def dark_current(self, T, DFM, pixel_area):
        """
          Parameters
          ----------
          T : float
              Detector Temperature
          DFM : float
              Dark current figure of merit
          pixel_area : float
              Area of pixel

          Returns
          -------
          DR : float
              Dark current rate e/s/pixels
        """
        Kb = 8.62e-5
        const = 2.55741439581387e15

        EgT = 1.1557 - (7.021e-4*T**2/(1108+T))
        DR = const*pixel_area*(T**1.5)*DFM*np.exp(-EgT/(2*Kb*T))
        return DR

    def generate_photons(self, image, patch_width, df):
        """
          This function creates sims based on ABmag  on a
          small patch (2D array) of size n_pix_s*n_pix_s.

          The patch with the sim is then added to the image array of size
          n_pix_m*n_pix_m using wcs object.

          Parameters
          ----------
          image       : numpy.ndarray
                        base image array for inserting star sims
          patch_width : int
                        number of pixels (length) in sim patch image
          df          : pandas.dataframe
                        Dataframe containing source list
          Returns
          -------
          image : numpy.ndarray
              Array with sims added based on df

        """
        patch_width_mid = patch_width//2

        x0, y0 = df['x'].astype(int), df['y'].astype(int)
        ABmag = df['mag'].values
        flux = self.zero_flux * 10**(-ABmag/2.5)
        patch = self.image_g_sub

        x1 = x0 - patch_width_mid
        x2 = x1 + patch_width
        y1 = y0 - patch_width_mid
        y2 = y1 + patch_width

        for x1_, x2_, y1_, y2_, flux_ in zip(x1, x2, y1, y2, flux):
            image[y1_:y2_, x1_:x2_] += flux_*patch

        image = image[patch_width-1:-patch_width+1,
                      patch_width-1:-patch_width+1]

        return image

    def __call__(self, det_params=None, n_stack=1, stack_type='median',
                 photometry='Aper', fwhm=3, sigma=3, detect_sources=False,
                 zero_point=None, **kwargs):
        """
          Parameters
          ----------
        det_params: dict, optional
          Dictionary contianing detector parameters. The default is None.
                    {     'shot_noise' :  str,
                          'M_sky'      :  float,
                          'qe_mean'    :  float,  photons to photoelectrons
                          'bias'       :  float,       electrons
                          'G1'         :  float,
                          'bit_res'    :  int,
                          'RN'         :  float,       elec/pix
                          'PRNU_frac'  :  float,       PRNU sigma
                          'T'          :  float,       K
                          'DFM'        :  float,       pA
                          'pixel_area' :  float,
                          'DCNU'       :  float        fraction
                          'DNFP'       :  float,       electrons
                          'NF'         :  float,       electrons
                          'FWC'        :  float,       electrons
                          'C_ray_r'    :  float        hits/second
                      }
          n_stack    : int, optional
                      Number of observations to be stacked. The default is 1.

          stack_type : str, optional
                      Stacking method. The default is 'median'.
        photometry : str,
                      Type of photometry to be employed
                      Choose from
                      'Aper' : Aperture photometry using Photutils
                      'PSF'  : PSF photometry using DAOPHOT
                      None   : Simulate without photometry
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
                Default None, zero point is calculated theoretically or using
                input catalog
        Simulates field by taking dataframe as input by inserting sim patches
        in image array using WCS
          Returns
          -------
          numpy.ndarray
          Final image array after adding all layers of simulation

        """

        if det_params is not None:

            self.det_params.update(det_params)
            A = pow(2, self.det_params['bit_res'])
            B = self.det_params['FWC']
            self.gain = A/B

            self.gain *= self.det_params['G1']

        digital_stack = []

        for i in range(n_stack):

            self.init_image_array()
            self.compute_coeff_arrays()

            # Source photons
            self.source_photons = self.generate_photons(self.image,
                                                        self.n_pix_sub,
                                                        self.sim_df)
            # Sky photons added to source
            self.light_array = (self.source_photons + self.sky_photons)

            # Source shot_noise
            if self.shot_noise:
                type_ = self.det_params['shot_noise']
                self.light_array = self.compute_shot_noise(self.light_array,
                                                           type_=type_)
            # QE pixel to pixel variation | Source photoelectrons
            self.light_array = self.light_array*self.det_params['qe_mean']

            # Photon Response (Quantum Efficiency) Non Uniformity
            if self.PRNU:
                self.light_array *= (1+self.PRNU_array)

            # Dark Current. Includes DCNU, DNFP and shot noise
            self.photoelec_array = self.light_array + self.DC_array

            # Addition of Quantization error, Bias and Noise floor
            self.charge = self.photoelec_array + self.QN_array \
                                               + self.det_params['NF'] \
                                               + self.bias_array
            # Photoelec to ADUs
            self.digital = (self.charge*self.gain).astype(int)

            # Full well condition
            br = pow(2, self.det_params['bit_res'])
            self.digital = np.where(self.digital >= br,
                                    pow(2, self.det_params['bit_res']),
                                    self.digital)
            digital_stack.append(self.digital)

        digital_stack = np.array(digital_stack)
        if n_stack > 1:
            if stack_type == 'median':
                self.digital = np.median(digital_stack, axis=0)
            elif stack_type == 'mean':
                self.digital = np.median(digital_stack, axis=0)
        if self.cosmic_rays:
            for i in range(self.n_cosmic_ray_hits):
                x = np.random.randint(0, self.n_x_main)
                y = np.random.randint(0, self.n_y_main)
                self.digital[x, y] = pow(2, self.det_params['bit_res'])

        self.wcs = self.create_wcs(self.n_x, self.n_y,
                                   self.ra, self.dec,
                                   self.pixel_scale, self.theta)

        self.sim_flag = True
        # Filtering out sources within Image
        x_left = 0
        x_right = self.n_x
        y_left = 0
        y_right = self.n_y

        self.img_df = self.init_df(df=self.sim_df.copy(),
                                   n_x=self.n_x, n_y=self.n_y,
                                   x_left=x_left, x_right=x_right,
                                   y_left=y_left, y_right=y_right)

        self.header = self.wcs.to_header()
        self.header['gain'] = self.det_params['G1']
        self.header['Temp'] = str(self.det_params['T']) + 'K'
        self.header['bias'] = self.det_params['bias']
        self.header['RN'] = self.det_params['RN']
        self.header['DR'] = np.mean(self.DR)
        self.header['NF'] = self.det_params['NF']

        super().__call__(df=self.img_df, wcs=self.wcs,
                         data=self.digital.astype(float),
                         photometry=photometry,
                         fwhm=fwhm, sigma=sigma,
                         detect_sources=detect_sources,
                         zero_point=zero_point)

    def add_distortion(self, xmap, ymap):
        self.x_map = xmap
        self.y_map = ymap
        # Interpolation to be added
        data = self.digital.astype(float).copy()
        self.org_digital = data
        distorted_img = cv2.remap(data, xmap.astype(np.float32),
                                  ymap.astype(np.float32), cv2.INTER_LANCZOS4)
        distorted_img = distorted_img.astype(int)
        self.digital = np.where(distorted_img > 0, distorted_img, 1)

    def remove_distortion(self):
        # undistort to be added
        self.digital = self.org_digital


class Mosaic(Imager):
    """
    A class to split bigger images to tiles and stitch them together

    """

    def __init__(self, df=None, coords=None, ras=None, decs=None,
                 tel_params=None, exp_time=100,
                 n_x=1000, n_y=1000, mos_n=1, mos_m=1, **kwargs):

        """
        Analyzer.init()

        Parameters
        ----------
        mos_n : int,
                number of tiles in RA direction

        mos_m : int,
                number of tiles in Dec direction
        """

        super().__init__(df=df, coords=coords, exp_time=exp_time, n_x=n_x,
                         n_y=n_y, tel_params=tel_params, **kwargs)
        self.n = mos_n
        self.m = mos_m

        if ras is None or decs is None:
            self.mosaic_ra = self.ra
            self.mosaic_dec = self.dec
            self.mosaic_n_x = n_x
            self.mosaic_n_y = n_y
            self.mosaic_df = df
            self.mosaic_wcs = self.create_wcs(self.mosaic_n_x, self.mosaic_n_y,
                                              self.mosaic_ra, self.mosaic_dec,
                                              self.pixel_scale)

            self.df_split(df, self.n, self.m, n_x, n_y, self.mosaic_wcs)
        else:
            self.ras = ras
            self.decs = decs

    def df_split(self, df, n, m, n_x, n_y, wcs):

        """
        Function to split dataframe based shape and number of tiles

        Parameters
        ----------
        n : int,
            number of tiles in RA direction

        m : int,
            number of tiles in Dec direction

        n_x : int,
            number of pixels in RA direction

        n_y : int,
            number of pixels in Dec direction
        """

        x_bins = np.linspace(0, n_x, n+1)
        y_bins = np.linspace(0, n_y, m+1)

        x_cens = 0.5*(x_bins[:-1] + x_bins[1:]) - 1
        y_cens = 0.5*(y_bins[:-1] + y_bins[1:]) - 1
        cens = wcs.array_index_to_world_values(y_cens, x_cens)
        ra_cens = cens[0]
        ra_cens = np.where(np.round(ra_cens, 1) == 360, ra_cens - 360, ra_cens)
        dec_cens = cens[1]

        self.ras = ra_cens
        self.decs = dec_cens
        self.filenames = []

    def __call__(self, det_params=None, n_stack=1, stack_type='median',
                 photometry=True, fwhm=None):
        """
          Imager.call()

          Calls the Imager class iteratively to generate image tiles
          and stitches the tiles together

        """
        # Flags
        n_x = self.n_x//self.n + 50
        n_y = self.n_y//self.m + 50

        for i in range(self.n):
            for j in range(self.m):

                df = self.mosaic_df
                coords = (self.ras[i], self.decs[j])
                exp_time = self.exp_time
                tel_params = self.tel_params
                det_params['T'] += np.random.randint(-3, 3)
                np.random.seed(i+2*j)
                super().__init__(df=df, coords=coords, exp_time=exp_time,
                                 n_x=n_x,
                                 n_y=n_y, tel_params=tel_params)
                np.random.seed(i+2*j)
                super().__call__(det_params=det_params, n_stack=n_stack,
                                 stack_type=stack_type, photometry=None)

                self.writeto(f'{i}{j}_mosaic.fits')
                self.filenames.append(f'{i}{j}_mosaic.fits')

    def make_mosaic(self):
        hdus = []
        for fil in self.filenames:
            hdu = fits.open(fil)[0]
            hdus.append(hdu)

        wcs_out, shape_out = find_optimal_celestial_wcs(hdus, frame='icrs',
                                                        auto_rotate=True)

        array, fp = reproject_and_coadd(hdus,
                                        wcs_out, shape_out=(shape_out),
                                        reproject_function=reproject_interp)

        self.wcs = wcs_out
        self.digital = array
        self.footprint = fp
