"""This modules contains classes for simulating Imaging, Mosaicing, and
    Spectroscopic modes"""
from pathlib import Path
from astropy.wcs import WCS
from astropy.io import fits

import numpy as np
import torch
import torch.nn.functional as F
from astropy.table import vstack

from tqdm.contrib import tzip

import cv2
import os
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
from .utils import bandpass, redshift_corr, select_mos, calc_mos_size
from .analysis import Analyzer, SpecAnalyzer

data_path = Path(__file__).parent.joinpath('data')


class Imager(Analyzer):
    """Imager class uses dataframe containing position and magntidue
    information to simulate image based on user defined telescope
    and detector characteristics
    """
    def __init__(self, df, coords=None, tel_params=None, n_x=1000,
                 n_y=1000, exp_time=100, plot=False, user_profiles=None):
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
        self.Bias = True

        # Detector
        self.Cal0 = False

        # Flux
        self.Cal1 = False

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
                           'theta': 0,
                           'M_sky': 27,
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

        self.det_params = {
                          'shot_noise': 'Gaussian',
                          'qe_response':  '',     # Wavelength dependence
                          'qe_mean': None,        # Effective QE
                          'bias': 35,             # electrons
                          'G1': 1,
                          'bit_res': 14,
                          'RN':  5,               # elec/pix
                          'PRNU_frac': 0.25/100,  # PRNU sigma
                          'T': 218,               # K
                          'DFM': 1.424e-2,        # 14.24 pA
                          'pixel_area': 1e-6,     # cm2
                          'DCNU': 0.1/100,        # percentage
                          'DNFP': 0.,             # electrons
                          'NF': 0.,               # electrons
                          'FWC': 1.4e5,           # electrons
                          'C_ray_r': 2/50         # hits/second
                          }


        self.df = df.copy()
        self.n_x = n_x
        self.n_y = n_y
        self.pixel_scale = self.tel_params['pixel_scale']
        self.theta = self.tel_params['theta']*np.pi/180

        self.response_funcs = self.tel_params['response_funcs']
        self.coeffs = self.tel_params['coeffs']

        # electrons/ADUs
        self.gain = self.det_params['FWC']/pow(2, self.det_params['bit_res'])
        self.gain /= self.det_params['G1']

        self.tel_area = np.pi*(self.tel_params['aperture']/2)**2

        self.exp_time = exp_time  # seconds

        self.psf_file = self.tel_params['psf_file']

        self.check_df()

        if coords is None:
            self.ra = np.median(self.df['ra'])
            self.dec = np.median(self.df['dec'])
        else:
            self.ra = coords[0]
            self.dec = coords[1]

        ra_n = np.round(self.ra, 3)
        dec_n = np.round(self.dec, 3)
        self.name = f" RA : {ra_n} degrees, Dec : {dec_n} degrees"

        self.generate_sim_field(plot)

    def generate_sim_field(self, plot):
        """This function creates array with FoV a bit wider
        than user defined size for flux conservation"""
        if self.df is not None:
            self.calc_zp(plot=plot)
            self.init_psf_patch()

            # Cropping df to sim_field
            x_left = self.n_pix_psf//2
            x_right = self.n_x_sim - self.n_pix_psf//2
            y_left = self.n_pix_psf//2
            y_right = self.n_y_sim - self.n_pix_psf//2

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
            self.fig_bp, self.ax_bp = fig, ax
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
        self.M_sky_p = self.tel_params['M_sky'] \
            - 2.5*np.log10(self.pixel_scale**2)
        self.sky_bag_flux = self.zero_flux*pow(10, -0.4*self.M_sky_p)

        if self.sky:
            if self.user_profiles['sky'] is not None:
                if self.user_profiles['sky'].shape == (self.n_x, self.n_y):
                    self.sky_photons = self.user_profiles['sky']
                else:
                    raise Exception(f"""User defined sky array shape: \
                    {self.user_profiles['sky'].shape} \
                    is not same as detector shape {(self.n_x, self.n_y)}""")
            else:
                self.sky_photons = self.compute_shot_noise(self.sky_bag_flux)
        else:
            self.sky_photons = 0

    def init_psf_patch(self, return_psf=False):
        """Creates PSF array from NPY or fits files"""
        ext = self.psf_file.split('.')[-1]

        if ext == 'npy':
            image = np.load(self.psf_file)
        elif ext == 'fits':
            image = fits.open(self.psf_file)[0].data

        if image.min()<0:
          print("PSF data has negative values!\n Modulus of PSF used instead")
          image =  abs(image)

        image /= image.sum()  # Flux normalized to 1

        self.psf = image

        self.n_pix_psf = self.psf.shape[0]

        # Defining shape of simulation field
        self.n_x_sim = self.n_x + 2*(self.n_pix_psf-1)
        self.n_y_sim = self.n_y + 2*(self.n_pix_psf-1)

        if return_psf:
            return image*self.zero_flux

    def init_df(self, df, n_x, n_y, x_left, x_right, y_left, y_right):
        """Bounds sources to boundary defined by x and y limits"""
        wcs = self.create_wcs(n_x, n_y, self.ra, self.dec,
                              self.pixel_scale, self.theta)

        coords = np.array([df['ra'], df['dec']])
        pix = np.array(wcs.world_to_array_index_values(coords.T))
        if len(df)<2:
          raise Exception("Input DataFrame must have atleast 2 sources")

        df['x'] = np.flip(pix[:, 0])
        df['y'] = np.flip(pix[:, 1])

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
        image = np.zeros((self.n_y_sim, self.n_x_sim))
        wcs = self.create_wcs(self.n_x_sim, self.n_y_sim,
                                   self.ra, self.dec, self.pixel_scale,
                                   self.theta)

        return image, wcs

    def xy_to_radec(self, df, n_x, n_y, pixel_scale):

        w = WCS(naxis=2)
        w.wcs.crpix = [n_x//2, n_y//2]
        w.wcs.cdelt = np.array([-pixel_scale/3600, pixel_scale/3600])
        w.wcs.crval = [10, 10]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        pos = np.array([df['x'], df['y']])
        coords = np.array(w.pixel_to_world_values(pos.T))
        df['ra'] = np.flip(coords[:, 0])
        df['dec'] = np.flip(coords[:, 1])

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
        w.wcs.crpix = [n_x//2, n_y//2]
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
        if isinstance(array, np.ndarray):
            n_x = array.shape[0]
            n_y = array.shape[1]
        else:
            n_x = self.n_y
            n_y = self.n_x

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

         # Quantum Efficiency
        if self.QE:
            if self.det_params['qe_mean'] is None and len(self.det_params['qe_response'])>2:
                if os.path.exists(self.det_params['qe_response']):
                    wav = np.linspace(1000, 10000, 10000)
                    flux = 3631/(3.34e4*wav**2)   # AB flux

                    resp = f"{self.det_params['qe_response']},1,100"

                    _, _, _, params = bandpass(wav, flux,
                                                [resp],
                                                plot=False)

                    _, _, _, _, flux_ratio = params
                    self.det_params['qe_mean'] = flux_ratio
                else:
                    raise Exception("Path does not exists!")
            else:
                self.det_params['qe_mean'] = 1
        else:
          self.det_params['qe_mean'] = 1

        # Bias
        if self.Bias:
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
        else:
          self.bias_array = 0

        # Photon Response Not Uniformity
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

        # Dark Current
        self.compute_DC()

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
    def compute_DC(self):
      n_x = self.n_y
      n_y = self.n_x
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
                sigma = self.det_params['DCNU']
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
          self.DR = 0
          self.DC_array = 0

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

    def generate_photons(self, image, patch_width, df, zero_flux):
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
        flux = zero_flux* 10**(-ABmag/2.5)
        patch = self.psf

        x1 = x0 - patch_width_mid
        x2 = x1 + patch_width
        y1 = y0 - patch_width_mid
        y2 = y1 + patch_width

        for x1_, x2_, y1_, y2_, flux_ in tzip(x1, x2, y1, y2, flux):
            image[y1_:y2_, x1_:x2_] += flux_*patch

        image = image[patch_width-1:-patch_width+1,
                      patch_width-1:-patch_width+1]

        return image

    def make_ccd_image(self, light_array):

        # QE pixel to pixel variation | Source photoelectrons
        self.source_photoelec = light_array*self.det_params['qe_mean']

        # Photon Response (Quantum Efficiency) Non Uniformity
        if self.PRNU:
            self.source_photoelec *= (1+self.PRNU_array)

        # Dark Current. Includes DCNU, DNFP and shot noise
        self.photoelec_array = self.source_photoelec + self.DC_array

        # Addition of Quantization error, Bias and Noise floor
        self.charge = self.photoelec_array + self.QN_array \
                                           + self.det_params['NF'] \
                                           + self.bias_array

        # Photoelec to ADUs
        digital = np.round(self.charge/self.gain).astype(int)

        # Full well condition
        digital = np.where(digital >= pow(2, self.det_params['bit_res']),
                           pow(2, self.det_params['bit_res']), digital)
        return digital

    @property
    def bias_frame(self):
      exp_time = self.exp_time
      self.exp_time = 0
      self.compute_DC()
      bias = self.make_ccd_image(0)
      self.exp_time = exp_time
      self.compute_DC()
      return bias

    @property
    def dark_frame(self):
        return (self.make_ccd_image(0) - self.bias_frame)

    @property
    def flat_frame(self):
        flat = self.make_ccd_image(100)
        flat = flat/flat.mean()
        return flat

    def __call__(self, det_params=None, photometry='Aper', fwhm=1.5, sigma=5,
                 detect_sources=False, ZP=None, **kwargs):
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

          Returns
          -------
          numpy.ndarray
          Final image array after adding all layers of simulation

        """

        if det_params is not None:

            self.det_params.update(det_params)
            A = pow(2, self.det_params['bit_res'])
            B = self.det_params['FWC']

            self.gain = B/A

            self.gain /= self.det_params['G1']

        image, _ = self.init_image_array()

        # Source photons
        self.source_photons = self.generate_photons(image,
                                                    self.n_pix_psf,
                                                    self.sim_df,
                                                    self.zero_flux)
        # Sky photons added to source
        self.light_array = (self.source_photons + self.sky_photons)

        # Source shot_noise
        if self.shot_noise:
            type_ = self.det_params['shot_noise']
            self.light_array = self.compute_shot_noise(self.light_array,
                                                        type_=type_)
        # Compute Coefficient arrays
        self.compute_coeff_arrays()

        self.digital = self.make_ccd_image(self.light_array)

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
        self.header['gain'] = (np.mean(self.det_params['G1']), 'e/ADU')
        self.header['Temp'] = (str(np.mean(self.det_params['T'])),'K')
        self.header['bias'] = (np.mean(self.det_params['bias']),'e')
        self.header['RN'] = (np.mean(self.det_params['RN']),'e')
        self.header['DR'] = (np.mean(self.DR),'e/s/pix')
        self.header['NF'] = (np.mean(self.det_params['NF']),'e')
        self.header['EXPTIME'] = (self.exp_time, 'seconds')
        self.header['BUNIT'] = 'ADU'
        self.header['CAL_LEV']= 0

        if ZP is None:
            QE = self.det_params['qe_mean']
            zero_p_flux = self.zero_flux*QE
            ZP = 2.5*np.log10(zero_p_flux)
            self.ZP = ZP

        if self.Cal0:
          # Bias and Dark Correction
          self.digital = self.digital - self.dark_frame - self.bias_frame
          # Flat field correction
          self.digital /= self.flat_frame

          self.header['CAL_LEV']= 1

        self.org_digital = self.digital.astype(float).copy()


        self.header['ZP'] = self.ZP

        super().__call__(df=self.img_df, wcs=self.wcs,
                         data=self.digital.astype(float),
                         photometry=photometry, fwhm=fwhm, sigma=sigma,
                         detect_sources=detect_sources, ZP=ZP)

    def add_stars(self, image_array, zero_flux, df):
        x_size = image_array.shape[1]
        y_size = image_array.shape[0]

        patch_width = self.n_pix_psf
        if self.n_x == x_size and self.n_y == y_size:
          image, _ = self.init_image_array()

          # Source photons
          source_photons = self.generate_photons(image,
                                                 self.n_pix_psf,
                                                 df,
                                                 zero_flux)
          if self.shot_noise:
                type_ = self.det_params['shot_noise']
                source_photons = self.compute_shot_noise(source_photons,
                                                         type_=type_)
        else:
         print(f"Expected shape: ({self.n_x}, {self.n_y}).\n Provided shape: ({x_size}, {y_size})")

        return source_photons + image_array

    def add_distortion(self, xmap, ymap):
        """Function for addition distortion using
        x and y mappings"""
        self.x_map = xmap
        self.y_map = ymap
        # Interpolation to be added
        data = self.digital.astype(float).copy()
        distorted_img = cv2.remap(data, xmap.astype(np.float32),
                                  ymap.astype(np.float32), cv2.INTER_LANCZOS4)
        distorted_img = distorted_img.astype(int)
        self.digital = np.where(distorted_img > 0, distorted_img, 1)

    def remove_distortion(self):
        """Function for returning the image to state
        before adding distortion"""
        # undistort to be added
        self.digital = self.org_digital

    def __del__(self):
      for i in self.__dict__:
        del i


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


class Spectrograph(Imager, Analyzer, SpecAnalyzer):
    """Class for simulating Spectra"""
    def __init__(self, df, coords=None, tel_params={}, n_x=1000,
                 n_y=1000, exp_time=100, plot=False, user_profiles={},
                 **kwargs):

        if torch.cuda.is_available():
            self.device = 'cuda'
        else:
            self.device = 'cpu'
        self.tel_params = tel_params
        bin_min = self.tel_params['lambda1'] - self.tel_params['dellambda']/2
        bin_max = self.tel_params['lambda2'] + self.tel_params['dellambda']
        step = self.tel_params['dellambda']
        self.bins = np.arange(bin_min, bin_max, step)

        self.spec_bins = 0.5*(self.bins[1:] + self.bins[:-1])
        if len(self.spec_bins) % 2 == 0:
            raise Exception("""Spectroscopy bins need to be odd in number
                                Increase lambda2""")

        self.interpolate_mult = 3

        super().__init__(df, coords, tel_params, n_x, n_y, exp_time, plot,
                         user_profiles, **kwargs)
        self.sim_df = redshift_corr(self.sim_df)

        self.QE = False

    def check_df(self):

        if 'flux' not in self.df.keys() or 'wav' not in self.df.keys():
            raise Exception("'flux' or 'wav' column not found input dataframe")

        if 'ra' not in self.df.keys() or 'dec' not in self.df.keys():
            if 'x' in self.df.keys() and 'y' in self.df.keys():
                print("Converting xy to ra-dec")
                self.df = self.xy_to_radec(self.df, self.n_x, self.n_y,
                                           self.pixel_scale)
            else:
                raise Exception("'ra','dec','x',or 'y', \
                 columns not found in input dataframe ")

    def generate_sim_field(self, plot):
        if self.df is not None:
            self.calc_zp(plot=plot)
            self.init_psf_patch()
            spec_corr = len(self.spec_bins)//2

            x_left = spec_corr + self.n_pix_psf
            x_right = self.n_x_sim - spec_corr - self.n_pix_psf
            y_left = self.n_pix_psf
            y_right = self.n_y_sim - self.n_pix_psf - 1

            self.sim_df = self.init_df(df=self.df,
                                       n_x=self.n_x_sim, n_y=self.n_y_sim,
                                       x_left=x_left, x_right=x_right,
                                       y_left=y_left, y_right=y_right)
        else:
            print("df cannot be None")

    def calc_zp(self, plot=False):
        if len(self.response_funcs) > 0:

            wav = np.linspace(1200, 3200, 20000)
            flux = 3631/(3.34e4*wav**2)   # AB flux

            fig, ax, data, params = bandpass(wav, flux, self.response_funcs,
                                             plot=plot)

            lambda_, _, Reff = data
            lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

            self.lambda_phot = lambda_phot
            self.int_flux = int_flux
            self.W_eff = W_eff
            self.int_flux_Jy = int_flux_Jy
            self.flux_ratio = flux_ratio

            self.Reff = self.bin_xy_norm(lambda_, Reff)

            self.flux_coeffs = self.Reff.copy()
            self.flux_coeffs *= self.exp_time*self.coeffs*self.tel_area

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

        self.sky_bag_flux = self.zero_flux*pow(10, -0.4*self.M_sky_p)

        # Unit Conversion
        self.UC = 1.51e3*(self.tel_params['dellambda']/self.spec_bins)
        self.UC *= self.flux_coeffs

        if self.sky:
            if self.user_profiles['sky'] is not None:
                if self.user_profiles['sky'].shape == (self.n_x, self.n_y):
                    self.sky_photons = self.user_profiles['sky']
                else:
                    raise Exception(f"""User defined sky array shape: \
                    {self.user_profiles['sky'].shape} \
                    is not same as detector shape {(self.n_x, self.n_y)}""")
            else:
                self.sky_photons = self.compute_shot_noise(self.sky_bag_flux)
        else:
            self.sky_photons = 0

    def init_MOS_df(self, L, B, PA, lw, df):
        delta = np.arctan(B/L)
        d = np.sqrt(L**2 + B**2)/2
        om = PA - delta

        x_corr = self.n_x_sim//2 - d*np.cos(om) - lw/2
        y_corr = self.n_y_sim//2 - d*np.sin(om) - B*np.cos(PA)

        t = df.copy()
        x0 = B*np.sin(PA) + lw/2 + x_corr
        y0 = y_corr
        t = t[t['y'] > np.tan(PA)*(t['x'] - x0) + y0]

        x0 = lw/2 + x_corr
        y0 = B*np.cos(PA) + y_corr
        t = t[t['y'] < np.tan(PA)*(t['x'] - x0) + y0]

        x0 = lw/2 + x_corr
        y0 = B*np.cos(PA) + y_corr
        t = t[t['x'] > -np.tan(PA)*(t['y'] - y0) + x0]

        x0 = L*np.cos(PA) + lw/2 + x_corr
        y0 = B*np.cos(PA) + L*np.sin(PA) + y_corr
        t = t[t['x'] < -np.tan(PA)*(t['y'] - y0) + x0]

        return t

    def init_MOS(self, L, B, PA):
        PA *= np.pi/180
        lw = len(self.spec_bins)
        x_size, y_size = calc_mos_size(L, B, PA, lw)
        if x_size < self.n_x and y_size < self.n_y:
            self.L = L
            self.B = B
            self.PA = PA

            self.sim_df = self.init_MOS_df(L, B, PA, lw, self.sim_df)
            self.sim_df['objid'] = np.arange(0, len(self.sim_df), 1)
        else:
            print(f"""n_x should be greater than {x_size} \n
                    n_y should be greater than {y_size}""")

    def select_MOS_sources(self, df=None, ids=[], radius=5, min_sep=8):
        if df is None:
            df = self.sim_df

        x = df['x'].value
        y = df['y'].value

        cen_x = x[np.argmin(abs(x-np.mean(x)))]
        cen_y = y[np.argmin(abs(y-np.mean(y)))]

        temp = df[['objid', 'x', 'y']].to_pandas()
        mos_df, res_df, cen_x, cen_y = select_mos(temp,
                                                  cen_x, cen_y, ids, radius,
                                                  min_sep, self.n_y_sim)
        t = []
        for id in mos_df['objid']:
            t.append(df[df['objid'] == id])

        mos_df = vstack(t)
        self.sim_df = mos_df

        t = []
        for id in res_df['objid']:
            t.append(df[df['objid'] == id])

        res_df = vstack(t)

        return mos_df, res_df, cen_x, cen_y

    def init_psf_patch(self, return_psf=False, plot=False):

        ext = self.psf_file.split('.')[-1]

        if ext == 'npy':
            image = np.load(self.psf_file)
        elif ext == 'fits':
            image = fits.open(self.psf_file)[0].data

        image /= image.sum()  # Flux normalized to 1

        self.psf = image

        self.n_pix_psf = self.psf.shape[0]

        spec_corr = len(self.spec_bins) + 1
        self.n_x_sim = self.n_x + spec_corr + 2*(self.n_pix_psf-1)
        self.n_y_sim = self.n_y + 2*(self.n_pix_psf-1)

        if return_psf:
            return image*self.zero_flux

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
        image = np.zeros((self.n_y_sim, self.n_x_sim))
        self.mask = np.zeros((self.n_y_sim, self.n_x_sim))
        wcs = self.create_wcs(self.n_x_sim, self.n_y_sim,
                              self.ra, self.dec, self.pixel_scale,
                              self.theta)

        return image, wcs

    def bin_xy(self, x, y):
        step = self.tel_params['dellambda']/self.interpolate_mult
        x_new = np.arange(x.min(), x.max(), step)
        y_new = np.interp(x_new, x, y)
        x, y = x_new, y_new

        x_new = self.bins

        bin_indices = np.digitize(x, x_new)
        binned_y = []
        for i in range(1, len(x_new)):
            if i in bin_indices:
                binned_y.append(y[bin_indices == i].mean())
            else:
                raise Exception("""Interpolation Error.
                                Increase interpolate_mult""")

        return np.array(binned_y)

    def bin_xy_norm(self, x, y):
        step = self.tel_params['dellambda']/self.interpolate_mult
        x_new = np.arange(x.min(), x.max(), step)
        y_new = np.interp(x_new, x, y)
        x, y = x_new, y_new

        x_new = self.bins

        bin_indices = np.digitize(x, x_new)
        binned_y = []
        for i in range(1, len(x_new)):
            if i in bin_indices:
                binned_y.append(np.median(y[bin_indices == i]))
            else:
                raise Exception("""Interpolation Error. \
                Increase interpolate_mult""")

        return np.array(binned_y)

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
        disp_width = len(self.spec_bins)
        disp_width_mid = disp_width//2

        patch = self.psf
        device = self.device

        # Convert to tensor once, before the loop
        image = torch.tensor(image, dtype=torch.float32, device=device)
        dl = self.tel_params['dellambda']
        lw = self.spec_bins

        x0 = df['x'].value.astype(int)
        y0 = df['y'].value.astype(int)

        x1 = x0 - patch_width_mid - disp_width_mid - 1
        x2 = x1 + patch_width + disp_width - 1
        y1 = y0 - patch_width_mid
        y2 = y1 + patch_width

        fluxes = []

        for row in df:
            x = row['wav']
            y = row['flux']
            fluxes.append(self.bin_xy(x, y))

        for x1_, x2_, y1_, y2_, flux in zip(x1, x2, y1, y2, fluxes):

            flux *= lw**2*3.34e4
            flux *= 1.51e3*self.flux_coeffs*(dl/lw)
            flux = torch.tensor(flux, dtype=torch.float32, device=device)

            # Generating star using sim

            data = torch.tensor(flux, dtype=torch.float32, device=device)
            data = data.view(1, 1, 1, -1)
            kernel = torch.tensor(patch, dtype=torch.float32, device=device)
            kernel = kernel.view(1, 1, patch_width, patch_width)
            spectrum = F.conv2d(data, kernel, padding=patch_width-1)

            image[y1_:y2_, x1_:x2_] += spectrum.squeeze()
            y0 = y1_ + patch_width_mid
            self.mask[y0 - 3:y0 + 4, x1_:x2_] += 0.5

            # Empty GPU cache after each iteration
            torch.cuda.empty_cache()

        # Convert to numpy array once, outside of the loop
        image = image.detach().cpu().numpy()
        image = image[patch_width-1:-patch_width+1,  # Cropping Image
                      patch_width+disp_width_mid-1:
                      -disp_width_mid-patch_width-1]
        return image

    def __call__(self, det_params=None, n_stack=1, stack_type='median',
                 spectroscopy=True, fwhm=None, detect_sources=False,
                 spec_width=50, sky_sep=5):

        zero_p_flux = self.zero_flux
        zero_p_flux *= self.gain
        ZP = 2.5*np.log10(zero_p_flux)
        self.ZP = ZP
        self.UC *= self.gain
        self.UC *= 3.34e4*self.spec_bins**2
        super().__call__(det_params=det_params, n_stack=n_stack,
                         photometry=None, stack_type=stack_type, ZP=1)

        if spectroscopy:
            self.extract(data=self.digital.copy(), wcs=self.wcs,
                         df=self.img_df, width=spec_width, sep=sky_sep)
