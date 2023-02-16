from astropy.wcs import WCS
from astropy.io import fits
from scipy.constants import c
import os
from pathlib import Path
import json
import requests
from .utils import bandpass
from .analysis import *


data_path = Path(__file__).parent.joinpath()

class Imager(Analyzer):
  def __init__(self, df, coords = None, tel_params = None, n_x = 1000, 
               n_y = 1000, exp_time = 100, plot = False,**kwargs):  
    super().__init__() 
    """Image simulation using Source catalog. Class which simulates the field 
    
    and creates image

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
                                          number of times to multiply filter
                                          profile
                  'coeffs'         : float, filter coefficients if not
                                          response_funcs
                  } 
    n_x      : int, 
               number of pixels along RA direction
    n_y      : int, 
               number of pixels along Dec direction
    exp_time : float, 
               Exposure time in seconds 

    """

    # Flags
    self.shot_noise  = True
    self.QE          = True
    self.sky         = True
    self.PRNU        = True
    self.DC          = True
    self.DNFP        = True
    self.QN          = True
    self.cosmic_rays = False

    # Telescope and Detector Parameters
    self.tel_params = {'aperture'       : 100, # cm
                       'pixel_scale'    : 0.1,
                       'psf_file'       :f'{data_path}/data/PSF/INSIST/on_axis_hcipy.npy',
                       'response_funcs' : [],
                       'coeffs'         : 1,
                       'theta'          : 0
                       }
    
    if tel_params is not None:
      self.tel_params.update(tel_params)
    self.tel_params.update(kwargs)

    self.det_params = {'shot_noise': 'Gaussian',
                      'M_sky'      :  27,
                      'qe_sigma'   :  0.01,       # Pixel to pixel fluctuation
                      'bias'       :  35,         # electrons
                      'G1'         :  1,
                      'bit_res'    :  14,         
                      'RN'         :  5,          # elec/pix
                      'PRNU_frac'  :  0.25/100,   # PRNU sigma
                      'T'          :  223,        # K
                      'DFM'        :  1.424e-2,   # 14.24 pA
                      'pixel_area' :  1e-6,       # 
                      'DN'         :  0.1/100,    # 
                      'NF'         :  0,          # electrons 
                      'FWC'        :  1.4e5,      # electrons 
                      'C_ray_r'    :  2/50        # hits/second
                   }
    self.det_params.update(kwargs)

    self.df    = df.dropna().copy()
    self.n_x  = n_x
    self.n_y  = n_y
    self.pixel_scale    = self.tel_params['pixel_scale']
    self.theta          = self.tel_params['theta']*np.pi/180

    
    self.response_funcs = self.tel_params['response_funcs']
    self.coeffs         = self.tel_params['coeffs']

    self.gain           = self.det_params['G1']*pow(2,
                            self.det_params['bit_res'])/self.det_params['FWC'] # DN/electrons
    self.tel_area       = np.pi*(self.tel_params['aperture']/2)**2
                            
    self.exp_time       = exp_time  # seconds
    self.sim_run        = False 
    
    self.psf_file       = self.tel_params['psf_file']

    # Input Dataframe
    if 'mag' not in self.df.keys():
      raise Exception("'mag' column not found input dataframe")
    
    if 'ra' not in self.df or 'dec' not in self.df.keys():
      if 'x' in self.df.keys() and 'y' in self.df.keys():
        self.ra  = 10
        self.dec = 10
      else:
        raise Exception("'ra','dec','x',or 'y', columns not found in input dataframe ")

    elif coords is None:
      self.ra  = np.median(self.df.ra)
      self.dec = np.median(self.df.dec)
    else:
      self.ra  = coords[0]
      self.dec = coords[1]

    ra_n  = np.round(self.ra,3)
    dec_n = np.round(self.dec,3)
    self.name = f" RA : {ra_n} degrees, Dec : {dec_n} degrees"

    if self.df is not None:
        self.calc_zp(plot = plot)
        self.init_psf_patch() 

        # Cropping df to sim_field
        x_left  = self.n_pix_sub - 1
        x_right = self.n_x_sim  - self.n_pix_sub - 1
        y_left  = self.n_pix_sub - 1
        y_right = self.n_y_sim  - self.n_pix_sub - 1

        self.sim_df = self.init_df(df = self.df, 
                                   n_x = self.n_x_sim, n_y = self.n_y_sim,
                                   x_left = x_left   , x_right= x_right,
                                   y_left = y_left   , y_right = y_right)
        if len(self.sim_df)<1:
          raise Exception("Not Enough sources inside FoV. Increase n_x and n_y")
    else:
        print("df cannot be None")

    # Init Cosmic Rays
    #area = ((n_x*self.pixel_scale)/3600)*((n_y*self.pixel_scale)/3600)
    
    #eff_area = (0.17/area)*self.exp_time
    
    #self.n_cosmic_ray_hits = int(eff_area*self.det_params['C_ray_r'])

  def calc_zp(self, plot = False):
    if len(self.response_funcs)>0:

        wav  = np.linspace(1000,10000,10000)
        flux = 3631/(3.34e4*wav**2)   # AB flux
        
        fig, ax, _ , params = bandpass(wav,flux,self.response_funcs
                                         , plot = plot)
        
        lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params
    
        self.lambda_phot = lambda_phot
        self.int_flux    = int_flux
        self.W_eff       = W_eff
        self.int_flux_Jy = int_flux_Jy
        self.flux_ratio  = flux_ratio
        
        self.photons     = 1.51e3*self.int_flux_Jy*(W_eff/lambda_phot)*flux_ratio
        self.zero_flux   = self.exp_time*self.tel_area*self.photons*self.coeffs
     
        filt_dat  = np.loadtxt(f'{data_path}/data/Sky_mag.dat')
        wav  = filt_dat[:,0]
        flux = filt_dat[:,1]
    
        _, _, _ , params = bandpass(wav,flux,self.response_funcs
                                         , plot = False)
        
        int_flux = params[1]
        self.det_params['M_sky'] = int_flux
        self.M_sky_p   = self.det_params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)

    else :

      print("Response functions not provided. Using default values")
      self.int_flux_Jy = 3631
      self.W_eff       = 1000
      self.lambda_phot = 2250

      self.photons   = 1.51e3*self.int_flux_Jy*(self.W_eff/self.lambda_phot)
      self.zero_flux = self.exp_time*self.tel_area*self.photons*self.coeffs
      self.M_sky_p   = self.det_params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)

  def init_psf_patch(self, return_psf = False): 

    ext = self.psf_file.split('.')[-1]
    
    if ext=='npy':
        image =  np.load(self.psf_file)
    elif ext=='fits':
        image = fits.open(self.psf_file)[0].data
         
    image            /= image.sum()  # Flux normalized to 1
    self.image_g_sub  = image
    self.sky_bag_flux = self.zero_flux*pow(10,-0.4*self.M_sky_p)  

    self.n_pix_sub    = self.image_g_sub.shape[0]

    # Defining shape of simulation field
    self.n_x_sim = self.n_x + 2*(self.n_pix_sub-1)
    self.n_y_sim = self.n_y + 2*(self.n_pix_sub-1)

    if return_psf:
      return image*self.zero_flux

  def init_df(self, df , n_x, n_y, x_left, x_right, y_left, y_right):  

      wcs = self.create_wcs(n_x, n_y, self.ra,self.dec, 
                           self.pixel_scale, self.theta)


      if 'x' not in self.df.keys():
        c     = SkyCoord(df['ra'],df['dec'],unit=u.deg)
        pix   = wcs.world_to_array_index(c)
        
        df['x'] = pix[1]
        df['y'] = pix[0]

      else:
        df['x'] = df['x'] + x_left
        df['y'] = df['y'] + y_left
      
      # Cropping Dataframe based on FoV
      x_min_cut  = (df['x']>x_left) 
      x_max_cut  = (df['x']<x_right)

      df = df[x_min_cut & x_max_cut]

      y_min_cut  = (df['y']>y_left) 
      y_max_cut  = (df['y']<y_right)

      df = df[y_min_cut & y_max_cut]

      return df

  def init_image_array(self, return_img = False):
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
    self.image    = np.zeros((self.n_y_sim, self.n_x_sim)) #  
    self.wcs      = self.create_wcs(self.n_x_sim, self.n_y_sim,
                                    self.ra, self.dec, self.pixel_scale,
                                    self.theta)
    if return_img:
      return self.image, self.wcs


  def create_wcs(self,n_x, n_y, ra,dec,pixel_scale, theta = 0 ):
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
    w.wcs.pc    = np.array([[ np.cos(theta),-np.sin(theta)],
                            [ np.sin(theta), np.cos(theta)]])
    return w

  def compute_coeff_arrays(self):
    """
    
    Computed coefficients based on input parameters
    Returns
    -------
    None.
    """
    n_x = self.n_y 
    n_y = self.n_x

    if self.QE:
        self.qe_array =  np.random.normal(loc=0, 
                                          scale=self.det_params['qe_sigma'],
                                          size=(n_x, n_y))
    else:
      self.qe_array = 0

    self.bias_array =  np.random.normal(loc=self.det_params['bias'], 
                                      scale=self.det_params['RN'],
                                      size=(n_x, n_y))
    if self.PRNU:

      self.PRNU_array =  np.random.normal(loc=0, 
                                        scale = self.det_params['PRNU_frac'],
                                        size=(n_x, n_y))
    else:
       self.PRNU_array = 0
       
    if self.DC:
      self.DR = self.dark_current(self.det_params['T'], self.det_params['DFM'], 
                                  self.det_params['pixel_area'])

      self.DC_array = np.random.normal(loc = self.DR*self.exp_time, 
                                          scale = np.sqrt(self.DR*self.exp_time),
                                          size=(n_x, n_y))
    else: 
        self.DC_array = 0
    if self.DNFP and self.DC:
      self.DNFP_array =  np.random.lognormal(mean= 0, 
                                    sigma = self.exp_time*self.DR*self.det_params['DN'],
                                    size=(n_x, n_y))
      self.DC_array*=(1 + self.DNFP_array)
    else:
        self.DNFP_array = 0

    if self.QN:
      self.QN_value = (self.det_params['FWC']/(pow(2,self.det_params['bit_res']) #electrons
                                            *np.sqrt(12)))

      self.QN_array = self.QN_value*np.random.randint(-1,2,size = (n_x,n_y))
    else :
      self.QN_array = 0

  def dark_current(self,T, DFM, pixel_area):
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
    Kb    = 8.62e-5
    const	= 2.55741439581387e15

    EgT	= 1.1557 - (7.021e-4*T**2/(1108+T))
    DR	= const*pixel_area*(T**1.5)*DFM*np.exp(-EgT/(2*Kb*T))
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
    ABmag  = df['mag'].values
    flux   = self.zero_flux * 10**(-ABmag/2.5)
    patch  = self.image_g_sub
    
    x1     = x0 - patch_width_mid
    x2     = x1 + patch_width
    y1     = y0 - patch_width_mid
    y2     = y1 + patch_width
    
    for x1_,x2_,y1_,y2_,flux_ in zip(x1,x2,y1,y2,flux):
      image[y1_:y2_, x1_:x2_] += flux_*patch

    image = image[patch_width-1:-patch_width+1, patch_width-1:-patch_width+1]

    return image

  def compute_shot_noise(self,array,type_ = 'Gaussian'):
    """

    Parameters
    ----------
    array : numpy.ndarray
        input array
    type_ : str, optional
         The default is 'gaussian'.

    Returns
    -------
    shot_noise : numpy.ndarray
         Return array with shot noise
    """

    if type(array) == np.float64 :
      n_x = self.n_y
      n_y = self.n_x

    else :
      n_x = array.shape[0]
      n_y = array.shape[1]

    if type_ == 'Gaussian':
        shot_noise = np.random.normal(loc=array, scale=np.sqrt(array), 
                                      size = (n_x, n_y))
    elif type_ =='Poisson':
        shot_noise = np.random.poisson(lam=array,
                                       size = (n_x, n_y)).astype(array.dtype)
    else:
      print('Invalid type')

    return shot_noise  


  def __call__(self,det_params = None,n_stack = 1,stack_type = 'median',
               photometry = 'Aper', fwhm = None, detect_sources = False,
               ZP = None, **kwargs):
    """
     Parameters
     ----------
    det_params: dict, optional  
     Dictionary contianing detector parameters. The default is None.
                {     'shot_noise' :  str,
                      'M_sky'      :  float,
                      'qe_sigma'   :  float,       Pixel to pixel fluctuation
                      'bias'       :  float,       electrons
                      'G1'         :  float,
                      'bit_res'    :  int,         
                      'RN'         :  float,       elec/pix
                      'PRNU_frac'  :  float,       PRNU sigma
                      'T'          :  float,       K
                      'DFM'        :  float,       pA
                      'pixel_area' :  float,       
                      'DN'         :  float,     
                      'NF'         :  float,       electrons 
                      'FWC'        :  float,       electrons 
                      'C_ray_r'    :  float        hits/second
                  }
                   
     n_stack    : int, optional
                  Number of observations to be stacked. The default is 1.

     stack_type : str, optional
                  Stacking method. The default is 'median'.
     
     
    Simulates field by taking dataframe as input by inserting sim patches in 
    image array using WCS
     Returns
     -------
     numpy.ndarray
     Final image array after adding all layers of simulation
    
    """
    self.sim_flag = True
    if det_params is not None:
      
      self.det_params.update(det_params)
      self.gain        = self.det_params['G1']*pow(2,
                         self.det_params['bit_res'])/self.det_params['FWC']

    digital_stack = []
    self.compute_coeff_arrays()
    
    for i in range(n_stack):

      self.init_image_array()
      
      self.source_photons   = self.generate_photons(self.image, self.n_pix_sub,
                                                    self.sim_df)
      
      if self.shot_noise:
         self.source_photons  = self.compute_shot_noise(self.source_photons 
                              ,type_ = self.det_params['shot_noise'])
 
      self.source_photoelec = self.source_photons*(1+self.qe_array)

      if self.sky:
        self.sky_photoelec = self.compute_shot_noise(self.sky_bag_flux,
                                                     'Gaussian')*(1+ self.qe_array)
        
        self.light_array = self.source_photoelec +  self.sky_photoelec
      else:
        self.light_array = self.source_photoelec 

      if self.PRNU:
        self.light_array*=(1+self.PRNU_array)

      if self.DC:
        self.photoelec_array = self.light_array + self.DC_array
      else:
        self.photoelec_array  = self.light_array

      self.charge = self.photoelec_array + self.QN_array + self.det_params['NF'] + self.bias_array 
      
      self.digital = (self.charge*self.gain).astype(int)

      # Full well condition
      self.digital = np.where(self.digital>=pow(2, self.det_params['bit_res']),pow(2, self.det_params['bit_res']), self.digital)
      
      digital_stack.append(self.digital)
    
    digital_stack = np.array(digital_stack)
    if n_stack>1:
      if stack_type   == 'median':
        self.digital  = np.median(digital_stack, axis = 0)
      elif stack_type == 'mean':
        self.digital  = np.median(digital_stack, axis = 0)
        
    if self.cosmic_rays:
        for i in range(self.n_cosmic_ray_hits):
          x = np.random.randint(0,self.n_x_main)
          y = np.random.randint(0,self.n_y_main)
          self.digital[x,y] = pow(2, self.det_params['bit_res'])


    self.wcs = self.create_wcs(self.n_x,self.n_y, 
                               self.ra,self.dec, self.pixel_scale, self.theta)
    
    # Filtering out sources within Image
    x_left  = 0
    x_right = self.n_x
    y_left  = 0
    y_right = self.n_y

    self.img_df = self.init_df(df = self.sim_df.copy(), 
                                n_x = self.n_x    , n_y = self.n_y,
                                x_left = x_left   , x_right= x_right,
                                y_left = y_left   , y_right = y_right)

    self.header = self.wcs.to_header()
    self.header['gain'] = self.det_params['G1']
    self.header['Temp'] = str(self.det_params['T']) + 'K' 
    self.header['bias'] = self.det_params['bias']
    self.header['RN']   = self.det_params['RN']
    self.header['DR']   = self.DR
    self.header['NF']   = self.det_params['NF']

    super().__call__(df = self.img_df, wcs = self.wcs, 
                     data = self.digital.astype(float),
                     photometry = photometry, fwhm = fwhm,
                     detect_sources = detect_sources, ZP = ZP)




