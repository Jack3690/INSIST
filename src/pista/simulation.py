import numpy as np
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from scipy.constants import c
import os
from pathlib import Path
import json
import requests
from .utils import bandpass
data_path = Path(__file__).parent.joinpath()

class Imager():
  def __init__(self,df, cols , tel_params, n_pix, exp_time,**kwargs):   
    """Image simulation using Source catalog. Class which simulates the field 
    and creates image

    Parameters
    ----------
    name  : str, Catalog name or coordinates of the source

    df   : pd.DataFrame, Pandas dataframe with source catalog

    cols  :,dict dict object with column name conversions for ra,dec,mag_nuv. 
    Eg {'RA': 'ra','Dec' : 'dec', 'ABmag' : 'mag_nuv'}

      'ra' (degrees)

      'dec' (degrees)

      'mag_nuv' (ABmag)
    psf_file : fits, 
    exp_time : float, Exposure time in seconds 
    n_pix    : int,number of pixels along one axis in the detector
    """
    self.cols   = cols


    # Flags
    self.shot_noise = True
    self.QE         = True
    self.sky        = True
    self.PRNU       = True
    self.DC         = True
    self.DNFP       = True
    self.QN         = True
    self.cosmic_rays= True

    # Parameters
    self.tel_params = {'aperture'   : 100, # cm
                       'pixel_scale': 0.1,
                       'psf_file'   :f'{data_path}/data/off_axis_hcipy.npy',
                       'response_funcs' : [],
                       'coeffs'         : 1
                       }
    
    if tel_params is not None:
        self.tel_params.update(tel_params)
    self.tel_params.update(kwargs)

    self.det_params = {'shot_noise' : 'Gaussian',
                       'M_sky' : 27,
                      'qe'         :  0.5,
                      'qe_sigma'   :  0.01,
                      'bias'       :  35,          # electrons
                      'G1'         :  1,
                      'bit_res'    :  14,        
                      'RN'         :  5,          # elec/pix
                      'PRNU_frac'  :  0.25/100,
                      'T'          :  223,        # K
                      'DFM'        :  1.424e-2,   # 14.24pA
                      'pixel_area' :  1e-6,       # 
                      'DN'         :  0.1/100,
                      'NF'         :  0,     # electrons 
                      'FWC'        : 1.4e5,
                      'C_ray_r'    :  2/50
                   }
    self.det_params.update(kwargs)
    
    self.n_pix_main  = n_pix
    self.pixel_scale = self.tel_params['pixel_scale']
    self.fov         = (n_pix*self.pixel_scale)/3600 # Degrees 
    self.radius      = self.fov/2

    self.response_funcs = self.tel_params['response_funcs']
    self.coeffs         = self.tel_params['coeffs']

    self.gain           = self.det_params['G1']*pow(2,
                            self.det_params['bit_res'])/self.det_params['FWC']
    self.tel_area = np.pi*(self.tel_params['aperture']/2)**2
                            
    self.exp_time    = exp_time  # seconds
    self.df          = df
    self.sim_run = False 
    
    self.psf_file = self.tel_params['psf_file']
    
    self.n_cosmic_ray_hits = int(self.exp_time*self.det_params['C_ray_r'])

    if self.df is not None:
        self.init_df()
        self.init_psf_patch() 
    else:
        print("df cannot be None")

  def init_psf_patch(self, return_psf = False):
      
    if len(self.response_funcs)>1:

        wav = np.linspace(1000,25000,10000)
        flux = (c*1e2*3.631e-20)/(wav**2*1e-8) 
        
        fig, ax, data, params= bandpass(wav,flux,self.response_funcs
                                         , plot = False)
        
        lambda_phot, int_flux, W_eff = params
    
        self.lambda_phot = lambda_phot
        self.int_flux    = int_flux
        self.W_eff       = W_eff
        self.int_flux_Jy =  (int_flux*lambda_phot**2*1e-8)/(c*1e2*1e-23)
        
        self.photons     =  1.51e3*self.int_flux_Jy*(W_eff/lambda_phot)
        
        self.zero_mag =  self.exp_time*self.tel_area*self.photons*self.coeffs
     
        filt_dat  = np.loadtxt(f'{data_path}/data/Sky_mag.dat')
        wav  = filt_dat[:,0]
        flux = filt_dat[:,1]
    
        fig, ax, data, params= bandpass(wav,flux,self.response_funcs
                                         , plot = False)
        
        lambda_eff, int_flux, W_eff = params
        self.det_params['M_sky'] = int_flux
        self.M_sky_p         = self.det_params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)
    else :
      self.zero_mag =  self.exp_time*(1.51e3*1000/2250)*self.tel_area*3631*self.coeffs
      self.M_sky_p         = self.det_params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)

    ext = self.psf_file.split('.')[-1]
    
    if ext=='npy':
        image =  np.load(self.psf_file)
    elif ext=='fits':
        image = fits.open(self.psf_file)[0].data
         
    image /= image.sum()  # Flux normalized to 1
    self.image_g_sub = image
    F_sky_p           = self.zero_mag*pow(10,-0.4*self.M_sky_p)
    self.sky_bag_flux = F_sky_p    
    self.zero_flux    = self.zero_mag

    if return_psf:
      return image*self.zero_flux

  def init_df(self):     

      print("Input : Dataframe.")

      if self.cols is not None:
          self.df = self.df.rename(columns = self.cols) 
      self.ra   = (self.df['ra'].max()+self.df['ra'].min())/2
      self.dec  = (self.df['dec'].max()+self.df['dec'].min())/2
      
      ra_min = (self.df['ra']>self.ra - self.radius) 
      ra_max = (self.df['ra']<self.ra + self.radius)
      self.df = self.df[ ra_min & ra_max ]
      dec_min = (self.df['dec']>self.dec - self.radius)
      dec_max = (self.df['dec']<self.dec + self.radius)
      self.df = self.df[dec_min & dec_max]
    
      self.name = f" RA : {np.round(self.ra,3)} degrees, Dec : {np.round(self.dec,3)} degrees"

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
          
    self.n_pix_sub  = self.image_g_sub.shape[0]
    if self.n_pix_sub%2!=0:
        self.n_pix_main += self.n_pix_sub -1
    else:
        self.n_pix_main += self.n_pix_sub

    self.image    = np.zeros((self.n_pix_main, self.n_pix_main))
    self.wcs      = self.create_wcs(self.n_pix_main,self.ra, self.dec, self.pixel_scale)
    if return_img:
      return self.image, self.wcs


  def create_wcs(self,npix,ra,dec,pixel_scale):
    """
    Parameters
    ----------
    npix : int
        number of pixels in image
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
    w.wcs.crpix = [(npix-1)//2, (npix-1)//2]
    w.wcs.cdelt = np.array([-pixel_scale/3600, self.pixel_scale/3600])
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w

  def compute_coeff_arrays(self):
    """
    
    Computed coefficients based on input parameters
    Returns
    -------
    None.

    """

    n_pix = self.n_pix_main #- self.n_pix_sub-1

    if self.QE:
      if self.det_params['qe']>0 and self.det_params['qe']<1:

        self.qe_array =  np.random.normal(loc=self.det_params['qe'], 
                                          scale=self.det_params['qe_sigma'],
                                          size=(n_pix, n_pix))
      else:
        print('QE should in the range (0,1]')
    else:
      self.qe_array = 1

    self.bias_array =  np.random.normal(loc=self.det_params['bias'], 
                                      scale=self.det_params['RN'],
                                      size=(n_pix, n_pix))
    if self.PRNU:

      self.PRNU_array =  np.random.normal(loc=0, 
                                        scale = self.det_params['PRNU_frac'],
                                        size=(n_pix, n_pix))
    else:
       self.PRNU_array = 0
       
    if self.DC:
      self.DR = self.dark_current(self.det_params['T'], self.det_params['DFM'], 
                                  self.det_params['pixel_area'])

      self.DC_array = np.random.normal(loc = self.DR*self.exp_time, 
                                          scale = np.sqrt(self.DR*self.exp_time),
                                          size=(n_pix, n_pix))
    else: 
        self.DC_array = 0
    if self.DNFP and self.DC:
      self.DNFP_array =  np.random.lognormal(mean= 0, 
                                    sigma = self.exp_time*self.DR*self.det_params['DN'],
                                    size=(n_pix, n_pix))
      self.DC_array*=(1 + self.DNFP_array)
    else:
        self.DNFP_array = 0

    if self.QN:
      self.QN_value = (self.det_params['FWC']/(pow(2,self.det_params['bit_res'])
                                            *np.sqrt(12)))

      self.QN_array = self.QN_value*np.random.randint(-1,2,size = (n_pix,n_pix))
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
    Kb  = 8.62e-5
    const	= 2.55741439581387e15

    EgT	= 1.1557 - (7.021e-4*T**2/(1108+T))
    DR	= const*pixel_area*(T**1.5)*DFM*np.exp(-EgT/(2*Kb*T))
    return DR
    
  def generate_photons(self, image,npix_m, npix_s,df):
    """
      This function creates PSFs based on ABmag  on a 
      small patch (2D array) of size n_pix_s*n_pix_s. 
      
      The patch with the PSF is then added to the image array of size 
      n_pix_m*n_pix_m using wcs object.

      Parameters
      ----------
      image : numpy.ndarray
         base image array for inserting star PSFs
      npix_m : int
          number of pixels (length) in base image
      npix_s : int
          number of pixels (length) in psf patch image
      df : pandas.dataframe
          Dataframe containing source list
          

      Returns
      -------
      image : numpy.ndarray
          Array with PSFs added based on df

      """
    patch_width = npix_s//2
    for i, row in df.iterrows():

        c = SkyCoord(row['ra'],row['dec'],unit=u.deg)
        pix = self.wcs.world_to_array_index(c)
        ABmag = row['mag']

        flux  = self.zero_flux*10**(-ABmag/2.5)  # Photo-elec per second

        patch =  flux*self.image_g_sub
    
        x1 = abs(pix[0]) - patch_width
        x2 = x1 + npix_s
        y1 = pix[1] - patch_width
        y2 = y1 + npix_s

        image[ x1: x2, y1:y2 ] += patch

    if npix_s%2!=0:
        image = image[patch_width:-patch_width,patch_width:-patch_width]
    else:
        image = image[patch_width:-patch_width,patch_width:-patch_width]

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
      n_pix = self.n_pix_main
    else :
      n_pix = array.shape[0]

    if type_ == 'Gaussian':
        shot_noise = np.random.normal(loc=array, scale=np.sqrt(array), size = (n_pix, n_pix))
    elif type_ =='Poisson':
        shot_noise = np.random.poisson(lam=array, size = (n_pix, n_pix)).astype(array.dtype)
    else:
      print('Invalid type')

    return shot_noise  


  def __call__(self,det_params = None,n_stack =1,stack_type = 'median'):
    """
    
     Parameters
     ----------
    det_params: dict, optional
     Dictionary contianing detector parameters. The default is None.
     n_stack : int, optional
     Number of observations to be stacked. The default is 1.
     stack_type : str, optional
     Stacking method. The default is 'median'.
     
     
    Simulates field by taking dataframe as input by inserting PSF patches in 
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
      self.M_sky_p     = self.det_params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)
      self.init_psf_patch() 

    digital_stack = []

    self.compute_coeff_arrays()
    
    for i in range(n_stack):

      self.init_image_array()
      
      self.source_photons   = self.generate_photons(self.image,self.n_pix_main,
                                                    self.n_pix_sub, self.df)
      
      if self.shot_noise:
         self.source_photons  = self.compute_shot_noise(self.source_photons 
                              ,type_ = self.det_params['shot_noise'])
         
      self.source_photoelec = self.source_photons*self.qe_array
         
      self.n_pix_main = self.source_photoelec.shape[0]

      if self.sky:
        self.sky_photoelec = self.compute_shot_noise(self.sky_bag_flux,
                                                     'Gaussian')*self.qe_array
        
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
          x = np.random.randint(0,self.n_pix_main)
          y = np.random.randint(0,self.n_pix_main)
          self.digital[x,y] = pow(2, self.det_params['bit_res'])

    self.wcs = self.create_wcs(self.n_pix_main,self.ra,self.dec, self.pixel_scale)

    self.header = self.wcs.to_header()
    
    return self.digital