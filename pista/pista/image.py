import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
from pathlib import Path
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from matplotlib import colors as col
from astropy.wcs import WCS

from photutils import aperture as aper
from photutils.aperture import aperture_photometry


sb.set_style('darkgrid')
matplotlib.rcParams['font.size']=12
matplotlib.rcParams['figure.figsize']=(10,10)

data_path = Path(__file__).parent.joinpath()

class PSF_gen():

  def __init__(self, df, cols, axis, mode,exp_time,fov):   
    """PSF Generator using Source catalog. Class which simulates the field 
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
    axis     : str, 'Off' or "On", Secondary of the telescope.
    mode     : str, 'HCIPy' or 'Zeemax', PSF patch generation method
    exp_time : float, Exposure time in seconds 
    fov      : float, Field of View in degrees
    """
    self.axis = axis
    self.mode = mode
    self.cols = cols
    self.fov  = fov
    self.radius = fov/2
    self.shot_noise = True
    self.sky  = True
    self.PRNU = True
    self.DC   = True
    self.DNFP = False
    self.pixel_scale = 0.1
    
    self.params = {'shot_noise' : 'Gaussian',
                   'M_sky'      :  27.5,
                   'qe'         :  0.5,
                   'qe_sigma'   :  0.01,
                   'bias'       :  733,
                   'gain'       :  1.85,
                   'bit_res'    :  16,
                   'RN'         :  3,
                   'PRNU_frac'  :  0.02,
                   'T'          :  218,        # K
                   'DFM'        :  1.424e-2,   # 14.24pA
                   'pixel_area' :  1e-6,       #
                   'DN'         :  0.4}
    self.full_well_capacity = pow(2, self.params['bit_res'])

    self.M_sky_p     = self.params['M_sky'] - 2.5*np.log10(self.pixel_scale**2)
    self.exp_time    = exp_time  # seconds
    self.df          = df
    if self.df is not None:
        self.init_df()
        self.init_psf_patch() 
    else:
        raise print("df cannot be None")

  def init_psf_patch(self, return_psf = False):
    """
    Parameters
    ----------
    return_psf : bool, 
        DESCRIPTION.

    Returns
    -------
    numpy.ndarray
        psf patch array scaled with zero_point based on axis and mode

    """
    self.zero_mag_s_on =  self.exp_time*1.51e3*3631*np.pi*(100/2)**2*(1500/2250)*0.8**6*0.95**2*0.68*0.83 # Photons

    self.zero_mag_s_off = self.exp_time*1.51e3*3631*np.pi*(100/2)**2*(1500/2250)*0.8**5*0.95**2*0.83      # Photons

    if self.mode == 'Zeemax':
      if self.axis =='On':
        image =  np.load(f'{data_path}/data/On_PSF_Zmax.npy')
        image /= image.sum()
        self.image_g_sub = image
        F_sky_p           = self.zero_mag_s_on*pow(10,-0.4*self.M_sky_p)
        self.sky_bag_flux = F_sky_p    
        self.zero_flux    = self.zero_mag_s_on 

      elif self.axis=='Off':
        image  = np.load(f'{data_path}/data/Off_PSF_Zmax.npy')
        image /= image.sum()
        self.image_g_sub  = image
        F_sky_p           = self.zero_mag_s_off*pow(10,-0.4*self.M_sky_p)
        self.sky_bag_flux = F_sky_p    
        self.zero_flux    = self.zero_mag_s_off 

    elif self.mode =='HCIPy':
      if self.axis =='On':
        image  = np.load(f'{data_path}/data/on_axis_hcipy.npy')
        image /= image.sum()
        self.image_g_sub =  image
        F_sky_p           = self.zero_mag_s_on*pow(10,-0.4*self.M_sky_p)
        self.sky_bag_flux = F_sky_p    
        self.zero_flux    = self.zero_mag_s_on  

      elif self.axis=='Off':
        image  = np.load(f'{data_path}/data/off_axis_hcipy.npy')
        image /= image.sum()
        self.image_g_sub  = image
        F_sky_p           = self.zero_mag_s_off*pow(10,-0.4*self.M_sky_p)
        self.sky_bag_flux = F_sky_p    
        self.zero_flux    = self.zero_mag_s_off 
    else:
        print("Invalid mode!")
    if return_psf:
      return image*self.zero_flux

  def init_df(self):     

    if self.cols is not None:
      self.df = self.df.rename(columns = self.cols) 
      
    self.ra   = (self.df['ra'].max()+self.df['ra'].min())/2
    self.dec  = (self.df['dec'].max()+self.df['dec'].min())/2
    
    self.df = self.df[ (self.df['ra']>self.ra - self.radius) & (self.df['ra']<self.ra + self.radius)]
    self.df = self.df[ (self.df['dec']>self.dec - self.radius) & (self.df['dec']<self.dec + self.radius)]
    
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
    del_ra  = self.df.ra.max()  - self.df.ra.min()
    del_dec = self.df.dec.max() - self.df.dec.min()
    
    n_pix_main = del_ra*3600/self.pixel_scale if del_ra>=del_dec else del_dec*3600/self.pixel_scale
  
    self.n_pix_main = int(n_pix_main) + 2*self.n_pix_sub

    if self.n_pix_main <=10000:
      self.image    = np.zeros((self.n_pix_main, self.n_pix_main))
      self.wcs      = self.create_wcs(self.n_pix_main,self.ra, self.dec, self.pixel_scale)
      
    else:
      print("FoV is too big.")

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
    if self.params['qe']>0 and self.params['qe']<1:

      n_pix = self.n_pix_main - self.n_pix_sub

      if n_pix > 205 :
        n_pix -= self.n_pix_sub

      self.qe_array =  np.random.normal(loc=self.params['qe'], 
                                        scale=self.params['qe_sigma'],
                                        size=(n_pix, n_pix))
    else:
      print('QE should in the range (0,1]')


    self.bias_array =  np.random.normal(loc=self.params['bias'], 
                                      scale=self.params['RN'],
                                      size=(n_pix, n_pix))
    if self.PRNU:

      self.PRNU_array =  np.random.normal(loc=0, 
                                        scale = self.params['PRNU_frac'],
                                        size=(n_pix, n_pix))
    if self.DC:
      self.DR = self.dark_current(self.params['T'], self.params['DFM'], 
                                  self.params['pixel_area'])

      self.DC_array = np.random.normal(loc = self.DR*self.exp_time, 
                                          scale = np.sqrt(self.DR*self.exp_time),
                                          size=(n_pix, n_pix))


    if self.DNFP and self.DC:
      self.DNFP_array =  np.random.lognormal(mean= 0, 
                                    sigma = self.exp_time*self.DR*self.params['DN'],
                                    size=(n_pix, n_pix))
      self.DC_array*=(1 + self.DNFP_array)

  
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
    if npix_s%2 ==0:
      patch_width_l = npix_s//2
      patch_width_r = npix_s//2

    else:
      patch_width_l = npix_s//2 
      patch_width_r = npix_s//2 +1

    for i, row in df.iterrows():

        c = SkyCoord(row['ra'],row['dec'],unit=u.deg)
        pix = self.wcs.world_to_array_index(c)
        ABmag = row['mag_nuv']

        flux  = self.zero_flux*10**(-ABmag/2.5)  # Photo-elec per second

        patch =  flux*self.image_g_sub

        x1 = pix[0] - patch_width_l
        x2 = pix[0] + patch_width_r
        y1 = pix[1] - patch_width_l
        y2 = pix[1] + patch_width_r

        image[ x1: x2, y1:y2 ] += patch

    image   = image[patch_width_l-1:-patch_width_r-1,patch_width_l-1:-patch_width_r-1]
    if image.shape[0]>205:
      image = image[patch_width_l-1:-patch_width_r-1,patch_width_l-1:-patch_width_r-1]

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


  def __call__(self,params = None,n_stack =1,stack_type = 'median'):
    """
    
     Parameters
     ----------
     params : dict, optional
     Dictionary contianing simulation parametes. The default is None.
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
    if params is not None:
      self.params.update(params)
    readout_stack = []

    self.init_image_array()
    self.compute_coeff_arrays()
    
    for i in range(n_stack):

      self.init_image_array()

      self.photon_array    = self.generate_photons(self.image,self.n_pix_main,self.n_pix_sub, self.df)
      self.photoelec_array = self.photon_array*self.qe_array

      if self.shot_noise:
          self.photoelec_array     = self.compute_shot_noise(self.photoelec_array,type_ = self.params['shot_noise'])

      self.n_pix_main = self.photoelec_array.shape[0]

      if self.sky:
        self.sky_array = self.compute_shot_noise(self.sky_bag_flux,'Gaussian')
        self.sky_array *= self.qe_array
        
        self.light_array = self.photoelec_array +  self.sky_array
      else:
        self.light_array = self.photoelec_array

      if self.PRNU:
        self.light_array*=(1+self.PRNU_array)

      if self.DC:
        self.charge = self.light_array + self.DC_array
      else:
        self.charge  = self.light_array

      readout = (self.charge + self.bias_array)/self.params['gain']   # ADU

      # Full well 

      self.readout = np.where(readout>self.full_well_capacity, self.full_well_capacity, readout)
      
      readout_stack.append(self.readout)
    
    readout_stack = np.array(readout_stack)
    if n_stack>1:
      if stack_type == 'median':
        self.readout = np.median(readout_stack, axis = 0)
      elif stack_type == 'mean':
        self.readout = np.median(readout_stack, axis = 0)

    self.wcs = self.create_wcs(self.n_pix_main,self.ra,self.dec, self.pixel_scale)

    self.header = self.wcs.to_header()
    
    return self.readout

class PSF(PSF_gen):
  def __init__(self, df = None, cols = None, axis = 'On',mode = 'Zeemax',
               exp_time = 100, fov = 0.02):
      super().__init__(df=df, cols=cols, axis=axis, mode=mode,
                       exp_time = exp_time, fov = fov)
      """
      
      A class to visualize and analyze the simulated image
      
      Parameters
      ----------
      name  : str, Catalog name or coordinates of the source

      df   : pd.DataFrame, Pandas dataframe with source catalog

      cols  :,dict dict object with column name conversions for ra,dec,mag_nuv. 
      Eg {'RA': 'ra','Dec' : 'dec', 'ABmag' : 'mag_nuv'}

        'ra' (degrees)

        'dec' (degrees)

        'mag_nuv' (ABmag)
      axis     : str, 'Off' or "On", Secondary of the telescope.
      mode     : str, 'HCIPy' or 'Zeemax', PSF patch generation method
      exp_time : float, Exposure time in seconds 
      fov      : float, Field of View in degrees
      
     

      Returns
      -------
      None.

      """

  def __call__(self,params= None, n_stack =1, stack_type ='median'):
    super().__call__(params = params, n_stack = n_stack, stack_type = stack_type)
    """
    
    Performs PSF simulation and PSF Photometry
    
    TBW
    """
    
    
    if self.axis == 'off':
      zero_p_flux = (3010478142.88666  + self.params['bias'])*self.params['gain']
    else:
      zero_p_flux = (1516505736.205873 + self.params['bias'])*self.params['gain']

    c          = SkyCoord(self.df['ra'], self.df['dec'],unit=u.deg)
    data       = self.readout
    wcs        = self.wcs
    pix        = wcs.world_to_array_index(c)

    position        = [(i,j) for i,j in zip(pix[1],pix[0])]

    aperture        = aper.CircularAperture(position, r=0.3/0.1)
    ap_pix          = np.count_nonzero(aperture.to_mask()[0])
    aperture_bag    = aper.CircularAnnulus(position, r_in = 0.3/0.1, r_out = 1/0.1)
    bag_mask        = aperture_bag.to_mask()[0]

    bag_flux        = bag_mask.get_values(data)

    # Median bag flux
    bag_flux_med    = np.sort(bag_flux)[len(bag_flux)//2]

    phot_table      = aperture_photometry(data, [aperture, aperture_bag])

    phot_table['sky_flux'] = ap_pix*bag_flux_med
    phot_table['flux']     = phot_table['aperture_sum_0'].value - phot_table['sky_flux'].value

    phot_table['flux_err'] = np.sqrt( phot_table['flux'].value  + phot_table['sky_flux'].value )
 
    phot_table['SNR']      = phot_table['flux'].value/ phot_table['flux_err'].value

    phot_table['mag_in']   = self.df.mag_nuv
    phot_table['mag_0.3']  = -2.5*np.log10(phot_table['flux'].value/(zero_p_flux*self.exp_time))
    phot_table['mag_err']  = 1.087*phot_table['flux_err'].value/phot_table['flux'].value
    self.phot_table = phot_table


  def show_field(self,figsize=(10,10)):
    """
    Function for creating a scatter plot of sources within the FoV

    Returns
    -------
    fig, ax
    """
    if self.wcs is None :
      self.init_image_array()

    fig, ax = plt.subplots(1,1,figsize=figsize)
    ax.scatter(self.df['ra'],self.df['dec'],marker='.',color='black')
    ax.set_title(f" Requested Center : {self.name} \n FoV : {np.round(self.pixel_scale*(self.n_pix_main-2*self.n_pix_sub )/3600,3)} degrees | {len(self.df)} sources")
    ax.invert_xaxis()
    ax.set_xlabel('RA (Degrees)')
    ax.set_ylabel('Dec (Degrees)')
    return fig,ax

  def show_image(self, source = 'Readout', figsize = (15,10), download = False):
    """
    Function for plotting the simulated field with PSFs

    Returns
    -------
    fig, ax
    """

    fig = plt.figure(figsize = figsize)
    norm = None

    if source =='Readout':
      data  = self.readout
      norm = col.LogNorm()
    elif source == 'Sky':
      data = self.sky_array
    elif source == 'DC':
      data = self.DC_array
    elif source == 'QE':
      data = self.qe_array
    elif source =='Bias':
      data = self.bias_array + self.DC_array
    elif source == 'PRNU':
      data = self.PRNU_array
    elif source == 'DNFP':
      norm = col.LogNorm()
      data = self.DNFP_array

    ax = fig.add_subplot(projection=self.wcs)
    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth('3') 
    img = ax.imshow(data,cmap='jet' , norm = norm)
    plt.colorbar(img)
    ax.set_title(f'{source} \nRequested center : {self.name}')
    ax.grid(False)
    if download:
        fig.savefig(f"{source}.png", format = 'png')
    return fig,ax

  def show_hist(self, source = 'Readout',bins = None,figsize=(15,8)):
    fig, ax = plt.subplots(1,1,figsize=figsize)

    if source =='Readout':
      data  = self.readout.ravel()
    elif source == 'Sky':
      data = self.sky_array.ravel()
    elif source == 'DC':
      data = self.DC_array.ravel()
    elif source == 'QE':
      data = self.qe_array.ravel()
    elif source =='Bias':
      data = (self.bias_array + self.DC_array).ravel()
    elif source == 'PRNU':
      data = self.PRNU_array.ravel()
    elif source == 'DNFP':
      data = self.DNFP_array.ravel()

    if bins is None:
      bins  = np.linspace(data.min(), data.max(), 20)
    ax.hist(data, bins = bins)
    ax.set_title(f'{source} histogram')
    ax.set_ylabel('Count')
    ax.set_yscale('log')
    return fig, ax

  def writeto(self,name):
    """
    Function for downloading a fits file of simulated field
    """
    if np.all(self.image) !=None:
      hdu = fits.PrimaryHDU(self.charge, header = self.header)
      hdu.wcs = self.wcs
      hdul = fits.HDUList([hdu])
      hdul.writeto(f'{name}',overwrite= True)
    else:
      print("Generate PSF")
      
def test():

    ra      = [0]
    dec     = [0]
    mag_nuv = [9]
    df      = pd.DataFrame(zip(ra,dec,mag_nuv), columns= ['ra', 'dec', 'mag_nuv'])
    cols = {'RAJ2000' : 'ra', 'DEJ2000' :'dec', 'F275W' : 'mag_nuv'}
    psf = PSF( df = df,  axis='Off',cols = cols ,mode='HCIPy',exp_time = 1.)
    psf()
    fig, ax = psf.show_image(download=True)
    return fig, ax

