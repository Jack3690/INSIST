from photutils import aperture as aper
from photutils.aperture import aperture_photometry
from matplotlib import colors as col
import matplotlib.pyplot as plt
from astropy.io import fits
from .simulation import *

data_path = Path(__file__).parent.joinpath()

class Analyzer(Imager):
  def __init__(self, df = None, coords = None,tel_params=None,exp_time = 100,
               n_x = 1000, n_y = 1000, plot = False,  **kwargs):
      
      super().__init__(df=df, coords=coords, tel_params=tel_params,
                       exp_time=exp_time, n_x = n_x,  n_y = n_y,
                       plot = plot, **kwargs)
      """
      A class to visualize and analyze the simulated image
      
      Parameters
      ----------
      Imager.init() 

      Returns
      -------
      None.

      """
  def __call__(self,det_params= None, n_stack =1, stack_type ='median', 
               do_photometry =True, fwhm = None):
    super().__call__(det_params = det_params, n_stack = n_stack, 
                     stack_type = stack_type)
    """
    
    Performs sim simulation and sim Photometry
    
    Imager.call()

    do_photometry : Bool, Default : True
                    Do Aperture Photometry 
    """
    
    if do_photometry and len(self.df)>1:
      self.data_jy, self.phot_table = self.photometry(self.digital.astype(float),
                                                   self.wcs,self.df, fwhm)
  def photometry(self,data,wcs,df, fwhm):
      if fwhm is None:
          fwhm = self.pixel_scale*3

      c          = SkyCoord(df['ra'], df['dec'],unit=u.deg)
      pix        = wcs.world_to_array_index(c)

      position   = [(i,j) for i,j in zip(pix[1],pix[0])]

      self.aps   = aper.CircularAperture(position, r= fwhm/self.pixel_scale)
      ap_pix     = np.count_nonzero(self.aps.to_mask()[0])
      self.bags  = aper.CircularAnnulus(position, r_in = fwhm/self.pixel_scale, 
                                        r_out = 3*fwhm/self.pixel_scale)
      bag_pix    = np.count_nonzero(self.bags.to_mask()[0])
     
      phot_table      = aperture_photometry(data, [self.aps, self.bags])

      phot_table['sky_flux'] = phot_table['aperture_sum_1']*(ap_pix/bag_pix)
      phot_table['flux']     = phot_table['aperture_sum_0'].value - phot_table['sky_flux'].value
      phot_table['flux_err'] = np.sqrt( phot_table['flux'].value  + phot_table['sky_flux'].value )
  
      phot_table['SNR']      = phot_table['flux'].value/ phot_table['flux_err'].value
      phot_table['mag_in']   = df['mag'].values
      
      if len(phot_table)>3:
        if phot_table['SNR'].max()>=5:
            temp_table = phot_table[phot_table['SNR']>5]
        else:
            temp_table = phot_table.sort('SNR', ascending = False)
      
        zero_p_flux = 0
        for i in range(3):
          zero_p_flux += temp_table['flux'].value[i]/pow(10,
                                        -0.4*temp_table['mag_in'].value[i])
        zero_p_flux/=3
        print('Estimated Zero point using 3 stars')
        
      elif len(phot_table)>0:
          zero_p_flux = 0
          for i in range(len(phot_table)):
            zero_p_flux += phot_table['flux'].value[i]/pow(10,
                                         -0.4*phot_table['mag_in'].value[i])
            zero_p_flux/=len(phot_table)
      else:
        zero_p_flux = self.zero_flux*self.exp_time
        
      data_jy= data*(3631/zero_p_flux)
      
      self.header['EXPTIME'] = self.exp_time
      self.header['ZPT']     = zero_p_flux
      self.header['BUNIT']   = 'DN'
      phot_table['mag_out']  = -2.5*np.log10(phot_table['flux']/zero_p_flux)
      phot_table['mag_err']  = 1.082/phot_table['SNR']

      self.phot_table = phot_table

      return data_jy, phot_table

  def show_field(self,figsize=(10,10)):
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

    wcs = self.create_wcs(self.n_x, self.n_y, self.ra,self.dec, self.pixel_scale)
    
    ra_max, dec_max = wcs.array_index_to_world_values(self.n_y,0)
    ra_min, dec_min = wcs.array_index_to_world_values(0,self.n_x)
    
    # Cropping Dataframe based on FoV
    ra_min_cut  = (self.df['ra']>ra_min) 
    ra_max_cut  = (self.df['ra']<ra_max)

    df = self.df[ ra_min_cut & ra_max_cut ]

    dec_min_cut = (df['dec']>dec_min)
    dec_max_cut = (df['dec']<dec_max)
    df      = df[ dec_min_cut & dec_max_cut ]

    fov_x  = ra_max - ra_min
    fov_y  = abs(dec_max - dec_min)

    fig, ax = plt.subplots(1,1,figsize=figsize)
    ax.scatter(df['ra'],df['dec'],marker='.',color='black')
    ax.set_title(f"""Requested Center : {self.name} | {len(df)} sources
    Fov(RA) : {fov_x} | Fov(Dec) : {fov_y} """)
    ax.invert_xaxis()
    ax.set_xlabel('RA (Degrees)')
    ax.set_ylabel('Dec (Degrees)')
    return fig,ax

  def show_image(self, source = 'Digital', fig = None, ax = None,cmap = 'jet', 
                 figsize = (15,10), download = False, show_wcs = True):
    """
    Function for plotting the simulated field image

    Source: str,
            Choose from
                         'Digital' : Final digial image
                         'Charge'  : electrons, Light(Source + sky) + Dark Current + Noises
                         'Source'  : Source + Sky + Noises
                         'Sky'     : Sky + shot_noise
                         'DC'      : Dark Current + DNFP
                         'QE'      : Quantum efficiency fluctuation across detector
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
    if np.all(self.image) !=None :
        if fig is None or ax is None:
            fig = plt.figure(figsize = figsize)
        norm = None
        
        if source == 'Digital':
          data  = self.digital
          norm = col.LogNorm()
        elif source =='Charge':
          data  = self.charge
          norm = col.LogNorm()
        elif source =='Source':
          data  = self.light_array
          norm = col.LogNorm()
        elif source == 'Sky':
          data = self.sky_photoelec
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
        elif source == 'QN':
          data = self.QN_array
        else:
          print("Invalid Input")
          return None
    
        if show_wcs:
            ax = fig.add_subplot(projection=self.wcs)
        else:
            ax = fig.add_subplot()
        ax.patch.set_edgecolor('black')  
        ax.patch.set_linewidth('3') 
        img = ax.imshow(data,cmap=cmap , norm = norm)
        plt.colorbar(img,ax = ax)
        ax.set_title(f'{source} \nRequested center : {self.name}')
        ax.grid(False)
        if download:
            fig.savefig(f"{source}.png", format = 'png')
        return fig,ax
    else:
        print("Run Simulation")
        
  def show_hist(self, source = 'Digital',bins = None,
                 fig = None, ax = None,figsize=(15,8)):
    """
    Function for plotting histogram of various stages of simulation

    Parameters
    ----------

    Source: str,
            Choose from
                      'Digital' : Final digial image
                      'Charge'  : electrons, Light(Source + sky) + Dark Current + Noises
                      'Source'  : Source + Sky + Noises
                      'Sky'     : Sky + shot_noise
                      'DC'      : Dark Current + DNFP
                      'QE'      : Quantum efficiency fluctuation across detector
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
   
    if np.all(self.image) !=None :
        if fig is None or ax is None: 
            fig, ax = plt.subplots(1,1,figsize=figsize)
    
        if source == 'Digital':
          data  = self.digital.ravel()
        elif source =='Charge':
          data  = self.charge.ravel()
          norm = col.LogNorm()
        elif source =='Source':
          data  = self.light_array
        elif source == 'Sky':
          data = self.sky_photoelec.ravel()
        elif source == 'DC':
          data = self.DC_array.ravel()
        elif source == 'QE':
          data = self.qe_array
        elif source =='Bias':
          data = (self.bias_array + self.DC_array).ravel()
        elif source == 'PRNU':
          data = self.PRNU_array.ravel()
        elif source == 'DNFP':
          data = self.DNFP_array.ravel()
        elif source == 'QN':
          data = self.QN_array.ravel()
    
        if bins is None:
          bins  = np.linspace(data.min(), data.max(), 20)
        ax.hist(data, bins = bins)
        ax.set_title(f'{source} histogram')
        ax.set_ylabel('Count')
        ax.set_yscale('log')
        return fig, ax
    else:
        print("Run Simulation")

  def getImage(self,source = 'Digital'):
    """
    Function of retrieving image array at different stages of simulation.

    Parameters
    ----------

    Source: str,
        Choose from
                      'Digital' : Final digial image
                      'Charge'  : electrons, Light(Source + sky) + Dark Current + Noises
                      'Source'  : Source + Sky + Noises
                      'Sky'     : Sky + shot_noise
                      'DC'      : Dark Current + DNFP
                      'QE'      : Quantum efficiency fluctuation across detector
                      'Bias'    : Charge offset
                      'PRNU'    : Photon Response Non-Uniformity
                      'DNFP'    : Dark Noise Fixed Pattern
                      'QN'      : Quantization Noise
    """
      
    if np.all(self.image) !=None :
        if source == 'Digital':
            data  = self.digital
        elif source =='Charge':
            data  = self.charge
        elif source == 'Sky':
            data = self.sky_photoelec
        elif source == 'DC':
            data = self.DC_array
        elif source == 'QE':
            data = self.qe_array
        elif source =='Bias':
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

  def writeto(self,name,source = 'Digital', user_source = None):

    """
    Function for downloading a fits file of simulated field image

    Parameters
    ----------
    name : str
           filename, Example : simulation.fits

    Source: str,
        Choose from
                      'Digital' : Final digial image
                      'Charge'  : electrons, Light(Source + sky) + Dark Current + Noises
                      'Source'  : Source + Sky + Noises
                      'Sky'     : Sky + shot_noise
                      'DC'      : Dark Current + DNFP
                      'QE'      : Quantum efficiency fluctuation across detector
                      'Bias'    : Charge offset
                      'PRNU'    : Photon Response Non-Uniformity
                      'DNFP'    : Dark Noise Fixed Pattern
                      'QN'      : Quantization Noise

    user_source : numpy.ndarray
                  2D numpy array user wants to save as FITS
    """
    if np.all(self.image) !=None :
      if user_source is not None and type(user_source)==np.ndarray:
        data = user_source
      elif source == 'Digital':
        data  = self.digital
        norm = col.LogNorm()
      elif source =='Charge':
        data  = self.charge
        norm = col.LogNorm()
      elif source =='Source':
          data  = self.light_array
      elif source == 'Sky':
        data = self.sky_photoelec
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
      elif source == 'QN':
        data = self.QN_array
    
      else:
        print(f"{source} is not a valid source")

      hdu = fits.PrimaryHDU(data, header = self.header)
      hdu.wcs = self.wcs
      hdul = fits.HDUList([hdu])
      hdul.writeto(f'{name}',overwrite= True)
    else:
      print("Run Simulation")