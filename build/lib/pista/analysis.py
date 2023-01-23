from photutils import aperture as aper
from photutils.aperture import aperture_photometry
from matplotlib import colors as col
import matplotlib.pyplot as plt
from astropy.io import fits

from photutils import aperture as aper
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
from photutils.isophote import Ellipse

from photutils.detection import IRAFStarFinder, DAOStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup, DAOPhotPSFPhotometry

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.coordinates import SkyCoord
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy import units as u

import numpy as np

class Analyzer(object):
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
  def __call__(self, df = None, wcs = None, data = None,
               photometry = None, detect_sources = False, fwhm = None):
    """
    
    Performs sim simulation and sim Photometry
    
    Imager.call()

    do_photometry : Bool, Default : True
                    Do Aperture Photometry 
    """
    
    if photometry == 'Aper':
      self.aper_photometry(data, wcs, df, fwhm,
                           detect_sources)
    elif photometry == 'PSF':
      self.psf_photometry(data, wcs, df, fwhm,
                          detect_sources)

  def aper_photometry(self,data,wcs,df, fwhm, detect):
      if fwhm is None:
          fwhm = 3
      else:
        fwhm/=self.pixel_scale

      if detect:

        mean, median, std = sigma_clipped_stats(data, sigma=3.0) 
        daofind = DAOStarFinder(fwhm=fwhm, threshold = 5.*std)
        sources = daofind(data - median)
        positions = np.transpose((sources['xcentroid'], sources['ycentroid'])) 

      else:
        c          = SkyCoord(df['ra'], df['dec'],unit=u.deg)
        pix        = wcs.world_to_array_index(c)
        positions   = [(i,j) for i,j in zip(pix[1],pix[0])]
      
      self.aps   = aper.CircularAperture(positions, r= fwhm)
      ap_pix     = np.count_nonzero(self.aps.to_mask()[0])
      self.bags  = aper.CircularAnnulus(positions, r_in = fwhm, 
                                        r_out = 3*fwhm)
      bag_pix    = np.count_nonzero(self.bags.to_mask()[0])
     
      phot_table      = aperture_photometry(data, [self.aps, self.bags])

      phot_table['sky_flux'] = phot_table['aperture_sum_1']*(ap_pix/bag_pix)
      phot_table['flux']     = phot_table['aperture_sum_0'].value - phot_table['sky_flux'].value
      phot_table['flux_err'] = np.sqrt( phot_table['flux'].value  + phot_table['sky_flux'].value )
  
      phot_table['SNR']      = phot_table['flux'].value/ phot_table['flux_err'].value

      if not detect:
        phot_table['mag_in']   = df['mag'].values 
        phot_table['ra']       = df['ra'].values
        phot_table['dec']      = df['dec'].values 
        if len(phot_table)>3: 
          zero_p_flux = 0

          for i in range(3):
            zero_p_flux += phot_table['flux'].value[i]/pow(10,
                                          -0.4*phot_table['mag_in'].value[i])
          zero_p_flux/=3
          phot_table['mag_out']  = -2.5*np.log10(phot_table['flux']/zero_p_flux)
          phot_table['mag_err']  = 1.082/phot_table['SNR']
      
      self.header['EXPTIME'] = self.exp_time
      self.header['BUNIT']   = 'DN'
      self.phot_table = phot_table
      

  def psf_photometry(self,data,wcs,df, fwhm, detect):
      if fwhm is None:
          fwhm = 3
      else:
        fwhm/=self.pixel_scale

      mean, median, std = sigma_clipped_stats(data, sigma=3.0) 
      sigma_psf = fwhm*gaussian_fwhm_to_sigma
      fitter = LevMarLSQFitter()
      psf_model = IntegratedGaussianPRF(sigma=sigma_psf)

      photometry = DAOPhotPSFPhotometry(  crit_separation = 3,
                                          threshold = mean + 5*std,
                                          fwhm = fwhm,
                                          psf_model=psf_model,
                                          fitter=LevMarLSQFitter(),
                                          fitshape=(5, 5))

      result_tab = photometry(image=data)

      self.phot_table = result_tab

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
 
    # Cropping Dataframe based on FoV
    x_min_cut = (self.df['x'] > self.n_pix_sub - 1) 
    x_max_cut = (self.df['x'] < self.n_y_sim   - self.n_pix_sub - 1)

    df        = self.df[ x_min_cut & x_max_cut ]

    y_min_cut = (self.df['y']>self.n_pix_sub - 1) 
    y_max_cut = (self.df['y']<self.n_x_sim - self.n_pix_sub -1)

    df = df[y_min_cut & y_max_cut]

    ra_max  = df['ra'].max()
    ra_min  = df['ra'].min()
    dec_max = df['dec'].max()
    dec_min = df['dec'].min()

    fov_x  = np.round(ra_max - ra_min,5)
    fov_y  = np.round(abs(dec_max - dec_min),5)

    fig, ax = plt.subplots(1,1,figsize=figsize)
    ax.scatter(df['ra'],df['dec'],marker='.',color='black')
    ax.set_title(f"""Requested Center : {self.name} | {len(df)} sources
    Fov(RA) : {fov_x} | Fov(Dec) : {fov_y} """)
    ax.invert_xaxis()
    ax.set_xlabel('RA (Degrees)')
    ax.set_ylabel('Dec (Degrees)')
    return fig,ax

  def show_image(self, source = 'Digital', fig = None, ax = None, cmap = 'jet', 
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
            if show_wcs:
              ax = fig.add_subplot(projection=self.wcs)
            else:
              ax = fig.add_subplot()

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
