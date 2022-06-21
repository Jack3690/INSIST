from photutils import aperture as aper
from photutils.aperture import aperture_photometry
from matplotlib import colors as col
from .simulation import *

data_path = Path(__file__).parent.joinpath()

class Analyzer(Imager):
  def __init__(self, df = None, cols = None,
               psf_file = f'{data_path}/data/off_axis_hcipy.npy',
               exp_time = 100,n_pix = 2000, response_funcs = None,band = 'G'):
      super().__init__(df=df, cols=cols, exp_time = exp_time,
                       psf_file = psf_file, n_pix = n_pix, 
                       response_funcs = response_funcs, band = band )
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
      mode     : str, 'HCIPy' or 'Zemax', PSF patch generation method
      exp_time : float, Exposure time in seconds 
      fov      : float, Field of View in degrees
      
     

      Returns
      -------
      None.

      """

  def __call__(self,params= None, n_stack =1, stack_type ='median', do_photometry =True):
    super().__call__(params = params, n_stack = n_stack, stack_type = stack_type)
    """
    
    Performs PSF simulation and PSF Photometry
    
    TBW
    """
    
    if do_photometry: 
      c          = SkyCoord(self.df['ra'], self.df['dec'],unit=u.deg)
      data       = self.digital.astype(float)
      wcs        = self.wcs
      pix        = wcs.world_to_array_index(c)

      position   = [(i,j) for i,j in zip(pix[1],pix[0])]

      self.aps   = aper.CircularAperture(position, r=0.3/0.1)
      ap_pix     = np.count_nonzero(self.aps.to_mask()[0])
      self.bags  = aper.CircularAnnulus(position, r_in = 0.6/0.1, r_out = 1/0.1)
      bag_pix    = np.count_nonzero(self.bags.to_mask()[0])

      phot_table      = aperture_photometry(data, [self.aps, self.bags])

      phot_table['sky_flux'] = phot_table['aperture_sum_1']*(ap_pix/bag_pix)
      phot_table['flux']     = phot_table['aperture_sum_0'].value - phot_table['sky_flux'].value
      phot_table['flux_err'] = np.sqrt( phot_table['flux'].value  + phot_table['sky_flux'].value )
  
      phot_table['SNR']      = phot_table['flux'].value/ phot_table['flux_err'].value
      phot_table['mag_in']   = self.df['mag_nuv'].values
      
      zero_p_flux = 0
      for i in range(3):
        zero_p_flux += phot_table['flux'].value[i]/pow(10,-0.4*phot_table['mag_in'].value[i])
      zero_p_flux/=3

      self.header['EXPTIME'] = self.exp_time
      self.header['ZPT']     = zero_p_flux
      self.header['BUNIT']   = 'DN'
      phot_table['mag_out']  = -2.5*np.log10(phot_table['flux']/zero_p_flux)
      phot_table['mag_err']  =1.087/phot_table['SNR']

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

  def show_image(self, source = 'Digital', figsize = (15,10), download = False):
    """
    Function for plotting the simulated field with PSFs

    Returns
    -------
    fig, ax
    """

    fig = plt.figure(figsize = figsize)
    norm = None
    
    if source == 'Digital':
      data  = self.digital
      norm = col.LogNorm()
    elif source =='Charge':
      data  = self.charge
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

  def show_hist(self, source = 'Digital',bins = None,figsize=(15,8)):
    fig, ax = plt.subplots(1,1,figsize=figsize)

    if source == 'Digital':
      data  = self.digital.ravel()
    elif source =='Charge':
      data  = self.charge.ravel()
      norm = col.LogNorm()
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

  def writeto(self,name,source = 'Digital', user_source = None):
    """
    Function for downloading a fits file of simulated field
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
      print("Generate PSF")