from .analysis import *

class Mosaic(Analyzer):
  """
  A class to split bigger images to tiles and stitch them together

  """

  def __init__(self, df = None, coords = None,tel_params=None,exp_time = 100,
               n_x = 1000, n_y  = 1000, mos_n = 1, mos_m = 1, **kwargs):
    
      """
      Analyzer.init()

      Parameters
      ----------
      mos_n : int,
              number of tiles in RA direction

      mos_m : int,
              number of tiles in Dec direction
      """

      super().__init__(df=df, coords=coords, exp_time = exp_time, n_x = n_x, 
                    n_y = n_y, tel_params=tel_params, **kwargs)
      
      self.n = mos_n
      self.m = mos_m

      self.mosaic_img, _ = self.init_image_array(return_img = True)
      patch_width_mid    = self.n_pix_sub//2

      self.mosaic_img=self.mosaic_img[2*patch_width_mid:-2*patch_width_mid,
                  2*patch_width_mid:-2*patch_width_mid]

      self.mosaic_name = self.name    
      self.mosaic_ra   = self.ra
      self.mosaic_dec  = self.dec  
      self.mosaic_n_x  = n_x
      self.mosaic_n_y  = n_y
      self.mosaic_df   = df
      self.mosaic_wcs  = self.create_wcs(self.mosaic_n_x, self.mosaic_n_y, 
                                      self.mosaic_ra,self.mosaic_dec,
                                      self.pixel_scale)

      self.df_split(df, self.n, self.m, n_x, n_y, self.mosaic_wcs)

  def df_split(self,df,n,m, n_x, n_y, wcs):

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

    x_bins = np.linspace(0,n_x,n+1)
    y_bins = np.linspace(0,n_y,m+1)

    x_cens = 0.5*(x_bins[:-1] + x_bins[1:]) - 1
    y_cens = 0.5*(y_bins[:-1] + y_bins[1:]) - 1
    
    cens = wcs.array_index_to_world_values(y_cens,x_cens)
    ra_cens  = cens[0]
    ra_cens  = np.where(np.round(ra_cens,1)==360,ra_cens-360,ra_cens)
    dec_cens = cens[1]

    self.ras  = ra_cens
    self.decs = dec_cens
                   
  def __call__(self,det_params = None, n_stack = 1, stack_type = 'median', 
               do_photometry = True, fwhm = None):
    
    """
      Analyzer.call()

      Calls the Analyzer clssiteratively to generate image tiles and stitches
      tiles together

    """

    # Flags
    self.mosaic_shot_noise = self.shot_noise
    self.mosaic_QE         = self.QE         
    self.mosaic_sky        = self.sky
    self.mosaic_PRNU       = self.PRNU
    self.mosaic_DC         = self.DC 
    self.mosaic_DNFP       = self.DNFP
    self.mosaic_QN         = self.QN
    self.mosaic_cosmic_rays= self.cosmic_rays

    n_x = self.n_x//self.n
    n_y = self.n_y//self.m
    
    for i in range(self.n):
      for j in range(self.m):

        df         = self.mosaic_df
        coords     = (self.ras[i],self.decs[j])
        exp_time   = self.exp_time
        tel_params = self.tel_params
        
        super().__init__(df = df, coords=coords, exp_time = exp_time, n_x = n_x, 
                    n_y = n_y, tel_params=tel_params)
        
        self.shot_noise = self.mosaic_shot_noise
        self.QE         = self.mosaic_QE         
        self.sky        = self.mosaic_sky
        self.PRNU       = self.mosaic_PRNU
        self.DC         = self.mosaic_DC 
        self.DNFP       = self.mosaic_DNFP
        self.QN         = self.mosaic_QN
        self.cosmic_rays= self.mosaic_cosmic_rays
        super().__call__(det_params = det_params, n_stack = n_stack, 
                          stack_type = stack_type)
        self.mosaic_img[n_y*j:n_y*(j+1):,n_x*i:n_x*(i+1)] = self.digital

    self.init_df(self.mosaic_n_x, self.mosaic_n_y)
    
    if do_photometry:
      self.data_jy, self.phot_table = self.photometry(self.mosaic_img.astype(float),
                                                   self.mosaic_wcs,self.df, fwhm)


  def show_mosaic(self, fig = None, ax = None,cmap = 'jet', 
                 figsize = (15,10), download = False, show_wcs = True):
    """
    Function for plotting the simulated Mosaic

    fig : matplotlib.pyplot.figure
          User defined figure
    ax  : matplotlib.pyplot.axes
          User defined axes
    cmap : str,
           matplotlib.pyplot colormap
    figsize  : tuple
    download : bool
    show_wcs : bool
               If true adds WCS projection to the image
    Returns
    -------

    fig, ax
    """
    if np.all(self.image) !=None :
        if fig is None or ax is None:
            fig = plt.figure(figsize = figsize)
        norm = None
        
        data  = self.mosaic_img
        norm = col.LogNorm()
    
        if show_wcs:
            ax = fig.add_subplot(projection=self.mosaic_wcs)
        else:
            ax = fig.add_subplot()
        ax.patch.set_edgecolor('black')  
        ax.patch.set_linewidth('3') 
        img = ax.imshow(data,cmap=cmap , norm = norm)
        plt.colorbar(img,ax = ax)
        ax.set_title(f'Mosaic \nRequested center : {self.mosaic_name}')
        ax.grid(False)
        if download:
            fig.savefig(f"mosaic.png", format = 'png')
        return fig,ax
    else:
        print("Run Simulation")

  def writemosaic(self,name):
    """
    Function for downloading a fits file of simulated image

    Parameters
    ----------
    name : str
           filename, Example : simulation.fits
    """
    if np.all(self.image) !=None :
      data  = self.mosaic_img
      hdu = fits.PrimaryHDU(data, header = self.mosaic_wcs.to_header())
      hdu.wcs = self.wcs
      hdul = fits.HDUList([hdu])
      hdul.writeto(f'{name}',overwrite= True)
    else:
      print("Run Simulation")

