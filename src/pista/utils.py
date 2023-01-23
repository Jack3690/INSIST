import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb

from scipy.interpolate import interp1d, interp2d
from scipy.integrate import trapz
from pathlib import Path
from astropy.modeling import models
from scipy.constants import c

sb.set_style('darkgrid')
matplotlib.rcParams['font.size']=12
matplotlib.rcParams['figure.figsize']=(10,10)

data_path = Path(__file__).parent.joinpath()

def bandpass(wav, flux, inputs, plot = True, fig = None, ax = None):
  """
  Function to convolve response functions

  Parameters
  ----------
  wav : numpy.ndarray 
        wavelenth in angstrom
  flux: numpy.ndarray
        flux normalized to [0,1]
  plot: bool,
        If true shows plots with input and convolved response functions

  fig : matplotlib.pyplot.figure
        User defined figure
  ax  : matplotlib.pyplot.axes
        User defined axes

  Returns
  -------

  fig, ax, data, params

  data : tuple,
        (wavelenth array, flux_array, convolved flux array)
  
  params: tuple,
          (effective wavelength, integrated flux, Effective Width)
  
  """
  lambda_   =  wav 
  flux_AB   = flux 

  if plot:
    if fig is None or ax is None:
        fig, ax = plt.subplots(1,1, figsize = (12,8))
    ax.plot(lambda_ ,flux_AB/flux_AB.max(),label = r'$F(\lambda)$', alpha = 0.7 )

  R_eff = 1

  for i in inputs:
    file_name = i.split(',')[0]
    n         = float(i.split(',')[1])
    f_max     = float(i.split(',')[2])
    
    filt_dat  = np.loadtxt(file_name)
    wav  = filt_dat[:,0]
    flux = filt_dat[:,1]
    
    if flux.max()>1:
      flux/=f_max

    indices  = np.where( (wav>lambda_ [0]) & (wav<lambda_[-1]))
    wav_new  = wav[indices]
    flux_new = flux[indices]
    
    wav_new  =  [lambda_ [0]] + [wav_new[0]- 1] + list(wav_new) + [wav_new[-1]+ 1] + [lambda_[-1]]
    flux_new =  [0]           + [0]             + list(flux_new) +   [0]           + [0]
    
    wav_new  = np.array(wav_new)
    flux_new = np.array(flux_new)

    filter_func = interp1d(wav_new,flux_new)

    flux_out    = filter_func(lambda_)

    R_eff      *= pow(flux_out,n)
    
    if plot:
      ax.plot(lambda_ ,flux_out/flux_out.max(),label=f"{file_name.split('/')[-1][:-4]}x{n}", alpha = 0.7)

  conv_flux     = R_eff*flux_AB
  conv_flux_Jy  = R_eff*flux_AB*lambda_**2*3.34e4
 
  lambda_phot = trapz(lambda_**2*conv_flux,lambda_)/trapz(lambda_*conv_flux,lambda_)

  W_eff       = trapz(R_eff,lambda_)/R_eff.max()

  int_flux    = trapz(lambda_*conv_flux,lambda_)/trapz(lambda_*R_eff, lambda_)
  int_flux_Jy = trapz(lambda_*conv_flux_Jy,lambda_)/trapz(lambda_*R_eff, lambda_)

  # Comparing to a square filter with same width
  R_sq = np.where(R_eff>0,1,0)

  flux_ratio = trapz(R_eff,lambda_)/trapz(R_sq,lambda_)

  data       =  lambda_, conv_flux, R_eff
  params     =  lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio

  if plot:
    ax.plot(lambda_,conv_flux/conv_flux.max(),label = 'Convloved Flux', linewidth = 5)
    y = np.linspace(0,1)
    x = y*0 + lambda_phot
    label = r'$\lambda_{phot} = $' + f'{round(lambda_phot,3)}' + r' $\AA$'
    ax.plot(x,y,'--', color = 'black',label = label )
   
    ax.set_xlabel(r'$\AA$')
    ax.set_ylabel(r'Normalized Flux')
    fig.suptitle('Bandpass', fontsize = 20, y = 0.95)
    ax.legend()
    
  return fig, ax, data, params

def generate_psf(npix,sigma, function = 'Gaussian'):

    """
    Function for generating user defined PSF

    npix : int,
           number of pixels along one axis for pixel array
    
    sigma: float,
           standard deviation of the PSF in pixels

    function: str,
               type of PSF function
    
    Returns
    -------

    numpy.ndarray
    
    """
    x       = np.linspace(0,1000,npix)
    y       = x
    yy,xx   = np.meshgrid(x,y)
    if function == 'Gaussian':
        psf     = models.Gaussian2D(1,500,500,sigma,sigma)(xx,yy)
        psf/=psf.sum()
    np.save('user_defined_psf.npy',psf)
    return psf  

def spectra_to_mags(df,inputs):
  mags = []
  for row in df:
    wav  = row['wav']
    flux = row['flux']
    out           = bandpass(wav, flux,inputs = inputs, 
                                          plot = False)
    params         = out[3]
    int_flux_Jy    = params[2]
    ABmag          = -2.5*np.log10(int_flux_Jy/3631)

    if ABmag == np.nan:
      ABmag  = 100
    mags.append(ABmag)
  df['mag']  = mags

  return df
    