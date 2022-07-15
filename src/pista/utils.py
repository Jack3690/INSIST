
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

def bandpass(wav, flux , inputs, plot = True, fig = None, ax = None):
  # AB system
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
    
    filt_dat  = np.loadtxt(f'{file_name}')
    wav  = filt_dat[:,0]
    flux = filt_dat[:,1:]
    
    if flux.max()>50:
      flux/=100
    flux = flux.prod(axis = 1)

    indices  = np.where( (wav>lambda_ [0]) & (wav<lambda_[-1]))
    wav_new  = wav[indices]
    flux_new = flux[indices]
    
    wav_new  =  [lambda_ [0]] + [wav_new[0]- 1] + list(wav_new) + [wav_new[-1]+ 1] + [lambda_[-1]]
    flux_new =  [0]           + [0]             + list(flux_new) +   [0]           + [0]

    wav_new  = np.array(wav_new)
    flux_new = np.array(flux_new)

    filter_func = interp1d(wav_new,flux_new)

    flux_out    = filter_func(lambda_)

    R_eff      *= flux_out**n
    
    if plot:
      ax.plot(lambda_ ,flux_out/flux_out.max(),label=f"{file_name.split('/')[-1][:-4]}x{n}", alpha = 0.7)

  conv_flux = R_eff*flux_AB
 
  lambda_phot = trapz(lambda_**2*conv_flux,lambda_)/trapz(lambda_*conv_flux,lambda_)

  W_eff      = trapz(R_eff,lambda_)/R_eff.max()

  int_flux   = trapz(lambda_*conv_flux,lambda_)/trapz(lambda_*R_eff, lambda_)

  data       =  lambda_, flux_AB, conv_flux
  params     =  lambda_phot, int_flux, W_eff

  if plot:
    ax.plot(lambda_,conv_flux/conv_flux.max(),label = 'Convloved Flux')
    y = np.linspace(0,1)
    x = y*0 + lambda_phot
    label = r'$\lambda_{phot} = $' + f'{round(lambda_phot,3)}' + r' $\AA$'
    ax.plot(x,y,'--', color = 'black',label = label )
   
    ax.set_xlabel(r'$\AA$')
    ax.set_ylabel(r'Normalized Flux')
    
  return fig, ax, data, params

def generate_psf(npix,sigma, function = 'Gaussian'):
    x       = np.linspace(0,1000,npix)
    y       = x
    yy,xx   = np.meshgrid(x,y)
    if function == 'Gaussian':
        psf     = models.Gaussian2D(1,500,500,sigma,sigma)(xx,yy)
        psf/=psf.sum()
    return psf

    
    