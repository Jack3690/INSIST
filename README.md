# Indian Spectroscopic and Imaging Space Telescope [![Documentation Status](https://readthedocs.org/projects/insist/badge/?version=latest)](https://insist.readthedocs.io/en/latest/?badge=latest)[![PyPI - Python Version](https://img.shields.io/pypi/v/insist-pista.svg)](https://pypi.org/project/insist-pista/) [![PyPI - Python Version](https://img.shields.io/pypi/pyversions/insist-pista)](https://pypi.org/project/insist-pista/)  ![PyPI - Downloads](https://img.shields.io/pypi/dm/insist-pista)

## Introduction

This is a repository dedicated to data simulations for the Indian Spectroscopic and Imaging Space Telescope (INSIST) project. It contains Jupyter Notebooks which can be hosted on cloud platforms such as [Google Colab](https://colab.research.google.com/notebooks/intro.ipynb?utm_source=scs-index), [Binder](https://mybinder.org/), and [Gradient](https://gradient.run/notebooks), and webtools which can be accessed using [Binder](https://mybinder.org/).

# Packages
## PISTA : Python Image Simulation and Testing Application
A python package aimed at simulating astronomical images. The routine simulates individual stars and adds different noises. The input parameter space is designed to inculcate observational parameters, telescope parameters and detector characteristics.


### Installation
```
pip install insist-pista
```

# Web Tools
* PISTA Webtool                 : [Streamlit](https://jack3690-insist-webtools-pista-webtool-d79yxm.streamlitapp.com/)
* Exposure Time Calculator Tool : [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jack3690/INSIST/main?urlpath=%2Fvoila%2Frender%2FExposure_Time_Calculator_Tool.ipynb%3Fvoila-theme%3Ddark) (Work in progress)

 
## Notebooks

* Image Simulation             : [Image_Simulation.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/Image_Simulation.ipynb)
* PISTA                        : [PISTA.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PISTA.ipynb)
* PSF Simulator                : [PSF_Simulator.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PSF_Simulator.ipynb)
* PSF Analysis                 : [PSF_Analysis.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PSF_Analysis.ipynb)
* UV Stellar Catalog Generator : [UV_Stellar_Catalog.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/UV_Stellar_Catalog.ipynb)
* Exposure Time Calculator     : [Exposure_Time_Calculator.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/Exposure_Time_Calculator)

# Usage

## [Image_Simulation.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/Image_Simulation.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/Image_Simulation.ipynb)

This notebook provides step by step instructions on how to use PISTA package for resolved stellar population simulation.

## [PISTA.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PISTA.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/PISTA.ipynb)

![](https://github.com/Jack3690/INSIST/blob/main/docs/pista_flow.png?raw=True) 
This notebook  expands on the low-level framework of PISTA package. 

## [PSF_Simulator.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PSF_Simulator.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/PSF_Simulator.ipynb)

In this notebook, we explore how to simulate Points Spread Functions (PSF) with various parameter inputs. The objective is to able to generate a simulated image of a field which can be compared with actual data. This will be helpful at various stages of INSIST developement such as Pipiline validation, science cases pre-observation simulation etc.

### Single PSF Generator

In this section we use Astropy to generate Gaussian and Airy disk PSFs based on FWHM and pixel scale, normalized based on ABmag

![](https://github.com/Jack3690/INSIST/blob/main/docs/SPG.png?raw=True) 

## PSF Simulation using HCIPy

In this section we explore how to use HCIPy for generating PSF for different telescopes
![](https://github.com/Jack3690/INSIST/blob/main/docs/psf.png?raw=True) 

## [PSF_Analysis.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/PSF_Analysis.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/PSF_Analysis.ipynb)

This notebook contains Python routines used for comparing off axis and on axis PSFs generated using Zemax and HCIPy. The aim is to quantitatively study how presence of an on-axis secondary modifies the PSF, and how it would affect the expected science cases.
![](https://github.com/Jack3690/INSIST/blob/main/docs/off_axis_vs_on_axis.png?raw=True) 


## [UV_Stellar_Catalog.ipynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/UV_Stellar_Catalog.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/UV_Stellar_Catalog.ipynb)

This notebook contains routines which allows users to predict UV fluxes of sources using their PanSTARRS optical band magntidues through SED Fitting. The objective is to able to generate a catalog of sources in UV band using optical band magnitudes as input. For SED fitting we utilize Kurucz models http://kurucz.harvard.edu/grids.html

### Convolving Filters with Stellar models
![](https://github.com/Jack3690/INSIST/blob/main/docs/filter_conv.png?raw=True) 

### SED fitting validation using HST PHAT M31 Survey
![](https://github.com/Jack3690/INSIST/blob/main/docs/sed_fit.png?raw=True)

## [Exposure_Time_Calculatoripynb](https://github.com/Jack3690/INSIST/blob/main/notebooks/Exposure_Time_Calculator.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/notebooks/Exposure_Time_Calculator)

This notebook contains routines which allows users to calculate exposure time for a range of telescope surveys such as GALEX, UVIT, LSST etc with magnitude and SNR as inputs.

The package includes a GUI designed using PyQT5. (Work in Progress)

![](https://github.com/Jack3690/INSIST/blob/main/docs/PISTA.png?raw=True) 
# Conclusion/Disclaimer

Please add the following acknowledgment if you use our package in your work.

"This work has made use of Python Image Simulation and Testing Application (PISTA) developed as part of the INdian Spectroscopic and Imaging Space Telescope (INSIST) project."

If you have any questions or suggestions for improvements to this repo,
please contact the owners of the repository.


## References
* https://docss.hcipy.org/stable/
* https://esa.gitlab.io/pyxel/
