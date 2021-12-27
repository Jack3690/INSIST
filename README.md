# Indian Spectroscopic and Imaging Space Telescope 

## Introduction

This is a repository dedicated to data simulations for the Indian Spectroscopic and Imaging Space Telescope (INSIST) project. It contains Jupyter Notebooks which can be hosted on cloud platforms such as [Google Colab](https://colab.research.google.com/notebooks/intro.ipynb?utm_source=scs-index), [Binder](https://mybinder.org/), and [Gradient](https://gradient.run/notebooks), and webtools which can be accessed using [Binder](https://mybinder.org/).

This repository contains codes for the following:

*  PSFs using catalogs and Casjobs
*  PSFs based on different Telescope apertures and obstructions

# Web Tools
* PSF Simulation Tool : [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jack3690/INSIST/main?urlpath=%2Fvoila%2Frender%2FPSF_Simulator_Tool.ipynb%3Fvoila-theme%3Ddark) (Work in progress)
 
## Notebooks

* PSF Simulator : [PSF_Simulator.ipynb](https://github.com/Jack3690/INSIST/blob/main/PSF_Simulator.ipynb)
* UV Stellar Catalog Generator : [UV_Stellar_Catalog.ipynb](https://github.com/Jack3690/INSIST/blob/main/UV_Stellar_Catalog.ipynb)


# Usage

## [PSF_Simulator.ipynb](https://github.com/Jack3690/INSIST/blob/main/PSF_Simulator.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/PSF_Simulator.ipynb)
This notebook can be used to understand how to generate Point Spread Functions of sources in a given field using user defined catalogs and CasJobs.

### Sections

#### Single PSF Generator

In this section we use Astropy to generate Gaussian and Airy disk PSFs based on FWHM and pixel scale, normalized based on ABmag
##PSF Generator using Source Catalog

![plot](./doc/SPG.png) 

#### PSF using source Catalog and CasJobs

In this section we explore how to use a catalog which is either directly acquired or acquired from MAST Casjobs to generate PSF distribuition of a field.


## [UV_Stellar_Catalog.ipynb](https://github.com/Jack3690/INSIST/blob/main/UV_Stellar_Catalog.ipynb) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Jack3690/INSIST/blob/main/UV_Stellar_Catalog.ipynb)

# Conclusion/Disclaimer

If you have any questions or suggestions for improvements to this repo,
please contact the owners of the repository.

This is not an official product.


## References
