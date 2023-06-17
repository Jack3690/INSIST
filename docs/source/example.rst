**********
Quickstart
**********

Basic simulation using package data. PISTA package data includes filter profiles and PSF profiles for INSIST. 
This can be utilized for simulating images using INSIST specifications.

Initialization
==============

.. jupyter-execute::

  import pista as pt
  from pathlib import Path
  from astropy.table import Table
  import matplotlib.pyplot as plt
  
  data_path = Path(pt.__file__).parent.joinpath()
  %matplotlib inline

Source Catalog
==============
.. jupyter-execute::

  tab = Table.read(f'{data_path}/data/sample.fits')  # FITS Table
  df  = tab.to_pandas()             # PISTA requires pandas DataFrame

Telescope parameters
====================

.. jupyter-execute::

  tel_params = {
                  'aperture'       : 100,
                  'pixel_scale'    : 0.1,
                  'psf_file'       : f'{data_path}/data/PSF/INSIST/off_axis_poppy.npy',
                  'response_funcs' :  [ f'{data_path}/data/INSIST/UV/Filter.dat,1,100',    
                                        f'{data_path}/data/INSIST/UV/Coating.dat,5,100',   # 5 mirrors
                                        f'{data_path}/data/INSIST/UV/Dichroic.dat,2,100',  # 2 dichroics
                                      ],                                
                } 

Initialize Imager object
==========================

.. jupyter-execute::
  
  sim = pt.Imager(df = df,tel_params = tel_params, n_x = 500, n_y = 500, exp_time = 2400)
  sim.show_field(figsize=(12, 10);
  
Detector parameters
===================

.. jupyter-execute::

  det_params = {  'shot_noise' :  'Gaussian',
                  'G1'         :  1,
                  'PRNU_frac'  :  0.25/100,
                  'qe_response': [f'{data_path}/data/INSIST/UV/QE.dat,1,100'],
                  'RN'         :  3,
                  'T'          :  218,        
                  'DN'         :  0.01/100     
                }

Simulate Image
==============
.. jupyter-execute::

  sim(det_params = det_params, photometry = 'Aper')
  sim.show_image()
