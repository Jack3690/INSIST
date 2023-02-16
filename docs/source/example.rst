******
Usage
******

Basic simulation using package data.

Initialization
==============
.. code-block:: python

  import pista as pis
  from pathlib import Path
  from astropy.table import Table
  data_path = Path(pis.__file__).parent.joinpath()

Source
======

.. code-block:: python

  tab = Table.read(f'{data_path}/data/sample.fits')  # FITS Table
  df  = tab.to_pandas()             # PISTA requires pandas DataFrame

Telescope Parameters
====================

.. code-block:: python

  tel_params = {
                  'aperture'       : 100,
                  'pixel_scale'    : 0.1,
                  'psf_file'       : f'{data_path}/data/PSF/INSIST/on_axis_hcipy.npy',
                  'response_funcs' :  [ f'{data_path}/data/INSIST/UV/Filter.dat,1,100',    
                                        f'{data_path}/data/INSIST/UV/Coating.dat,5,100',   # 5 mirrors
                                        f'{data_path}/data/INSIST/UV/Dichroic.dat,2,100',  # 2 dichroics
                                        f'{data_path}/data/INSIST/UV/QE.dat,1,100'
                                      ],                                
                } 

Initializing Imager Object
==========================

.. code-block:: python

  sim = pis.Imager(df = df,tel_params = tel_params, n_x = 500, n_y = 500, exp_time = 2400)
  sim.show_field()
  plt.show()

.. jupyter-execute::