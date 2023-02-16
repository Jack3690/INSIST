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
 

.. code-block:: python
Source 
======
.. code-block:: python
    # Example FITS table
    tab = Table.read(f'{data_path}/data/sample.fits')  # FITS Table
    df  = tab.to_pandas()    # PISTA requires pandas DataFrame

.. code-block:: python
Telescope parameters
=====================
.. code-block:: python
  # Sample inputs from package data
  tel_params ={
            'aperture'       : 100,
            'pixel_scale'    : 0.1,
            'psf_file'       : f'{data_path}/data/PSF/INSIST/on_axis_hcipy.npy',
            'response_funcs' :  [ f'{data_path}/data/INSIST/UV/Filter.dat,1,100',    
                                  f'{data_path}/data/INSIST/UV/Coating.dat,5,100',   # 5 mirrors
                                  f'{data_path}/data/INSIST/UV/Dichroic.dat,2,100',   # 2 dichroics
                                  f'{data_path}/data/INSIST/UV/QE.dat,1,100'
                                ],                                
            } 
.. code-block:: python

