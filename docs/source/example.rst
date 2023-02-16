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

.. text:: print(data_path)
  :alt: output text


.. code-block:: python

  tab = Table.read(f'{data_path}/data/sample.fits')  # FITS Table
  df  = tab.to_pandas()             # PISTA requires pandas DataFrame
