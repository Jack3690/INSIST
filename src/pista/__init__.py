__version__ = "1.0.27"
__author__ = 'Avinash CK'
__credits__ = 'Indian Institute of Astrophysics'

from .simulation import Imager
from .simulation import Mosaic
from .analysis import Analyzer
from pathlib import Path

data_dir = Path(__file__).parent.joinpath('data')
