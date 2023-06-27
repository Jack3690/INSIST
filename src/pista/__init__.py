__version__ = "2.0.0"
__author__ = 'Avinash CK'
__credits__ = 'Indian Institute of Astrophysics'

from .simulation import Imager, Mosaic
from .analysis import Analyzer, SpecAnalyzer
from pathlib import Path

data_dir = Path(__file__).parent.joinpath('data')
