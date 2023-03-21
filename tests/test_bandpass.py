"""This modules contains several test routines for
different modules in this python package"""
import numpy as np
from pista.utils import bandpass

data_path = 'src/pista/data'


def test_zeropoint():
    """Testing if zero point in AB system is 3631"""
    wav = np.arange(1000, 10000, 1)
    flux = 3631/(wav**2*3.34e4)
    inputs = [
                f'{data_path}/INSIST/UV/Coating.dat,5,100',
                f'{data_path}/INSIST/UV/Filter.dat,1,100',
                f'{data_path}/INSIST/UV/Dichroic.dat,2,100',
              ]

    _, _, data, params = bandpass(wav, flux,  inputs, False)

    lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

    assert round(int_flux_Jy, 0) == 3631
