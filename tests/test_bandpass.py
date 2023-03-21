"""This modules contains several test routines for
different modules in this python package"""

from src.pista.utils import bandpass
import numpy as np


def test_zeropoint():
    wav = np.arange(1000, 10000, 1)
    flux = 3631/(wav**2*3.34e4)
    data_path = 'src/pista/data'
    inputs = [
                f'{data_path}/INSIST/UV/Coating.dat,5,100',
                f'{data_path}/INSIST/UV/Filter.dat,1,100',
                f'{data_path}/INSIST/UV/Dichroic.dat,2,100',
              ]

    _, _, data, params = bandpass(wav, flux,  inputs, False)

    lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

    assert round(int_flux_Jy, 0) == 3631
