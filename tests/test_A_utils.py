"""This modules contains several test routines for
different modules in this python package"""
import numpy as np
import pandas as pd
import pytest
from pista.utils import bandpass, generate_psf, Xmatch
from pathlib import Path

# Zero Magnitude star flux in AB system
wav = np.arange(1000, 10000, 1)
flux = 3631/(wav**2*3.34e4)

data_dir = Path(__file__).parent.joinpath('data')
print(data_dir)


@pytest.fixture
def hst_bandpass():

    inputs = [
                f'{data_dir}/HST_WFC3_UVIS1.F275W.dat,1,1',
              ]

    _, _, data, params = bandpass(wav, flux,  inputs, False)

    lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

    tel_area = np.pi*(240/2)**2
    photons = 1.51e3*int_flux_Jy*(W_eff/lambda_phot)*tel_area*flux_ratio

    zp = 2.5*np.log10(photons)

    return zp


@pytest.fixture
def uvit_bandpass():
    inputs = [
                f'{data_dir}/Astrosat_UVIT.F148W.dat,1,1',
              ]

    _, _, data, params = bandpass(wav, flux,  inputs, False)

    lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

    photons = 1.51e3*int_flux_Jy*(W_eff/lambda_phot)*flux_ratio

    zp = 2.5*np.log10(photons)

    return zp


def test_bandpass(hst_bandpass, uvit_bandpass):
    """Testing zero points in AB system"""

    # HST
    zero_point = hst_bandpass

    assert np.round(zero_point, 2) == 24.16

    # UVIT
    zero_point = uvit_bandpass

    assert np.round(zero_point, 2) == 18.06


def test_generate_psf():
    params = {'sigma_x': 0.8,
              'sigma_y': 0.8,
              'gamma': 1,
              'alpha': 1}

    psf = generate_psf(101, params, 'Gaussian')
    assert len(psf) == 101
    psf = generate_psf(101, params, 'Moffat')
    assert len(psf) == 101


def test_Xmatch():
    ra = [10/3600, 15/3600]
    dec = [10/3600, 10/3600]
    df1 = pd.DataFrame(zip(ra, dec), columns=['ra', 'dec'])

    ra = [15/3600, 10/3600]
    dec = [10/3600, 10/3600]
    df2 = pd.DataFrame(zip(ra, dec), columns=['ra', 'dec'])

    matched = Xmatch(df1, df2, r=5)

    assert matched[0][1] == 1 and matched[1][0] == 1 and matched[0][2] == 0
