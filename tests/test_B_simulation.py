"""Tests for all Simulation Classes"""
from pista import Imager
from pista import data_dir
import pytest
import pandas as pd


@pytest.fixture
def gen_imager_input():
    ra = [0]
    dec = [0]
    mag = [18]

    df = pd.DataFrame(zip(ra, dec, mag), columns=['ra', 'dec', 'mag'])
    response_funcs = [f'{data_dir}/INSIST/UV/Filter.dat, 1, 100',
                      f'{data_dir}/INSIST/UV/Coating.dat, 5, 100',
                      f'{data_dir}/INSIST/UV/Dichroic.dat, 2, 100',
                      ]
    tel_params = {
                'aperture': 100,
                'pixel_scale': 0.1,
                'psf_file': f'{data_dir}/PSF/INSIST/off_axis_poppy.npy',
                'response_funcs': response_funcs,
                'coeffs': 1,
                'theta': 0
                }
    return df, tel_params


def test_Imager_init(gen_imager_input):
    df, tel_params = gen_imager_input

    sim = Imager(df=df, tel_params=tel_params, n_x=200, n_y=500,
                 exp_time=1)

    assert hasattr(sim, 'zero_flux')
    assert hasattr(sim, 'sky_bag_flux')
    assert hasattr(sim, 'psf')
    assert hasattr(sim, 'sim_df')
    assert 'x' in sim.sim_df.keys() and 'y' in sim.sim_df.keys()


def test_Imager_call(gen_imager_input):
    df, tel_params = gen_imager_input

    sim = Imager(df=df, tel_params=tel_params, n_x=200, n_y=300,
                 exp_time=1)

    sim(photometry=None)

    assert hasattr(sim, 'wcs')
    assert hasattr(sim, 'digital')
    assert hasattr(sim, 'img_df')
    assert hasattr(sim, 'light_array')
    assert hasattr(sim, 'DC_array')
    assert sim.digital.shape == (300, 200)
