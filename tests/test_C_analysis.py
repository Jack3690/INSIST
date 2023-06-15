"""Tests for all Analysis Classes"""
from pista import Analyzer
from pista import Imager
from pista import data_dir
import os
import pytest
import pandas as pd


@pytest.fixture
def init_Imager():
    ra = [0, 1/3600, 2/3600]
    dec = [0, 0, 0]
    mag = [18, 19, 20]

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

    sim = Imager(df=df, tel_params=tel_params, n_x=200, n_y=500,
                 exp_time=600)
    sim()

    return sim


def test_Analyzer(init_Imager):
    sim = init_Imager
    an = Analyzer()

    df = sim.df
    wcs = sim.wcs
    data = sim.digital
    ZP = sim.ZP
    an(df=df, wcs=wcs, data=data, photometry='Aper', ZP=ZP)

    assert hasattr(an, 'phot_table')
    assert len(an.phot_table) == 3
    assert 'mag_in' in an.phot_table.keys()
    assert 'mag_out' in an.phot_table.keys()
    assert 'SNR' in an.phot_table.keys()


def test_plotting(init_Imager):
    sim = init_Imager
    an = Analyzer()

    df = sim.df
    wcs = sim.wcs
    data = sim.digital
    ZP = sim.ZP
    an(df=df, wcs=wcs, data=data, photometry='Aper', ZP=ZP)
    an.n_x = sim.n_x
    an.n_y = sim.n_y

    an.n_x_sim = sim.n_x_sim
    an.n_y_sim = sim.n_y_sim

    an.sim_df = sim.sim_df
    an.pixel_scale = sim.pixel_scale
    an.ra = sim.ra
    an.dec = sim.dec

    an.name = ''

    an.show_field()
    an.show_image()
    an.writeto('test.fits')
    if os.path.exists('test.fits'):
        os.remove('test.fits')
