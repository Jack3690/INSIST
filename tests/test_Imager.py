from pista import Imager
from pista import data_dir
import pandas as pd
import numpy as np
from pathlib import Path

data_path = Path(__file__).parent.joinpath()

ra = [0]
dec = [0]
mag = [18]

df = pd.DataFrame(zip(ra, dec, mag), columns=['ra', 'dec', 'mag'])

tel_params = {
            'aperture': 100,
            'pixel_scale': 0.1,
            'psf_file': f'{data_dir}/PSF/INSIST/off_axis_poppy.npy',
            'response_funcs': [f'{data_dir}/INSIST/UV/Filter.dat, 1, 100',
                               f'{data_dir}/INSIST/UV/Coating.dat, 5, 100',
                               f'{data_dir}/INSIST/UV/Dichroic.dat, 2, 100',
                               ],
            'coeffs': 1,
            'theta': 0
            }


def test_Imager():
    sim = Imager(df=df, tel_params=tel_params, n_x=200, n_y=200,
                 exp_time=1)
    sim()
    assert np.round(sim.phot_table['mag_out'][0], 0) == 18
