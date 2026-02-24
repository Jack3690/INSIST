"""This modules contains classes for simulating Spectroscopic modes"""
import torch
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import vstack
from tqdm.contrib import tzip
from .analysis import SpecAnalyzer, Analyzer
from .utils import bandpass, redshift_corr, select_mos, calc_mos_size
from .imaging import Imager
data_path = Path(__file__).parent.joinpath('data')

class Spectrograph(Imager, Analyzer, SpecAnalyzer):
    """Class for simulating Spectra"""
    def __init__(self, df, coords=None, tel_params={}, n_x=1000,
                 n_y=1000, exp_time=100, plot=False, user_profiles={},
                 **kwargs):

        if torch.cuda.is_available():
            self.device = 'cuda'
        else:
            self.device = 'cpu'
        self.tel_params = tel_params
        bin_min = self.tel_params['lambda1'] - self.tel_params['dellambda']/2
        bin_max = self.tel_params['lambda2'] + self.tel_params['dellambda']
        step = self.tel_params['dellambda']
        self.bins = np.arange(bin_min, bin_max, step)

        self.spec_bins = 0.5*(self.bins[1:] + self.bins[:-1])
        if len(self.spec_bins) % 2 == 0:
            raise Exception("""Spectroscopy bins need to be odd in number
                                Increase lambda2""")

        self.interpolate_mult = 3

        super().__init__(df, coords, tel_params, n_x, n_y, exp_time, plot,
                         user_profiles, **kwargs)
        self.sim_df = redshift_corr(self.sim_df)

        self.QE = False

    def check_df(self):

        if 'flux' not in self.df.keys() or 'wav' not in self.df.keys():
            raise Exception("'flux' or 'wav' column not found input dataframe")

        if 'ra' not in self.df.keys() or 'dec' not in self.df.keys():
            if 'x' in self.df.keys() and 'y' in self.df.keys():
                print("Converting xy to ra-dec")
                self.df = self.xy_to_radec(self.df, self.n_x, self.n_y,
                                           self.pixel_scale)
            else:
                raise Exception("'ra','dec','x',or 'y', \
                 columns not found in input dataframe ")

    def generate_sim_field(self, plot):
        if self.df is not None:
            self.calc_zp(plot=plot)
            self.init_psf_patch()
            spec_corr = len(self.spec_bins)//2

            x_left = spec_corr + self.n_pix_psf
            x_right = self.n_x_sim - spec_corr - self.n_pix_psf
            y_left = self.n_pix_psf
            y_right = self.n_y_sim - self.n_pix_psf - 1

            self.sim_df = self.init_df(df=self.df,
                                       n_x=self.n_x_sim, n_y=self.n_y_sim,
                                       x_left=x_left, x_right=x_right,
                                       y_left=y_left, y_right=y_right)
        else:
            print("df cannot be None")

    def calc_zp(self, plot=False):
        if len(self.response_funcs) > 0:

            wav = np.linspace(1200, 3200, 20000)
            flux = 3631/(3.34e4*wav**2)   # AB flux

            fig, ax, data, params = bandpass(wav, flux, self.response_funcs,
                                             plot=plot)

            lambda_, _, Reff = data
            lambda_phot, int_flux, int_flux_Jy, W_eff, flux_ratio = params

            self.lambda_phot = lambda_phot
            self.int_flux = int_flux
            self.W_eff = W_eff
            self.int_flux_Jy = int_flux_Jy
            self.flux_ratio = flux_ratio

            self.Reff = self.bin_xy_norm(lambda_, Reff)

            self.flux_coeffs = self.Reff.copy()
            self.flux_coeffs *= self.exp_time*self.coeffs*self.tel_area

            filt_dat = np.loadtxt(self.tel_params['sky_resp'])

            wav = filt_dat[:, 0]
            flux = filt_dat[:, 1]

            _, _, _, params = bandpass(wav, flux, self.response_funcs,
                                       plot=False)

            int_flux = params[1]
            self.det_params['M_sky'] = int_flux
        else:

            print("Response functions not provided. Using default values")
            self.int_flux_Jy = 3631
            self.W_eff = 1000
            self.lambda_phot = 2250
            self.flux_ratio = 1

        self.photons = 1.51e3*self.int_flux_Jy*(self.W_eff/self.lambda_phot)
        self.photons *= self.flux_ratio

        self.zero_flux = self.exp_time*self.tel_area*self.photons
        self.zero_flux *= self.coeffs
        self.M_sky_p = self.det_params['M_sky'] \
            - 2.5*np.log10(self.pixel_scale**2)

        self.sky_bag_flux = self.zero_flux*pow(10, -0.4*self.M_sky_p)

        # Unit Conversion
        self.UC = 1.51e3*(self.tel_params['dellambda']/self.spec_bins)
        self.UC *= self.flux_coeffs

        if self.sky:
            if self.user_profiles['sky'] is not None:
                if self.user_profiles['sky'].shape == (self.n_x, self.n_y):
                    self.sky_photons = self.user_profiles['sky']
                else:
                    raise Exception(f"""User defined sky array shape: \
                    {self.user_profiles['sky'].shape} \
                    is not same as detector shape {(self.n_x, self.n_y)}""")
            else:
                self.sky_photons = self.compute_shot_noise(self.sky_bag_flux)
        else:
            self.sky_photons = 0

    def init_MOS_df(self, L, B, PA, lw, df):
        delta = np.arctan(B/L)
        d = np.sqrt(L**2 + B**2)/2
        om = PA - delta

        x_corr = self.n_x_sim//2 - d*np.cos(om) - lw/2
        y_corr = self.n_y_sim//2 - d*np.sin(om) - B*np.cos(PA)

        t = df.copy()
        x0 = B*np.sin(PA) + lw/2 + x_corr
        y0 = y_corr
        t = t[t['y'] > np.tan(PA)*(t['x'] - x0) + y0]

        x0 = lw/2 + x_corr
        y0 = B*np.cos(PA) + y_corr
        t = t[t['y'] < np.tan(PA)*(t['x'] - x0) + y0]

        x0 = lw/2 + x_corr
        y0 = B*np.cos(PA) + y_corr
        t = t[t['x'] > -np.tan(PA)*(t['y'] - y0) + x0]

        x0 = L*np.cos(PA) + lw/2 + x_corr
        y0 = B*np.cos(PA) + L*np.sin(PA) + y_corr
        t = t[t['x'] < -np.tan(PA)*(t['y'] - y0) + x0]

        return t

    def init_MOS(self, L, B, PA):
        PA *= np.pi/180
        lw = len(self.spec_bins)
        x_size, y_size = calc_mos_size(L, B, PA, lw)
        if x_size < self.n_x and y_size < self.n_y:
            self.L = L
            self.B = B
            self.PA = PA

            self.sim_df = self.init_MOS_df(L, B, PA, lw, self.sim_df)
            self.sim_df['objid'] = np.arange(0, len(self.sim_df), 1)
        else:
            print(f"""n_x should be greater than {x_size} \n
                    n_y should be greater than {y_size}""")

    def select_MOS_sources(self, df=None, ids=[], radius=5, min_sep=8):
        if df is None:
            df = self.sim_df

        x = df['x'].value
        y = df['y'].value

        cen_x = x[np.argmin(abs(x-np.mean(x)))]
        cen_y = y[np.argmin(abs(y-np.mean(y)))]

        temp = df[['objid', 'x', 'y']].to_pandas()
        mos_df, res_df, cen_x, cen_y = select_mos(temp,
                                                  cen_x, cen_y, ids, radius,
                                                  min_sep, self.n_y_sim)
        t = []
        for id in mos_df['objid']:
            t.append(df[df['objid'] == id])

        mos_df = vstack(t)
        self.sim_df = mos_df

        t = []
        for id in res_df['objid']:
            t.append(df[df['objid'] == id])

        res_df = vstack(t)

        return mos_df, res_df, cen_x, cen_y

    def init_psf_patch(self, return_psf=False, plot=False):

        ext = self.psf_file.split('.')[-1]

        if ext == 'npy':
            image = np.load(self.psf_file)
        elif ext == 'fits':
            image = fits.open(self.psf_file)[0].data

        image /= image.sum()  # Flux normalized to 1

        self.psf = image

        self.n_pix_psf = self.psf.shape[0]

        spec_corr = len(self.spec_bins) + 1
        self.n_x_sim = self.n_x + spec_corr + 2*(self.n_pix_psf-1)
        self.n_y_sim = self.n_y + 2*(self.n_pix_psf-1)

        if return_psf:
            return image*self.zero_flux

    def init_image_array(self, return_img=False):
        """
        Creates a base image array for adding photons

        Parameters
        ----------
        return_img : bool, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        numpy.ndarray
            if return_img is true return base image array

        """
        image = np.zeros((self.n_y_sim, self.n_x_sim))
        self.mask = np.zeros((self.n_y_sim, self.n_x_sim))
        wcs = self.create_wcs(self.n_x_sim, self.n_y_sim,
                              self.ra, self.dec, self.pixel_scale,
                              self.theta)

        return image, wcs

    def bin_xy(self, x, y):
        step = self.tel_params['dellambda']/self.interpolate_mult
        x_new = np.arange(x.min(), x.max(), step)
        y_new = np.interp(x_new, x, y)
        x, y = x_new, y_new

        x_new = self.bins

        bin_indices = np.digitize(x, x_new)
        binned_y = []
        for i in range(1, len(x_new)):
            if i in bin_indices:
                binned_y.append(y[bin_indices == i].mean())
            else:
                raise Exception("""Interpolation Error.
                                Increase interpolate_mult""")

        return np.array(binned_y)

    def bin_xy_norm(self, x, y):
        step = self.tel_params['dellambda']/self.interpolate_mult
        x_new = np.arange(x.min(), x.max(), step)
        y_new = np.interp(x_new, x, y)
        x, y = x_new, y_new

        x_new = self.bins

        bin_indices = np.digitize(x, x_new)
        binned_y = []
        for i in range(1, len(x_new)):
            if i in bin_indices:
                binned_y.append(np.median(y[bin_indices == i]))
            else:
                raise Exception("""Interpolation Error. \
                Increase interpolate_mult""")

        return np.array(binned_y)

    def generate_photons(self, image, patch_width, df):
        """
        This function creates sims based on ABmag  on a
        small patch (2D array) of size n_pix_s*n_pix_s.

        The patch with the sim is then added to the image array of size
        n_pix_m*n_pix_m using wcs object.

        Parameters
        ----------
        image       : numpy.ndarray
                    base image array for inserting star sims
        patch_width : int
                    number of pixels (length) in sim patch image
        df          : pandas.dataframe
                    Dataframe containing source list


        Returns
        -------
        image : numpy.ndarray
            Array with sims added based on df

        """
        patch_width_mid = patch_width//2
        disp_width = len(self.spec_bins)
        disp_width_mid = disp_width//2

        patch = self.psf
        device = self.device

        # Convert to tensor once, before the loop
        image = torch.tensor(image, dtype=torch.float32, device=device)
        dl = self.tel_params['dellambda']
        lw = self.spec_bins

        x0 = df['x'].value.astype(int)
        y0 = df['y'].value.astype(int)

        x1 = x0 - patch_width_mid - disp_width_mid - 1
        x2 = x1 + patch_width + disp_width - 1
        y1 = y0 - patch_width_mid
        y2 = y1 + patch_width

        fluxes = []

        for row in df:
            x = row['wav']
            y = row['flux']
            fluxes.append(self.bin_xy(x, y))

        for x1_, x2_, y1_, y2_, flux in tzip(x1, x2, y1, y2, fluxes):

            flux *= lw**2*3.34e4
            flux *= 1.51e3*self.flux_coeffs*(dl/lw)
            flux = torch.tensor(flux, dtype=torch.float32, device=device)

            # Generating star using sim

            data = torch.tensor(flux, dtype=torch.float32, device=device)
            data = data.view(1, 1, 1, -1)
            kernel = torch.tensor(patch, dtype=torch.float32, device=device)
            kernel = kernel.view(1, 1, patch_width, patch_width)
            spectrum = F.conv2d(data, kernel, padding=patch_width-1)

            image[y1_:y2_, x1_:x2_] += spectrum.squeeze()
            y0 = y1_ + patch_width_mid
            self.mask[y0 - 3:y0 + 4, x1_:x2_] += 0.5

            # Empty GPU cache after each iteration
            torch.cuda.empty_cache()

        # Convert to numpy array once, outside of the loop
        image = image.detach().cpu().numpy()
        image = image[patch_width-1:-patch_width+1,  # Cropping Image
                      patch_width+disp_width_mid-1:
                      -disp_width_mid-patch_width-1]
        return image

    def __call__(self, det_params=None, n_stack=1, stack_type='median',
                 spectroscopy=True, fwhm=None, detect_sources=False,
                 spec_width=50, sky_sep=5):

        zero_p_flux = self.zero_flux
        zero_p_flux *= self.gain
        ZP = 2.5*np.log10(zero_p_flux)
        self.ZP = ZP
        self.UC *= self.gain
        self.UC *= 3.34e4*self.spec_bins**2
        super().__call__(det_params=det_params, n_stack=n_stack,
                         photometry=None, stack_type=stack_type, ZP=1)

        if spectroscopy:
            self.extract(data=self.digital.copy(), wcs=self.wcs,
                         df=self.img_df, width=spec_width, sep=sky_sep)
