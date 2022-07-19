from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets

import pista as Sim
from pista.analysis import *
from pista.utils import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits
from astropy.table import Table
from scipy.constants import c
from matplotlib.ticker import LogLocator
from matplotlib import gridspec as gs
import matplotlib
from pathlib import Path

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

data_path = Path(Sim.__file__).parent.joinpath()
matplotlib.rcParams['font.size']=7
class HoverTracker(QObject):
    positionChanged = pyqtSignal(QPoint)
    clicked         = pyqtSignal(QPoint)

    def __init__(self, widget):
        super().__init__(widget)
        self._widget = widget
        self.widget.setMouseTracking(True)
        self.widget.installEventFilter(self)

    @property
    def widget(self):
        return self._widget

    def eventFilter(self, obj, event):
        if obj is self.widget and event.type() == QEvent.MouseMove:
            self.positionChanged.emit(event.pos())
        elif obj is self.widget and event.type()== QEvent.MouseButtonPress:
            self.clicked.emit(event.pos())
        return super().eventFilter(obj, event)
    
class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        self.df             = None
        self.psf_file       = None
        self.sim            = None
        self.count_sims = 1
        self.response_funcs = []
        
        super(Ui, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi(f'{data_path}/data/pista.ui', self) # Load the .ui file
        self.upload.clicked.connect(self.upload_source)
        self.set_default.clicked.connect(self.setdefault)
        self.upload_psf.clicked.connect(self.upload_psf_btn)
        self.upload_coating.clicked.connect(self.upload_coating_btn)
        self.upload_filters.clicked.connect(self.upload_filters_btn)
        self.simulate.clicked.connect(self.simulate_btn)
        self.upload_QE_wav.clicked.connect(self.upload_QE_wav_btn)
        self.generate.clicked.connect(self.generate_1d)
        self.download_image.clicked.connect(self.download_img)
    
 
    def on_position_changed(self, p):
        x = (p.x()-96)/(303-96)
        y = 1-(p.y()-37)/(248-37)
        
        x = int(x*self.n_pix_value)
        y = int(y*self.n_pix_value)
        self.x = x
        self.y = y
        
        ra,dec = self.sim.wcs.all_pix2world(x,y,1)
        self.ra = ra
        self.dec = dec
        
        if x>self.n_pix_value or x<0 or y>self.n_pix_value or y<0:
            ra = ' '
            dec = ' '   
        else:
            ra = np.round(ra,4)
            dec = np.round(dec,4)
            
        self.statusBar.showMessage(f'RA: {ra,p.x()}, Dec: {dec,p.y()}')
    
    def on_position_clicked(self,p):
        x = (p.x()-96)/(303-96)
        y = 1-(p.y()-37)/(248-37)
        x = x*self.n_pix_value
        y = y*self.n_pix_value
        self.x = x
        self.y = y
        self.x_c.setText(f'{int(x)}')
        self.y_c.setText(f'{self.n_pix_value-int(y)}')
        self.canvas_img.figure.clear()
        self.ax_img.cla()
        self.ax_img.set_title('')
        self.canvas_img.figure, self.ax_img = self.sim.show_image(self.output_select.currentText(),
                           fig = self.canvas_img.figure, ax = self.ax_img,
                           cmap = 'gray')
        w = int(self.n_pix_value*0.035)
        self.ax_img.annotate('+',(x-w,y-w), color = 'red')
        self.canvas_img.draw_idle()
        self.output_select.raise_() 
        self.generate.click()
        
    def check_input(self,var):
        flag = True
        for i in var:
            if not i.isnumeric() and i not in['.','e','-','+'] :
                flag = False
        if len(var)<1:
            flag = False
        return flag
            
    def check_df_params(self,df):
        flag = True
        
        for i in ['ra','dec','mag']:
            if i not in df.keys():
                flag = False
                break
        return flag
                
    def upload_source(self):
        filters = "FITS (*.fit*);;CSV (*.csv);;DAT (*.dat)"
        self.source_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.source.setText(self.source_file.split('/')[-1])
        if len(self.source.text())>4:
            ext = self.source_file.split('.')[-1]
            if ext =='fits' or ext =='fit':
                tab = Table.read(self.source_file)
                df  = tab.to_pandas()
                df  = df.dropna()
                if self.check_df_params(df):
                    self.df = df
                else:
                    self.error_box.setText('Error: Invalid input')
                    self.abmag.setChecked(True)                 
            elif ext =='csv':
                df = pd.read_csv(self.source_file)
                df = df.dropna()
                if self.check_df_params(df):
                    self.df = df
                else:
                    self.error_box.setText('Error: Invalid input')
                    self.abmag.setChecked(True)       
    
    def upload_psf_btn(self):
        filters = "FITS (*.fits);;NPY (*.npy)"
        psf_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.psf_filename.setText(psf_file.split('/')[-1])
        if os.path.exists(psf_file):
            self.psf_file = psf_file
        else :
            self.psf_file = None
        
    def upload_coating_btn(self):
        filters = "DAT (*.dat)"
        coating_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.coating_filename.setText(coating_file.split('/')[-1])
        self.n_mirrors = self.mirrors.value()
        self.coating_file = coating_file
      
            
    def upload_filters_btn(self):
        filters = "DAT (*.dat)"
        filter_files = QFileDialog.getOpenFileNames(filter = filters)[0]
        filter_names = [i.split('/')[-1] for i in filter_files]
        self.filter_filenames.setText(str(filter_names))
        
        self.filter_files = filter_files
        
    def upload_QE_wav_btn(self):
        filters = "DAT (*.dat)"
        self.QE_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.QE_filename.setText(self.QE_file.split('/')[-1])
  
        
    
    def setdefault(self):
        self.response_funcs = []
        self.abmag.setText('20')
        self.n_pix.setText("1000")
        psf_file = f'{data_path}/data/off_axis_hcipy.npy'
        self.psf_file = psf_file

        self.psf_filename.setText(psf_file.split('/')[-1])
        if self.df is None:
            ra = [0]
            dec= [0]
            mag_nuv = [10]
            self.df = pd.DataFrame(zip(ra,dec,mag_nuv), columns = ['ra','dec','mag'])
           
        sim = Analyzer(df = self.df, exp_time = 60)
        
        self.abmag.setText('10')
        self.check_abmag.setChecked(True)
        self.exp_time.setText('60')
        self.shot_noise_type.setCurrentText('Gaussian')
        self.sky_shot_noise_type.setCurrentText('Gaussian')
        
        self.qe_mean.setText(str(sim.params['qe']))
        self.qe_sigma.setText(str(sim.params['qe_sigma']))
        self.bias.setText(str(sim.params['bias']))
        self.gain.setText(str(sim.params['G1']))
        self.RN.setText(str(sim.params['RN']))
        self.NF.setText(str(sim.params['NF']))
        self.bit_res.setText(str(sim.params['bit_res']))
        self.FWC.setText(str(sim.params['FWC']))
     
        self.PRNU.setChecked(True)
        self.PRNU_frac.setText(str(sim.params['PRNU_frac']))
        
        self.DC.setChecked(True)
        self.DNFP.setChecked(True)
        self.det_T.setText(str(sim.params['T']))
        self.DFM.setText(str(sim.params['DFM']))
        self.pix_area.setText(str(sim.params['pixel_area']))
        self.DN.setText(str(sim.params['DN']))
        
        coating = f'{data_path}/data/UV/Coating.dat'
        self.coating_file = coating
        self.coating_filename.setText(coating.split('/')[-1])
        self.mirrors.setValue(5)
        self.n_mirrors = 5
        self.n_stack.setValue(1)
        self.stack_type.setCurrentText('Mean')
        filters = [f'{data_path}/data/UV/Filter.dat', 
                   f'{data_path}/data/UV/Dichroic.dat',
                   f'{data_path}/data/UV/Dichroic.dat']
        
        self.pixel_
        self.filter_files =  [f'{i},1' for i in filters]
        self.filter_filenames.setText(str([i.split('/')[-1] for i in filters]))
        self.pixel_scale.setText("0.1")
        self.QE_file = f'{data_path}/data/QE.dat'
        self.QE_filename.setText(self.QE_file.split('/')[-1])   
        self.cosmic_rays.setChecked(True)
        self.QE.setChecked(True)
        self.QN.setChecked(True)

    def check_params(self):    
        try:
            if not os.path.exists(self.source_file):
                self.source_file = None
        except:
            self.source_file = None
            
        if len(self.source.text())<1:
            self.check_abmag.setChecked(True)

        elif self.source_file is not None:
            ext = self.source_file.split('.')[-1]
            if ext =='fits' or ext =='fit':
                tab = Table.read(self.source_file)
                df  = tab.to_pandas()
                if self.check_df_params(df):
                    self.df = df
                else:
                    self.error_box.setText('Error: Invalid input')
                    self.abmag.setChecked(True)
                    
            elif ext =='csv':
                df = pd.read_csv(self.source_file)
                if self.check_df_params(df):
                    self.df = df
                else:
                    self.error_box.setText('Error: Invalid input')
                    self.abmag.setChecked(True)
    
            
        if self.check_abmag.isChecked():
            abmag = self.abmag.text()
            if len(abmag)>0 and self.check_input(abmag):
                abmag = float(abmag)
            else : 
                self.error_box.setText("Error: ABmag should be a real number."
                                        + " Default ABmag set")
                abmag = 10
                self.abmag.setText('10')
            ra = [0]
            dec= [0]
            mag = [abmag]
            self.df = pd.DataFrame(zip(ra,dec,mag), 
                                   columns = ['ra','dec','mag'])
            
        exp_time = self.exp_time.text()
        if len(exp_time)>0 and self.check_input(exp_time):
            self.exp_time_value = float(exp_time)
            if self.exp_time_value<0:
                self.error_box.setText("Error: Exposure should be a positive real number."
                                        + " Default Exposure time set")
                self.exp_time.value = 60         
        else :
            self.exp_time.setText('60')
            self.exp_time_value = 60
            self.error_box.setText("Error: Exposure should be a positive real number."
                                    + " Default Exposure time set")
        
        pixel_scale = self.pixel_scale.text()
        if len(pixel_scale)> 0 and self.check_input(pixel_scale):
            pixel_scale = abs(float(pixel_scale))
            self.ps_value = pixel_scale
        else:
            pixel_scale = 0.1
            self.pixel_scale.setText('0.1')
            self.ps_value = pixel_scale
            
        psf_filename = self.psf_filename.text()
        if self.psf_file is None or len(psf_filename)<1 :
            fwhm = self.fwhm.text()
            if len(fwhm)>0:
                if self.check_input(fwhm):
                    fwhm   = float(fwhm)
                else:
                    self.fwhm.setText('0.1')
                    fwhm = 0.1
                          
                fwhm_in = fwhm/self.ps_value
                sigma   = fwhm_in*gaussian_fwhm_to_sigma
                
                psf = generate_psf(1001,sigma)
                np.save(f'{data_path}/data/user_defined_psf.npy', psf)
                self.psf_file =f'{data_path}/data/user_defined_psf.npy'
                self.psf_filename.setText(self.psf_file.split('/')[-1])
            else:
                self.psf_file = f'{data_path}/data/off_axis_hcipy.npy'
                self.psf_filename.setText(self.psf_file.split('/')[-1])
           
        self.response_funcs= []
        # Coating
        if len(self.coating_filename.text())<1:
            coating = f'{data_path}/data/UV/Coating.dat'
            self.coating_file = coating
            self.coating_filename.setText(coating.split('/')[-1])
            self.mirrors.setValue(5)
            self.n_mirrors = 5
            self.response_funcs.append(f'{coating},5')
        else:
            if f'{self.coating_file},{self.n_mirrors}' not in self.response_funcs:
                if os.path.exists(self.coating_file):
                    self.response_funcs.append(f'{self.coating_file},{self.n_mirrors}')
            
        # Filters
        if len(self.filter_filenames.text() )<1:
            filters = [f'{data_path}/data/UV/Filter.dat', 
                       f'{data_path}/data/UV/Dichroic.dat',
                       f'{data_path}/data/UV/Dichroic.dat']
            self.filter_files = filters
            self.response_funcs+= [f'{i},1' for i in filters]
            self.filter_filenames.setText(str([i.split('/')[-1] for i in filters]))
        else:
            for filt in self.filter_files:
                if f'{filt},{1}' not in self.response_funcs:
                    if os.path.exists(filt):
                        self.response_funcs.append(f'{filt},{1}')
        # QE
        if self.QE.isChecked():
            if len(self.QE_filename.text())<1:
                self.QE_file = f'{data_path}/data/QE.dat'
                self.QE_filename.setText('QE.dat')
                self.response_funcs.append(f'{self.QE_file},1')
            else:
                if f'{self.QE_file},{1}' not in self.response_funcs:
                    if os.path.exists(self.QE_file):
                        self.response_funcs.append(f'{self.QE_file},1')
               
            qe_mean = self.qe_mean.text() 
            if len(qe_mean)<1:
                self.qe_mean.setText('0.5')
                self.params['qe_mean'] = 0.5
            elif self.check_input(qe_mean):
                qe_mean = float(qe_mean)     
                if qe_mean>0 and qe_mean <= 1:
                    self.params['qe_mean'] = qe_mean
                else:
                    self.qe_mean.setText('0.5')
                    self.params['qe_mean'] = 0.5
            
            qe_sigma = self.qe_sigma.text()
            if len(qe_sigma)<1:
                self.qe_sigma.setText('0.01')
                self.params['qe_sigma'] = 0.01
            elif self.check_input(qe_sigma):
                qe_sigma  = float(qe_sigma)  
                
                if qe_sigma <0.6*self.params['qe_mean']:
                    self.params['qe_sigma'] = qe_sigma
                else:
                    self.qe_sigma.setText('0.01')
                    self.params['qe_sigma'] = 0.01
                
        bias = self.bias.text()
        if len(bias)<1:
            self.bias.setText('35')
            self.params['bias'] = 35
        elif self.check_input(bias):
            bias  = float(bias)  
            
            if bias >= 0:
                self.params['bias'] = bias
            else:
                self.bias.setText('35')
                self.params['bias'] = 35
                
        gain = self.gain.text()
        if len(gain)<1:
            self.gain.setText('1')
            self.params['gain'] = 1
            
        elif self.check_input(gain):
            gain  = float(gain)    
            if gain >= 0:
                self.params['gain'] = gain
            else:
                self.gain.setText('1')
                self.params['gain'] = 1
                
        bit_res = self.bit_res.text()
        if len(bit_res)<1:
            self.bit_res.setText('14')
            self.params['bit_res'] = 14
            
        elif bit_res.isnumeric():
            bit_res  = int(bit_res)    
            if bit_res >= 0:
                self.params['bit_res'] = bit_res
            else:
                self.bit_res.setText('14')
                self.params['bit_res'] = 14
                
        RN = self.RN.text()
        if len(RN)<1:
            self.RN.setText('5')
            self.params['RN'] = 5
            
        elif self.check_input(RN):
            RN  = float(RN)    
            if RN >= 0:
                self.params['RN'] = RN
            else:
                self.RN.setText('5')
                self.params['RN'] = 5
                
        NF = self.NF.text()
        if len(NF)<1:
            self.NF.setText('0')
            self.params['NF'] = 0
            
        elif self.check_input(NF):
            NF  = float(NF)    
            if NF >= 0:
                self.params['NF'] = NF
            else:
                self.NF.setText('0')
                self.params['NF'] = 0
                
        FWC = self.FWC.text()
        if len(FWC)<1:
            self.FWC.setText('140000')
            self.params['FWC'] = 1.4e5
            
        elif self.check_input(FWC):
            FWC  = float(FWC)    
            if FWC >= 0:
                self.params['FWC'] = FWC
            else:
                self.FWC.setText('140000')
                self.params['FWC'] = 1.4e5
                
        n_pix = self.n_pix.text()
        if len(n_pix)<1:
            self.n_pix.setText('1000')
            self.n_pix_value = int(1000)
            
        elif n_pix.isnumeric():
            n_pix  = int (n_pix)    
            if n_pix > 1:
                self.n_pix_value  = n_pix
            else:
                self.n_pix.setText('1000')
                self.n_pix_value  = int(1000)
        
        if self.PRNU.isChecked():
            PRNU_frac = self.PRNU_frac.text()
            if len(PRNU_frac)<1:
                self.PRNU_frac.setText('0.0025')
                self.params['PRNU_frac'] = 0.25/100
                
            elif self.check_input(PRNU_frac):
                PRNU_frac  = float(PRNU_frac)    
                if PRNU_frac >= 0:
                    self.params['PRNU_frac'] = PRNU_frac
                else:
                    self.PRNU_frac.setText('0.0025')
                    self.params['PRNU_frac'] = 0.25/100
                    
        if self.DNFP.isChecked():
            self.DC.setChecked(True)
            DN = self.DN.text()
            if len(DN)<1:
                self.DN.setText('0.001')
                self.params['DN'] = 0.1/100
                
            elif self.check_input(DN):
                DN  = float(DN)    
                if DN >= 0:
                    self.params['DN'] = DN
                else:
                    self.DN.setText('0.001')
                    self.params['DN'] =0.1/100
        if self.DC.isChecked():
          
            DFM = self.DFM.text()
            if len(DFM)<1:
                self.DFM.setText('0.01424')
                self.params['DFM'] =0.01424
                
            elif self.check_input(DFM):
                DFM  = float(DFM)    
                if DFM >= 0:
                    self.params['DFM'] = DFM
                else:
                    self.DFM.setText('0.01424')
                    self.params['DFM'] = 0.01424
                    
            det_T = self.det_T.text()
            if len(det_T)<1:
                self.det_T.setText('223')
                self.params['det_T'] = 223
                
            elif self.check_input(det_T):
                det_T  = float(det_T)    
                if det_T >= 0:
                    self.params['det_T'] = det_T
                else:
                    self.det_T.setText('223')
                    self.params['det_T'] = 223
            
            pix_area = self.pix_area.text()
            if len(pix_area)<1:
                self.pix_area.setText('1e-6')
                self.params['pix_area'] = 1e-6
                
            elif self.check_input(pix_area):
                pix_area  = float(pix_area)    
                if pix_area >= 0:
                    self.params['pix_area'] = pix_area
                else:
                    self.pix_area.setText('1e-6')
                    self.params['pix_area'] = 1e-6
            
    def simulate_btn(self):
        
        self.params = {}
        self.error_box.setText('')
        self.check_params()
        self.statusBar.showMessage('Simulating Image')
        sim = Analyzer(df = self.df, cols = {'mag_nuv' :'mag'}, psf_file = self.psf_file
                       , exp_time = self.exp_time_value,n_pix= self.n_pix_value,
                       response_funcs = self.response_funcs
                       , pixel_scale = self.ps_value)
        n = int(self.n_stack.value())
        stack_type = self.stack_type.currentText()
        
        sim.cosmic_rays = self.cosmic_rays.isChecked()
        sim.QE      = self.QE.isChecked()
        sim.QN      = self.QN.isChecked()
        sim.PRNU    = self.PRNU.isChecked()
        sim.DC      = self.DC.isChecked()
        sim.DNFP    = self.DNFP.isChecked()
        if self.shot_noise_type.currentText()=='None':
            sim.shot_noise = False
        if self.sky_shot_noise_type.currentText()=='None':
            sim.sky = False

        sim(params = self.params,n_stack=n,stack_type=stack_type)
        self.sim = sim
        self.make_plots()
        
    def make_plots(self):
        
        plt.ioff()
        sim = self.sim
        # Main Image    
        if self.count_sims>1:
            self.canvas_img.figure.clear()
            self.ax_img.cla()
            self.ax_img.set_title('')
            self.sim.show_image(self.output_select.currentText(),
                               fig = self.canvas_img.figure, ax = self.ax_img,
                               cmap = 'gray')
            self.canvas_img.draw_idle()
            self.output_select.raise_() 
            
        else:    
            self.canvas_img = FigureCanvas()  
            self.ax_img = self.canvas_img.figure.add_subplot()
            
            self.canvas_img.figure.clear()   
            self.sim.show_image(self.output_select.currentText(),
                               fig = self.canvas_img.figure, ax = self.ax_img,
                               cmap = 'gray')

            
            self.image_box.addWidget(self.canvas_img)
            self.output_select.raise_() 
            
        self.image_hover = HoverTracker(self.canvas_img)
        self.image_hover.positionChanged.connect(self.on_position_changed)
        self.image_hover.clicked.connect(self.on_position_clicked)
            
        # Panel
        if self.count_sims>1:
            
            self.figure_panel.clear()
            
            self.wav_zp = np.linspace(1000,9000,10000)
            self.flux_zp = (c*1e2*3.631e-20)/(self.wav_zp**2*1e-8) 
            
            fig = self.figure_panel
            g = gs.GridSpec(2, 2,wspace = 0.35,hspace = 0.7, 
                            width_ratios=(1,0.35), top =0.91, right=0.91, 
                            bottom = 0.15)   
            axes = []
            
            ax = self.canvas_panel.figure.add_subplot(g[0,0])
            
            fig, ax, data, params = bandpass(self.wav_zp,self.flux_zp,
                                             self.response_funcs, fig = fig,
                                             ax = ax)
            ax.set_title('Bandpass | Zero Point', fontsize = 10)
            ax.set_xlabel(r'$\AA$',fontsize = 10)
            
            
            # PSF
            ax = self.canvas_panel.figure.add_subplot(g[0,1])
            ext = self.psf_file.split('.')[-1]
            if ext =='npy':
                self.psf_data = np.load(self.psf_file)
                ax.set_title('Point Spread Function', fontsize = 10)
            elif ext =='fits':
                hdu = fits.open(self.psf_file)[0]
                self.psf_data = hdu.data
                ax.set_title(f"Point Spread Function\n Pixel scale : {hdu.header['PIXELSCL']}"
                         , fontsize = 7)
            if self.psf_data.min()<0:
                self.psf_data-= self.psf_data.min() + 1e-3
            ax = self.canvas_panel.figure.add_subplot(g[0,1])
            
            ax.imshow(self.psf_data, cmap = 'jet',
                               extent = (-1,1,-1,1), norm = col.LogNorm())
            ax.grid(False)

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            
            # Bandpass Sky
            filt_dat  = np.loadtxt(f'{data_path}/data/Sky_mag.dat')
            self.wav_sky   = filt_dat[:,0]
            self.flux_sky  = filt_dat[:,1]
  
            ax = self.canvas_panel.figure.add_subplot(g[1,0])
            fig, ax, data, params = bandpass(self.wav_sky,self.flux_sky
                                             ,self.response_funcs, fig = fig,
                                             ax = ax)
            ax.set_title('Bandpass | Sky', fontsize = 10)
            ax.set_xlabel(r'$\AA$',fontsize = 10)
            
            # Photon Response
            
            ax = self.canvas_panel.figure.add_subplot(g[1,1])
            
            x = sim.source_photons.ravel()
            y = x/np.sqrt(sim.source_photons).ravel()
            x = np.where(x==0,np.nan,x)
            ax.plot(x,y,'.', color ='red')
            
            if self.DC.isChecked():
                x = (sim.photoelec_array-sim.DC_array).ravel()
                y = x/np.sqrt(sim.photoelec_array + sim.DC_array**2
                              + (sim.photoelec_array*sim.PRNU_array)**2 + 
                              (sim.DC_array*sim.DNFP_array)**2).ravel()
                
                x = np.where(x==0,np.nan,x)
                y = np.where(y==0,np.nan,y)
                
                ax.plot(x,y,'.', color ='purple')
            
            x =  np.linspace(1,1e11)
            y = x*0 +3
            
            ax.plot(x,y)
           
            ax.set_xlabel('Counts')
            ax.set_ylabel('SNR')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1,1e10)
            
            ax.xaxis.set_major_locator(LogLocator(numticks=5))
            ax.xaxis.set_minor_locator(LogLocator(numticks=5,subs=np.arange(2,10)))
            ax.set_ylim(1e-2,1e6)
            ax.grid(True, which="both",axis = 'both', ls="-", color='0.65')
            ax.set_title('Photon Response', fontsize = 10)

            self.figure_panel.canvas.draw_idle()
            
        else:
            g = gs.GridSpec(2, 2,wspace = 0.35,hspace = 0.7, 
                            width_ratios=(1,0.35), top =0.91, right=0.91, 
                            bottom = 0.15)                      
            self.canvas_panel = FigureCanvas()
            fig = self.canvas_panel.figure
            # Bandpass Zero point
            self.wav_zp = np.linspace(1000,25000,10000)
            self.flux_zp = (c*1e2*3.631e-20)/(self.wav_zp**2*1e-8) 
            axes = []
            
            ax = self.canvas_panel.figure.add_subplot(g[0,0])
            
            fig, ax, data, params = bandpass(self.wav_zp,self.flux_zp,
                                             self.response_funcs, fig = fig,
                                             ax = ax)
            ax.set_title('Bandpass | Zero Point', fontsize = 10)
            ax.set_xlabel(r'$\AA$',fontsize = 10)
        
            
            
            # PSF
            ax = self.canvas_panel.figure.add_subplot(g[0,1])
            ext = self.psf_file.split('.')[-1]
            if ext =='npy':
                self.psf_data = np.load(self.psf_file)
                ax.set_title('Point Spread Function', fontsize = 10)
            elif ext =='fits':
                hdu = fits.open(self.psf_file)[0]
                self.psf_data = hdu.data
                ax.set_title(f"Point Spread Function\n Pixel scale : {hdu.header['PIXELSCL']}"
                             , fontsize = 7)
            if self.psf_data.min()<=0:
                self.psf_data+= -self.psf_data.min() + 1e-2

        
            ax.imshow(self.psf_data, cmap = 'jet',
                               extent = (-1,1,-1,1), norm = col.LogNorm())
            ax.grid(False)
           
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            
            
            
            # Bandpass Sky
            filt_dat  = np.loadtxt(f'{data_path}/data/Sky_mag.dat')
            self.wav_sky   = filt_dat[:,0]
            self.flux_sky  = filt_dat[:,1]
  
            ax = self.canvas_panel.figure.add_subplot(g[1,0])
            fig, ax, data, params = bandpass(self.wav_sky,self.flux_sky
                                             ,self.response_funcs, fig = fig,
                                             ax = ax)
            ax.set_title('Bandpass | Sky', fontsize = 10)
            ax.set_xlabel(r'$\AA$',fontsize = 10)
            
            # Photon Response
            
            ax = self.canvas_panel.figure.add_subplot(g[1,1])
            
            x = sim.source_photons.ravel()
            y = x/np.sqrt(sim.source_photons).ravel()
            x = np.where(x==0,np.nan,x)
            ax.plot(x,y,'.', color ='red')
            if self.DC.isChecked():
                x = (sim.photoelec_array-sim.DC_array).ravel()
                y = x/np.sqrt(sim.photoelec_array + sim.DC_array**2
                          + (sim.photoelec_array*sim.PRNU_array)**2 + 
                          (sim.DC_array*sim.DNFP_array)**2).ravel()
                
                ax.set_ylim(1e-2,1e11)
                x = np.where(x==0,np.nan,x)
                y = np.where(y==0,np.nan,y)
                
               
            ax.plot(x,y,'.', color ='purple')
            
            x =  np.linspace(1,1e11)
            y = x*0 +3
            
            ax.plot(x,y)
           
            ax.set_xlabel('Counts')
            ax.set_ylabel('SNR')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1,1e10)
            
            ax.xaxis.set_major_locator(LogLocator(numticks=5))
            ax.xaxis.set_minor_locator(LogLocator(numticks=5,subs=np.arange(2,10)))
            
            ax.grid(True, which="both",axis = 'both', ls="-", color='0.65')
            ax.set_title('Photon Response', fontsize = 10)
            
            
            self.figure_panel = self.canvas_panel.figure
            
            self.panel_box.addWidget(self.canvas_panel)
            
       
        # Main Hist
        if self.count_sims>1:
            self.ax_hist.cla()
            self.ax_hist.set_title('')
            self.sim.show_hist(self.output_select.currentText(),
                                fig = self.canvas_hist.figure
                                , ax = self.ax_hist)
            self.canvas_hist.draw_idle()
            
        else:         
            self.canvas_hist = FigureCanvas()
            self.ax_hist     = self.canvas_hist.figure.add_subplot()
            self.sim.show_hist(self.output_select.currentText(),
                                fig = self.canvas_hist.figure
                                , ax = self.ax_hist)           
            
            self.output_hist.addWidget(self.canvas_hist)
            
        # 1D slice
        if self.count_sims>1:
            self.ax_1d.cla()
            self.ax_1d.set_title('')
            data = sim.getImage(self.output_select.currentText())
            
            x = self.x_c.text()
            if self.check_input(x):
                if int(x)< self.n_pix_value:
                    x = int(x)
                else:
                    x = self.n_pix_value//2
            else:
                 x = self.n_pix_value//2
            
            self.x_c.setText(f'{x}')
            
            y = self.y_c.text()
            if self.check_input(y):
                if int(y)< self.n_pix_value:
                    y = int(y)
                else:
                    y = self.n_pix_value//2
            else:
                 y = self.n_pix_value//2
            
            self.y_c.setText(f'{y}')
            
            w = self.n_pix_value//2
            if self.direction.currentText()=="Horizontal":
                
                self.ax_1d.set_xlabel('x')
                self.ax_1d.set_ylabel('Flux')
                
                data = data[:,y]
            elif self.direction.currentText()=="Vertical":
                
                 self.ax_1d.set_xlabel('y')
                 self.ax_1d.set_ylabel('Flux')
                 data = data[x,:]
                 
            self.ax_1d.plot(data)
            if self.scale.currentText()=='Log':
                self.ax_1d.set_yscale('log')
            self.canvas_1d.draw_idle()
            
        else:                
            self.canvas_1d = FigureCanvas()
            self.slice_1d.addWidget(self.canvas_1d)
            self.ax_1d = self.canvas_1d.figure.add_subplot()
            self.ax_1d.grid(False)
            data = sim.getImage(self.output_select.currentText())
            
            x = self.x_c.text()
            if self.check_input(x):
                if int(x)< self.n_pix_value:
                    x = int(x)
                else:
                    x = self.n_pix_value//2
            else:
                 x = self.n_pix_value//2
            
            self.x_c.setText(f'{x}')
            
            y = self.y_c.text()
            if self.check_input(y):
                if int(y)< self.n_pix_value:
                    y = int(y)
                else:
                    y = self.n_pix_value//2
            else:
                 y = self.n_pix_value//2
            
            self.y_c.setText(f'{y}')
            
            w = self.n_pix_value//10
            if self.direction.currentText()=="Horizontal":
                
                self.ax_1d.set_xlabel('x')
                self.ax_1d.set_ylabel('Flux')
                
                data = data[:,y]
            elif self.direction.currentText()=="Vertical":
                
                 self.ax_1d.set_xlabel('y')
                 self.ax_1d.set_ylabel('Flux')
                 data = data[x,:]
                 
            self.ax_1d.plot(data)
            if self.scale.currentText()=='Log':
                self.ax_1d.set_yscale('log')
                
        self.count_sims+=1
        
    def generate_1d(self):
        if self.count_sims>1:
            self.ax_1d.cla()
            self.ax_1d.set_title('')
            data = self.sim.getImage(self.output_select.currentText())
            
            x = self.x_c.text()
            if self.check_input(x):
                if int(x)< self.n_pix_value:
                    x = int(x)
                else:
                    x = self.n_pix_value//2
            else:
                 x = self.n_pix_value//2
            
            self.x_c.setText(f'{x}')
            
            y = self.y_c.text()
            if self.check_input(y):
                if int(y)< self.n_pix_value:
                    y = int(y)
                else:
                    y = self.n_pix_value//2
            else:
                 y = self.n_pix_value//2
            
            self.y_c.setText(f'{y}')
            
            w = self.n_pix_value//10
            if self.direction.currentText()=="Horizontal":
                
                self.ax_1d.set_xlabel('x')
                self.ax_1d.set_ylabel('Flux')
                
                data = data[:,y]
            elif self.direction.currentText()=="Vertical":
                
                 self.ax_1d.set_xlabel('y')
                 self.ax_1d.set_ylabel('Flux')
                 data = data[x,:]
                 
            self.ax_1d.plot(data)
            if self.scale.currentText()=='Log':
                self.ax_1d.set_yscale('log')
                
   
            self.canvas_1d.draw_idle()
        else:
            self.error_box.setText("Error: Image not generated. Click Simulate First!")
         
    def download_img(self):
        source = self.output_select.currentText()
        if self.sim is not None:
            self.sim.writeto(name = f'{source}.fits',source = source)   
            self.statusBar.showMessage(f'Downloading Image to : {data_path}')
            
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv) # Create an instance of QtWidgets.QApplication
    window = Ui() # Create an instance of our class
    window.show() # Show the GUI
    app.exec_()   # Start the application
    