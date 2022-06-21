from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets

import pista as Sim
from pista.analysis import Analyzer
from pista.utils import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.stats import gaussian_fwhm_to_sigma

from pathlib import Path

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

data_path = Path(Sim.__file__).parent.joinpath()
class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        self.df             = None
        self.psf_file       = None
        self.response_funcs = []
        self.count_sims = 1
        
        super(Ui, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi(f'{data_path}/data/pista.ui', self) # Load the .ui file
        self.upload.clicked.connect(self.upload_source)
        self.set_default.clicked.connect(self.setdefault)
        self.upload_psf.clicked.connect(self.upload_psf_btn)
        self.upload_coating.clicked.connect(self.upload_coating_btn)
        self.upload_filters.clicked.connect(self.upload_filters_btn)
        self.simulate.clicked.connect(self.simulate_btn)
        
    def check_input(self,var):
        flag = True
        for i in var:
            if not i.isnumeric() and i not in['.','e','-','+'] :
                flag = False
        return flag
            
        
    def upload_source(self):
        filters = "CSV (*.csv);;DAT (*.dat)"
        source_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.source.setText(source_file.split('/')[-1])
        if len(self.source.text())>=1:
            self.df = pd.read_csv(source_file)
        
        
    def upload_psf_btn(self):
        filters = "FITS (*.fits);;NPY (*.npy)"
        psf_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.psf_filename.setText(psf_file.split('/')[-1])
        self.psf_file = psf_file
        
    def upload_coating_btn(self):
        filters = "DAT (*.dat)"
        coating_file = QFileDialog.getOpenFileName(filter = filters)[0]
        self.coating_filename.setText(coating_file.split('/')[-1])
        n = self.mirrors.value()
        if f'{coating_file},{n}' not in self.response_funcs:
            self.response_funcs.append(f'{coating_file},{n}')
            
    def upload_filters_btn(self):
        filters = "DAT (*.dat)"
        filter_files = QFileDialog.getOpenFileNames(filter = filters)[0]
        filter_names = [i.split('/')[-1] for i in filter_files]
        self.filter_filenames.setText(str(filter_names))
        
        for filt in filter_files:
            if f'{filt},{1}' not in self.response_funcs:
                self.response_funcs.append(f'{filt},{1}')
    
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
        self.coating_filename.setText(coating.split('/')[-1])
        self.response_funcs.append(f'{coating},5')
        self.mirrors.setValue(5)
        filters = [f'{data_path}/data/UV/Filter.dat', 
                   f'{data_path}/data/UV/Dichroic.dat',
                   f'{data_path}/data/UV/Dichroic.dat']
        
        self.response_funcs+= [f'{i},1' for i in filters]
        self.filter_filenames.setText(str([i.split('/')[-1] for i in filters]))
        self.pixel_scale.setText("0.1")
        qe_filename = f'{data_path}/data/UV/QE.dat'
        self.QE_filename.setText(qe_filename.split('/')[-1])        
        
    def check_params(self):
        
        if len(self.source.text())<1:
            self.check_abmag.setChecked(True)  
            
        if self.check_abmag.isChecked():
            abmag = self.abmag.text()
            if len(abmag)>=1 and self.check_input(abmag):
                abmag = float(abmag)
            else : 
                abmag = 10
                self.abmag.setText('10')
            ra = [0]
            dec= [0]
            mag = [abmag]
            self.df = pd.DataFrame(zip(ra,dec,mag), 
                                   columns = ['ra','dec','mag'])
            
        exp_time = self.exp_time.text()
        if len(exp_time)>=1 and self.check_input(exp_time):
            self.exp_time_value = float(exp_time)
        else :
            self.exp_time.setText('60')
            self.exp_time_value = 60
             
        if self.psf_file is None:
            fwhm = self.fwhm.text()
            if len(fwhm)>=1:
                if self.check_input(fwhm):
                    fwhm   = float(fwhm)
                else:
                    self.fwhm.setText('0.1')
                    fwhm = 0.1
                    
                pixel_scale = self.pixel_scale.text()
                if len(pixel_scale)> 1 and self.check_input(pixel_scale):
                    pixel_scale = float(pixel_scale)
                else:
                    pixel_scale = 0.1
                    self.pixel_scale.setText('0.1')
                
                          
                fwhm_in = fwhm/pixel_scale
                sigma   = fwhm_in*gaussian_fwhm_to_sigma
                
                psf = generate_psf(1001,sigma)
                np.save(f'{data_path}/data/user_defined_psf.npy', psf)
                self.psf_file =f'{data_path}/data/user_defined_psf.npy'
                self.psf_filename.setText(self.psf_file.split('/')[-1])
            else:
                self.psf_file = f'{data_path}/data/off_axis_hcipy.npy'
                self.psf_filename.setText(self.psf_file.split('/')[-1])
                
        if len(self.coating_filename.text() )<1:
            coating = f'{data_path}/data/UV/Coating.dat'
            self.coating_filename.setText(coating.split('/')[-1])
            self.mirrors.setValue(5)
            self.response_funcs.append(f'{coating},5')
            
        if len(self.filter_filenames.text() )<1:
            filters = [f'{data_path}/data/UV/Filter.dat', 
                       f'{data_path}/data/UV/Dichroic.dat',
                       f'{data_path}/data/UV/Dichroic.dat']
            
            self.response_funcs+= [f'{i},1' for i in filters]
            self.filter_filenames.setText(str([i.split('/')[-1] for i in filters]))
           
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
            if n_pix >= 0:
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
        self.check_params()
        sim = Analyzer(df = self.df, cols = {'mag_nuv' :'mag'}, psf_file = self.psf_file
                       , exp_time = self.exp_time_value,n_pix= self.n_pix_value,
                       response_funcs = self.response_funcs)
        sim(params = self.params)
        self.sim = sim
        
        try:
            fig_img, ax_img = self.sim.show_image(self.output_select.currentText())
        except:
            #self.display_out.setText(f"Error: {self.image_select.currentText()} data not generated")
            self.image_select.setCurrentText('Readout')
            fig_img, ax_img = self.psf.show_image(self.image_select.currentText())
    
        if self.count_sims>=1:
            self.image_box.removeWidget(self.canvas_img)
            self.canvas_img.deleteLater()      
            self.figure_img = fig_img
            self.canvas_img = FigureCanvas(self.figure_img)
            self.image_box.addWidget(self.canvas_img)
            self.output_select.raise_()
        else :
            self.figure_img = fig_img
            self.canvas_img = FigureCanvas(self.figure_img)
            self.image_box  = QtWidgets.QHBoxLayout(self.output_image)
            self.image_box.setObjectName('image_box')
            self.image_box.addWidget(self.canvas_img)
            self.output_select.raise_()
            
    def plot_image(self,fig,ax,out_var):
    # Image
       try:
           fig_img, ax_img = self.sim.show_image(self.output_select.currentText())
       except:
           #self.display_out.setText(f"Error: {self.image_select.currentText()} data not generated")
           self.image_select.setCurrentText('Readout')
           fig_img, ax_img = self.psf.show_image(self.image_select.currentText())

       if self.count_sims>=1:
           self.image_box.removeWidget(self.canvas_img)
           self.canvas_img.deleteLater()      
           self.figure_img = fig_img
           self.canvas_img = FigureCanvas(self.figure_img)
           self.image_box.addWidget(self.canvas_img)
           self.image_select.raise_()
       else :
           self.figure_img = fig_img
           self.canvas_img = FigureCanvas(self.figure_img)
           self.image_box  = QtWidgets.QHBoxLayout(self.output_image)
           self.image_box.setObjectName('image_box')
           self.image_box.addWidget(self.canvas_img)
           self.image_select.raise_()
                
  

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv) # Create an instance of QtWidgets.QApplication
    window = Ui() # Create an instance of our class
    window.show() # Show the GUI
    app.exec_()   # Start the application