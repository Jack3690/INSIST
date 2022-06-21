from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets
from PSF_Sim.image import PSF
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class Ui_Dialog(object):
    def __init__(self):
        self.df = None
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(823, 635)
        self.label_3 = QLabel(Dialog)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setGeometry(QRect(20, 479, 47, 13))
        self.label_4 = QLabel(Dialog)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setGeometry(QRect(20, 503, 47, 13))
        self.mode = QComboBox(Dialog)
        self.mode.addItem("")
        self.mode.addItem("")
        self.mode.setObjectName(u"mode")
        self.mode.setGeometry(QRect(90, 473, 111, 22))
        self.axis = QComboBox(Dialog)
        self.axis.addItem("")
        self.axis.addItem("")
        self.axis.setObjectName(u"axis")
        self.axis.setGeometry(QRect(90, 500, 111, 22))
        self.label_5 = QLabel(Dialog)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setGeometry(QRect(629, 17, 101, 21))
        self.label_5.setStyleSheet(u"font: 25 12pt \"Bookman Old Style\";")
        self.qe_mean = QLineEdit(Dialog)
        self.qe_mean.setObjectName(u"qe_mean")
        self.qe_mean.setGeometry(QRect(669, 53, 133, 16))
        self.label_6 = QLabel(Dialog)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setGeometry(QRect(521, 106, 121, 16))
        self.bias = QLineEdit(Dialog)
        self.bias.setObjectName(u"bias")
        self.bias.setGeometry(QRect(669, 100, 133, 20))
        self.label_8 = QLabel(Dialog)
        self.label_8.setObjectName(u"label_8")
        self.label_8.setGeometry(QRect(521, 55, 93, 16))
        self.bit_res = QLineEdit(Dialog)
        self.bit_res.setObjectName(u"bit_res")
        self.bit_res.setGeometry(QRect(669, 143, 133, 16))
        self.label_7 = QLabel(Dialog)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setGeometry(QRect(521, 145, 62, 16))
        self.label_9 = QLabel(Dialog)
        self.label_9.setObjectName(u"label_9")
        self.label_9.setGeometry(QRect(521, 126, 21, 16))
        self.gain = QLineEdit(Dialog)
        self.gain.setObjectName(u"gain")
        self.gain.setGeometry(QRect(669, 124, 133, 16))
        self.label_10 = QLabel(Dialog)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setGeometry(QRect(521, 164, 121, 16))
        self.RN = QLineEdit(Dialog)
        self.RN.setObjectName(u"RN")
        self.RN.setGeometry(QRect(669, 163, 133, 16))
        self.PRNU_cb = QCheckBox(Dialog)
        self.PRNU_cb.setObjectName(u"PRNU_cb")
        self.PRNU_cb.setGeometry(QRect(520, 185, 50, 17))
        self.label_11 = QLabel(Dialog)
        self.label_11.setObjectName(u"label_11")
        self.label_11.setGeometry(QRect(520, 205, 69, 16))
        self.PRNU_frac = QLineEdit(Dialog)
        self.PRNU_frac.setObjectName(u"PRNU_frac")
        self.PRNU_frac.setGeometry(QRect(669, 203, 133, 16))
        self.DC_cb = QCheckBox(Dialog)
        self.DC_cb.setObjectName(u"DC_cb")
        self.DC_cb.setGeometry(QRect(520, 226, 85, 17))
        self.DC_cb.setChecked(True)
        self.label_12 = QLabel(Dialog)
        self.label_12.setObjectName(u"label_12")
        self.label_12.setGeometry(QRect(520, 247, 85, 16))
        self.DFM = QLineEdit(Dialog)
        self.DFM.setObjectName(u"DFM")
        self.DFM.setGeometry(QRect(669, 242, 133, 17))
        self.DNFP_cb = QCheckBox(Dialog)
        self.DNFP_cb.setObjectName(u"DNFP_cb")
        self.DNFP_cb.setGeometry(QRect(520, 285, 142, 17))
        self.DN = QLineEdit(Dialog)
        self.DN.setObjectName(u"DN")
        self.DN.setGeometry(QRect(669, 305, 133, 17))
        self.label_13 = QLabel(Dialog)
        self.label_13.setObjectName(u"label_13")
        self.label_13.setGeometry(QRect(520, 307, 16, 16))
        self.Temp = QLineEdit(Dialog)
        self.Temp.setObjectName(u"Temp")
        self.Temp.setGeometry(QRect(669, 264, 133, 17))
        self.label_14 = QLabel(Dialog)
        self.label_14.setObjectName(u"label_14")
        self.label_14.setGeometry(QRect(520, 265, 81, 16))
        self.label_15 = QLabel(Dialog)
        self.label_15.setObjectName(u"label_15")
        self.label_15.setGeometry(QRect(310, 400, 101, 21))
        self.label_15.setStyleSheet(u"font: 25 12pt \"Bookman Old Style\";")
        self.sky_cb = QCheckBox(Dialog)
        self.sky_cb.setObjectName(u"sky_cb")
        self.sky_cb.setGeometry(QRect(230, 424, 85, 17))
        self.sky_cb.setChecked(True)
        self.sky_mag = QLineEdit(Dialog)
        self.sky_mag.setObjectName(u"sky_mag")
        self.sky_mag.setGeometry(QRect(379, 446, 133, 17))
        self.label_16 = QLabel(Dialog)
        self.label_16.setObjectName(u"label_16")
        self.label_16.setGeometry(QRect(230, 450, 85, 16))
        self.label_17 = QLabel(Dialog)
        self.label_17.setObjectName(u"label_17")
        self.label_17.setGeometry(QRect(80, 400, 101, 21))
        self.label_17.setStyleSheet(u"font: 25 12pt \"Bookman Old Style\";")
        self.label_17.setTextFormat(Qt.AutoText)
        self.shot_noise_cb = QCheckBox(Dialog)
        self.shot_noise_cb.setObjectName(u"shot_noise_cb")
        self.shot_noise_cb.setGeometry(QRect(230, 487, 85, 17))
        self.shot_noise_cb.setChecked(True)
        self.shot_noise_type = QComboBox(Dialog)
        self.shot_noise_type.addItem("")
        self.shot_noise_type.addItem("")
        self.shot_noise_type.setObjectName(u"shot_noise_type")
        self.shot_noise_type.setGeometry(QRect(379, 504, 133, 17))
        self.label_19 = QLabel(Dialog)
        self.label_19.setObjectName(u"label_19")
        self.label_19.setGeometry(QRect(231, 510, 47, 13))
        self.label_18 = QLabel(Dialog)
        self.label_18.setObjectName(u"label_18")
        self.label_18.setGeometry(QRect(20, 584, 85, 16))
        self.exp_time = QLineEdit(Dialog)
        self.exp_time.setObjectName(u"exp_time")
        self.exp_time.setGeometry(QRect(120, 580, 81, 20))
        self.output_image = QFrame(Dialog)
        self.output_image.setObjectName(u"output_image")
        self.output_image.setGeometry(QRect(20, 20, 491, 361))
        self.output_image.setFrameShape(QFrame.StyledPanel)
        #self.output_image.setFrameShadow(QFrame.Raised)
        self.image_select = QComboBox(self.output_image)
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.addItem("")
        self.image_select.setObjectName(u"image_select")
        self.image_select.setGeometry(QRect(0, 0, 69, 21))
        self.upload = QPushButton(Dialog)
        self.upload.setObjectName(u"upload")
        self.upload.setGeometry(QRect(60, 449, 111, 23))
        self.upload.clicked.connect(self.selectFile)
        self.lineEdit = QLineEdit(Dialog)
        self.lineEdit.setObjectName(u"lineEdit")
        self.lineEdit.setGeometry(QRect(22, 423, 181, 20))
        self.simulate = QPushButton(Dialog)
        self.simulate.setObjectName(u"simulate")
        self.simulate.setGeometry(QRect(114, 530, 81, 23))
        self.simulate.clicked.connect(self.run_simulation)
        self.n_stack = QLineEdit(Dialog)
        self.n_stack.setObjectName(u"n_stack")
        self.n_stack.setGeometry(QRect(379, 559, 133, 17))
        self.label_20 = QLabel(Dialog)
        self.label_20.setObjectName(u"label_20")
        self.label_20.setGeometry(QRect(230, 563, 85, 10))
        self.stack_method = QComboBox(Dialog)
        self.stack_method.addItem("")
        self.stack_method.addItem("")
        self.stack_method.setObjectName(u"stack_method")
        self.stack_method.setGeometry(QRect(379, 582, 133, 17))
        self.label_21 = QLabel(Dialog)
        self.label_21.setObjectName(u"label_21")
        self.label_21.setGeometry(QRect(230, 587, 47, 13))
        self.output_hist = QFrame(Dialog)
        self.output_hist.setObjectName(u"output_hist")
        self.output_hist.setGeometry(QRect(520, 380, 291, 221))
        self.output_hist.setFrameShape(QFrame.StyledPanel)
        #self.output_hist.setFrameShadow(QFrame.Raised)
        self.hist_select = QComboBox(self.output_hist)
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.addItem("")
        self.hist_select.setObjectName(u"hist_select")
        self.hist_select.setGeometry(QRect(0, 0, 69, 21))
        self.label_22 = QLabel(Dialog)
        self.label_22.setObjectName(u"label_22")
        self.label_22.setGeometry(QRect(320, 533, 101, 21))
        self.label_22.setStyleSheet(u"font: 25 12pt \"Bookman Old Style\";")
        self.set_default = QPushButton(Dialog)
        self.set_default.setObjectName(u"set_default")
        self.set_default.setGeometry(QRect(21, 530, 81, 23))
        self.set_default.clicked.connect(self.setdefault_btn)
        self.label_23 = QLabel(Dialog)
        self.label_23.setObjectName(u"label_23")
        self.label_23.setGeometry(QRect(632, 55, 31, 16))
        self.qe_sigma = QLineEdit(Dialog)
        self.qe_sigma.setObjectName(u"qe_sigma")
        self.qe_sigma.setGeometry(QRect(669, 74, 133, 16))
        self.label_24 = QLabel(Dialog)
        self.label_24.setObjectName(u"label_24")
        self.label_24.setGeometry(QRect(632, 74, 31, 16))
        self.fov = QLineEdit(Dialog)
        self.fov.setObjectName(u"fov")
        self.fov.setGeometry(QRect(120, 557, 81, 16))
        self.label_25 = QLabel(Dialog)
        self.label_25.setObjectName(u"label_25")
        self.label_25.setGeometry(QRect(20, 557, 91, 16))
        self.display_out = QLabel(Dialog)
        self.display_out.setObjectName(u"display_out")
        self.display_out.setGeometry(QRect(524, 350, 191, 16))
        self.display_out.setStyleSheet("QLabel {  color : red; }")
        self.retranslateUi(Dialog)
            
        self.count_sims = 0
        QMetaObject.connectSlotsByName(Dialog)
    # setupUi
    
    def selectFile(self):
        self.lineEdit.setText(QFileDialog.getOpenFileName()[0])
        if len(self.lineEdit.text())>1:
            self.df = pd.read_csv(self.lineEdit.text())
            
    def setdefault_btn(self):
        self.mode.setCurrentText('HCIPy')
        self.axis.setCurrentText('Off')
        self.exp_time.setText("60")
        
        if self.df is None:
            ra = [0]
            dec= [0]
            mag_nuv = [9]
            self.df = pd.DataFrame(zip(ra,dec,mag_nuv), columns = ['ra','dec','mag_nuv'])
        psf = PSF(df = self.df, mode ='HCIPy', axis = 'Off', exp_time = 60)
        
        self.qe_mean.setText(str(psf.params['qe']))
        self.qe_sigma.setText(str(psf.params['qe_sigma']))
        self.bias.setText(str(psf.params['bias']))
        self.gain.setText(str(psf.params['gain']))
        self.bit_res.setText(str(psf.params['bit_res']))
        self.shot_noise_type.setCurrentText('Gaussian')
        self.PRNU_frac.setText(str(psf.params['PRNU_frac']))
        self.sky_mag.setText(str(psf.params['M_sky']))
        self.DFM.setText(str(psf.params['DFM']))
        self.Temp.setText(str(psf.params['T']))
        self.RN.setText(str(psf.params['RN']))
        self.DN.setText(str(psf.params['DN']))
        
        
    def run_simulation(self):
        self.display_out.setText("")
        if self.df is not None:
            self.count_sims+=1
            mode = self.mode.currentText()
            axis = self.axis.currentText()
            exp_time = float(self.exp_time.text()) if np.all(np.char.isnumeric(self.exp_time.text().split('.'))) else 600 
            self.exp_time.setText(f"{exp_time}")
            
            # FoV
            if np.all(np.char.isnumeric(self.fov.text().split('.'))) :
                fov = float(self.fov.text())        
            else:
                fov = 0.02
                self.fov.setText(f"{fov}")
                
            self.psf = PSF(df = self.df, mode = mode, 
                      axis = axis, exp_time = exp_time, fov = fov)
            
            params = {}
            # Shot noise
            if self.shot_noise_cb.isChecked() :
                self.psf.shot_noise = True
                params['shot_noise'] = self.shot_noise_type.currentText()
            else:
                self.psf.shot_noise = False
            
            # Sky                
            if self.sky_cb.isChecked() :
                self.psf.sky = True
                if np.all(np.char.isnumeric(self.sky_mag.text().split('.'))):
                    params['M_sky'] = float(self.sky_mag.text())        
                else:
                    self.sky_mag.setText(f"{self.psf.params['M_sky']}")
            else:
                self.psf.sky = False
            
            # QE Mean
            if np.all(np.char.isnumeric(self.qe_mean.text().split('.'))) :
                params['qe'] = float(self.qe_mean.text())        
            else:
                self.qe_mean.setText(f"{self.psf.params['qe']}")
            # QE Sigma  
            if np.all(np.char.isnumeric(self.qe_sigma.text().split('.'))) :
                params['qe_sigma'] = float(self.qe_sigma.text())        
            else:
                self.qe_sigma.setText(f"{self.psf.params['qe_sigma']}")
            
            # Bias
            if np.all(np.char.isnumeric(self.bias.text().split('.'))) :
                params['bias'] = float(self.bias.text())        
            else:
                self.bias.setText(f"{self.psf.params['bias']}")
                
            # Gain
            if np.all(np.char.isnumeric(self.gain.text().split('.'))) :
                params['gain'] = float(self.gain.text())        
            else:
                self.gain.setText(f"{self.psf.params['gain']}")
                
            # Bit Resolution
            if np.all(np.char.isnumeric(self.bit_res.text().split('.'))) :
                params['bit_res'] = float(self.bit_res.text())        
            else:
                self.bit_res.setText(f"{self.psf.params['bit_res']}")
                
            # Read Noise
            if np.all(np.char.isnumeric(self.RN.text().split('.'))) :
                params['RN'] = float(self.RN.text())        
            else:
                self.RN.setText(f"{self.psf.params['RN']}")
                
            # PRNU
            if self.PRNU_cb.isChecked() :
                self.psf.PRNU = True
                if np.all(np.char.isnumeric(self.PRNU_frac.text().split('.'))) :
                    params['PRNU_frac'] = float(self.PRNU_frac.text())        
                else:
                    self.PRNU_frac.setText(f"{self.psf.params['PRNU_frac']}")
     
            else:
                self.psf.PRNU = False
                
            # Dark current
            if self.DC_cb.isChecked() :
                self.psf.DC = True
                if np.all(np.char.isnumeric(self.DFM.text().split('.'))) :
                    params['DFM'] = float(self.DFM.text())        
                else:
                    self.DFM.setText(f"{self.psf.params['DFM']}")
                    
                if np.all(np.char.isnumeric(self.Temp.text().split('.'))) :
                    params['T'] = float(self.Temp.text())        
                else:
                    self.Temp.setText(f"{self.psf.params['T']}")

            else:
                self.psf.DC = False
                
            # DNFP
            if self.DNFP_cb.isChecked() :
                self.psf.DNFP = True
                if np.all(np.char.isnumeric(self.DN.text().split('.'))) :
                    params['DN'] = float(self.DN.text())        
                else:
                    self.DN.setText(f"{self.psf.params['DN']}")

            else:
                self.psf.DNFP = False
                
            
            # Stacking
            stack_method = 'median'
            if self.n_stack.text().isnumeric() :
                n_stack = int(self.n_stack.text())    
                if n_stack>1:
                    stack_method = self.stack_method.currentText()
            else:
                self.n_stack.setText("1")
           
            self.psf(params = params, n_stack=n_stack, stack_type = stack_method)
            
            # Image
            try:
                fig_img, ax_img = self.psf.show_image(self.image_select.currentText())
            except:
                self.display_out.setText(f"Error: {self.image_select.currentText()} data not generated")
                self.image_select.setCurrentText('Readout')
                fig_img, ax_img = self.psf.show_image(self.image_select.currentText())

            if self.count_sims>1:
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
                
            # Hist
            try:
                fig_hist, ax_hist = self.psf.show_hist(self.hist_select.currentText())
            except:
                self.display_out.setText(f"Error: {self.hist_select.currentText()} data not generated")
                self.hist_select.setCurrentText('Readout')
                fig_hist, ax_hist = self.psf.show_hist(self.hist_select.currentText())
                
            if self.count_sims>1:
                plt.close()
                self.hist_box.removeWidget(self.canvas_hist)
                self.canvas_hist.deleteLater()                
                self.figure_hist = fig_hist
                self.canvas_hist = FigureCanvas(self.figure_hist)
                self.hist_box.addWidget(self.canvas_hist)
                self.hist_select.raise_()
            else :
                self.figure_hist = fig_hist
                self.canvas_hist = FigureCanvas(self.figure_hist)
                self.hist_box = QtWidgets.QHBoxLayout(self.output_hist)
                self.hist_box.setObjectName('hist_box')
                self.hist_box.addWidget(self.canvas_hist)
                self.hist_select.raise_()

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"PISTA", None))
        self.label_3.setText(QCoreApplication.translate("Dialog", u"Mode ", None))
        self.label_4.setText(QCoreApplication.translate("Dialog", u"Axis", None))
        self.mode.setItemText(0, QCoreApplication.translate("Dialog", u"HCIPy", None))
        self.mode.setItemText(1, QCoreApplication.translate("Dialog", u"Zeemax", None))

        self.axis.setItemText(0, QCoreApplication.translate("Dialog", u"Off", None))
        self.axis.setItemText(1, QCoreApplication.translate("Dialog", u"On", None))

        self.label_5.setText(QCoreApplication.translate("Dialog", u"Detector", None))
        self.qe_mean.setText(QCoreApplication.translate("Dialog", u"0.5", None))
        self.label_6.setText(QCoreApplication.translate("Dialog", u"Bias  (photo electrons)", None))
        self.bias.setText(QCoreApplication.translate("Dialog", u"100", None))
        self.label_8.setText(QCoreApplication.translate("Dialog", u"Quantum Efficiency", None))
        self.bit_res.setText(QCoreApplication.translate("Dialog", u"16", None))
        self.label_7.setText(QCoreApplication.translate("Dialog", u"Bit resolution", None))
        self.label_9.setText(QCoreApplication.translate("Dialog", u"Gain", None))
        self.gain.setText(QCoreApplication.translate("Dialog", u"1", None))
        self.label_10.setText(QCoreApplication.translate("Dialog", u"Read Noise (e/s/pix)", None))
        self.RN.setText(QCoreApplication.translate("Dialog", u"0", None))
        self.PRNU_cb.setText(QCoreApplication.translate("Dialog", u"PRNU", None))
        self.label_11.setText(QCoreApplication.translate("Dialog", u"PRNU Fraction", None))
        self.PRNU_frac.setText(QCoreApplication.translate("Dialog", u"0.5", None))
        self.DC_cb.setText(QCoreApplication.translate("Dialog", u"Dark Current", None))
        self.label_12.setText(QCoreApplication.translate("Dialog", u"Dark Current FoM", None))
        self.DFM.setText(QCoreApplication.translate("Dialog", u"0.01", None))
        self.DNFP_cb.setText(QCoreApplication.translate("Dialog", u"Dark Noise Fixed Pattern", None))
        self.DN.setText(QCoreApplication.translate("Dialog", u"0.5", None))
        self.label_13.setText(QCoreApplication.translate("Dialog", u"DN", None))
        self.Temp.setText(QCoreApplication.translate("Dialog", u"273", None))
        self.label_14.setText(QCoreApplication.translate("Dialog", u"Temperature (K)", None))
        self.label_15.setText(QCoreApplication.translate("Dialog", u"Observation", None))
        self.sky_cb.setText(QCoreApplication.translate("Dialog", u"Sky", None))
        self.sky_mag.setText(QCoreApplication.translate("Dialog", u"27.5", None))
        self.label_16.setText(QCoreApplication.translate("Dialog", u"Sky Magnitude", None))
        self.label_17.setText(QCoreApplication.translate("Dialog", u"Source", None))
        self.shot_noise_cb.setText(QCoreApplication.translate("Dialog", u"Shot Noise", None))
        self.shot_noise_type.setItemText(0, QCoreApplication.translate("Dialog", u"Gaussian", None))
        self.shot_noise_type.setItemText(1, QCoreApplication.translate("Dialog", u"Poisson", None))

        self.label_19.setText(QCoreApplication.translate("Dialog", u"Type", None))
        self.label_18.setText(QCoreApplication.translate("Dialog", u"Exposure time", None))
        self.exp_time.setText(QCoreApplication.translate("Dialog", u"600", None))
        self.image_select.setItemText(0, QCoreApplication.translate("Dialog", u"Readout", None))
        self.image_select.setItemText(1, QCoreApplication.translate("Dialog", u"Sky", None))
        self.image_select.setItemText(2, QCoreApplication.translate("Dialog", u"DC", None))
        self.image_select.setItemText(3, QCoreApplication.translate("Dialog", u"QE", None))
        self.image_select.setItemText(4, QCoreApplication.translate("Dialog", u"Bias", None))
        self.image_select.setItemText(5, QCoreApplication.translate("Dialog", u"PRNU", None))
        self.image_select.setItemText(6, QCoreApplication.translate("Dialog", u"DNFP", None))

        self.upload.setText(QCoreApplication.translate("Dialog", u"Upload Data Frame", None))
        self.simulate.setText(QCoreApplication.translate("Dialog", u"Simulate Image", None))
        self.n_stack.setText(QCoreApplication.translate("Dialog", u"1", None))
        self.label_20.setText(QCoreApplication.translate("Dialog", u"n_stack", None))
        self.stack_method.setItemText(0, QCoreApplication.translate("Dialog", u"Median", None))
        self.stack_method.setItemText(1, QCoreApplication.translate("Dialog", u"Mean", None))

        self.label_21.setText(QCoreApplication.translate("Dialog", u"Method", None))
        self.hist_select.setItemText(0, QCoreApplication.translate("Dialog", u"Readout", None))
        self.hist_select.setItemText(1, QCoreApplication.translate("Dialog", u"Sky", None))
        self.hist_select.setItemText(2, QCoreApplication.translate("Dialog", u"DC", None))
        self.hist_select.setItemText(3, QCoreApplication.translate("Dialog", u"QE", None))
        self.hist_select.setItemText(4, QCoreApplication.translate("Dialog", u"Bias", None))
        self.hist_select.setItemText(5, QCoreApplication.translate("Dialog", u"PRNU", None))
        self.hist_select.setItemText(6, QCoreApplication.translate("Dialog", u"DNFP", None))

        self.label_22.setText(QCoreApplication.translate("Dialog", u"Stacking", None))
        self.set_default.setText(QCoreApplication.translate("Dialog", u"Set default", None))
        self.label_23.setText(QCoreApplication.translate("Dialog", u"Mean", None))
        self.qe_sigma.setText(QCoreApplication.translate("Dialog", u"1", None))
        self.label_24.setText(QCoreApplication.translate("Dialog", u"Sigma", None))
        self.fov.setText(QCoreApplication.translate("Dialog", u"0.02", None))
        self.label_25.setText(QCoreApplication.translate("Dialog", u"Field of View (deg)", None))
        self.display_out.setText("")
    # retranslateUi

if __name__ == "__main__":
    import sys
    app        = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui         = Ui_Dialog()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
    
    




        