# should raise and exception# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 17:50:53 2017

@author: Rui
"""

from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import gui_mce.gui_mce
import numpy as np
import matplotlib.pyplot as plt
import pickle
import plots


# Boltzmann Constant
k_B = 8.6173324*(10**(-5))  # eV K^-1
# Bohr magneton
mu_B = 5.7883818066*(10**(-5))  # eV T^-1
# Convert units in entropy plots from eV/K to another, needs to be defined but
# will not be used in the GUI
Conv = 1


class GUI(QtWidgets.QMainWindow, gui_mce.gui_mce.Ui_MainWindow):
    def __init__(self, parent=None):
        super(GUI, self).__init__(parent)
        self.setupUi(self)

        # add actions to menu buttons
        self.actionExit.triggered.connect(QtWidgets.qApp.quit)
        self.actionAbout.triggered.connect(self.about_msg)
        self.actionSuggestion.triggered.connect(self.suggestion_msg)

        # actions for database
        self.database()  # create/load database
        self.pushButton_Save.clicked.connect(self.save_update_to_db)
        self.pushButton_Delete.clicked.connect(self.delete_entry_in_db)

        self.new_selection()
        self.comboBox_Database.activated.connect(self.new_selection)

        self.pushButton_Run.clicked.connect(self.run_plots)

        # for debugging
        # self.actionDebug.triggered.connect(self.print_vars)

        # Used mainly for testing, it's boring to always type the values
        self.lineEdit_Ti.setText('5')
        self.lineEdit_Tf.setText('250')
        self.lineEdit_DeltaT.setText('1')
        self.lineEdit_Bi.setText('0')
        self.lineEdit_Bf.setText('4')
        self.lineEdit_DeltaB.setText('1')
        self.lineEdit_Deltasigma.setText('0.001')

    def database(self):
        try:
            # reload dictionary from file
            database_file = open('DataBase_MCE', 'rb')
            self.db_dict = pickle.load(database_file)
            database_file.close()
        except IOError:
            print 'Creating new database'
            # create file with Gd and Tb compounds if database does not exist
            self.db_dict = {
            'Gd5Si2Ge2': {'J': 7/2., 'gJ': 2., 'TC1': 251., 'TC2': 308.,
                          'DeltaF': 0.36, 'ThetaD1': 250., 'ThetaD2': 278.,
                          'N': 36., 'Nm': 20.},
            'Tb5Si2Ge2': {'J': 6., 'gJ': 3/2., 'TC1': 112., 'TC2': 200.,
                          'DeltaF': 0.11, 'ThetaD1': 153., 'ThetaD2': 170.,
                          'N': 36., 'Nm': 20.}}
            database_file = open('DataBase_MCE', 'wb')
            pickle.dump(self.db_dict, database_file)
            database_file.close()

        self.comboBox_Database.clear()
        self.comboBox_Database.addItems(self.db_dict.keys())

    def print_vars(self):  # for debugging
        print "===================================================="
        print "J =", self.lineEdit_J.text(), type(self.lineEdit_J.text())
        print "gJ =", self.lineEdit_gJ.text()
        print "TC1 =", self.lineEdit_TC1.text()
        print "TC2 =", self.lineEdit_TC2.text()
        print "ThetaD1 =", self.lineEdit_ThetaD1.text()
        print "ThetaD2 =", self.lineEdit_ThetaD2.text()
        print "DeltaF =", self.lineEdit_DeltaF.text()
        print "N =", self.lineEdit_N.text()
        print "Nm =", self.lineEdit_Nm.text()
        print "----------------------------------------------------"
        print "Ti =", self.lineEdit_Ti.text()
        print "Tf =", self.lineEdit_Tf.text()
        print "DeltaT =", self.lineEdit_DeltaT.text()
        print "Bi =", self.lineEdit_Bi.text()
        print "Bf =", self.lineEdit_Bf.text()
        print "DeltaB =", self.lineEdit_DeltaB.text()
        print "Deltasigma =", self.lineEdit_Deltasigma.text()
        print "----------------------------------------------------"
        print "MvsT =", self.checkBox_M_vs_T.checkState()
        print "MvsB =", self.checkBox_M_vs_B.checkState()
        print "MvsTB =", self.checkBox_M_vs_TB.checkState()
        print "UvsT =", self.checkBox_U_vs_T.checkState()
        print "M_hys_vs_T =", self.checkBox_Mhys_vs_TB.checkState()
        print "FvsT =", self.checkBox_F_vs_T.checkState()
        print "FvsB =", self.checkBox_F_vs_B.checkState()
        print "trans_temp =", self.checkBox_Ts_vs_B.checkState()
        print "F_M_vs_T =", self.checkBox_FM_vs_T.checkState()
        print "F_M_vs_M =", self.checkBox_FM_vs_M.checkState()
        print "F_L_vs_T =", self.checkBox_FL_vs_T.checkState()
        print "FtotvsM =", self.checkBox_Ftot_vs_M.checkState()
        print "Ftot_heatcool =", self.checkBox_Ftotheatcool_vs_TB.checkState()
        print "S_M_vs_T =", self.checkBox_SM_vs_T.checkState()
        print "S_L_VS_T =", self.checkBox_SL_vs_T.checkState()
        print "S_tot_vs_T =", self.checkBox_Stot_vs_T.checkState()
        print "DeltaS_vs_T =", self.checkBox_DeltaS_vs_T.checkState()
        print "max_DeltaS_vs_B =", self.checkBox_maxDeltaS_vs_B.checkState()
        print "S_M_vs_M =", self.checkBox_SM_vs_M.checkState()
        print "save =", self.checkBox_Savetxt.checkState()

    def check_parameters(self):
        if (self.is_number(self.lineEdit_J.text()) and
            self.is_number(self.lineEdit_gJ.text()) and
            self.is_number(self.lineEdit_TC1.text()) and
            self.is_number(self.lineEdit_TC2.text()) and
            self.is_number(self.lineEdit_ThetaD1.text()) and
            self.is_number(self.lineEdit_ThetaD2.text()) and
            self.is_number(self.lineEdit_DeltaF.text()) and
            self.is_number(self.lineEdit_N.text()) and
            self.is_number(self.lineEdit_Nm.text())
            ):
            return True
        else:
            return False

    def save_update_to_db(self):
        if self.check_parameters:
            print "Saving to db"
            self.db_dict[str(self.lineEdit_Name.text())] = {
                'J': float(self.lineEdit_J.text()),
                'gJ': float(self.lineEdit_gJ.text()),
                'TC1': float(self.lineEdit_TC1.text()),
                'TC2': float(self.lineEdit_TC2.text()),
                'DeltaF': float(self.lineEdit_DeltaF.text()),
                'ThetaD1': float(self.lineEdit_ThetaD1.text()),
                'ThetaD2': float(self.lineEdit_ThetaD2.text()),
                'N': float(self.lineEdit_N.text()),
                'Nm': float(self.lineEdit_Nm.text())}
            # open database file and save the new entry to it
            database_file = open('DataBase_MCE', 'wb')
            pickle.dump(self.db_dict, database_file)
            database_file.close()

            self.comboBox_Database.clear()
            # update values in combobox
            self.comboBox_Database.addItems(self.db_dict.keys())
        else:
            QtWidgets.QMessageBox.about(self, 'Error', 'ValueError: could not convert string to float.\n\nCheck entries.')

    def delete_entry_in_db(self):
        name = str(self.comboBox_Database.currentText())
        del self.db_dict[name]

        # open database file and update it
        database_file = open('DataBase_MCE', 'wb')
        pickle.dump(self.db_dict, database_file)
        database_file.close()

        self.comboBox_Database.clear()
        # update values in combobox
        self.comboBox_Database.addItems(self.db_dict.keys())

    def new_selection(self):
        name = self.comboBox_Database.currentText()
        self.lineEdit_Name.setText(name)
        self.lineEdit_J.setText(str(self.db_dict[name]['J']))
        self.lineEdit_gJ.setText(str(self.db_dict[name]['gJ']))
        self.lineEdit_TC1.setText(str(self.db_dict[name]['TC1']))
        self.lineEdit_TC2.setText(str(self.db_dict[name]['TC2']))
        self.lineEdit_DeltaF.setText(str(self.db_dict[name]['DeltaF']))
        self.lineEdit_ThetaD1.setText(str(self.db_dict[name]['ThetaD1']))
        self.lineEdit_ThetaD2.setText(str(self.db_dict[name]['ThetaD2']))
        self.lineEdit_N.setText(str(self.db_dict[name]['N']))
        self.lineEdit_Nm.setText(str(self.db_dict[name]['Nm']))

    def run_plots(self):
        print "\n"
        # print self.db_dict

        if self.update_variables():
            plots.plt_M(MvsT=self.checkBox_M_vs_T.checkState(),
                        MvsB=self.checkBox_M_vs_B.checkState(),
                        MvsTB=self.checkBox_M_vs_TB.checkState(),
                        UvsT=self.checkBox_U_vs_T.checkState(),
                        M_hys_vs_T=self.checkBox_Mhys_vs_TB.checkState(),
                        save=self.checkBox_Savetxt.checkState(),
                        TT=self.TT,
                        BB=self.BB,
                        J1=float(self.lineEdit_J.text()),
                        J2=float(self.lineEdit_J.text()),
                        TC1=float(self.lineEdit_TC1.text()),
                        TC2=float(self.lineEdit_TC2.text()),
                        lamb1=self.lamb1,
                        lamb2=self.lamb2,
                        Delta_T=self.Delta_T,
                        Delta_B=self.Delta_B,
                        theta_D1=float(self.lineEdit_ThetaD1.text()),
                        theta_D2=float(self.lineEdit_ThetaD2.text()),
                        gJ=float(self.lineEdit_gJ.text()),
                        F01=float(self.lineEdit_DeltaF.text()),
                        F02=0.,
                        Nm=float(self.lineEdit_Nm.text()),
                        N=float(self.lineEdit_N.text())
                        )

            plots.plt_F(FvsT=self.checkBox_F_vs_T.checkState(),
                        FvsB=self.checkBox_F_vs_B.checkState(),
                        trans_temp=self.checkBox_Ts_vs_B.checkState(),
                        F_M_vs_T=self.checkBox_FM_vs_T.checkState(),
                        F_M_vs_M=self.checkBox_FM_vs_M.checkState(),
                        F_L_vs_T=self.checkBox_FL_vs_T.checkState(),
                        FtotvsM=self.checkBox_Ftot_vs_M.checkState(),
                        Ftot_heatcool=self.checkBox_Ftotheatcool_vs_TB.checkState(),
                        save=self.checkBox_Savetxt.checkState(),
                        TT=self.TT,
                        BB=self.BB,
                        J1=float(self.lineEdit_J.text()),
                        J2=float(self.lineEdit_J.text()),
                        TC1=float(self.lineEdit_TC1.text()),
                        TC2=float(self.lineEdit_TC2.text()),
                        lamb1=self.lamb1,
                        lamb2=self.lamb2,
                        Delta_T=self.Delta_T,
                        Delta_B=self.Delta_B,
                        theta_D1=float(self.lineEdit_ThetaD1.text()),
                        theta_D2=float(self.lineEdit_ThetaD2.text()),
                        gJ=float(self.lineEdit_gJ.text()),
                        F01=float(self.lineEdit_DeltaF.text()),
                        F02=0.,
                        Nm=float(self.lineEdit_Nm.text()),
                        N=float(self.lineEdit_N.text()),
                        sigma=self.sigma)

            plots.plt_S(S_M_vs_T=self.checkBox_SM_vs_T.checkState(),
                        S_L_VS_T=self.checkBox_SL_vs_T.checkState(),
                        S_tot_vs_T=self.checkBox_Stot_vs_T.checkState(),
                        DeltaS_vs_T=self.checkBox_DeltaS_vs_T.checkState(),
                        max_DeltaS_vs_B=self.checkBox_maxDeltaS_vs_B.checkState(),
                        S_M_vs_M=self.checkBox_SM_vs_M.checkState(),
                        save=self.checkBox_Savetxt.checkState(),
                        TT=self.TT,
                        BB=self.BB,
                        J1=float(self.lineEdit_J.text()),
                        J2=float(self.lineEdit_J.text()),
                        TC1=float(self.lineEdit_TC1.text()),
                        TC2=float(self.lineEdit_TC2.text()),
                        lamb1=self.lamb1,
                        lamb2=self.lamb2,
                        Delta_T=self.Delta_T,
                        Delta_B=self.Delta_B,
                        theta_D1=float(self.lineEdit_ThetaD1.text()),
                        theta_D2=float(self.lineEdit_ThetaD2.text()),
                        gJ=float(self.lineEdit_gJ.text()),
                        F01=float(self.lineEdit_DeltaF.text()),
                        F02=0.,
                        Nm=float(self.lineEdit_Nm.text()),
                        N=float(self.lineEdit_N.text()),
                        sigma=self.sigma)
            plt.show()
        print 'end'

    def check_ranges(self):
        if (self.is_number(self.lineEdit_Ti.text()) and
            self.is_number(self.lineEdit_Tf.text()) and
            self.is_number(self.lineEdit_DeltaT.text()) and
            self.is_number(self.lineEdit_Bi.text()) and
            self.is_number(self.lineEdit_Bf.text()) and
            self.is_number(self.lineEdit_DeltaB.text()) and
            self.is_number(self.lineEdit_Deltasigma.text())):
            return True
        else:
            return False

    def update_variables(self):
        if self.check_parameters() and self.check_ranges():
            print 'Updating variables.'
            # Temperature interval, in Kelvin
            self.Delta_T = np.arange(float(self.lineEdit_Ti.text()),
                                     float(self.lineEdit_Tf.text())+float(self.lineEdit_DeltaT.text()),
                                     float(self.lineEdit_DeltaT.text()))

            # Temperature interval, in Kelvin
            self.Delta_B = np.arange(float(self.lineEdit_Bi.text()),
                                     float(self.lineEdit_Bf.text())+float(self.lineEdit_DeltaB.text()),
                                     float(self.lineEdit_DeltaB.text()))

            # Domain, Grid with Temperature and Magnetic Field Values
            self.TT, self.BB = np.meshgrid(self.Delta_T,
                                           self.Delta_B)

            # Curie Constants divided by Vacuum Permeability (C/mu_0)
            self.C = (mu_B**2.)*float(self.lineEdit_Nm.text())*(float(self.lineEdit_gJ.text())**2.)*float(self.lineEdit_J.text())*(float(self.lineEdit_J.text()) + 1.)/(3.*k_B)

            # Parameter of the strength of the Molecular Field, lambda
            self.lamb1 = float(self.lineEdit_TC1.text())/self.C
            self.lamb2 = float(self.lineEdit_TC2.text())/self.C

            self.sigma = np.arange(-1.,
                                   1.+float(self.lineEdit_Deltasigma.text()),
                                   float(self.lineEdit_Deltasigma.text()))
            return True
        else:
            print 'Inputs missing or inputs type wrong.\nCheck inputs.'
            QtWidgets.QMessageBox.about(self, 'Error', ('Inputs missing or inputs type wrong.\nCheck inputs.'))
            return False

    def is_number(self, n):
        try:
            n = float(n)
            return True
        except Exception:
            # QtWidgets.QMessageBox.about(self, 'Error','Input can only be a number.')
            pass

    def about_msg(self):
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle("About")
        msg.setText(
        '''Magnetocaloric Effect v1.0\nCopyright (c) 2017 Rui M. Costa.
Licensed under the terms of GNU GPL v3.0

Created by Rui Costa, Bernardo Bordalo, João P. Araújo and André Pereira.

Python 2.7.13 64bits, Qt 5.6.2, PyQt5 6.0 on Windows''')
# This software was developed using PyQt, numpy,.. (acabar)\n

        msg.setIconPixmap(QtGui.QPixmap("IFIMUP_icon.png"))
        msg.exec_()

    def suggestion_msg(self):
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle("Suggestion")
        msg.setText('For bug reports and feature requests, please send to '
                    'up201100341@fc.up.pt with the subject "MCE Program - '
                    'Suggestion". Alternatively, the code is available at '
                    'https://github.com/RuiMarioCosta/MagCalc to anyone '
                    'that wants to make changes.')
        msg.exec_()


def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GUI()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
