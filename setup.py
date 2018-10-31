# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:44:59 2017

@author: Rui
"""

# setup.py
from distutils.core import setup
import py2exe
import matplotlib

import sys
sys.setrecursionlimit(5000)

DATA = [
        ("gui_mce", ["gui_mce/fig_J.png",
                     "gui_mce/fig_gJ.png",
                     "gui_mce/fig_TC1.png",
                     "gui_mce/fig_TC2.png",
                     "gui_mce/fig_ThetaD1.png",
                     "gui_mce/fig_ThetaD2.png",
                     "gui_mce/fig_DeltaF.png",
                     "gui_mce/fig_N.png",
                     "gui_mce/fig_Nm.png",
                     "gui_mce/fig_Ti.png",
                     "gui_mce/fig_Tf.png",
                     "gui_mce/fig_DeltaT.png",
                     "gui_mce/fig_Bi.png",
                     "gui_mce/fig_Bf.png",
                     "gui_mce/fig_DeltaB.png",
                     "gui_mce/fig_Deltasigma.png",
                     "gui_mce/fig_M_vs_T.png",
                     "gui_mce/fig_M_vs_B.png",
                     "gui_mce/fig_M_vs_TB.png",
                     "gui_mce/fig_U_vs_T.png",
                     "gui_mce/fig_Mhys_vs_T.png",
                     "gui_mce/fig_SM_vs_T.png",
                     "gui_mce/fig_SL_vs_T.png",
                     "gui_mce/fig_Stot_vs_T.png",
                     "gui_mce/fig_DeltaS_vs_T.png",
                     "gui_mce/fig_maxDeltaS_vs_T.png",
                     "gui_mce/fig_SM_vs_sigma.png",
                     "gui_mce/fig_F_vs_T.png",
                     "gui_mce/fig_F_vs_B.png",
                     "gui_mce/fig_Ts_vs_B.png",
                     "gui_mce/fig_FM_vs_T.png",
                     "gui_mce/fig_FM_vs_sigma.png",
                     "gui_mce/fig_FL_vs_T.png",
                     "gui_mce/fig_Ftot_vs_sigma.png",
                     "gui_mce/fig_Ftotheatcool_vs_TB.png",
                     ])
    ]
DATA += matplotlib.get_py2exe_datafiles()

#DATA=[matplotlib.get_py2exe_datafiles(),
#      ('imageformats',['C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qjpeg4.dll',
#                       'C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qgif4.dll',
#                       'C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qico4.dll',
#                       'C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qmng4.dll',
#                       'C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qsvg4.dll',
#                       'C:\\Program Files/Anaconda2/Lib/site-packages/PySide/plugins/imageformats/qtiff4.dll',
#                       'C:\\Program Files/Anaconda2/Library\plugins\imageformats\qjpeg.dll'
#                       ])]
setup(windows=[{"script":"main.py"}],
      data_files=DATA,
      options={"py2exe" : {"includes" : ["matplotlib.backends.backend_qt5agg",
                                         'scipy',
                                         'scipy.integrate',
                                         'scipy.special.*',
                                         'scipy.linalg.*',
                                         'scipy.sparse.csgraph._validation',
                                         'PyQt5.QtCore',
                                         'PyQt5.QtGui',
                                         'PyQt5.QtWidgets']}})

#setup(name = "spam",
#            version = "1.0",
#            py_modules = ['libspam'],
#            packages = ['spampkg'],
#            scripts = ['runspam.py'],
#            )