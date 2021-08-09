# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'JV_Analyser_GUI.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(1405, 962)
        MainWindow.setMinimumSize(QtCore.QSize(0, 0))
        MainWindow.setMaximumSize(QtCore.QSize(1600, 1100))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/icons/NOA1.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setAutoFillBackground(False)
        MainWindow.setDocumentMode(False)
        MainWindow.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.figures_of_merit = QtWidgets.QGroupBox(self.centralwidget)
        self.figures_of_merit.setGeometry(QtCore.QRect(1100, 80, 281, 281))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.figures_of_merit.setFont(font)
        self.figures_of_merit.setAlignment(QtCore.Qt.AlignCenter)
        self.figures_of_merit.setObjectName("figures_of_merit")
        self.label_voc = QtWidgets.QLabel(self.figures_of_merit)
        self.label_voc.setGeometry(QtCore.QRect(15, 25, 51, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_voc.setFont(font)
        self.label_voc.setObjectName("label_voc")
        self.result_voc = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_voc.setGeometry(QtCore.QRect(130, 25, 141, 21))
        self.result_voc.setDragEnabled(True)
        self.result_voc.setReadOnly(True)
        self.result_voc.setObjectName("result_voc")
        self.label_jsc = QtWidgets.QLabel(self.figures_of_merit)
        self.label_jsc.setGeometry(QtCore.QRect(15, 54, 91, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_jsc.setFont(font)
        self.label_jsc.setObjectName("label_jsc")
        self.result_jsc = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_jsc.setGeometry(QtCore.QRect(130, 54, 141, 20))
        self.result_jsc.setDragEnabled(True)
        self.result_jsc.setReadOnly(True)
        self.result_jsc.setObjectName("result_jsc")
        self.label_rp = QtWidgets.QLabel(self.figures_of_merit)
        self.label_rp.setGeometry(QtCore.QRect(15, 82, 81, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_rp.setFont(font)
        self.label_rp.setObjectName("label_rp")
        self.result_rp = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_rp.setGeometry(QtCore.QRect(130, 82, 141, 20))
        self.result_rp.setDragEnabled(True)
        self.result_rp.setReadOnly(True)
        self.result_rp.setObjectName("result_rp")
        self.label_rs = QtWidgets.QLabel(self.figures_of_merit)
        self.label_rs.setGeometry(QtCore.QRect(15, 110, 81, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_rs.setFont(font)
        self.label_rs.setObjectName("label_rs")
        self.result_rs = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_rs.setGeometry(QtCore.QRect(130, 110, 141, 21))
        self.result_rs.setDragEnabled(True)
        self.result_rs.setReadOnly(True)
        self.result_rs.setObjectName("result_rs")
        self.label_vmp = QtWidgets.QLabel(self.figures_of_merit)
        self.label_vmp.setGeometry(QtCore.QRect(15, 139, 51, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_vmp.setFont(font)
        self.label_vmp.setObjectName("label_vmp")
        self.result_vmp = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_vmp.setGeometry(QtCore.QRect(130, 139, 141, 20))
        self.result_vmp.setDragEnabled(True)
        self.result_vmp.setReadOnly(True)
        self.result_vmp.setObjectName("result_vmp")
        self.label_jmp = QtWidgets.QLabel(self.figures_of_merit)
        self.label_jmp.setGeometry(QtCore.QRect(15, 167, 91, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_jmp.setFont(font)
        self.label_jmp.setObjectName("label_jmp")
        self.result_jmp = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_jmp.setGeometry(QtCore.QRect(130, 167, 141, 20))
        self.result_jmp.setDragEnabled(True)
        self.result_jmp.setReadOnly(True)
        self.result_jmp.setObjectName("result_jmp")
        self.label_pmax = QtWidgets.QLabel(self.figures_of_merit)
        self.label_pmax.setGeometry(QtCore.QRect(15, 195, 71, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_pmax.setFont(font)
        self.label_pmax.setObjectName("label_pmax")
        self.result_pmax = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_pmax.setGeometry(QtCore.QRect(130, 195, 141, 21))
        self.result_pmax.setDragEnabled(True)
        self.result_pmax.setReadOnly(True)
        self.result_pmax.setObjectName("result_pmax")
        self.label_ff = QtWidgets.QLabel(self.figures_of_merit)
        self.label_ff.setGeometry(QtCore.QRect(15, 224, 91, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_ff.setFont(font)
        self.label_ff.setObjectName("label_ff")
        self.result_ff = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_ff.setGeometry(QtCore.QRect(130, 224, 141, 20))
        self.result_ff.setDragEnabled(True)
        self.result_ff.setReadOnly(True)
        self.result_ff.setObjectName("result_ff")
        self.label_eff = QtWidgets.QLabel(self.figures_of_merit)
        self.label_eff.setGeometry(QtCore.QRect(15, 252, 91, 18))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_eff.setFont(font)
        self.label_eff.setObjectName("label_eff")
        self.result_eff = QtWidgets.QLineEdit(self.figures_of_merit)
        self.result_eff.setGeometry(QtCore.QRect(130, 252, 141, 20))
        self.result_eff.setDragEnabled(True)
        self.result_eff.setReadOnly(True)
        self.result_eff.setObjectName("result_eff")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(1100, 360, 281, 491))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.groupBox.setFont(font)
        self.groupBox.setAlignment(QtCore.Qt.AlignCenter)
        self.groupBox.setObjectName("groupBox")
        self.analysis_results = QtWidgets.QGroupBox(self.groupBox)
        self.analysis_results.setGeometry(QtCore.QRect(10, 120, 261, 361))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.analysis_results.setFont(font)
        self.analysis_results.setAlignment(QtCore.Qt.AlignCenter)
        self.analysis_results.setObjectName("analysis_results")
        self.label_a1dark = QtWidgets.QLabel(self.analysis_results)
        self.label_a1dark.setGeometry(QtCore.QRect(12, 157, 51, 16))
        self.label_a1dark.setObjectName("label_a1dark")
        self.label_j0light = QtWidgets.QLabel(self.analysis_results)
        self.label_j0light.setGeometry(QtCore.QRect(12, 326, 121, 16))
        self.label_j0light.setObjectName("label_j0light")
        self.label_gsdark = QtWidgets.QLabel(self.analysis_results)
        self.label_gsdark.setGeometry(QtCore.QRect(12, 22, 121, 16))
        self.label_gsdark.setObjectName("label_gsdark")
        self.label_gslight = QtWidgets.QLabel(self.analysis_results)
        self.label_gslight.setGeometry(QtCore.QRect(12, 56, 121, 16))
        self.label_gslight.setObjectName("label_gslight")
        self.label_rshdark = QtWidgets.QLabel(self.analysis_results)
        self.label_rshdark.setGeometry(QtCore.QRect(12, 90, 131, 16))
        self.label_rshdark.setObjectName("label_rshdark")
        self.label_rshlight = QtWidgets.QLabel(self.analysis_results)
        self.label_rshlight.setGeometry(QtCore.QRect(12, 123, 131, 16))
        self.label_rshlight.setObjectName("label_rshlight")
        self.label_a1light = QtWidgets.QLabel(self.analysis_results)
        self.label_a1light.setGeometry(QtCore.QRect(12, 191, 51, 16))
        self.label_a1light.setObjectName("label_a1light")
        self.label_a2dark = QtWidgets.QLabel(self.analysis_results)
        self.label_a2dark.setGeometry(QtCore.QRect(12, 225, 61, 16))
        self.label_a2dark.setObjectName("label_a2dark")
        self.label_a2light = QtWidgets.QLabel(self.analysis_results)
        self.label_a2light.setGeometry(QtCore.QRect(12, 259, 51, 16))
        self.label_a2light.setObjectName("label_a2light")
        self.label_j0dark = QtWidgets.QLabel(self.analysis_results)
        self.label_j0dark.setGeometry(QtCore.QRect(12, 292, 121, 16))
        self.label_j0dark.setObjectName("label_j0dark")
        self.lineEdit_gsdark = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_gsdark.setGeometry(QtCore.QRect(160, 24, 91, 22))
        self.lineEdit_gsdark.setText("")
        self.lineEdit_gsdark.setDragEnabled(True)
        self.lineEdit_gsdark.setReadOnly(True)
        self.lineEdit_gsdark.setCursorMoveStyle(QtCore.Qt.LogicalMoveStyle)
        self.lineEdit_gsdark.setClearButtonEnabled(False)
        self.lineEdit_gsdark.setObjectName("lineEdit_gsdark")
        self.lineEdit_gslight = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_gslight.setGeometry(QtCore.QRect(160, 58, 91, 22))
        self.lineEdit_gslight.setDragEnabled(True)
        self.lineEdit_gslight.setReadOnly(True)
        self.lineEdit_gslight.setObjectName("lineEdit_gslight")
        self.lineEdit_rshdark = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_rshdark.setGeometry(QtCore.QRect(160, 91, 91, 22))
        self.lineEdit_rshdark.setDragEnabled(True)
        self.lineEdit_rshdark.setReadOnly(True)
        self.lineEdit_rshdark.setObjectName("lineEdit_rshdark")
        self.lineEdit_rshlight = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_rshlight.setGeometry(QtCore.QRect(160, 125, 91, 22))
        self.lineEdit_rshlight.setDragEnabled(True)
        self.lineEdit_rshlight.setReadOnly(True)
        self.lineEdit_rshlight.setObjectName("lineEdit_rshlight")
        self.lineEdit_a1dark = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_a1dark.setGeometry(QtCore.QRect(160, 157, 91, 22))
        self.lineEdit_a1dark.setDragEnabled(True)
        self.lineEdit_a1dark.setReadOnly(True)
        self.lineEdit_a1dark.setObjectName("lineEdit_a1dark")
        self.lineEdit_a1light = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_a1light.setGeometry(QtCore.QRect(160, 191, 91, 22))
        self.lineEdit_a1light.setDragEnabled(True)
        self.lineEdit_a1light.setReadOnly(True)
        self.lineEdit_a1light.setObjectName("lineEdit_a1light")
        self.lineEdit_a2dark = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_a2dark.setGeometry(QtCore.QRect(160, 225, 91, 22))
        self.lineEdit_a2dark.setDragEnabled(True)
        self.lineEdit_a2dark.setReadOnly(True)
        self.lineEdit_a2dark.setObjectName("lineEdit_a2dark")
        self.lineEdit_a2light = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_a2light.setGeometry(QtCore.QRect(160, 259, 91, 22))
        self.lineEdit_a2light.setDragEnabled(True)
        self.lineEdit_a2light.setReadOnly(True)
        self.lineEdit_a2light.setObjectName("lineEdit_a2light")
        self.lineEdit_j0dark = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_j0dark.setGeometry(QtCore.QRect(160, 294, 91, 22))
        self.lineEdit_j0dark.setDragEnabled(True)
        self.lineEdit_j0dark.setReadOnly(True)
        self.lineEdit_j0dark.setObjectName("lineEdit_j0dark")
        self.lineEdit_j0light = QtWidgets.QLineEdit(self.analysis_results)
        self.lineEdit_j0light.setGeometry(QtCore.QRect(160, 328, 91, 22))
        self.lineEdit_j0light.setDragEnabled(True)
        self.lineEdit_j0light.setReadOnly(True)
        self.lineEdit_j0light.setObjectName("lineEdit_j0light")
        self.analysis_range = QtWidgets.QGroupBox(self.groupBox)
        self.analysis_range.setGeometry(QtCore.QRect(10, 20, 261, 101))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.analysis_range.setFont(font)
        self.analysis_range.setAlignment(QtCore.Qt.AlignCenter)
        self.analysis_range.setObjectName("analysis_range")
        self.label_conductance = QtWidgets.QLabel(self.analysis_range)
        self.label_conductance.setGeometry(QtCore.QRect(13, 31, 80, 16))
        self.label_conductance.setObjectName("label_conductance")
        self.label_ideality_factor = QtWidgets.QLabel(self.analysis_range)
        self.label_ideality_factor.setGeometry(QtCore.QRect(12, 61, 84, 16))
        self.label_ideality_factor.setObjectName("label_ideality_factor")
        self.spin_conductance_min = QtWidgets.QDoubleSpinBox(self.analysis_range)
        self.spin_conductance_min.setGeometry(QtCore.QRect(121, 31, 61, 22))
        self.spin_conductance_min.setMinimum(-99.0)
        self.spin_conductance_min.setSingleStep(0.01)
        self.spin_conductance_min.setProperty("value", -0.5)
        self.spin_conductance_min.setObjectName("spin_conductance_min")
        self.spin_conductance_max = QtWidgets.QDoubleSpinBox(self.analysis_range)
        self.spin_conductance_max.setGeometry(QtCore.QRect(191, 31, 61, 22))
        self.spin_conductance_max.setMinimum(-99.0)
        self.spin_conductance_max.setSingleStep(0.01)
        self.spin_conductance_max.setObjectName("spin_conductance_max")
        self.spin_ideality_min = QtWidgets.QDoubleSpinBox(self.analysis_range)
        self.spin_ideality_min.setGeometry(QtCore.QRect(121, 61, 61, 22))
        self.spin_ideality_min.setMinimum(-99.0)
        self.spin_ideality_min.setSingleStep(0.01)
        self.spin_ideality_min.setProperty("value", 0.4)
        self.spin_ideality_min.setObjectName("spin_ideality_min")
        self.spin_ideality_max = QtWidgets.QDoubleSpinBox(self.analysis_range)
        self.spin_ideality_max.setGeometry(QtCore.QRect(191, 61, 61, 22))
        self.spin_ideality_max.setSingleStep(0.01)
        self.spin_ideality_max.setProperty("value", 0.6)
        self.spin_ideality_max.setObjectName("spin_ideality_max")
        self.logo = QtWidgets.QLabel(self.centralwidget)
        self.logo.setGeometry(QtCore.QRect(1160, -10, 211, 101))
        self.logo.setScaledContents(False)
        self.logo.setObjectName("logo")
        self.label_loaded_file = QtWidgets.QLabel(self.centralwidget)
        self.label_loaded_file.setGeometry(QtCore.QRect(10, 880, 49, 16))
        self.label_loaded_file.setObjectName("label_loaded_file")
        self.lineEdit_loaded_file = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_loaded_file.setGeometry(QtCore.QRect(70, 880, 561, 21))
        self.lineEdit_loaded_file.setObjectName("lineEdit_loaded_file")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 10, 1081, 861))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.tabWidget.setFont(font)
        self.tabWidget.setAutoFillBackground(False)
        self.tabWidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.tabWidget.setObjectName("tabWidget")
        self.plots_tab = QtWidgets.QWidget()
        self.plots_tab.setObjectName("plots_tab")
        self.layoutWidget = QtWidgets.QWidget(self.plots_tab)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 1041, 801))
        self.layoutWidget.setObjectName("layoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.layoutWidget)
        self.gridLayout.setContentsMargins(10, 10, 10, 10)
        self.gridLayout.setSpacing(0)
        self.gridLayout.setObjectName("gridLayout")
        self.djdv_plot = PlotWidget(self.layoutWidget)
        self.djdv_plot.setAutoFillBackground(False)
        self.djdv_plot.setObjectName("djdv_plot")
        self.gridLayout.addWidget(self.djdv_plot, 0, 1, 1, 1)
        self.dvdj_plot = PlotWidget(self.layoutWidget)
        self.dvdj_plot.setObjectName("dvdj_plot")
        self.gridLayout.addWidget(self.dvdj_plot, 1, 0, 1, 1)
        self.jjsc_plot = PlotWidget(self.layoutWidget)
        self.jjsc_plot.setObjectName("jjsc_plot")
        self.gridLayout.addWidget(self.jjsc_plot, 1, 1, 1, 1)
        self.jv_plot = PlotWidget(self.layoutWidget)
        self.jv_plot.setObjectName("jv_plot")
        self.gridLayout.addWidget(self.jv_plot, 0, 0, 1, 1)
        self.tabWidget.addTab(self.plots_tab, "")
        self.sheet_tab = QtWidgets.QWidget()
        self.sheet_tab.setObjectName("sheet_tab")
        self.tableView = QtWidgets.QTableView(self.sheet_tab)
        self.tableView.setGeometry(QtCore.QRect(10, 10, 1051, 811))
        self.tableView.setIconSize(QtCore.QSize(0, 0))
        self.tableView.setSortingEnabled(False)
        self.tableView.setObjectName("tableView")
        self.tabWidget.addTab(self.sheet_tab, "")
        self.options_tab = QtWidgets.QWidget()
        self.options_tab.setObjectName("options_tab")
        self.checkBox_autodark = QtWidgets.QCheckBox(self.options_tab)
        self.checkBox_autodark.setGeometry(QtCore.QRect(30, 60, 301, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.checkBox_autodark.setFont(font)
        self.checkBox_autodark.setChecked(True)
        self.checkBox_autodark.setTristate(False)
        self.checkBox_autodark.setObjectName("checkBox_autodark")
        self.checkBox_analyse_dark = QtWidgets.QCheckBox(self.options_tab)
        self.checkBox_analyse_dark.setGeometry(QtCore.QRect(10, 40, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.checkBox_analyse_dark.setFont(font)
        self.checkBox_analyse_dark.setChecked(True)
        self.checkBox_analyse_dark.setObjectName("checkBox_analyse_dark")
        self.checkBox_analyse_light = QtWidgets.QCheckBox(self.options_tab)
        self.checkBox_analyse_light.setGeometry(QtCore.QRect(10, 20, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.checkBox_analyse_light.setFont(font)
        self.checkBox_analyse_light.setChecked(True)
        self.checkBox_analyse_light.setObjectName("checkBox_analyse_light")
        self.tabWidget.addTab(self.options_tab, "")
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1405, 21))
        self.menubar.setNativeMenuBar(True)
        self.menubar.setObjectName("menubar")
        self.menuFiles = QtWidgets.QMenu(self.menubar)
        self.menuFiles.setObjectName("menuFiles")
        MainWindow.setMenuBar(self.menubar)
        self.actionLoad_Data = QtWidgets.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/images/icons/open.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionLoad_Data.setIcon(icon1)
        self.actionLoad_Data.setObjectName("actionLoad_Data")
        self.actionSave_Data = QtWidgets.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/images/icons/save.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave_Data.setIcon(icon2)
        self.actionSave_Data.setObjectName("actionSave_Data")
        self.actionOptions = QtWidgets.QAction(MainWindow)
        self.actionOptions.setObjectName("actionOptions")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.menuFiles.addAction(self.actionLoad_Data)
        self.menuFiles.addAction(self.actionSave_Data)
        self.menubar.addAction(self.menuFiles.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.checkBox_autodark.stateChanged['int'].connect(self.checkBox_analyse_dark.click)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "JV Analyser Mk 3"))
        self.figures_of_merit.setTitle(_translate("MainWindow", "Figures of Merit"))
        self.label_voc.setText(_translate("MainWindow", "V<sub>OC</sub> (V)"))
        self.result_voc.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_jsc.setText(_translate("MainWindow", "J<sub>SC</sub> (mA/cm<sup>2</sup>)"))
        self.result_jsc.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_rp.setText(_translate("MainWindow", "R<sub>P</sub> (mS/cm<sup>2</sup>)"))
        self.result_rp.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_rs.setText(_translate("MainWindow", "<html><head/><body><p>R<span style=\" vertical-align:sub;\">S</span> (mS/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.result_rs.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_vmp.setText(_translate("MainWindow", "<html><head/><body><p>V<span style=\" vertical-align:sub;\">MP</span> (V)</p></body></html>"))
        self.result_vmp.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_jmp.setText(_translate("MainWindow", "<html><head/><body><p>J<span style=\" vertical-align:sub;\">MP </span>(mA/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.result_jmp.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_pmax.setText(_translate("MainWindow", "<html><head/><body><p>P<span style=\" vertical-align:sub;\">MAX </span>(mW)</p></body></html>"))
        self.result_pmax.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_ff.setText(_translate("MainWindow", "<html><head/><body><p>Fill Factor (%)</p></body></html>"))
        self.result_ff.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.label_eff.setText(_translate("MainWindow", "<html><head/><body><p>Efficiency (%)</p></body></html>"))
        self.result_eff.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.groupBox.setTitle(_translate("MainWindow", "Analysis"))
        self.analysis_results.setTitle(_translate("MainWindow", "Results"))
        self.label_a1dark.setText(_translate("MainWindow", "<html><head/><body><p>A<span style=\" vertical-align:sub;\">1 </span>Dark</p></body></html>"))
        self.label_j0light.setText(_translate("MainWindow", "<html><head/><body><p>J<span style=\" vertical-align:sub;\">0</span> Light (mA/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.label_gsdark.setText(_translate("MainWindow", "<html><head/><body><p>G<span style=\" vertical-align:sub;\">S</span> Dark (mS/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.label_gslight.setText(_translate("MainWindow", "<html><head/><body><p>G<span style=\" vertical-align:sub;\">S</span> Light (mS/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.label_rshdark.setText(_translate("MainWindow", "<html><head/><body><p>R<span style=\" vertical-align:sub;\">SH</span> Dark (kΩ • cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.label_rshlight.setText(_translate("MainWindow", "<html><head/><body><p>R<span style=\" vertical-align:sub;\">SH</span> Light (kΩ • cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.label_a1light.setText(_translate("MainWindow", "<html><head/><body><p>A<span style=\" vertical-align:sub;\">1</span> Light</p></body></html>"))
        self.label_a2dark.setText(_translate("MainWindow", "<html><head/><body><p>A<span style=\" vertical-align:sub;\">2</span> Dark</p></body></html>"))
        self.label_a2light.setText(_translate("MainWindow", "<html><head/><body><p>A<span style=\" vertical-align:sub;\">2</span> Light</p></body></html>"))
        self.label_j0dark.setText(_translate("MainWindow", "<html><head/><body><p>J<span style=\" vertical-align:sub;\">0</span> Dark (mA/cm<span style=\" vertical-align:super;\">2</span>)</p></body></html>"))
        self.lineEdit_gsdark.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_gslight.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_rshdark.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_rshlight.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_a1dark.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_a1light.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_a2dark.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_a2light.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_j0dark.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.lineEdit_j0light.setPlaceholderText(_translate("MainWindow", "0.00"))
        self.analysis_range.setTitle(_translate("MainWindow", "Input Parameters (V)"))
        self.label_conductance.setText(_translate("MainWindow", "Conductance: "))
        self.label_ideality_factor.setText(_translate("MainWindow", "Ideality Factor:"))
        self.logo.setText(_translate("MainWindow", "<html><head/><body><p><img src=\":/images/icons/noa2.png\"/></p></body></html>"))
        self.label_loaded_file.setText(_translate("MainWindow", "File:"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.plots_tab), _translate("MainWindow", "Resulting Plots"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.sheet_tab), _translate("MainWindow", "Raw Data"))
        self.checkBox_autodark.setToolTip(_translate("MainWindow", "<html><head/><body><p>The program needs to load the data related to the &quot;dark JV curve&quot;. It can either assume the name, by adding &quot;_dark&quot; to the end of the loaded file, or the user can select it manually.</p></body></html>"))
        self.checkBox_autodark.setText(_translate("MainWindow", "Let JV Analyser assume \"dark\" data filename?"))
        self.checkBox_analyse_dark.setText(_translate("MainWindow", "Analyse Dark JV data?"))
        self.checkBox_analyse_light.setText(_translate("MainWindow", "Analyse Light JV data?"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.options_tab), _translate("MainWindow", "Options"))
        self.menuFiles.setTitle(_translate("MainWindow", "Files"))
        self.actionLoad_Data.setText(_translate("MainWindow", "Load Data"))
        self.actionSave_Data.setText(_translate("MainWindow", "Save Data"))
        self.actionOptions.setText(_translate("MainWindow", "Options"))
        self.actionAbout.setText(_translate("MainWindow", "About"))
from pyqtgraph import PlotWidget
import images_rc