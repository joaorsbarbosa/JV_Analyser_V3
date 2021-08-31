from PyQt5 import QtWidgets, QtCore, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox  # necessary import of the file dialog window
#from pyqtgraph import PlotWidget
import pyqtgraph as pg
from JV_Analyser_GUI import Ui_MainWindow  # this import will load the .py file with the GUI

import pandas as pd  # Pandas dataframes are a super useful tool to handle and process data.
import os, sys   # Necessary modules to allow the python code to interact with Windows.
from scipy import interpolate
from scipy import stats
from scipy import constants
import numpy as np

# TODO: Add standard errors of slope and intercepts
# TODO: Add saving feature


##########################################
############ CONFIGURATIONS ##############
##########################################
# ---------------------------------------#
# ---------------- PLOT -----------------#
# ---------------------------------------#

#QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, False)

pg.setConfigOption('background', 'w')  # Sets the background of the plots as WHITE
pg.setConfigOption('foreground', 'k')  # Sets the lines of the plots as BLACK
pg.setConfigOptions(antialias=True)  # Makes the lines smoother (removes "stair casing" effect). Purely visual setting. Does not affect the data.

# ---------------------------------------#

class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, obj=None, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)

        # ---------------- PLOT SETTINGS -----------------#

        #In order for the program to open with the plots with the correct axis names, etc, it needs to be set before the signals part of the code.
        self.jv_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J (mA/cm<sup>2</sup>)</span>")
        self.jv_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.jv_plot.showGrid(x=True, y=True)


        self.djdv_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">dJ/dV (mS cm <sup>-2</sup>)</span>")
        self.djdv_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.djdv_plot.showGrid(x=True, y=True)

        self.dvdj_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">dV/dJ (Ω cm <sup>2</sup>)</span>")
        self.dvdj_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">(J+J<sub>SC</sub>-GV)<sup>-1</sup> (mA<sup>-1</sup> cm<sup>2</sup>)</span>")
        self.dvdj_plot.showGrid(x=True, y=True)

        self.jjsc_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J+J<sub>SC</sub>-GV (mA/cm<sup>2</sup>)</span>")
        self.jjsc_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V-RJ (V)</span>")
        self.jjsc_plot.showGrid(x=True, y=True)
        self.jjsc_plot.setLogMode(False, True)

        # Some additional plots will be added, so the user can check the quality of the fits.
        # The settings of the plots will be defined

        self.interpolation_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J (mA/cm<sup>2</sup>)</span>")
        self.interpolation_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.interpolation_plot.showGrid(x=True, y=True)
        self.interpolation_plot.setTitle("Interpolated J-V")

        self.linear_reg_rs_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J (mA/cm<sup>2</sup>)</span>")
        self.linear_reg_rs_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.linear_reg_rs_plot.showGrid(x=True, y=True)
        self.linear_reg_rs_plot.setTitle("Series Resistance Linear Regression")

        self.linear_reg_rsh_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J (mA/cm<sup>2</sup>)</span>")
        self.linear_reg_rsh_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.linear_reg_rsh_plot.showGrid(x=True, y=True)
        self.linear_reg_rsh_plot.setTitle("Shunt Conductance Linear Regression")
        # ---------------------------------------#

        ################## SIGNALS ##################
        # A signal is the connection between what happens in the GUI and the python code.

        self.actionLoad_Data.triggered.connect(self.execute_processing)  # when the button "Load Data" is pressed, it will run the execute_processing function
        self.spin_conductance_min.valueChanged.connect(self.refresh_analysis)
        self.spin_conductance_max.valueChanged.connect(self.refresh_analysis)
        self.spin_ideality_min.valueChanged.connect(self.refresh_analysis)
        self.spin_ideality_max.valueChanged.connect(self.refresh_analysis)
        self.radioButton_djdv_gsh.toggled.connect(self.refresh_analysis)
        self.radioButton_jv_gsh.toggled.connect(self.refresh_analysis)
        self.actionSave_Data.triggered.connect(self.save_results)
    def execute_processing(self):
        # As far as I know, Qt slots cannot return any data. Thus, I created this function to run all the necessary code to process the data when a file is loaded.
        # The data processing starts by loading the .csv data of both light and dark curves

        # Python will run the load_data() function. If the user quits the file dialog without selecting a file, it would crash the program, if it was not for the
        # error handling "try" and "except" functions. If the user selects a file, the code inside "try:" will run. If not, an error will be raised, and the "except:"
        # code will just ignore it and do nothing.

        self.plot_light_JV = True  # For now, this feature is not implemented. JV Light data will always be required.
        self.plot_dark_JV = self.checkBox_analyse_dark.isChecked()
        try: # A janky way to handler errors but...it works 🤷‍
            self.data = self.load_data(self.plot_light_JV, self.plot_dark_JV)
            # the load_data() function will return the light JV and dark JV data in dataframes, with said dataframes inside a tuple
            dataframe_light_JV = self.data[0]  # the light JV data sits in the 1st position of the tuple
            dataframe_dark_JV = self.data[1]  # the dark JV data sits in the 2nd position of the tuple

            # The "self.figures_of_merit" will make this dataframe "global" inside the MainWindow class.
            self.figures_of_merit_data = self.compute_fom(dataframe_light_JV)
            self.update_jv_plot(dataframe_light_JV, dataframe_dark_JV, self.plot_light_JV, self.plot_dark_JV)
            self.refresh_analysis()

        except:
            pass

    def load_data(self, plot_light_JV, plot_dark_JV ):
        # the line below will open a file dialog window, and ask for the user to select a .csv file. The third part of the .getOpenFileName()
        # with "CSV Files (*.csv)" will apply a filter so that only .csv files show up. This will return a tuple with the position
        # zero containing the path and file name, and the second containing the data filter mentioned previously. Since the program
        # only requires the path and filename the ", _" is used to discard the second tuple value.
        if plot_light_JV:
            light_data_work_path, _ = QFileDialog.getOpenFileName(self, "Select an input file with Light JV data", "", "CSV Files (*.csv)")
            # with the file selected, the program needs to change it's work path. The previous function will return a string with the entire path and the name of the file included.
            # However, the os.chdir function will not work if give the filename. It needs to be removed. To do so, .rpartition is used to partition the previous string in the last "/"
            # found. This will return a tuple with the previous string split. Since the program only requires the 1st value of the tuple, the index [0] is used to specify that.
            os.chdir(light_data_work_path.rpartition('/')[0])

            # for the data ingestion, Pandas only requires the name of the file, assuming that the work path as been changed
            light_data_file = os.path.basename(light_data_work_path)

            if light_data_file:
                self.lineEdit_loaded_file.setText(light_data_file)  # with the file selected, the work path and file name will be displayed in the lineEdit "loaded_file"

                light_data = self.parse_file(light_data_file)
            if plot_dark_JV:
                if self.checkBox_autodark.isChecked():
                    try:
                        file_name_list = list(light_data_file)
                        dot_position = light_data_file.find('.')
                        file_name_list.insert(dot_position, '_dark')
                        dark_data_file = ''.join(file_name_list)
                        dark_data = self.parse_file(dark_data_file)

                    except:
                        # If the "_dark" file does not exist, an error message will be raised, and the user will have to select the dark file manually.
                        error_dialog = QMessageBox()
                        error_dialog.setIcon(QMessageBox.Critical)
                        error_dialog.setText("The automatic \'dark\' file detection failed!")
                        error_dialog.setInformativeText("Please select the dark JV data manually.")
                        error_dialog.setWindowTitle("Error!")
                        error_dialog.exec_()
                        try:
                            # Will try to get the user to select the dark JV file.
                            dark_data_work_path, _ = QFileDialog.getOpenFileName(self, "Select an input file with Dark JV data", "", "CSV Files (*.csv)")
                            dark_data_file = os.path.basename(dark_data_work_path)
                            dark_data = self.parse_file(dark_data_file)
                        except:
                            error_dialog = QMessageBox()
                            error_dialog.setIcon(QMessageBox.Critical)
                            error_dialog.setText("The user failed to select the Dark JV data!")
                            error_dialog.setWindowTitle("Error!")
                            error_dialog.exec_()
                else:
                    # if the user does not want the automatic dark file detection to occur, the user will have to select the dark file manually.
                    dark_data_work_path, _ = QFileDialog.getOpenFileName(self, "Select an input file with Dark JV data", "", "CSV Files (*.csv)")
                    dark_data_file = os.path.basename(dark_data_work_path)
                    dark_data = self.parse_file(dark_data_file)
            else:
                # if the user does not want to analyse the dark JV data, a dataframe similar to the light JV one will be created, but filled with zeros.
                dark_data = pd.DataFrame(0, columns=light_data.columns, index=light_data.index)
        else:
            try:
                dark_data_work_path, _ = QFileDialog.getOpenFileName(self, "Select an input file with Dark JV data", "", "CSV Files (*.csv)")
                os.chdir(dark_data_work_path.rpartition('/')[0])
                dark_data_file = os.path.basename(dark_data_work_path)
                dark_data = self.parse_file(dark_data_file)
                # if the user does not want to analyse the light JV data, a dataframe similar to the dark JV one will be created, but filled with zeros.
                light_data = pd.DataFrame(0, columns=dark_data.columns, index=dark_data.index)
            except:
                # if the user does not select anything, the except function will just ignore the error and end the load_data function without any data in return
                # this will generate an error which will be handled in the execute_processing function.
                pass

        light_data.I = light_data.I * 1000  # The light current is being multiplied by 1000 in order to convert it from A/cm² to mA/cm²
        dark_data.I = dark_data.I * 1000  # The dark current is being multiplied by 1000 in order to convert it from A/cm² to mA/cm²
        return light_data, dark_data

    def parse_file(self, file_name):
        # This function will take the file from the load_data function and load/process the csv data so that it can be used later
        data = pd.read_csv(file_name, header=15, sep=';')
        return data

    def update_jv_plot(self, dataframe_light_JV, dataframe_dark_JV, plot_light_JV, plot_dark_JV):
        # the purpose of this function is to update the JV plot. Basing off the Hegedus' paper (10.1002/pip.518), the JV will be the top left plot.
        # The light JV data will be plotted in red. The dark JV data will be plotted in dark, if the user wants to plot it.
        self.jv_plot.clear()
       #self.jv_plot.setBackground('w')
        self.jv_plot.addLegend()
        if plot_light_JV:  # If the user selected to process Light JV data, plot the light JV data
            pen = pg.mkPen(color=(255, 0, 0), width=3)  # To change the color of the plot, you need assign a color to the "pen". In this case (RGB) 255, 0, 0 is RED.
            # The width is also changed to 3 px wide
            self.jv_plot.plot(dataframe_light_JV.V, dataframe_light_JV.I, name="Light J-V", pen=pen, symbol='o', symbolSize=5, symbolBrush='r')  # if plot_dark_JV else None

        if plot_dark_JV: # If the user selected to process Dark JV data, plot the Dark JV data
            pen = pg.mkPen(color=(0, 0, 0), width=3)  # Changing the pen color back to black, in order to plot the Dark JV data.
            self.jv_plot.plot(dataframe_dark_JV.V, dataframe_dark_JV.I, name = "Dark J-V", pen=pen, symbol='o', symbolSize=5, symbolBrush='k')

    def update_djdv_plot(self,dataframe_light_JV, dataframe_dark_JV, plot_light_JV, plot_dark_JV):

        self.djdv_plot.clear()
        self.djdv_plot.addLegend()

        # Getting the voltage values for the dJdV analysis
        conductance_shunt_x_min = self.spin_conductance_min.value()
        conductance_shunt_x_max = self.spin_conductance_max.value()
        # Making the interpolation functions of the light and dark curves so that the 1st derivative can be calculated
        if plot_light_JV:
            # Creating an interpolation of the light JV
            djdv_light_spline = interpolate.InterpolatedUnivariateSpline(dataframe_light_JV.V, dataframe_light_JV.I)
            # Calculating the derivative to the n=1 degree (1st degree)
            djdv_light_spline_derivative = djdv_light_spline.derivative(n=1)
            # This will create a dataframe with the voltages between the selected range
            voltage_interval_light = dataframe_light_JV.V[(dataframe_light_JV["V"] >= conductance_shunt_x_min) & (dataframe_light_JV["V"] <= conductance_shunt_x_max)]
            voltage_interval_light = voltage_interval_light.reset_index(drop=True)
            # Linear regression to calculate the shunt conductance

            djdv_light_lin_regress = stats.linregress(voltage_interval_light, djdv_light_spline_derivative(voltage_interval_light)) 
            djdv_light_slope = djdv_light_lin_regress.slope  # The slope will be the shunt conductance
            djdv_light_rsquared = djdv_light_lin_regress.rvalue**2  # Calculates the r-squared
            pen = pg.mkPen(color=(255, 0, 0), width=3)  # To change the color of the plot, you need assign a color to the "pen". In this case (RGB) 255, 0, 0 is RED.
            # The width is also changed to 3 px wide
            self.djdv_plot.plot(voltage_interval_light, djdv_light_spline_derivative(voltage_interval_light), name="Light dJdV", pen=pen, symbol='o', symbolSize=5, symbolBrush='r')
            pen = pg.mkPen(color=(0, 0, 255), width=3)
            self.djdv_plot.plot(voltage_interval_light, djdv_light_lin_regress.intercept + djdv_light_lin_regress.slope * voltage_interval_light, name="Light Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="b")
        else:
            # If the user does not want to analyse the light JV data, these variables will be filled with 0
            djdv_light_slope = 0
            djdv_light_rsquared = 0
            
        if plot_dark_JV:
            # Creating an interpolation of the dark JV
            djdv_dark_spline = interpolate.InterpolatedUnivariateSpline(dataframe_dark_JV.V, dataframe_dark_JV.I)
            # Calculating the derivative to the n=1 degree (1st degree)
            djdv_dark_spline_derivative = djdv_dark_spline.derivative(n=1)
            # This will create a dataframe with the voltages between the selected range
            voltage_interval_dark = dataframe_dark_JV.V[(dataframe_dark_JV["V"] >= conductance_shunt_x_min) & (dataframe_dark_JV["V"] <= conductance_shunt_x_max)]
            voltage_interval_dark = voltage_interval_dark.reset_index(drop=True)
            djdv_dark_lin_regress = stats.linregress(voltage_interval_dark, djdv_dark_spline_derivative(voltage_interval_dark))
            djdv_dark_slope = djdv_dark_lin_regress.slope
            djdv_dark_rsquared = djdv_dark_lin_regress.rvalue**2
            pen = pg.mkPen(color=(0, 0, 0), width=3)  # Changing the pen color back to black, in order to plot the Dark JV data.
            self.djdv_plot.plot(voltage_interval_dark, djdv_dark_spline_derivative(voltage_interval_dark), name="Dark dJdV", pen=pen, symbol='o', symbolSize=5, symbolBrush='k')
            pen = pg.mkPen(color=(0, 255, 0), width=3)
            self.djdv_plot.plot(voltage_interval_dark, djdv_dark_lin_regress.intercept + djdv_dark_lin_regress.slope * voltage_interval_dark, name="Dark Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="g")
           
        else:
            # If the user does not want to analyse the dark JV data, these variables will be filled with 0
            djdv_dark_slope = 0
            djdv_dark_rsquared = 0

        djdv_results = dict({"djdv_light_condutance": djdv_light_slope, "djdv_light_rsquared": djdv_light_rsquared, "djdv_dark_condutance": djdv_dark_slope, "djdv_dark_rsquared": djdv_dark_rsquared})
        # For the values calculated by the advanced analysis
        self.result_rshlight.setText(str(round(djdv_results["djdv_light_condutance"], 2)))
        self.result_rshlight_rsquared.setText(str(round(djdv_results["djdv_light_rsquared"], 3)))
        self.result_rshdark.setText(str(round(djdv_results["djdv_dark_condutance"], 2)))
        self.result_rshdark_rsquared.setText(str(round(djdv_results["djdv_dark_rsquared"], 3)))

        return djdv_results

    def update_dvdj_plot(self,dataframe_light_JV, dataframe_dark_JV,shunt_conductance_light, shunt_conductance_dark, plot_light_JV, plot_dark_JV):

        # Lab's temperature for ideality factor calculation:
        T = 293.15  # 20℃
        q = constants.elementary_charge
        k = constants.k

        self.dvdj_plot.clear()
        self.dvdj_plot.addLegend()
        dvdj_x_min = self.spin_ideality_min.value()
        dvdj_x_max = self.spin_ideality_max.value()

        if plot_light_JV:
            # The voltage range is adjusted in accordance to the user's input
            voltage_interval_light = dataframe_light_JV.V[(dataframe_light_JV["V"] >= dvdj_x_min) & (dataframe_light_JV["V"] <= dvdj_x_max)]
            voltage_interval_light = voltage_interval_light.reset_index(drop=True)
            current_interval_light = dataframe_light_JV.I[(dataframe_light_JV["V"] >= dvdj_x_min) & (dataframe_light_JV["V"] <= dvdj_x_max)]
            current_interval_light = current_interval_light.reset_index(drop=True)
            # Now the short circuit current will be calculated
            current_spline_light = interpolate.InterpolatedUnivariateSpline(dataframe_light_JV.V, dataframe_light_JV.I)
            current_short_circuit_light = abs(current_spline_light(0))
            # The (J+Jsc)^-1 is calculated
            # The X axis is calculated with the -GV correction, per Hegedus et. al. paper
            jjsc_light = (current_spline_light(voltage_interval_light) + current_short_circuit_light - shunt_conductance_light * voltage_interval_light) ** (-1)
            # The 1st value of this array needs to be dropped, as the np.diff does not calculate the difference of the 1st value
            jjsc_light = jjsc_light.to_numpy()  # The -GV operation will transform the numpy array back into a pandas series. Therefore, it needs to be transformed into a ndarray
            jjsc_light = np.delete(jjsc_light, 0) # Eliminates the 1st value in order to correct the array's lengths
            # Now the derivative dV/dJ is calculated. numpy diff is used instead of pandas diff, since the later will return a NaN in the first row
            dVdJ_light = np.diff(voltage_interval_light)/np.diff(current_interval_light)*1000  # multiplying by a thousand to go from kΩ to Ω

            dvdj_light_lin_regress = stats.linregress(jjsc_light, dVdJ_light)  # Calculating the linear regression
            dvdj_light_slope = dvdj_light_lin_regress.slope  # The slope will be used to calculate the ideality factor
            dvdj_light_intercept = dvdj_light_lin_regress.intercept  # The intercept will be the shunt conductance
            # dvdj_light_intercept_stderr = dvdj_light_lin_regress.intercept_stderr # Standard error
            dvdj_light_rsquared = dvdj_light_lin_regress.rvalue ** 2  # Calculates the r-squared
            # For the ideality factor calculation
            A1_light = ((dvdj_light_slope/ 1000) * q) / (k * T)
            A1_light_rsquared = dvdj_light_rsquared
            # Plot
            pen = pg.mkPen(color=(255, 0, 0), width=3)  # To change the color of the plot, you need assign a color to the "pen". In this case (RGB) 255, 0, 0 is RED.
            # The width is also changed to 3 px wide
            self.dvdj_plot.plot(jjsc_light, dVdJ_light, name="Light", pen=pen, symbol='o', symbolSize=5, symbolBrush='r')
            pen = pg.mkPen(color=(0, 0, 255), width=3)
            self.dvdj_plot.plot(jjsc_light, dvdj_light_lin_regress.intercept + dvdj_light_lin_regress.slope * jjsc_light, name="Light Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="b")


        else:
            dvdj_light_slope = 0
            dvdj_light_intercept = 0
            dvdj_light_rsquared = 0
            A1_light = 0
            A1_light_rsquared = dvdj_light_rsquared
            
        if plot_dark_JV:

            # The voltage range is adjusted in accordance to the user's input
            voltage_interval_dark = dataframe_dark_JV.V[(dataframe_dark_JV["V"] >= dvdj_x_min) & (dataframe_dark_JV["V"] <= dvdj_x_max)]
            voltage_interval_dark = voltage_interval_dark.reset_index(drop=True)
            current_interval_dark = dataframe_dark_JV.I[(dataframe_dark_JV["V"] >= dvdj_x_min) & (dataframe_dark_JV["V"] <= dvdj_x_max)]
            current_interval_dark = current_interval_dark.reset_index(drop=True)

            # Now the short circuit current will be calculated
            current_spline_dark = interpolate.InterpolatedUnivariateSpline(dataframe_dark_JV.V, dataframe_dark_JV.I)
            current_short_circuit_dark = abs(current_spline_dark(0))
            # The (J+Jsc)^-1 is calculated
            jjsc_dark = (current_spline_dark(voltage_interval_dark) + current_short_circuit_dark - shunt_conductance_dark * voltage_interval_dark) ** (-1)
            jjsc_dark = jjsc_dark.to_numpy()  # The -GV operation will transform the numpy array back into a pandas series. Therefore, it needs to be transformed into a ndarray
            # The 1st value of this array needs to be dropped, as the np.diff does not calculate the difference of the 1st value
            jjsc_dark = np.delete(jjsc_dark, 0)
            # Now the derivative dV/dJ is calculated. numpy diff is used instead of pandas diff, since the later will return a NaN in the first row
            dVdJ_dark = np.diff(voltage_interval_dark) / np.diff(current_interval_dark)*1000 # multiplying by a thousand to go from kΩ to Ω

            dvdj_dark_lin_regress = stats.linregress(jjsc_dark, dVdJ_dark) # Calculating the linear regression
            dvdj_dark_slope = dvdj_dark_lin_regress.slope  # The slope will be used to calculate the ideality factor
            dvdj_dark_intercept = dvdj_dark_lin_regress.intercept  # The intercept will be the shunt conductance
            # dvdj_dark_intercept_stderr = dvdj_dark_lin_regress.intercept_stderr # Standard error
            dvdj_dark_rsquared = dvdj_dark_lin_regress.rvalue ** 2  # Calculates the r-squared
            # For the ideality factor calculation
            A1_dark = dvdj_dark_slope * q / (k * T)/1000
            A1_dark_rsquared = dvdj_dark_rsquared
            # Plot
            pen = pg.mkPen(color=(0, 0, 0), width=3)
            self.dvdj_plot.plot(jjsc_dark, dVdJ_dark, name="Dark", pen=pen, symbol='o', symbolSize=5, symbolBrush='k')
            pen = pg.mkPen(color=(0, 255, 0), width=3)
            self.dvdj_plot.plot(jjsc_dark, dvdj_dark_lin_regress.intercept + dvdj_dark_lin_regress.slope * jjsc_dark, name="Dark Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="g")
        else:
            
            dvdj_dark_slope = 0
            dvdj_dark_intercept = 0
            dvdj_dark_rsquared = 0
            A1_dark = 0
            A1_dark_rsquared = dvdj_dark_rsquared

        dvdj_results = dict({"dvdj_light_series_resistance": dvdj_light_intercept, "dvdj_light_rsquared": dvdj_light_rsquared, "A1_light": A1_light, "A1_light_rsquared": A1_light_rsquared,
                             "dvdj_dark_series_resistance": dvdj_dark_intercept, "djdv_dark_rsquared": dvdj_dark_rsquared, "A1_dark": A1_dark, "A1_dark_rsquared": A1_dark_rsquared})

        self.result_gslight.setText(str(round(dvdj_results["dvdj_light_series_resistance"], 2)))
        self.result_gslight_rsquared.setText(str(round(dvdj_results["dvdj_light_rsquared"], 3)))
        self.result_a1light.setText(str(round(A1_light, 2)))
        self.result_a1light_rsquared.setText(str(round(dvdj_results["A1_light_rsquared"], 3)))

        self.result_gsdark.setText(str(round(dvdj_results["dvdj_dark_series_resistance"], 2)))
        self.result_gsdark_rsquared.setText(str(round(dvdj_results["djdv_dark_rsquared"], 3)))
        self.result_a1dark.setText(str(round(A1_dark, 2)))
        self.result_a1dark_rsquared.setText(str(round(dvdj_results["A1_dark_rsquared"], 3)))

        return dvdj_results

    def update_jjsc_plot(self, dataframe_light_JV, dataframe_dark_JV, light_series_resistance, dark_series_resistance, shunt_conductance_light, shunt_conductance_dark, plot_light_JV, plot_dark_JV):

        # Lab's temperature for ideality factor calculation:
        T = 293.15  # 20℃
        q = constants.elementary_charge
        k = constants.k

        self.jjsc_plot.clear()
        self.jjsc_plot.addLegend()
        
        # Getting the voltage values to calculate the ideality factor from the interface
        jjsc_x_min = self.spin_ideality_min.value()
        jjsc_x_max = self.spin_ideality_max.value()
       
        if plot_light_JV:
            # The voltage interval to use for the calculations will be the one stipulated previously by the user.
            voltage_interval_light = dataframe_light_JV.V[(dataframe_light_JV["V"] >= jjsc_x_min) & (dataframe_light_JV["V"] <= jjsc_x_max)]
            # The Y axis needs the short circuit current, so it's calculated here:
            current_spline_light = interpolate.InterpolatedUnivariateSpline(dataframe_light_JV.V, dataframe_light_JV.I)
            current_short_circuit_light = abs(current_spline_light(0))
            # The X axis will the the voltage minus the series resistance multiplied by the short circuit current
            V_RJ_light = dataframe_light_JV.V - light_series_resistance * dataframe_light_JV.I / 1000  # The division by a thousand is to convert mA to Amps
            # The Y axis will be the sum of the current with the short circuit current, minus the shunt conductance multiplied by the voltage
            jjsc_gv_light = (current_spline_light(dataframe_light_JV.V) + current_short_circuit_light - shunt_conductance_light * dataframe_light_JV.V)
            # The plot will not be restricted by the voltage interval previously defined.
            # In order to plot the semi-log Y axis, the negative values need to be removed.
            # Therefore, the voltage displayed will be the voltages corresponding to positive current values
            V_RJ_light = V_RJ_light[(jjsc_gv_light > 0) & (V_RJ_light > 0)]
            # Same thing is done to the current.
            jjsc_gv_light = jjsc_gv_light[(jjsc_gv_light > 0) & (V_RJ_light > 0)]

            V_RJ_light_analysis = V_RJ_light  # This is done in order to keep the indexes, so that the data inside the selected voltage can be attained.
            V_RJ_light = V_RJ_light.reset_index(drop=True)  # The series is re-indexed back to zero

            jjsc_gv_light_analysis = jjsc_gv_light
            jjsc_gv_light = jjsc_gv_light.reset_index(drop=True)

            # Matching the current values to the selected voltage interval for analysis
            jjsc_gv_light_analysis = jjsc_gv_light_analysis.loc[voltage_interval_light.index]
            V_RJ_light_analysis = V_RJ_light_analysis.loc[voltage_interval_light.index]
            jjsc_gv_light_analysis = jjsc_gv_light_analysis.reset_index(drop=True)
            V_RJ_light_analysis = V_RJ_light_analysis.reset_index(drop=True)
            # Since it is not possible to calculate the linear regression of a logarithmic function, the Y variable goes through log()
            linear_regression_jjsc = stats.linregress(V_RJ_light_analysis, np.log(jjsc_gv_light_analysis))

            A2_light = q / (linear_regression_jjsc.slope * k * T)  # Calculating the A2 using the variables previously defined
            A2_light_rsquared = linear_regression_jjsc.rvalue**2

            # However, the J0 is calculated in the semilog plot, thus, the interception value goes trough exp()
            J0_light = np.exp(linear_regression_jjsc.intercept)
            J0_light_rsquared = linear_regression_jjsc.rvalue**2

            # Plotting in red
            pen = pg.mkPen(color=(255, 0, 0), width=3)
            self.jjsc_plot.plot(V_RJ_light, jjsc_gv_light, name="Light", pen=pen, symbol="o", symbolSize=5, symbolBrush='r')
            # Plotting the selected data
            pen = pg.mkPen(color=(0, 0, 255), width=3)
            self.jjsc_plot.plot(V_RJ_light_analysis, jjsc_gv_light_analysis, name="Selected Light", pen=pen, symbol="o", symbolSize=5, symbolBrush='b')
        else:
            A2_light = 0
            A2_light_rsquared = 0
            J0_light = 0
            J0_light_rsquared = 0

        if plot_dark_JV:

            # The voltage interval to use for the calculations will be the one stipulated previously by the user.
            voltage_interval_dark = dataframe_dark_JV.V[(dataframe_dark_JV["V"] >= jjsc_x_min) & (dataframe_dark_JV["V"] <= jjsc_x_max)]
            # The Y axis needs the short circuit current, so it's calculated here:
            current_spline_dark = interpolate.InterpolatedUnivariateSpline(dataframe_dark_JV.V, dataframe_dark_JV.I)
            current_short_circuit_dark = 0  # Per Hegedus suggestion
            # The X axis will the the voltage minus the series resistance multiplied by the short circuit current
            V_RJ_dark = dataframe_dark_JV.V - dark_series_resistance * dataframe_dark_JV.I / 1000  # The division by a thousand is to convert mA to Amps
            # The Y axis will be the sum of the current with the short circuit current, minus the shunt conductance multiplied by the voltage
            jjsc_gv_dark = (current_spline_dark(dataframe_dark_JV.V) + current_short_circuit_dark - shunt_conductance_dark * dataframe_dark_JV.V)
            # The plot will not be restricted by the voltage interval previously defined.
            # In order to plot the semi-log Y axis, the negative values need to be removed.
            # Therefore, the voltage displayed will be the voltages corresponding to positive current and voltage values
            V_RJ_dark = V_RJ_dark[(jjsc_gv_dark > 0) & (V_RJ_dark > 0)]
            V_RJ_dark_analysis = V_RJ_dark # This is done in order to keep the indexes, so that the data inside the selected voltage can be attained.
            # Same thing is done to the current.
            jjsc_gv_dark = jjsc_gv_dark[(jjsc_gv_dark > 0) & (V_RJ_dark > 0)]
            jjsc_gv_dark_analysis = jjsc_gv_dark

            # Since there are two conditions, mixing both data sets in the same condition, the re-indexing of both axis data needs to occur after the conditions, or the resulting data will be mismatched
            V_RJ_dark = V_RJ_dark.reset_index(drop=True)  # The series is re-indexed back to zero
            jjsc_gv_dark = jjsc_gv_dark.reset_index(drop=True)

            # Matching the current values to the selected voltage interval for analysis
            jjsc_gv_dark_analysis = jjsc_gv_dark_analysis.loc[voltage_interval_dark.index]
            V_RJ_dark_analysis = V_RJ_dark_analysis.loc[voltage_interval_dark.index]
            jjsc_gv_dark_analysis = jjsc_gv_dark_analysis.reset_index(drop=True)
            V_RJ_dark_analysis = V_RJ_dark_analysis.reset_index(drop=True)

            linear_regression_jjsc = stats.linregress(V_RJ_dark_analysis, np.log(jjsc_gv_dark_analysis))

            A2_dark = q / (linear_regression_jjsc.slope * k * T)
            A2_dark_rsquared = linear_regression_jjsc.rvalue ** 2

            J0_dark = np.exp(linear_regression_jjsc.intercept)
            J0_dark_rsquared = linear_regression_jjsc.rvalue ** 2

            # Plotting in black
            pen = pg.mkPen(color=(0, 0, 0), width=3)
            self.jjsc_plot.plot(V_RJ_dark, jjsc_gv_dark, name="Dark", pen=pen, symbol="o", symbolSize=5, symbolBrush='k')
            # Plotting the selected data
            pen = pg.mkPen(color=(0, 255, 0), width=3)
            self.jjsc_plot.plot(V_RJ_dark_analysis, jjsc_gv_dark_analysis, name="Selected Dark", pen=pen, symbol="o", symbolSize=5, symbolBrush='g')

        else:
            A2_dark = 0
            A2_dark_rsquared = 0
            J0_dark = 0
            J0_dark_rsquared = 0
        jjsc_results = dict({"A2_light": A2_light, "A2_light_rsquared": A2_light_rsquared, "J0_light": J0_light, "J0_light_rsquared": J0_light_rsquared,
                             "A2_dark": A2_dark, "A2_dark_rsquared": A2_dark_rsquared, "J0_dark": J0_dark, "J0_dark_rsquared": J0_dark_rsquared})

        self.result_a2light.setText(str(round(A2_light, 2)))
        self.result_a2light_rsquared.setText(str(round(A2_light_rsquared, 3)))
        self.result_j0light.setText('{:0.2e}'.format(J0_light))
        self.result_j0light_rsquared.setText(str(round(J0_light_rsquared, 3)))

        self.result_a2dark.setText(str(round(A2_dark, 2)))
        self.result_a2dark_rsquared.setText(str(round(A2_dark_rsquared, 3)))
        self.result_j0dark.setText('{:0.2e}'.format(J0_dark))
        self.result_j0dark_rsquared.setText(str(round(J0_dark_rsquared, 3)))

        return jjsc_results
    def compute_fom(self, dataframe_light_JV):
        # The aim of this function is to calculate the Figures-Of-Merit of the solar cell from the input data. Jsc, Voc, FF, efficiency, etc.
        # It used data from the LIGHT J-V CURVE!!
        voltage_sweep = dataframe_light_JV.V
        current_sweep = dataframe_light_JV.I

        current_spline = interpolate.InterpolatedUnivariateSpline(voltage_sweep, current_sweep)  # The JV data is interpolated using a univariate spline.
        # Check for more info: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html

        voltage_open_circuit = current_spline.roots()[0]  # The open circuit voltage will be the "X" (V) where the interpolation function finds the "Y" (J) at zero.
        current_short_circuit = abs(current_spline(0))  # The short circuit current will be the value of "Y" (J) calculated by the interpolation function at an "X" (V) of zero.

        # Calculate the shunt conductance:
        voltage_range_rsh = voltage_sweep.loc[voltage_sweep < 0]  # For the shunt conductance calculations, the negative voltage values will be used
        current_range_rsh = current_sweep.loc[voltage_sweep < 0]   # Same thing for the current. The current values for negative voltages were copied onto current_range_rsh
        conductance_shunt_lin_regress = stats.linregress(voltage_range_rsh,current_range_rsh)  # Calculates the linear regression of the previously set data range
        conductance_shunt = conductance_shunt_lin_regress.slope  # The shunt conductance value will be the slope of the linear regression of the negative voltages current curve
        conductance_shunt_rsquared = conductance_shunt_lin_regress.rvalue**2  # Calculates the R-Squared value of the previous linear regression fit

        # With the shunt conductance calculated, it's time to calculate the series conductance
        voltage_range_rs = voltage_sweep.loc[current_sweep > 0]  # For the series conductance, the voltage values where the current becomes positive (after Voc) will be used
        current_range_rs = current_sweep.loc[current_sweep > 0]  # For the series conductance, the positive current values will be used.
        # NOW THE CODE WILL FILTER THE CURRENT RANGE TO BE USED!!!
        # For the values to be used, it will be those with positive current value and up to 5 values after (and including) the 1st positive current value
        current_range_rs_diff = current_range_rs.diff()  # This will calculate the difference between the current values (1st derivative).
        # The aim is to detect when the current hits the current compliance of the source meter.
        current_range_rs_diff[current_range_rs_diff.index[0]] = current_range_rs_diff[current_range_rs_diff.index[0]+1]  # Because the 1st value cannot have its 1st derivative calculated, it will return a NaN.
        # The line below does the following: The values that are below 10% of the average value of the 1st derivative of the current will be dropped from the current_range_rs_diff variable
        current_range_rs_diff = current_range_rs_diff.drop(current_range_rs_diff[current_range_rs_diff.mean()*0.1 > current_range_rs_diff].index)
        # The line below does the following: It keeps only the first 5 values of the positive current data.
        current_range_rs_diff = current_range_rs_diff.drop(current_range_rs_diff[current_range_rs_diff.index > current_range_rs_diff.index[0] + 4].index)
        # Now that we have a series with the 1st derivative values higher than 10% of the average, the indexes can be used to "clean" the current_range_rs series.
        current_range_rs = current_range_rs.loc[current_range_rs_diff.index]
        voltage_range_rs = voltage_range_rs.loc[current_range_rs_diff.index]  # The voltage series needs to have the same length as the current one.
        resistance_series_lin_regress = stats.linregress(voltage_range_rs, current_range_rs)  # Linear regression of the prepared data
        resistance_series = 1000 / resistance_series_lin_regress.slope  # The series resistance will be the slope of the linear regression. Divided by a thousand to convert from mS to Ω
        resistance_series_rsquared = resistance_series_lin_regress.rvalue**2  # Calculates the r-squared, in order to determine how good the fit is.

        # Calculate the power, efficiency, fill factor
        power = voltage_sweep*current_sweep  # Since the polarity of the current was not changed, the "useful" power will be negative
        power_maximum = abs(power.min())  # As explained previously, the "useful" power is negative. Thus, the maximum power value is negative. Its converted into a positive value by the absolute function (abs)
        current_max_power = abs(current_sweep[power.idxmin()])  # From the index of the minimum of the power series, get the corresponding value of the current
        voltage_max_power = voltage_sweep[power.idxmin()]  # From the index of the minimum of the power series, get the corresponding value of the voltage
        fill_factor = (voltage_max_power * current_max_power) / (voltage_open_circuit * current_short_circuit) * 100
        power_incident = 100  # mW/cm². Later, if needed, this will be controlled by the user
        efficiency = (voltage_open_circuit * abs(current_short_circuit) * fill_factor) / power_incident  # Could have simply used the value of the maximum power, assuming 1 Sol, but this way it's more rigorous.
        # Will now create a list with the names of the columns that will be used to create a dataframe to hold and return the FOM values.
        figures_of_merit = dict({"voltage_open_circuit": voltage_open_circuit, "current_short_circuit": current_short_circuit, "fill_factor": fill_factor, "voltage_max_power": voltage_max_power,
                                 "current_max_power": current_max_power, "power_maximum": power_maximum, "efficiency": efficiency, "conductance_shunt": conductance_shunt,
                                 "conductance_shunt_rsquared": conductance_shunt_rsquared, "resistance_series": resistance_series, "resistance_series_rsquared": resistance_series_rsquared })

        # ---------------- #
        # ADDITIONAL PLOTS #
        # ---------------- #

        # In order for the user to determine the quality of the interpolation, the interpolated approximation will be plotted
        self.interpolation_plot.clear()
        self.interpolation_plot.addLegend()

        # The X axis needs to be defined. The X axis will have the same start and end as the voltage_sweep.
        interpolation_x = np.linspace(dataframe_light_JV['V'].iloc[0], dataframe_light_JV['V'].iloc[-1], len(dataframe_light_JV) * 5)  # The number of points plotted will be 5 times greater than those in the original data
        pen = pg.mkPen(color=(255, 0, 0), width=3)  # First, we will define the "red pen", to plot the original JV data.
        self.interpolation_plot.plot(dataframe_light_JV.V, dataframe_light_JV.I, name="Original Data", pen=pen, symbol="o", symbolSize=8, symbolBrush='r')  # This will plot the original JV data
        pen = pg.mkPen(color=(0, 0, 255), width=3)  # With the original data plotted, the interpolated data will be plotted in blue.
        self.interpolation_plot.plot(interpolation_x, current_spline(interpolation_x), name="Interpolated Data", pen=pen, symbol="o", symbolSize=5, symbolBrush='b')

        # Now the linear regression of the shunt conductance will be plotted
        self.linear_reg_rsh_plot.clear()
        self.linear_reg_rsh_plot.addLegend()
        pen = pg.mkPen(color=(255, 0, 0), width=3)  # First, we will define the "red pen", to plot the original JV data used for the JV calculation.
        # This line will plot the original data used for the linear regression
        self.linear_reg_rsh_plot.plot(voltage_range_rsh, current_range_rsh, name="J-V Data", pen=pen, symbol="o", symbolSize=8, symbolBrush="r")
        pen = pg.mkPen(color=(0, 0, 255), width=3)  # Setting the pen to blue to plot the linear regression.
        # This line will plot the linear regression
        self.linear_reg_rsh_plot.plot(voltage_range_rsh, conductance_shunt_lin_regress.intercept + conductance_shunt_lin_regress.slope*voltage_range_rsh, name="Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="b")

        # Time to plot the linear regression of the series resistance
        self.linear_reg_rs_plot.clear()
        self.linear_reg_rs_plot.addLegend()
        pen = pg.mkPen(color=(255, 0, 0), width=3)  # First, we will define the "red pen", to plot the original JV data used for the JV calculation.
        voltage_range_rs = voltage_range_rs.reset_index(drop=True)
        current_range_rs = current_range_rs.reset_index(drop=True)
        # This line will plot the original data used for the linear regression
        self.linear_reg_rs_plot.plot(voltage_range_rs, current_range_rs, name="JV Data", pen=pen, symbol="o", symbolSize=8, symbolBrush="r")
        pen = pg.mkPen(color=(0, 0, 255), width=3)  # Setting the pen to blue to plot the linear regression.
        # This line will plot the linear regression
        self.linear_reg_rs_plot.plot(voltage_range_rs, resistance_series_lin_regress.intercept + resistance_series_lin_regress.slope*voltage_range_rs, name="Linear Regression", pen=pen, symbol="o", symbolSize=5, symboBrush="b")

        # ---------------- #
        #    UPDATE GUI    #
        # ---------------- #

        self.result_voc.setText(str(round(voltage_open_circuit, 3)))
        self.result_jsc.setText(str(round(current_short_circuit, 2)))
        self.result_ff.setText(str(round(fill_factor, 2)))
        self.result_vmp.setText(str(round(voltage_max_power, 3)))
        self.result_jmp.setText(str(round(current_max_power, 2)))
        self.result_pmax.setText(str(round(power_maximum, 3)))
        self.result_eff.setText(str(round(efficiency, 3)))
        self.result_gsh.setText(str(round(conductance_shunt, 3)))
        self.result_gsh_squared.setText(str(round(conductance_shunt_rsquared, 3)))
        self.result_rs.setText(str(round(resistance_series, 2)))
        self.result_rs_squared.setText(str(round(resistance_series_rsquared, 3)))

        return figures_of_merit



    def save_results(self):

        # This is messy, but it will add all the results from the various dicts into a single dataframe with the names and data in two collumns
        save_data_dataframe = pd.DataFrame.from_dict(self.figures_of_merit_data, orient="index")
        save_data_dataframe = save_data_dataframe.append(pd.DataFrame.from_dict(self.djdv_results, orient="index"))
        save_data_dataframe = save_data_dataframe.append(pd.DataFrame.from_dict(self.dvdj_results, orient="index"))
        save_data_dataframe = save_data_dataframe.append(pd.DataFrame.from_dict(self.jjsc_results, orient="index"))

        save_file_name = self.lineEdit_loaded_file.text()
        save_file_name = save_file_name.split(".")[0]
        save_file_name = save_file_name + "_results.xlsx"

        # Get a directory to save the data
        save_directory = QFileDialog.getSaveFileName(self, "Choose a directory to save the resulting data", save_file_name, "Excel Files (*.xlsx)")
        # with the directory obtained, change the work path
        path_filename = save_directory[0]
        os.chdir(path_filename.rpartition("/")[0])
        save_data_dataframe.to_excel(path_filename.rsplit("/")[-1])

    def refresh_analysis(self):
        # if self.spin_ideality_min.value() == self.spin_ideality_max.value():
        #     self.spin_ideality_min.setValue(self.spin_ideality_max.value() - 0.01)
        # elif self.spin_ideality_min.value() > self.spin_ideality_max.value():
        #     self.spin_ideality_min.setValue(self.spin_ideality_max.value() - 0.01)

        try:
            dataframe_light_JV = self.data[0]  # the light JV data sits in the 1st position of the tuple
            dataframe_dark_JV = self.data[1]  # the dark JV data sits in the 2nd position of the tuple
            self.djdv_results = self.update_djdv_plot(dataframe_light_JV, dataframe_dark_JV, self.plot_light_JV, self.plot_dark_JV)
            # Feeds the Light and Dark J-V data, as well as the shunt conductance in order for the GV correction
            # Due to the sensitivity of the Hegedus method to calculate the shunt conductance, an option was added to use either the "Hegedus" conductance,
            # or the "simple" conductance, calculated directly from the light JV curve
            # This if will check if the dJdV radio button is checked. If it is, it will use the "Hegedus conductance"
            if self.radioButton_djdv_gsh.isChecked():
                conductance_variable_light = self.djdv_results["djdv_light_condutance"]
                conductance_variable_dark = self.djdv_results["djdv_dark_condutance"]
            else:
                # If it's not checked, then that must mean that the J-V curve button is pressed. If so, then the shunt conductance to be used will the LIGHT J-V one.
                conductance_variable_light = self.figures_of_merit_data["conductance_shunt"]
                conductance_variable_dark = self.figures_of_merit_data["conductance_shunt"]
            self.dvdj_results = self.update_dvdj_plot(dataframe_light_JV, dataframe_dark_JV, conductance_variable_light, conductance_variable_dark, self.plot_light_JV, self.plot_dark_JV)

            self.jjsc_results = self.update_jjsc_plot(dataframe_light_JV, dataframe_dark_JV, self.dvdj_results["dvdj_light_series_resistance"], self.dvdj_results["dvdj_dark_series_resistance"],
            conductance_variable_light, conductance_variable_dark, self.plot_light_JV, self.plot_dark_JV)
        except:
            pass
def main():


    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
