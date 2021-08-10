from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox  # necessary import of the file dialog window
from pyqtgraph import PlotWidget
import pyqtgraph as pg

from JV_Analyser_GUI import Ui_MainWindow  # this import will load the .py file with the GUI

import pandas as pd  # Pandas dataframes are a super useful tool to handle and process data.
import os, sys  # Necessary modules to allow the python code to interact with Windows.

import random
import numpy as np

##########################################
############ CONFIGURATIONS ##############
##########################################
# ---------------------------------------#
# ---------------- PLOT -----------------#
# ---------------------------------------#

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
        self.jv_plot.setLabel("bottom","<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.jv_plot.showGrid(x=True, y=True)


        self.djdv_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">dJ/dV (mS cm <sup>-2</sup>)</span>")
        self.djdv_plot.setLabel("bottom","<span style=\"color:black;font-size:18px\">V (V)</span>")
        self.djdv_plot.showGrid(x=True, y=True)

        self.dvdj_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">dV/dJ (Ω cm <sup>2</sup>)</span>")
        self.dvdj_plot.setLabel("bottom", "<span style=\"color:black;font-size:18px\">(J+J<sub>SC</sub>)<sup>-1</sup> (mA<sup>-1</sup> cm<sup>2</sup>)</span>")
        self.dvdj_plot.showGrid(x=True, y=True)

        self.jjsc_plot.setLabel("left", "<span style=\"color:black;font-size:18px\">J+J<sub>SC</sub>-GV (mA/cm<sup>2</sup>)</span>")
        self.jjsc_plot.setLabel("bottom","<span style=\"color:black;font-size:18px\">V-RJ (V)</span>")
        self.jjsc_plot.showGrid(x=True, y=True)
        # ---------------------------------------#

        ################## SIGNALS ##################
        # A signal is the connection between what happens in the GUI and the python code.

        self.actionLoad_Data.triggered.connect(self.execute_processing)  # when the button "Load Data" is pressed, it will run the execute_processing function


    def execute_processing(self):
        # As far as I know, Qt slots cannot return any data. Thus, I created this function to run all the necessary code to process the data when a file is loaded.
        # The data processing starts by loading the .csv data of both light and dark curves

        # Python will run the load_data() function. If the user quits the file dialog without selecting a file, it would crash the program, if it was not for the
        # error handling "try" and "except" functions. If the user selects a file, the code inside "try:" will run. If not, an error will be raised, and the "except:"
        # code will just ignore it and do nothing.


        try:
            data = self.load_data()
             # the load_data() function will return the light JV and dark JV data in dataframes, with said dataframes inside a tuple
            dataframe_light_JV = data[0]  # the light JV data sits in the 1st position of the tuple
            dataframe_dark_JV = data[1]  # the dark JV data sits in the 2nd position of the tuple

        except:
            pass
        plot_light_JV = True
        plot_dark_JV = True
        self.update_jv_plot(dataframe_light_JV, dataframe_dark_JV, plot_light_JV, plot_dark_JV)

    def load_data(self):
        # the line below will open a file dialog window, and ask for the user to select a .csv file. The third part of the .getOpenFileName()
        # with "CSV Files (*.csv)" will apply a filter so that only .csv files show up. This will return a tuple with the position
        # zero containing the path and file name, and the second containing the data filter mentioned previously. Since the program
        # only requires the path and filename the ", _" is used to discard the second tuple value.
        if self.checkBox_analyse_light.isChecked():
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
            if self.checkBox_analyse_dark.isChecked():
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
        return light_data, dark_data

    def parse_file(self, file_name):
        # This function will take the file from the load_data function and load/process the csv data so that it can be used later
        data = pd.read_csv(file_name, header=15, sep=';')
        return data

    def update_jv_plot(self, dataframe_light_JV, dataframe_dark_JV, plot_light_JV, plot_dark_JV):
        # the purpose of this function is to update the JV plot. Basing off the Hegedus' paper (10.1002/pip.518), the JV will be the top left plot.
        # The light JV data will be plotted in red. The dark JV data will be plotted in dark, if the user wants to plot it.
        self.jv_plot.setBackground('w')
        if plot_light_JV:  # If the user selected to process Light JV data, plot the light JV data
            pen = pg.mkPen(color=(255, 0, 0), width=3)  # To change the color of the plot, you need assign a color to the "pen". In this case (RGB) 255, 0, 0 is RED.
            # The width is also changed to 3 px wide
            self.jv_plot.plot(dataframe_light_JV.V, dataframe_light_JV.I, pen=pen)
        if plot_dark_JV: # If the user selected to process Dark JV data, plot the Dark JV data
            pen = pg.mkPen(color=(0, 0, 0), width=3)  # Changing the pen color back to black, in order to plot the Dark JV data.
            self.jv_plot.plot(dataframe_dark_JV.V, dataframe_dark_JV.I, pen=pen)



def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
