# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

# PyQt requires that a module is created to create the link between the UI widget and the matplotlib plot.
# This module does that function. It basically adds the matplotlib backend and necessary code for the plots to appear.
# More information in the book ISBN 978-1-847197-90-0 "Matplotlib for Python Developers". Its for PyQt4 but a lot of stuff still applies to PyQt5/6


class MplCanvas(FigureCanvas):

    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)  # 111 means that it will add a grid of 1x1 plots, and will add 1 plot to said grid
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.updateGeometry(self)


class MplWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
        # QWidget.__init__(self, parent)
        # self.canvas = FigureCanvas(Figure())
        # vertical_layout = QVBoxLayout()
        # vertical_layout.addWidget(self.canvas)
        # self.canvas.axes = self.canvas.figure.add_subplot(111)  # 111 means that it will add a grid of 1x1 plots, and will add 1 plot to said grid
        # self.canvas.figure.add_subplot()
        # self.setLayout(vertical_layout)
