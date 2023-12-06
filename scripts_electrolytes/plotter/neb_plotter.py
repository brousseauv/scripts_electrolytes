import matplotlib.pyplot as plt
from .plotter import Plotter
from matplotlib.ticker import StrMethodFormatter

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')

class NebPathPlotter(Plotter):

    def __init__(self, **kwargs):

        super(NebPathPlotter, self).__init__(**kwargs)


    def set_labels(self):                                                                                                                                                                 
        self.ax.set_xlabel(r'Normalized reaction coordinate', fontsize=self.labelsize)
        self.ax.set_ylabel(r'Energy (eV)', fontsize=self.labelsize)

        self.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
        self.ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))

        self.ax.tick_params(axis='both', labelsize=self.labelsize-2)
