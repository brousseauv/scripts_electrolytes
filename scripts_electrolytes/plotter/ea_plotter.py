import matplotlib.pyplot as plt
from .plotter import Plotter

plt.rc('text', usetex=True)

class EaPlotter(Plotter):

    def __init__(self, **kwargs):

        super(EaPlotter, self).__init__(**kwargs)


    def set_labels(self):                                                                                                                                                                 
        self.ax.set_xlabel(r'1000/T (K$^{-1}$)', fontsize=self.labelsize)
        self.ax.set_ylabel(r'D (cm$^2$/s)', fontsize=self.labelsize)
