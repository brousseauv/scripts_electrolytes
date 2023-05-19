import matplotlib.pyplot as plt
from .plotter import Plotter

plt.rc('text', usetex=True)

class MsdPlotter(Plotter):

    def __init__(self, **kwargs):

        super(MsdPlotter, self).__init__(**kwargs)


    def set_labels(self):                                                                                                                                                                 
        self.ax.set_xlabel(r'time (ps)', fontsize=self.labelsize)
        self.ax.set_ylabel(r'MSD (\AA$^2$)', fontsize=self.labelsize)
