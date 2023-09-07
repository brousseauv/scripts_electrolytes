import matplotlib.pyplot as plt
from .plotter import Plotter

plt.rc('text', usetex=True)

class MsdPlotter(Plotter):

    def __init__(self, **kwargs):

        super(MsdPlotter, self).__init__(**kwargs)


    def set_labels(self):                                                                                                                                                                 
        self.ax.set_xlabel(r'time (ps)', fontsize=self.labelsize)
        self.ax.set_ylabel(r'MSD (\AA$^2$)', fontsize=self.labelsize)

    def set_limits(self):
        
        if not self.ylim:
            self.ylim = (0, self.ax.get_ylim()[1])
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)


