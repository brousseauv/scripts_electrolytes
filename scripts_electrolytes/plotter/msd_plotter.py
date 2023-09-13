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


class DCPlotter(Plotter):

    def __init__(self, **kwargs):

        super(DCPlotter, self).__init__(**kwargs)

    
    def set_labels(self):

        self.ax.set_xlabel(r'$\Delta$ t (ps)', fontsize=self.labelsize)
        self.ax.set_ylabel(r'Diffusion coefficient (cm$^2$/s)', fontsize=self.labelsize)

    def set_limits(self):
        
        print('current ylim:', self.ylim)
        if not self.ylim:
            self.ylim = self.ax.get_ylim()
        if not self.xlim:
            self.xlim = (0, self.ax.get_xlim()[1])
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

