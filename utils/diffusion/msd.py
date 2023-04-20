import numpy as np
from ..interfaces.lammps_interface import extract_thermo, create_thermo_dataframe, read_msd
from colorpalettes import bright
from ..plotter import Plotter

class MsdData:

    ''' Base class for calculating diffusion coefficient using MSD'''

    def __init__(self, fname):

        self.fname = fname

    def extract_diffusion_coefficient(self):

        # from here, classes should already have a time and msd property

        # Assume units of anstrom^2/ps
        self.slope = np.polyfit(self.time, self.msd, 1)
        self.diffusion = 1E-4*self.slope[0]/6

        return self.diffusion


    
    def plot_diffusion_coefficient(self, **kwargs):
        # optional plot of average msd vs time and fit
        
        myplot = Plotter(**kwargs) 
        myplot.set_line2d_params(**kwargs)

        y = self.slope[0]*self.time + self.slope[1]
        myplot.ax.plot(self.time, self.msd, color=bright['blue'], linewidth=myplot.linewidth, linestyle='solid')
        myplot.ax.plot(self.time, y, color=bright['red'], linewidth=myplot.linewidth, linestyle='dashed')

        myplot.set_labels()
        try:
            myplot.set_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()
