import numpy as np
from .msd_output import MsdOutput
from ..utils.constants import boltzmann_evK
from ..plotter.ea_plotter import EaPlotter
from ..plotter.colorpalettes import bright

class ActivationEnergyData:

    def __init__(self, flist):

        if flist is not None:
            self.flist = flist
        else:
            raise Exception('flist should not be empty. Please provide a list of MsdData netCDF output files')

        self.extract_data()


    def extract_data(self):

        # define a flist of eventual MsdOutput objects, read them and add them to the relevant variables
        nfile = len(self.flist)
        self.temperature = np.zeros((nfile))
        self.diffusion_coefficient = np.zeros((nfile))

        for i, myfile in enumerate(self.flist):
            data = MsdOutput(fname=myfile)
            data.read_data()
            self.temperature[i] = data.temp
            self.diffusion_coefficient[i] = data.coeff

        self.inverse_temperature = 1./self.temperature

    def compute_activation_energy(self, plot=False, plot_verbose=True, **kwargs):

        # Linear fit of ln(D(T)) vs 1/T. 
        fit, cov = np.polyfit(self.inverse_temperature, np.log(self.diffusion_coefficient), 1, cov=True)
        self.activation_energy = -fit[0]*boltzmann_evK
        self.d0 = np.exp(fit[1])
        self.ea_std = np.sqrt(np.diag(cov))[0]*boltzmann_evK
        print('Activation energy = {:.3f}+-{:.3f} eV'.format(self.activation_energy, self.ea_std))
        
        if plot:
            self.plot_activation_energy(defname='activation_energy.png', verbose=plot_verbose, **kwargs)


    def plot_activation_energy(self, verbose=True, **kwargs):
        # this should be somewhere else!
        myplot = EaPlotter(**kwargs)
        myplot.set_line2d_params(**kwargs)

        myplot.ax.semilogy(1000*self.inverse_temperature, self.diffusion_coefficient, marker='o', color=bright['blue'], linestyle='None')
        y = self.d0*np.exp(-self.activation_energy/boltzmann_evK*self.inverse_temperature)
        myplot.ax.semilogy(1000*self.inverse_temperature, y, color=bright['red'])

        if verbose:
            myplot.ax.text(0.55, 0.90, r'E$_A$={:.3f}$\pm${:.3f} eV'.format(self.activation_energy, self.ea_std), fontsize=myplot.labelsize+2, transform=myplot.ax.transAxes)
        myplot.set_labels()
        try:
            myplot.add_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()
