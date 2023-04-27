import numpy as np
from colorpalettes import bright
from ..plotter import Plotter
import tqdm

class MsdData:

    ''' Base class for calculating diffusion coefficient using MSD'''

    def __init__(self, fname):

        self.fname = fname


    def compute_atoms_msd(self, displacements, msd_type):
        # Compute MSD for target atoms
        msd_atoms = np.zeros((self.nframes, self.natoms))

        if msd_type == 'timesliced':
            # Average all atom MSD on equivalent timeslices
            for t in tqdm.tqdm(range(self.nframes)):
                arr = displacements[t:, :, :] - displacements[:(self.nframes-t), :, :]
                msd_atoms[t, :] = np.mean(np.einsum('fad, fad -> fa', arr, arr), axis=0)

        elif msd_type == 'bare':
            # No timeslice averaging
            msd_atoms = np.einsum('fad, fad -> fa', displacements, displacements)

        return msd_atoms



    def extract_diffusion_coefficient(self):

        # from here, classes should already have a time and msd property
        # Also, one could decide to elimitate some initial and final part of the trajectory when computing the fit

        # Assume units of angstrom^2/ps
        self.slope = np.polyfit(self.time, self.msd, 1)
        # Assume 3D diffusion, for which the slope of MSD vs t is 6D
        self.diffusion = 1E-4*self.slope[0]/6

        return self.diffusion


    def extract_msd_errors(self):

        if self.msd_atoms is not None:
            return np.std(self.msd_atoms, axis=1)
        else:
            return None

    
    def plot_diffusion_coefficient(self, **kwargs):
        
        myplot = Plotter(**kwargs) 
        myplot.set_line2d_params(**kwargs)

        y = self.slope[0]*self.time + self.slope[1]
        myplot.ax.plot(self.time, self.msd, color=bright['blue'], linewidth=myplot.linewidth, linestyle='solid')
        myplot.ax.plot(self.time, y, color=bright['red'], linewidth=myplot.linewidth, linestyle='dashed')

        # if we have std on atoms at each time step, add it
        if self.msd_std is not None:
            myplot.ax.fill_between(self.time, self.msd-self.msd_std, self.msd+self.msd_std, color='gray', zorder=-1, alpha=0.3)

        myplot.ax.text(0.10, 0.90, r'D={:.3e} cm$^2$/s'.format(self.diffusion), fontsize=myplot.labelsize+2, transform=myplot.ax.transAxes)
        myplot.set_labels()
        try:
            myplot.set_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()
