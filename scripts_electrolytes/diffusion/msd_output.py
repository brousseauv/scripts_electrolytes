import numpy as np
import netCDF4 as nc
import os
from pymatgen.io.abinit.netcdf import NetcdfReader
from ..plotter.msd_plotter import MsdPlotter, DCPlotter
from ..plotter.colorpalettes import bright
from ..utils.functions import sort_consecutive_groups

class MsdOutput:

    ''' Base class to postprocess netCDF output files from MsdData class '''

    def __init__(self, fname=None):

        # what do I need to do here?
        try:
            os.path.exists(fname)
        except OSError:
            raise OSError('File {} does not exist!'.format(fname))
        else:
            if not fname.endswith('.nc'):
                raise Exception('files should be in netCDF format')
            else:
                self.fname = fname

    
    def read_data(self, mode='diffusion'):

        ''' Mode: defines the amount of data that will be read. 
                "diffusion": read only temperature, diffusion coefficient and msd_type
                "msd": read the above, plus time, timestep and averaged MSD(t)
                "msd_atoms": reads all the above, plus individual atoms msd
                Default: "diffusion"
        '''

        reader = NetcdfReader(self.fname)
        self.temp = reader.read_value('temperature')
        self.coeff = reader.read_value('diffusion_coefficient')
        self.msd_type = reader.rootgrp.getncattr('msd_type')

        if mode == 'msd' or mode == 'msd_atoms':
            self.atom_type = reader.rootgrp.getncattr('diffusing_atom_type')

            self.time = reader.read_value('time')
            self.timestep = reader.read_value('timestep')
            self.msd = reader.read_value('mean_squared_displacement')

            if mode == 'msd_atoms':
                self.natoms = reader.read_dimvalue('number_of_diffusing_atoms')
                self.msd_atoms = reader.read_value('mean_squared_displacement_individual_atoms')


    def extract_atomic_jumps(self, threshold=4.0, plot=False, **kwargs):

        ''' Locate possible atomic jumps in individual atoms MSD(t), 
            i.e. values of MSD(t) that are larger than the threshold value. 

            threshold: minimal squared displacement for detection (in Angstrom^2)
            Default: 4.0 \AA^2

            plot: Plot MSD(t) for the detected jumps for selected atoms and add visual guides as to where the jump occurs
            Default: False

            I will also need to add a diff_threshold so that I can compare the jump with the average MSD over the last N steps (to locate single jumps, not all timesteps where MSD
            is larger than threshold)
        '''

        jumping_atoms_list = []

        for a in range(self.natoms):
            mymsd = self.msd_atoms[:, a]
            if any(mymsd > threshold):
                # or a condition with a mean?
                where = [idx for idx in range(len(mymsd)) if mymsd[idx]>threshold 
                        and np.abs(mymsd[idx]-np.mean(mymsd[idx-200:idx]))>threshold/2]

#                        and np.abs(mymsd[idx]-mymsd[idx-100])>threshold/2 
                # Keep only the first of each block of consecutive timesteps
                # TO BE TESTED on multijump trajectories!!!!!
                start, end = sort_consecutive_groups(where)

                #print('\nFound jumps in atom {}, at timestep(s) {}'.format(a, start))
                print('\nFound jumps in atom {}, between timesteps {}-{}'.format(a, start, end))

                jumping_atoms_list.append(a)
                if plot:
                    self.plot_msd_atom(a, threshold, start, **kwargs)
        print('Found jumps in atoms:{}'.format(jumping_atoms_list))
        print('If jumps were found, remember that the timesteps are shifted (smaller) by discard_init_steps vs the dump file.')


    def plot_msd_atom(self, index, href, vref, **kwargs):

        myplot = MsdPlotter(figsize=(12,4), **kwargs)
        myplot.set_line2d_params(defname='msd_atom{}.png'.format(index),title='Atom {}'.format(index), xlim=(0, self.time[-1]), **kwargs)
        myplot.ax.plot(self.time, self.msd_atoms[:,index], color='k')
        myplot.ax.axhline(href, linestyle='dashed', color='blue')
        if len(vref) != 0:
            for val in vref:
                myplot.ax.axvline(val*self.timestep, linestyle='dashed', color='red')


        myplot.fig.subplots_adjust(bottom=0.15)
        myplot.set_labels()
        myplot.set_limits()
        try:
            myplot.add_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()


    def diffusion_coefficient_from_slices(self, plot=True, **kwargs):

        ''' Computes the diffusion coefficient as an average on non-overlapping slices of size delta t \in [1, time[-1]] '''

        self.read_data(mode='msd')
        nsteps = len(self.time)

        self.deltat = []
        self.diffusion_from_slices = []
        self.diffusion_std = []

        print('nsteps:',nsteps, self.time[-1])
        cut = 0
        for nslice in np.arange(1, nsteps+1):
            newcut = int(np.ceil((nsteps-1)/nslice))
            if newcut != cut:
                cut = newcut
                self.deltat.append(cut*self.timestep)

                print('for nslice {}, ncut={}'.format(nslice, cut))
                n = 0
                diff = []
                while n*cut+1<nsteps:
                    # For now they have overlapping start/endpoints
                    coeff = self.compute_diffusion_coefficient(self.time[n*cut:(n+1)*cut+1], self.msd[n*cut:(n+1)*cut+1])
                    diff.append(coeff)
                    n += 1

                self.diffusion_from_slices.append(np.mean(diff))
                self.diffusion_std.append(np.std(diff))

        if plot:
            self.plot_diffusion_from_slices(**kwargs)


    def compute_diffusion_coefficient(self, x, y):

        # Assume units of angstrom^2/ps
        slope = np.polyfit(x, y, 1)
        # Assume 3D diffusion, for which the slope of MSD vs t is 6D
        coefficient = 1E-4*slope[0]/6

        return coefficient


    def plot_diffusion_from_slices(self, **kwargs):

        myplot = DCPlotter(**kwargs) 
        myplot.set_line2d_params(defname='diffusion_deltat.png', **kwargs)

        myplot.ax.errorbar(self.deltat, self.diffusion_from_slices, yerr=self.diffusion_std, color=bright['blue'])
        myplot.ax.plot(self.deltat, self.diffusion_from_slices, color='black', linewidth=1.5)

        myplot.ax.semilogy()
#        if verbose:
#            myplot.ax.text(0.10, 0.90, r'D={:.3e} cm$^2$/s'.format(self.diffusion), fontsize=myplot.labelsize+2, transform=myplot.ax.transAxes)
        myplot.set_labels()
        mean = np.mean(self.diffusion_from_slices)
        myplot.ylim=(mean*1E-2, mean*3)
        myplot.set_limits()

        try:
            myplot.add_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()


