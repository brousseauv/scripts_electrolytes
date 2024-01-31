import numpy as np
import netCDF4 as nc
import os
from pymatgen.io.abinit.netcdf import NetcdfReader
from ..plotter.msd_plotter import MsdPlotter, DCPlotter
from ..plotter.colorpalettes import bright
from ..utils.functions import sort_consecutive_groups
from matplotlib.pyplot import subplots


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
        self.temp = reader.read_value('temperature')[0]
        self.coeff = reader.read_value('diffusion_coefficient')[0]
        self.msd_type = reader.rootgrp.getncattr('msd_type')

        if mode == 'msd' or mode == 'msd_atoms':
            self.atom_type = reader.rootgrp.getncattr('diffusing_atom_type')

            self.time = reader.read_value('time')
            self.timestep = reader.read_value('timestep')[0]
            self.msd = reader.read_value('mean_squared_displacement')

            if mode == 'msd_atoms':
                self.natoms = reader.read_dimvalue('number_of_diffusing_atoms')
                self.msd_atoms = reader.read_value('mean_squared_displacement_individual_atoms')

    def extract_atomic_jumps(self, threshold=4.0, plot=False, dist2=2.0, **kwargs):

        ''' Locate possible atomic jumps in individual atoms MSD(t),
            i.e. values of MSD(t) that are larger than the threshold value.

            threshold: minimal squared displacement for detection (in Angstrom^2)
            Default: 4.0 angstrom^2

            plot: Plot MSD(t) for the detected jumps for selected atoms and add visual guides as to
                  where the jump occurs
            Default: False

            dist2: Minimal difference between the current MSD and the average MSD over the last 200 steps so that the
                   jump is detected.
                   Default: 2.0 angstrom^2

            FIX ME: I will also need to add a diff_threshold so that I can compare the jump with the average MSD
            over the last N steps (to locate single jumps, not all timesteps where MSD
            is larger than threshold)
        '''

        jumping_atoms_list = []

        for a in range(self.natoms):
            mymsd = self.msd_atoms[:, a]
            if any(mymsd > threshold):
                # or a condition with a mean?
#                where = [idx for idx in range(len(mymsd)) if mymsd[idx]>threshold
#                        and np.abs(mymsd[idx]-np.mean(mymsd[idx-200:idx]))>threshold/2]
                where = [idx for idx in range(len(mymsd)) if mymsd[idx]>threshold
                         and np.abs(mymsd[idx]-np.mean(mymsd[idx-200:idx]))>dist2]

#                        and np.abs(mymsd[idx]-mymsd[idx-100])>threshold/2
                # Keep only the first of each block of consecutive timesteps
                # TO BE TESTED on multijump trajectories!!!!!
                # FIX ME: for longer trajectories, the printed output is indigest
                start, end = sort_consecutive_groups(where)

                print('\nFound jumps in atom {}, between timesteps {}-{}'.format(a, start, end))

                jumping_atoms_list.append(a)
                if plot:
                    self.plot_msd_atom(a, threshold, start, **kwargs)
        print('Found jumps in atoms:{}'.format(jumping_atoms_list))

    def plot_msd_atom(self, index, href, vref, **kwargs):

        myplot = MsdPlotter(figsize=(12, 4), **kwargs)
        myplot.set_line2d_params(defname='msd_atom{}.png'.format(index), title='Atom {}'.format(index),
                                 xlim=(0, self.time[-1]), **kwargs)
        myplot.ax.plot(self.time, self.msd_atoms[:, index], color='k')
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

        ''' Computes the diffusion coefficient as an average on non-overlapping slices
            of size delta t in [1, time[-1]] '''

        self.read_data(mode='msd')
        nsteps = len(self.time)

        self.deltat = []
        self.diffusion_from_slices = []
        self.diffusion_std = []

        print('nsteps:', nsteps, self.time[-1])
        cut = 0
        for nslice in np.arange(1, nsteps+1):
            newcut = int(np.ceil((nsteps-1)/nslice))
            if newcut != cut:
                cut = newcut
                self.deltat.append(cut*self.timestep)

                n = 0
                diff = []
                while n*cut+1<nsteps:
                    # For now they have overlapping start/endpoints
                    coeff = self.compute_diffusion_coefficient(self.time[n*cut:(n+1)*cut+1], self.msd[n*cut:(n+1)*cut+1])
                    diff.append(coeff)
                    n += 1

                self.diffusion_from_slices.append(np.mean(diff))
                self.diffusion_std.append(np.std(diff))

                print('for nslice {}, ncut={}, D={:.3e}cm^2/s, nintervals={}'.format(nslice, cut, np.mean(diff), len(diff)))
                if cut<10: 
                    print(np.asarray(diff)[::500])

        if plot:
            self.plot_diffusion_from_slices(**kwargs)


    def compute_diffusion_coefficient(self, x, y):

        # Assume units of angstrom^2/ps
        slope, cov = np.polyfit(x, y, 1, cov=True)
        # Assume 3D diffusion, for which the slope of MSD vs t is 6D
        coefficient = 1E-4*slope[0]/6
        std = 1E-4*np.sqrt(np.diag(cov)[0])/6

        return coefficient, slope, std


    def plot_diffusion_from_slices(self, **kwargs):

        fig, ax = subplots(1, 2, figsize=(12, 6))
        for j in range(2):
            myplot = DCPlotter(ax=ax[j], **kwargs) 
            myplot.set_line2d_params(defname='diffusion_deltat.png', **kwargs)

            myplot.ax.errorbar(self.deltat, self.diffusion_from_slices, yerr=self.diffusion_std, color=bright['blue'], alpha=0.4, zorder=1)
            myplot.ax.plot(self.deltat, self.diffusion_from_slices, color='black', linewidth=3, zorder=2)

            if j == 0:
                myplot.ax.set_yscale('log')
                myplot.ax.set_xscale('log')
    #        if verbose:
    #            myplot.ax.text(0.10, 0.90, r'D={:.3e} cm$^2$/s'.format(self.diffusion), fontsize=myplot.labelsize+2, transform=myplot.ax.transAxes)
            myplot.set_labels()
            myplot.ax.axhline(self.diffusion_from_slices[0], color='r', linestyle='dashed')

            mean = np.mean(self.diffusion_from_slices)
            myplot.ylim=(mean*1E-1, mean*5)
            myplot.xlim=(self.deltat[-1][0], self.deltat[0][0])
            myplot.set_limits()

            try:
                myplot.add_title()
            except:
                pass

        myplot.save_figure()
        myplot.show_figure()


    def recompute_diffusion_coefficient(self, discard_init_steps=0, rootname=None, plot=False, fill=True, verbose=True, **kwargs):
        ''' Recompute diffusion coefficient with a new value for discard_init_steps, i.e. 
            change the portion of the MD run used for slope calculation. 

            Input: 
            discard_init_steps: number of timesteps to discard from slope evaluation
                                default=0
            rootname: rootname for new output files (.dat and .nc formats)

            plot: should the new MSD(t) and diffusion slope be plotted
                  default: False

            fill: should the discarded portion of the MSD be shaded
                  default: True

            verbose: should the new diffusion coefficient be printed to the screen
                     default: True
        '''

        if not isinstance(discard_init_steps, int) and not isinstance(discard_init_steps, np.integer):
            raise TypeError('discard_init_steps should be an integer, but I got {} which is a {}'.format(discard_init_steps, type(discard_init_steps)))
        self.discard_init_steps = discard_init_steps

        if not rootname:
            rootname = os.path.splitext(os.path.basename(self.fname))[0] + f'_discard{discard_init_steps}'
        rootdir = os.path.dirname(self.fname)

        self.nc_output = os.path.join(rootdir, str(rootname+'.nc'))
        self.output = os.path.join(rootdir, str(rootname+'.dat'))

        self.read_data(mode='msd')

        self.coeff, self.slope, self.coeff_std = self.compute_diffusion_coefficient(self.time[discard_init_steps:], self.msd[discard_init_steps:])

        if verbose:
            print(f'D={self.coeff:.3e}+-{self.coeff_std:.3e} cm^2/s')
        self.write_data()        

        if plot:
            self.plot_diffusion_coefficient(fill, **kwargs)


    def plot_diffusion_coefficient(self, fill, **kwargs):

        myplot = MsdPlotter(**kwargs)
        myplot.set_line2d_params(defname='diffusion.png', **kwargs)

        y = self.slope[0]*self.time + self.slope[1]

        if myplot.linecolor:
            color = myplot.linecolor
        else:
            color = bright['red']

        myplot.ax.plot(self.time, self.msd, color=color, linewidth=1.0, linestyle='solid')
        myplot.ax.plot(self.time, y, color=color, linewidth=myplot.linewidth, linestyle='dashed')

        if self.discard_init_steps != 0:
            myplot.ax.axvline(self.time[self.discard_init_steps], color='black', linestyle='dashed', linewidth=0.5)
            ylims = myplot.ax.get_ylim()
            xlims = myplot.ax.get_xlim()
            if fill:
                myplot.ax.fill_between(np.linspace(0, self.time[self.discard_init_steps], 10), ylims[0], ylims[1], color='gray', alpha=0.4, zorder=-1)

        myplot.set_labels()
        try:
            myplot.add_title()
        except:
            pass

        myplot.save_figure()
        myplot.show_figure()

    def write_data(self):

        self.write_output()
        self.write_netcdf()


    def write_netcdf(self):

        with nc.Dataset(self.nc_output, 'w') as dts:

            dts.createDimension('one', 1)
            dts.createDimension('string100', 100)

            dts.setncattr('diffusing_atom_type', self.atom_type)
            dts.setncattr('msd_type', self.msd_type)
            dts.setncattr('msd_data_source_file', os.path.abspath(self.fname))

            data = dts.createVariable(
                    'temperature', 'd', ('one'))
            data.units = 'Kelvin'
            data[:] = self.temp

            data = dts.createVariable(
                    'diffusion_coefficient', 'd', ('one'))
            data.units = 'cm^2/s'
            data[:] = self.coeff

            data = dts.createVariable(
                    'total_runtime', 'd', ('one'))
            data.units = 'picosecond'
            data[:] = self.time[-1]

            data = dts.createVariable(
                    'timestep', 'd', ('one'))
            data.units = 'picosecond'
            data[:] = self.timestep

            data = dts.createVariable(
                    'discard_initial_timesteps', 'd',  ('one'))
            data.units = 'picosecond'
            data[:] = self.timestep * self.discard_init_steps

    def write_output(self):

        with open(self.output, 'w') as f:

            f.write('MSD data source file: {}\n'.format(self.fname))
            f.write('Temperature: {:.0f}K\n'.format(self.temp))
            f.write('Total runtime: {:.5f} ps\n'.format(self.time[-1]))
            f.write('Timestep: {:.5f} ps\n'.format(self.timestep))
            f.write('Initial {} ps has been discarded\n'.format(self.discard_init_steps*self.timestep))

            f.write('Diffusing atoms type: {}\n'.format(self.atom_type))
            f.write('MSD type: {}\n'.format(self.msd_type))
            f.write('Diffusion coefficient: {:.5e} cm^2/s\n'.format(self.coeff))
        f.close()
