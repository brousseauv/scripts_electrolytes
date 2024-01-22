import numpy as np
from ..plotter.colorpalettes import bright
from ..plotter.msd_plotter import MsdPlotter
import netCDF4 as nc
import os

class MsdData:

    ''' Base class for calculating diffusion coefficient using MSD'''

    def __init__(self, fname, rootname):

        self.fname = fname
        self.nc_output = str('OUT/'+rootname+'.nc')
        self.output = str('OUT/'+rootname+'.dat')
        try:
            os.mkdir('OUT/')
        except OSError:
            pass

        self.my_atoms = []


    def compute_atoms_msd(self, displacements):
        '''
            Compute MSD for targeted atoms.

            Options:
                timesliced: this will average each atom's MSD at timestep t on all equivalent timeslices equal to t
                bare: no timeslice averaging is done.
        '''
        msd_atoms = np.zeros((self.nframes, self.natoms))
        if self.msd_type == 'timesliced':
            for t in range(self.nframes):
                if t%1000 == 0:
                    print('Treating time interval {}'.format(t))
                arr = displacements[t:, :, :] - displacements[:(self.nframes-t), :, :]
                msd_atoms[t, :] = np.mean(np.einsum('fad, fad -> fa', arr, arr), axis=0)

        elif self.msd_type == 'bare':
            msd_atoms = np.einsum('fad, fad -> fa', displacements, displacements)

        return msd_atoms


    def extract_diffusion_coefficient(self):

        # FIX ME: from here, classes should already have a time and msd property
        # Also, one could decide to elimitate some initial and final part of the trajectory when computing the fit

        # Assume units of angstrom^2/ps
        self.slope = np.polyfit(self.time[self.discard_init_steps:], self.msd[self.discard_init_steps:], 1)

        # Assume 3D diffusion, for which the slope of MSD vs t is 6D
        self.diffusion = 1E-4*self.slope[0]/6
        # FIX ME: add standard deviation of diffusion coefficient fit

        return self.diffusion


    def extract_msd_errors(self):

        if self.msd_atoms is not None:
            return np.std(self.msd_atoms, axis=1)
        else:
            return None

    
    def plot_diffusion_coefficient(self, verbose=True, plot_all_atoms=False, **kwargs):
        
        myplot = MsdPlotter(**kwargs) 
        myplot.set_line2d_params(**kwargs)

        y = self.slope[0]*self.time + self.slope[1]

        if plot_all_atoms:
            for a in range(self.natoms):
                myplot.ax.plot(self.time, self.msd_atoms[:, a], linewidth=0.5*myplot.linewidth, linestyle='solid', alpha=0.6)
            myplot.ax.plot(self.time, self.msd, color='black', linewidth=1.5*myplot.linewidth, linestyle='solid')
            myplot.ax.plot(self.time, y, color=bright['red'], linewidth=1.5*myplot.linewidth, linestyle='dashed')
        else:
            myplot.ax.plot(self.time, self.msd, color=bright['blue'], linewidth=myplot.linewidth, linestyle='solid')
            myplot.ax.plot(self.time, y, color=bright['red'], linewidth=myplot.linewidth, linestyle='dashed')

            # FIX ME: this is a little crude, as msd-msd_std can be negative
            # Plotting the MSD-std makes only sense if we do not plot individual trajectories :P
            if self.msd_std is not None and self.plot_errors:
                myplot.ax.fill_between(self.time, self.msd-self.msd_std, self.msd+self.msd_std, color='gray', zorder=-1, alpha=0.3)

        if verbose:
            myplot.ax.text(0.10, 0.90, r'D={:.3e} cm$^2$/s'.format(self.diffusion), fontsize=myplot.labelsize+2, transform=myplot.ax.transAxes)
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

            dts.createDimension('number_of_frames', self.nframes)
            dts.createDimension('number_of_diffusing_atoms', self.natoms)
            dts.createDimension('one', 1)

            dts.setncattr('diffusing_atom_type', self.atom_type)
            dts.setncattr('msd_type', self.msd_type)
            dts.setncattr('data_source', self.data_source)

            data = dts.createVariable(
                    'temperature', 'd', ('one'))
            data.units = 'Kelvin'
            data[:] = self.temperature

            data = dts.createVariable(
                    'diffusion_coefficient', 'd', ('one'))
            data.units = 'cm^2/s'
            data[:] = self.diffusion

            data = dts.createVariable(
                    'time', 'd', ('number_of_frames'))
            data.units = 'picosecond'
            data[:] = self.time

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
            try:
                data[:] = self.thermo_step * self.discard_init_steps
            except AttributeError:
                data[:] = self.timestep * self.discard_init_steps

            data = dts.createVariable(
                    'mean_squared_displacement', 'd',
                    ('number_of_frames'))
            data.units = 'Angstrom^2'
            data[:] = self.msd

            data = dts.createVariable(
                    'mean_squared_displacement_individual_atoms', 'd',
                    ('number_of_frames', 'number_of_diffusing_atoms'))
            data.units = 'Angstrom^2'
            if self.msd_atoms is not None:
                data[:, :] = self.msd_atoms

            data = dts.createVariable(
                    'standard_deviation_mean_squared_displacement', 'd',
                    ('number_of_frames'))
            data.units = 'Angstrom^2'
            if self.msd_std is not None:
                data[:] = self.msd_std


    def write_output(self):

        with open(self.output, 'w') as f:

            f.write('Data source: {}\n'.format(self.data_source))
            f.write('Temperature: {:.0f}K\n'.format(self.temperature))
            f.write('Total runtime: {:.5f} ps\n'.format(self.time[-1]))
            f.write('Timestep: {:.5f} ps\n'.format(self.timestep))
            try:
                f.write('Initial {} ps has been discarded\n'.format(self.discard_init_steps*self.thermo_step))
            except AttributeError:
                f.write('Initial {} ps has been discarded\n'.format(self.discard_init_steps*self.timestep))

            f.write('Diffusing atoms type: {}\n'.format(self.atom_type))
            f.write('MSD type: {}\n'.format(self.msd_type))
            f.write('Diffusion coefficient: {:.5e} cm^2/s\n'.format(self.diffusion))
        f.close()
