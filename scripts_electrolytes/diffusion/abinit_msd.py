from .msd import MsdData
from abipy.dynamics.hist import HistFile
from ..utils.constants import atomic_timeunit, bohr_to_ang
import numpy as np
import logging
import os

class HistMsdData(MsdData):

    def __init__(self, fname, rootname='MsdData'):

        '''
            Input:
                fname: name of the netCDF HIST.nc file containing atomic positions information

                rootname: rootname for the .dat and .nc output files containing the computed information
        '''

        super(HistMsdData, self).__init__(fname, rootname)
        self.data = HistFile(fname)
        self.data_source = 'Abinit HIST file'
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

    @property
    def read_timestep(self):
        timestep_ps = self.data.reader.read_value('dtion')*atomic_timeunit/1E-12
        return timestep_ps

    @property
    def read_temperature(self):
        self.temperature = self.data.reader.read_value('mdtemp')[0]

    @property
    def read_positions(self):
        return self.data.reader.read_value('xcart')[:, self.atom_indices, :]*bohr_to_ang

    def get_atoms_for_diffusion(self):
        if self.atom_type == 'all':
            self.atom_indices = list(range(self.data.initial_structure.num_sites))
        else:
            if str(self.data.initial_structure.formula).find(self.atom_type) == -1:
                raise ValueError('Did not find atom_type {} in symbols {}'.format(self.atom_type, self.data.initial_structure.formula))
            self.atom_indices = self.data.initial_structure.indices_from_symbol(self.atom_type)


    def compute_msd_from_positions(self):

        displacements = np.zeros((self.nframes, self.natoms, 3))
        for i in range(self.nframes):
            displacements[i, :, :] = self.traj[i, :, :] - self.traj[0, :, :]

        self.msd_atoms = self.compute_atoms_msd(displacements)
        self.msd = np.mean(self.msd_atoms, axis=1)

    def compute_diffusion_coefficient(self, atom_type='all', msd_type='bare', discard_init_steps=0, discard_init_time_ps=None,
                                      discard_final_steps=None,
                                      plot=False, plot_errors=False, plot_verbose=True, plot_all_atoms=False, **kwargs):

        '''
            atom_type: for which atoms the MSD must be computed. Currently, possible options 
                       are "<atom symbol>" for a single atomic specie and  an "all" for all species.

            msd_type: how the MSD should be processed. 
                        "bare": no treatment of raw MSD
                        "timesliced": for each time interval t, the MSD of each individual atom is averaged over all possible
                        time slices equivalent to t (ex.: if t=2, average 2-0, 3-1, 4-2, etc.)

            discard_init_steps: Do not take the N first steps of the trajectory into account when computing the
                                diffusion coefficient.
                                Default: 0

            discard_init_time_ps: time interval (in ps) to discard from slope evaluation
                                default=None (not considered)

            discard_final_steps: Do not take the N last steps of the trajectory into account when computing the
                                diffusion coefficient.
                                Default: None

            plot: activate plotting of MSD vs t

            plot_errors: plot MSD(T) +- standard deviation on all atoms at each timestep, if available

            plot_verbose: print the numerical value of the diffuson coefficient on the plot.
                          Default: True

            plot_all_atoms: plots individual atom MSD in addition to the mean MSD and linear fit
                            Default: False

            **kwargs: optional arguments that will be passes to the Plotter object (see plotter/plotter.py)
        '''

        self.timestep = self.read_timestep
        self.atom_type = atom_type
        logging.info('Will average MSD(T) on {} atoms'.format(self.atom_type))

        self.msd_type = msd_type

        if not isinstance(discard_init_steps, int):
            raise TypeError('discard_init_steps should be an integer, but I got {} which is a {}'.format(discard_init_steps, type(discard_init_steps)))

        if discard_final_steps is not None:
            raise NotImplementedError('discard_final_steps was not tested for Abinit MD.')
            if not isinstance(discard_final_steps, int):
                raise TypeError('discard_final_steps should be an integer, but I got {} which is a {}'.format(discard_final_steps, type(discard_final_steps)))

        self.read_temperature
        self.get_atoms_for_diffusion()
        logging.info('Extracting trajectories...')
        self.traj = self.read_positions
        if discard_final_steps is not None:
            self.traj = self.traj[:-discard_final_steps]
        self.nframes, self.natoms = np.shape(self.traj)[:2]
        # check if the positions are wrapped or unwrapped, with condition like dx larger than half the unit cell?
#        for i, frame in enumerate(self.traj):
#            if np.any(frame > 12):
#                #print('found atoms outside the box in frame {}'.format(i))
#                #print(frame)
#                here = np.where(frame>12)[:][0]
#                if np.any(here>71):
#                    print('found non-Li diffusing atoms in frame {}!!!'.format(i))
        # From this test, positions seem to be unwrapped as at some points some atoms move outside the unit cell
    
        self.time = self.timestep*np.arange(self.nframes)

        if discard_init_time_ps:
            self.discard_init_steps = np.asarray([self.time>=discard_init_time_ps]).nonzero()[1][0]
        else:
            self.discard_init_steps = discard_init_steps

        logging.info('Computing MSD from atomic positions...') 
        self.compute_msd_from_positions()
        logging.info('... done!')

        self.diffusion = self.extract_diffusion_coefficient()
        self.msd_std = self.extract_msd_errors()
        logging.info(f'Diffusion coefficient: {self.diffusion:.3e}+-{self.diffusion_std:.3e} cm^2/s')

        self.write_data()

        if plot:
            self.plot_errors = plot_errors
            self.plot_diffusion_coefficient(defname='diffusion.png', verbose=plot_verbose, plot_all_atoms=plot_all_atoms, **kwargs)
