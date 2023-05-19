from .msd import MsdData
from abipy.dynamics.hist import HistFile
from ..utils.constants import atomic_timeunit, bohr_to_ang
import numpy as np
import logging
import os

class HistMsdData(MsdData):

    def __init__(self, fname, rootname='MsdData'):

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

    def compute_diffusion_coefficient(self, atom_type='all', msd_type='bare', plot=False, plot_errors=False, **kwargs):

        timestep = self.read_timestep
        self.atom_type = atom_type
        logging.info('Will average MSD(T) on {} atoms'.format(self.atom_type))

        self.msd_type = msd_type
        self.read_temperature
        self.get_atoms_for_diffusion()
        logging.info('Extracting trajectories...')
        self.traj = self.read_positions
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
    
        # need : self.time, self.msd_atoms, self.msd, self.msd_std
        self.time = timestep*np.arange(self.nframes)
        logging.info('Computing MSD from atomic positions...') 
        self.compute_msd_from_positions()
        logging.info('... done!')

        self.diffusion = self.extract_diffusion_coefficient()
        self.msd_std = self.extract_msd_errors()
        logging.info('Diffusion coefficient: {:.3e} cm^2/s'.format(self.diffusion))

        self.write_data()

        if plot:
            self.plot_errors = plot_errors
            self.plot_diffusion_coefficient(defname='diffusion.png', **kwargs)
