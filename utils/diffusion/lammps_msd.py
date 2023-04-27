import numpy as np
from ..interfaces.lammps_interface import extract_thermo, create_thermo_dataframe, read_msd_from_thermo, read_traj_from_dump
from .msd import MsdData
from ase.md.analysis import DiffusionCoefficient
import logging
import os

class LammpsMsdData(MsdData):

    def __init__(self, fname, filetype):

        if filetype not in ['thermo', 'dump', 'dump-netcdf']:
            raise Exception('''filetype should be one of the following: "thermo", "dump",
                               "dump-netcdf", but I got {}'''.format(filetype))

        self.filetype = filetype
        super(LammpsMsdData, self).__init__(fname)
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

    
    def compute_diffusion_coefficient(self, thermo_fname=None, timestep=None, atom_type='all', atomic_numbers=None, 
                                      msd_type='bare', plot=False, **kwargs):

        if self.filetype == 'thermo':
            if not thermo_fname:
                thermo_fname = extract_thermo(self.fname)
            data = create_thermo_dataframe(thermo_fname)
            logging.info('Extracting MSD from thermo file...')
            self.time, self.msd = read_msd_from_thermo(data)
            self.msd_atoms = None

        elif self.filetype == 'dump':
            if not timestep:
                raise ValueError('Missing value for timestep (in ps)')
            if not atomic_numbers:
                raise ValueError('Must define a list for atomic_numbers')
            if msd_type not in ['bare', 'timesliced']:
                raise ValueError('msd_type should be either "bare" or "timesliced" but I got {}'.format(msd_type))
            print('Will average MSD(T) on {} atoms'.format(atom_type))

            logging.info('Extracting trajectories...')
            self.traj = read_traj_from_dump(self.fname, atomic_numbers)
            self.get_atoms_for_diffusion(atom_type)
            self.nframes = len(self.traj)
            self.natoms = len(self.atom_indices)

            logging.info('Computing MSD from atomic positions...')
            self.compute_msd_from_positions(msd_type)
            logging.info('... done!')
            self.time = timestep*np.arange(len(self.msd)) # check the final shape, it should be just len=N timesteps

            # Just curious, does this work?!?
#            self.coeff = self.get_diffusion_ase(timestep)

        elif self.filetype == 'dump-netcdf':
            raise NotImplementedError('"dump-netcdf" filetype not yet implemented')

        self.diffusion = self.extract_diffusion_coefficient()
        self.msd_std = self.extract_msd_errors()
        print('Diffusion coefficient: {:.3e} cm^2/s'.format(self.diffusion))
        logging.info('Diffusion coefficient: {:.3e} cm^2/s'.format(self.diffusion))

        if plot:
            self.plot_diffusion_coefficient(**kwargs)


    def get_atoms_for_diffusion(self, atom_type):
        
        if atom_type == 'all':
            self.atom_indices = [i for i in range(len(self.traj[0].get_chemical_symbols()))]
        else:
            if str(self.traj[0].symbols).find(atom_type) == -1:
                raise ValueError('Did not find atom_type {} in symbols {}'.format(atom_type, self.traj[0].symbols))
            self.atom_indices = [i for i, symbol in enumerate(self.traj[0].get_chemical_symbols()) if "Li" in str(symbol)]

    def get_diffusion_ase(self, timestep):

        # Did not check if this works...
        coeff = DiffusionCoefficient(self.traj, timestep)

    def compute_msd_from_positions(self, msd_type):

        # Could add a cutoff for thermalization here, discarding the first N% steps and defining the t=0 at a later point in the MD run
        displacements = np.zeros((self.nframes, self.natoms, 3))
        for i, frame in enumerate(self.traj):
            displacements[i, :, :] = self.traj[i].get_positions()[self.atom_indices] - self.traj[0].get_positions()[self.atom_indices]

        self.msd_atoms = self.compute_atoms_msd(displacements, msd_type)
        self.msd = np.mean(self.msd_atoms, axis=1)



