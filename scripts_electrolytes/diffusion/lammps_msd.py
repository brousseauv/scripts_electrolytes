import numpy as np
from ..interfaces.lammps_interface import extract_thermo, create_thermo_dataframe, read_msd_from_thermo, read_traj_from_dump
from .msd import MsdData
from ase.md.analysis import DiffusionCoefficient
import logging
import os

class LammpsMsdData(MsdData):

    def __init__(self, fname, filetype, rootname='MsdData'):

        '''
            Input:
                fname: name of the file containing MSD/atomic positions information

                filetype: type of file from which the MSD has to be extracted.
                            "thermo": output of LAMMPS thermo command which includes MSD
                            "dump": UNWRAPPED LAMMPS trajectory file
                            "dump-netcdf": similar to "dump" but in netCDF format

                rootname: rootname for the .dat and .nc output files containing the computed information

        '''

        if filetype not in ['thermo', 'dump', 'dump-netcdf']:
            raise Exception('''filetype should be one of the following: "thermo", "dump",
                               "dump-netcdf", but I got {}'''.format(filetype))

        self.filetype = filetype
        super(LammpsMsdData, self).__init__(fname, rootname)
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

    
    def compute_diffusion_coefficient(self, thermo_fname=None, timestep=None, atom_type='all', atomic_numbers=None, 
                                      msd_type='bare', input_temperature=None, discard_init_steps=0, plot=False, 
                                      plot_errors=False, plot_verbose=True, plot_all_atoms=False, **kwargs):

        '''
            thermo_fname: path to the thermo.dat file is already extracted from lammps.log

            timestep: MD timestep, in picosecond

            atom_type: for which atoms the MSD must be computed. Currently, possible options 
                       are "<atom symbol>" for a single atomic specie and  an "all" for all species.

            atomic_numbers: list of atomic numbers with the same ordering as the dump file, i.e. [<atom type 0>, <atom type 1> ...]
                            for "dump" filetype only.

            msd_type: how the MSD should be processed. 
                        "bare": no treatment of raw MSD
                        "timesliced": for each time interval t, the MSD of each individual atom is averaged over all possible
                        time slices equivalent to t (ex.: if t=2, average 2-0, 3-1, 4-2, etc.)

            input_temperature: Running temperature of the MD run. For "dump" filetype only.

            discard_init_steps: Do not take the N first steps of the trajectory into account when computing MSD.
                                Default: 0
                                Note: only the default value can be used with "thermo" files.
                                (since the average MSD on atoms is by default taken from the initial step)

            plot: activate plotting of MSD vs t

            plot_errors: plot MSD(T) +- standard deviation on all atoms at each timestep, if available

            plot_verbose: print the numerical value of the diffuson coefficient on the plot.
                          Default: True

            plot_all_atoms: plots individual atom MSD in addition to the mean MSD and linear fit
                            Default: False

            **kwargs: optional arguments that will be passes to the Plotter object (see plotter/plotter.py)
        '''

        self.msd_type = msd_type

        if not isinstance(discard_init_steps, int):
            raise TypeError('discard_init_steps should be an integer, but I got {} which is a {}'.format(discard_init_steps, type(discard_init_steps)))

        if self.filetype == 'thermo':
            if not thermo_fname:
                thermo_fname = extract_thermo(self.fname)
            self.data_source = 'LAMMPS thermo data'
            data = create_thermo_dataframe(thermo_fname)
            logging.info('Extracting MSD from thermo file...')
            self.time, self.msd, self.temperature = read_msd_from_thermo(data)
            self.msd_atoms = None
            self.atom_type = 'See lammps input file'
            self.nframes = len(self.msd)
            self.natoms = None
            self.atom_type = 'See lammps input file'

            if discard_init_steps != 0:
                raise ValueError('With "thermo" files, discard_init_steps should be 0 (default value)')

        elif self.filetype == 'dump':
            if not timestep:
                raise ValueError('Missing value for timestep (in ps)')
            if not atomic_numbers:
                raise ValueError('Must define a list for atomic_numbers')
            if self.msd_type not in ['bare', 'timesliced']:
                raise ValueError('msd_type should be either "bare" or "timesliced" but I got {}'.format(self.msd_type))
            if not input_temperature:
                raise ValueError('Must specify MD running temperature as input_temperature (int)')
            elif not isinstance(input_temperature, int):
                raise ValueError('input_temperature must be an integer, but I got a {}'.format(type(input_temperature)))
            else:
                self.temperature = input_temperature

            self.data_source = 'LAMMPS .dump file'
            self.atom_type = atom_type
            logging.info('Will average MSD(T) on {} atoms'.format(self.atom_type))

            logging.info('Extracting trajectories...')
            self.traj = read_traj_from_dump(self.fname, atomic_numbers)
            # Discard some initial timesteps
            self.traj = self.traj[discard_init_steps:]
            self.get_atoms_for_diffusion()
            self.nframes = len(self.traj)
            self.natoms = len(self.atom_indices)
            logging.info('Computing MSD from atomic positions...')
            self.compute_msd_from_positions()
            logging.info('... done!')
            self.time = timestep*np.arange(len(self.msd))
            # Just curious, does this work?!?
#            self.coeff = self.get_diffusion_ase(timestep)

        elif self.filetype == 'dump-netcdf':
            self.data_source = 'LAMMPS .dump netCDF file'
            raise NotImplementedError('"dump-netcdf" filetype not yet implemented')

        self.diffusion = self.extract_diffusion_coefficient()
        self.msd_std = self.extract_msd_errors()
        logging.info('Diffusion coefficient: {:.3e} cm^2/s'.format(self.diffusion))

        if plot:
            self.plot_errors = plot_errors
            self.plot_diffusion_coefficient(defname='diffusion.png', verbose=plot_verbose, plot_all_atoms=plot_all_atoms, **kwargs)

        self.write_data()

    def get_atoms_for_diffusion(self):
        
        if self.atom_type == 'all':
            self.atom_indices = [i for i in range(len(self.traj[0].get_chemical_symbols()))]
        else:
            if str(self.traj[0].symbols).find(self.atom_type) == -1:
                raise ValueError('Did not find atom_type {} in symbols {}'.format(self.atom_type, self.traj[0].symbols))
            self.atom_indices = [i for i, symbol in enumerate(self.traj[0].get_chemical_symbols()) if "Li" in str(symbol)]

    def get_diffusion_ase(self, timestep):

        # Did not check if this works...
        # It's mostly for sanity check
        coeff = DiffusionCoefficient(self.traj, timestep)

    def compute_msd_from_positions(self):

        # Could add a cutoff for thermalization here, discarding the first N% steps and defining the t=0 at a later point in the MD run
        displacements = np.zeros((self.nframes, self.natoms, 3))
        for i, frame in enumerate(self.traj):
            displacements[i, :, :] = self.traj[i].get_positions()[self.atom_indices] - self.traj[0].get_positions()[self.atom_indices]

        self.msd_atoms = self.compute_atoms_msd(displacements)
        self.msd = np.mean(self.msd_atoms, axis=1)



