import numpy as np
from ase.io import read as ase_read
from ase import units
from .msd import MsdData
from ase.md.analysis import DiffusionCoefficient
import logging
import os
import warnings

class AseMsdData(MsdData):

    def __init__(self, fname, rootname='MsdData'):

        '''
            Input:
                fname: name of the ASE trajectory file containing MSD/atomic positions information

                rootname: rootname for the .dat and .nc output files containing the computed information

        '''

        fileext = os.path.splitext(fname)[-1]
        if fileext not in ['.traj', '.trj']:
            raise Exception('''ASE trajectory file should have one of the following extensions: ".traj", ".trj",
                               but I got "{}"'''.format(fileext))

        super(AseMsdData, self).__init__(fname, rootname)
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

    
    def compute_diffusion_coefficient(self, timestep=None, atom_type='all', msd_type='bare', input_temperature=None,
                                      discard_init_steps=0, discard_init_time_ps=None, plot=False, plot_errors=False,
                                      plot_verbose=True, plot_all_atoms=False, **kwargs):

        '''
            timestep: MD timestep, in picosecond. DO NOT USE ASE UNITS MODULE!


            atom_type: for which atoms the MSD must be computed. Currently, possible options 
                       are "<atom symbol>" for a single atomic specie and  an "all" for all species.

            msd_type: how the MSD should be processed. 
                        "bare": no treatment of raw MSD
                        "timesliced": for each time interval t, the MSD of each individual atom is averaged over all possible
                        time slices equivalent to t (ex.: if t=2, average 2-0, 3-1, 4-2, etc.)

            input_temperature: Running temperature of the MD run.

            discard_init_steps: Do not take the N first steps of the trajectory into account when fitting the diffusion
                                coefficient fromthe MSD.
                                Default: 0

            discard_init_time_ps: time interval (in ps) to discard from slope evaluation
                                default=None (not considered)

            plot: activate plotting of MSD vs t

            plot_errors: plot MSD(t) +- standard deviation on all atoms at each timestep, if available

            plot_verbose: print the numerical value of the diffuson coefficient on the plot.
                          Default: True

            plot_all_atoms: plots individual atom MSD in addition to the mean MSD and linear fit
                            Default: False

            **kwargs: optional arguments that will be passes to the Plotter object (see plotter/plotter.py)
        '''

        self.msd_type = msd_type

        if not isinstance(discard_init_steps, int):
            raise TypeError('discard_init_steps should be an integer, but I got {} which is a {}'.format(discard_init_steps, type(discard_init_steps)))

        if not timestep:
            raise ValueError('Missing value for timestep (in ps)')
        else:
            if not isinstance(timestep, float):
                raise ValueError('Timestep shoulf be a float, but I got {}'.format(type(timestep)))
            else:
                logging.warning(f'Timestep is set to {timestep:.2f}ps. It should NOT have been defined using ase.units module.')
                self.timestep = timestep
        if self.msd_type not in ['bare', 'timesliced']:
            raise ValueError('msd_type should be either "bare" or "timesliced" but I got {}'.format(self.msd_type))
        if not input_temperature:
            raise ValueError('Must specify MD running temperature as input_temperature (int)')
        elif not isinstance(input_temperature, int):
            raise ValueError('input_temperature must be an integer, but I got a {}'.format(type(input_temperature)))
        else:
            self.temperature = input_temperature


        self.atom_type = atom_type
        logging.info('Will average MSD(T) on {} atoms'.format(self.atom_type))

        logging.info('Extracting trajectories...')

        self.data_source = 'ASE trajectory file'
        ### check if I can have unwrapped positions!!!
        ## Idea: do it assuming they are wrapped, then use ase's built-in DiffusionCoefficient class and compare.
        self.traj = ase_read(self.fname, index=':')

        self.get_atoms_for_diffusion()
        self.nframes = len(self.traj)
        self.natoms = len(self.atom_indices)
        logging.info('Computing MSD from atomic positions...')
        self.compute_msd_from_positions()
        logging.info('... done!')

        self.time = timestep*np.arange(len(self.msd))

        if discard_init_time_ps:
            self.discard_init_steps = np.asarray([self.time>=discard_init_time_ps]).nonzero()[1][0]
        else:
            self.discard_init_steps = discard_init_steps

        # For sanity check, compare with ASE class
        self.coeff_from_ase = self.get_diffusion_ase(timestep)

        self.diffusion = self.extract_diffusion_coefficient()
        self.msd_std = self.extract_msd_errors()
        logging.info(f'Diffusion coefficient: {self.diffusion:.3e}+-{self.diffusion_std:.3e} cm^2/s')

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

        # convert timesteps in ASE units
        ase_timestep = timestep*1E3*units.fs
        # This is mostly for sanity check
        coeff = DiffusionCoefficient(self.traj, ase_timestep, atom_indices=self.atom_indices)
        coeff.calculate(ignore_n_images=self.discard_init_steps)
        slopes, std = coeff.get_diffusion_coefficients()

        conversion_slope = units.fs*1e-1
        for sym_index in range(coeff.no_of_types_of_atoms):
            print('Mean Diffusion Coefficient (from ASE) : %s = %.3e Å^2/cm; Std. Dev. = %.3e Å^2/cm' %
                  (coeff.types_of_atoms[sym_index], slopes[sym_index] * conversion_slope, std[sym_index] * conversion_slope))

    def compute_msd_from_positions(self):

        displacements = np.zeros((self.nframes, self.natoms, 3))
        for i, frame in enumerate(self.traj):
            displacements[i, :, :] = self.traj[i].get_positions()[self.atom_indices] - self.traj[0].get_positions()[self.atom_indices]

        self.msd_atoms = self.compute_atoms_msd(displacements)
        self.msd = np.mean(self.msd_atoms, axis=1)



