from ..interfaces.lammps_interface import read_traj_from_dump
from ..interfaces.ase_interface import ase_to_abistruct
from ..interfaces.mtp_interface import split_cfg_configs, convert_chunk_to_abivars
from ..interfaces.abinit_interface import abivars_to_abistruct
from ..interfaces.partn_interface import fix_species_in_xyz_mlip
import numpy as np
from ase.io import read as ase_read
from ase.db import connect as ase_connect

class DbReader:

    def __init__(self, fname):

        self.fname = fname


    def initialize_arrays(self, nconfig, natom):

        if self.has_energy:
            self.energy = np.zeros((nconfig))
        else:
            self.energy = [None] * nconfig

        if self.has_forces:
            self.forces = np.zeros((nconfig, natom, 3))
        else:
            self.forces = [None] * nconfig

        if self.has_stress:
            self.stresses = np.zeros((nconfig, 6))
        else:
            self.stresses = [None] * nconfig

        self.structures = []


class AseDbReader(DbReader):

    def __init__(self, fname):

        super(AseDbReader, self).__init__(fname)


    def load_database(self):

        db = ase_connect(self.fname)
        traj = [row for row in db.select()]  # This creates a list of AtomsRow objects
        
        self.set_properties(traj)
        nconfig = len(traj)
        natom = traj[0].natoms

        self.initialize_arrays(nconfig, natom)

        # they should already be in the correct units (eV, eV/ang, eV/ang^3)
        for i in range(nconfig):
            struct = traj[i].toatoms() # Convert AtomsRow object to Atoms object
            self.structures.append(ase_to_abistruct(struct))

            if self.has_energy:
                self.energy[i] = traj[i].data['energy']
            if self.has_forces:
                self.forces[i, :, :] = traj[i].data['forces']
            if self.has_stress:
                self.stresses[i, :] = traj[i].data['stresses']

#    
    def set_properties(self, traj):

        info = list(traj[0].data.keys())

        # Here we assume that all database entries contain the same properties. 
        # check if traj contains energies
        if 'energy' in info and traj[0].data['energy'] is not None:
            self.has_energy = True
        else:
            self.has_energy = False

        # check if traj contains forces
        if 'forces' in info and traj[0].data['forces'] is not None:
            self.has_forces = True
        else:
            self.has_forces = False

        # check if traj contains stresses
        if 'stresses' in info and traj[0].data['stresses'] is not None:
            self.has_stress = True
        else:
            self.has_stress = False


class MtpDbReader(DbReader):

    def __init__(self, fname, atomic_numbers = None):

        super(MtpDbReader, self).__init__(fname)

        if not atomic_numbers:
            raise ValueError('Must define a list for atomic_numbers')
        self.atomic_numbers = atomic_numbers


    def load_database(self):

        configs = split_cfg_configs(self.fname)
        nconfig = len(configs)

        #now, convert each chunk to abivars and define a Structure object
        for i, chunk in enumerate(configs):
            abivars, energy, forces, stresses = convert_chunk_to_abivars(chunk, self.atomic_numbers)

            if i == 0:
                self.set_properties(energy, forces, stresses)
                self.initialize_arrays(nconfig, abivars['natom'])

            self.structures.append(abivars_to_abistruct(abivars))

            if self.has_energy:
                self.energy[i] = energy
            if self.has_forces:
                self.forces[i, :, :] = forces
            if self.has_stress:
                self.stresses[i, :] = stresses


    def set_properties(self, e, f, s):

        # Here we assume that all database entries contain the same properties. 
        # check if data contains energies
        if e is not None:
            self.has_energy = True
        else:
            self.has_energy = False

        # check if data contains forces
        if f is not None:
            self.has_forces = True
        else:
            self.has_forces = False

        # check if data contains stresses
        if s is not None:
            self.has_stress = True
        else:
            self.has_stress = False


class XyzDbReader(DbReader):

    def __init__(self, fname, fix_species=False, symbols=None):

        super(XyzDbReader, self).__init__(fname)
        if fix_species:
            if not symbols:
                raise ValueError('when fix_species is activated in XyzReader, a list of atomic symbols as "symbols" should be provided')

            fix_species_in_xyz_mlip(self.fname, symbols)


    def load_database(self):

        traj = ase_read(filename=self.fname, index=':')
        
        self.set_properties(traj)
        nconfig = len(traj)
        natom = np.shape(traj[0].get_positions())[0]

        self.initialize_arrays(nconfig, natom)

        # they should already be in the correct units (eV, eV/ang, eV/ang^3)
        for i in range(nconfig):
            self.structures.append(ase_to_abistruct(traj[i]))

            if self.has_energy:
                self.energy[i] = traj[i].info['energy']
            if self.has_forces:
                self.forces[i, :, :] = traj[i].calc.get_forces()
            if self.has_stress:
                strs = traj[i].info['stress']
                self.stresses[i, :] = strs[0,0], strs[1,1], strs[2,2], strs[1,2], strs[0,2], strs[0,1]

    
    def set_properties(self, traj):

        info = list(traj[0].info.keys())

        # check if traj contains energies
        if 'energy' in info:
            self.has_energy = True
        else:
            self.has_energy = False

        # check if traj contains forces
        try:
            natom = np.shape(traj[0].calc.get_forces())[0]
            self.has_forces = True
        except:
            self.has_forces = False

        # check if traj contains stresses
        if 'stress' in info:
            self.has_stress = True
        else:
            self.has_stress = False




class DumpDbReader(DbReader):

    def __init__(self, fname, atomic_numbers = None):

        super(DumpDbReader, self).__init__(fname)

        if not atomic_numbers:
            raise ValueError('Must define a list for atomic_numbers')
        self.atomic_numbers = atomic_numbers


    def load_database(self):

        traj = read_traj_from_dump(self.fname, self.atomic_numbers)
        # now traj is a sequence of Atoms objects, must be converted to abistruct

        self.set_properties(traj)
        nconfig = len(traj)
        natom = np.shape(traj[0].get_positions())[0]

        self.initialize_arrays(nconfig, natom)

        # they should already be in the correct units (eV, eV/ang, eV/ang^3)
        for i in range(nconfig):
            self.structures.append(ase_to_abistruct(traj[i]))

            # There should NOT be total energy and stress data in the dump file. 
            # Treating only the forces. 
            if self.has_forces:
                self.forces[i, :, :] = traj[i].calc.get_forces()


    def set_properties(self, traj):

        info = list(traj[0].info.keys())

        # dump files are intended for atom properties. 
        # so they can contain forces but no energy/stress
        # for now, disabling has_energy and has_stress

        self.has_energy = False
        self.has_stress = False

        # check if traj contains forces
        try:
            natom = np.shape(traj[0].calc.get_forces())[0]
            self.has_forces = True
        except:
            self.has_forces = False
