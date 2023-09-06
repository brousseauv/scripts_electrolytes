import numpy as np
import os
from .neb import NebData, NebTraj
from ..interfaces.lammps_interface import read_neb_logfile, read_config_from_dump
import logging

class LammpsNebData(NebData):

    def __init__(self, fname='log.lammps', rootname='neb_from_lammps', rescale_energy=True):


        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
        
        try:
            os.path.exists(fname)
        except:
            raise NameError('File {} does not exists. Please provide correct path to the "log.lammps" output file containing NEB results.'.format(fname))

        super(LammpsNebData, self).__init__(fname, rootname)

        self.read_data(rescale_energy)
        self.print_barriers()

    def read_data(self, rescale_energy):
        self.forward_barrier, self.backward_barrier, self.reaction_coordinate_length, self.reaction_coordinate, self.potential_energy = read_neb_logfile(self.fname, rescale_energy)


class LammpsNebTraj(NebTraj):

    def __init__(self, in_rootname=None, nreplica=None,  atomic_numbers=None):

        self.check_input(in_rootname, nreplica, atomic_numbers)

        super(LammpsNebTraj, self).__init__(nreplica)

        self.read_neb_trajectory()


    def check_input(self, root, nimg, zval):

        ''' Check if the provided input contains all the required information'''
        if not os.path.exists('{}.1'.format(root)):
            raise ValueError('Could not find the first NEB dump file {}.1'.format(root))
        self.in_rootname = root

        if not nimg:
            raise Exception('Must specify the number of replicas (i.e. the total number of dump.X files)')

        if not os.path.exists('{}.{}'.format(root, nimg)):
            raise ValueError('Could not find the last NEB dump file {}.{} . Perhaps check nreplica?'.format(root, nimg))

        if not zval:
            raise Exception('Must provide a list of atomic numbers')
        self.atomic_numbers = zval


    def read_neb_trajectory(self):

        # This will be a list of ASE Atoms objects
        self.trajectory = []

        # LAMMPS NEB dump files for replicas start at 1
        for irep in np.arange(1, self.nreplica+1):
            fname = '{}.{}'.format(self.in_rootname, irep)
            replica = read_config_from_dump(fname, self.atomic_numbers, which=-1)
            self.trajectory.append(replica)
