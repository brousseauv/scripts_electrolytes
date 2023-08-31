import numpy as np
import os
from .neb import NebData, NebTraj
from ..interfaces.lammps_interface import read_neb_logfile
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

    def __init__(self, fname, rootname):

        # check if fname exists
        # super parent class

        # read the structures : they are loaded as ASE trajectory objects (one per replica)
        # see read_traj_from_dump
        # output them in a given format : lammps dump file, or ASE trajectory files? (NO you need the pro version for this)
        # or XYZ or CIF

        # or do I need classes for these? this is most likely just a conversion scripts... 
        # no need for a class, just a function? 
        # input: rootname for the dump.X files, rootname for output file, format of output file (xyz, dump or cif)
        x=1
