import numpy as np
import os
from .neb import NebData
from ..interfaces.lammps_interface import read_neb_logfile
import logging

class LammpsNebData(NebData):

    def __init__(self, fname='log.lammps', rootname='neb_from_lammps'):


        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
        
        try:
            os.path.exists(fname)
        except:
            raise NameError('File {} does not exists. Please provide correct path to the "log.lammps" output file containing NEB results.'.format(fname))

        super(LammpsNebData, self).__init__(fname, rootname)
        self.read_data()

        self.print_barriers()


    def read_data(self):
        self.forward_barrier, self.backward_barrier, self.reaction_coordinate_length, self.reaction_coordinate, self.potential_energy = read_neb_logfile(self.fname)
