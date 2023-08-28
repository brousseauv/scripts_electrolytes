import numpy as np
import os
from .neb import NebData
from ..interfaces.lammps_interface import read_neb_logfile
import logging

class LammpsNebData(NebData):

    def __init__(self, fname='log.lammps', rootname='neb'):


        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

        super(LammpsNebData, self).__init__(fname, rootname)
        self.read_data()

    def read_data(self):
        self.forward_barrier, self.backward_barrier, self.reaction_coordinate_length, self.reaction_coordinate, self.potential_energy = read_neb_logfile(self.fname)
