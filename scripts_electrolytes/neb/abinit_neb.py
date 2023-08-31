from .neb import NebData
from ..interfaces.abinit_histfile import read_neb_efs
import logging
import os
import numpy as np
from abipy.dynamics.hist import HistFile
from abipy.core.structure import Structure
from ..utils.constants import ha_to_ev, bohr_to_ang

class HistNebData(NebData):

    def __init__(self, fname=None, rootname='neb_from_hist', rescale_energy=True):
        
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

        try:
            os.path.exists(fname)
        except:
            raise NameError('File {} does not exist. Please provide the correct path to the Abinit HIST.nc file containing NEB results.'.format(fname))

        super(HistNebData, self).__init__(fname, rootname)
        self.read_data(rescale_energy)

        self.print_barriers()


    def read_data(self, rescale_energy):
        
        hist = HistFile(self.fname)
        self.potential_energy, forces, stresses = read_neb_efs(hist)
        if rescale_energy:
            self.potential_energy -= self.potential_energy[0]

        nimage = hist.reader.read_dimvalue('nimage')
        self.reaction_coordinate = np.linspace(0, 1, num=nimage, endpoint=True)

        saddle = np.argmax(self.potential_energy)
        self.forward_barrier = np.max(self.potential_energy) - np.min(self.potential_energy[:saddle])
        self.backward_barrier = np.max(self.potential_energy) - np.min(self.potential_energy[saddle:])
