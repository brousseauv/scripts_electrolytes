import numpy as np
import netCDF4 as nc
import os
from pymatgen.io.abinit.netcdf import NetcdfReader

class MsdOutput:

    ''' Base class to postprocess netCDF output files from MsdData class '''

    def __init__(self, fname=None):

        # what do I need to do here?
        try:
            os.path.exists(fname)
        except OSError:
            raise OSError('File {} does not exist!'.format(fname))
        else:
            if not fname.endswith('.nc'):
                raise Exception('files should be in netCDF format')
            else:
                self.fname = fname

    
    def read_data(self):

        reader = NetcdfReader(self.fname)
        self.temp = reader.read_value('temperature')
        self.coeff = reader.read_value('diffusion_coefficient')
        self.msd_type = reader.rootgrp.getncattr('msd_type')

        # I could add the other variables...
