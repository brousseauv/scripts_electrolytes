#!/usr/bin/env python
from scripts_electrolytes.diffusion.lammps_msd import LammpsMsdData
import numpy as np
import os
import sys

atomic_numbers = [3, 16, 7]
temp = 800 

## Using netCDF dump file:
filetype = 'dump-netcdf' 
output_filename = 'traj.nc'
timestep = 0.001 # in netCDF format, timestep is read from the file, but we need to provide a value 
data = LammpsMsdData(output_filename, filetype=filetype, rootname=f'mtp{temp}K_100ps_2x2x2_netcdf')
data.compute_diffusion_coefficient(timestep=timestep, atomic_numbers=atomic_numbers, input_temperature=temp,
                                   msd_type='bare', plot=True, atom_type='Li', discard_init_steps=0,
                                   savefig=False, showfig=True)

### Using text dump file:
filetype = 'dump' 
output_filename = 'unwrapped_traj.dump'
timestep = 50*0.001 # timestep between dumped configurations
data = LammpsMsdData(output_filename, filetype=filetype, rootname=f'mtp{temp}K_100ps_2x2x2_dump')
data.compute_diffusion_coefficient(timestep=timestep, atomic_numbers=atomic_numbers, input_temperature=temp,
                                   msd_type='bare', plot=True, atom_type='Li', discard_init_steps=0,
                                   savefig=False, showfig=True)

## Using "thermo" data stored in lammps.log file
filetype = 'thermo' 
output_filename = 'lammps.log'
timestep = 0.001 # simulation timestep 
data = LammpsMsdData(output_filename, filetype=filetype, rootname=f'mtp{temp}K_100ps_2x2x2_thermo')
data.compute_diffusion_coefficient(timestep=timestep, atomic_numbers=atomic_numbers, input_temperature=temp,
                                   msd_type='bare', plot=True, atom_type='Li', discard_init_steps=0,
                                   savefig=False, showfig=True)
