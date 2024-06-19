#!/usr/bin/env python
from scripts_electrolytes.diffusion.msd_output import MsdOutput
import numpy as np
import os
import sys

temp = 800
atomic_numbers = [3, 16, 7]

flush_nsteps = 50  # flush first 50 stepsa (those are the saved steps, not the actual timesteps)
flush_dt = 10  # Flush first 10 ps

mydir = 'OUT'
os.chdir(mydir)    
data = MsdOutput(f'mtp{temp}K_100ps_2x2x2_netcdf.nc')

# Discarding based on number of MD steps
data.recompute_diffusion_coefficient(discard_init_steps=flush_nsteps)
print(f'When discarding {flush_nsteps} MD steps, diffusion coefficient={data.coeff:.3e} cm^2/s')

# Discarding based on time
data.recompute_diffusion_coefficient(discard_init_time_ps=flush_dt)
print(f'When discarding {flush_dt} ps, diffusion coefficient={data.coeff:.3e} cm^2/s')
