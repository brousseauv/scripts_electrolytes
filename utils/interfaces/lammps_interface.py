import numpy as np
import pandas as pd
from ase.io import read
from ase.md.analysis import DiffusionCoefficient

''' Some functions to treat the outputs from a LAMMPS run'''

def extract_thermo(fname, out='thermo.dat'):

    f = open(fname, 'r')
    g = open(out, 'w')
    start = False
    
    for line in f.readlines():
        if line.find('Loop time of') != -1:
            break
        if start:
            g.write(line)
        if line.find('Per MPI rank') != -1:
            start = True

    return out

def create_thermo_dataframe(fname):

    cols = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, nrows=1).columns
    df = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, usecols=cols[:-1])

    return df


def read_thermo(data, key):

    data[key] = pd.to_numeric(data[key])
    return data[key].values


def read_msd_from_thermo(data):

    time = read_thermo(data, 'Time')
    msd = read_thermo(data, 'c_msd[4]')
    temp = read_thermo(data, 'Temp')[0]

    return time, msd, temp


def read_traj_from_dump(fname, atomic_numbers):

    traj= read(fname, format='lammps-dump-text', index=':')

    # as the lammps dump outputs only atom id (1,2,3...) and not type, ASE sees H, He, Li...
    # So, convert atom id to atomic masses
    for frame in traj:
        for a in range(len(frame.numbers)):
            frame.numbers[a] = atomic_numbers[frame.numbers[a]-1]
    return traj
