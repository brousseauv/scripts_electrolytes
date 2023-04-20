import numpy as np
import pandas as pd

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

    return out  # Return output fname??

def create_thermo_dataframe(fname):

    cols = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, nrows=1).columns
    df = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, usecols=cols[:-1])

    return df

def read_thermo(data, key):

    data[key] = pd.to_numeric(data[key])
    return data[key].values

def read_msd(data):

    time = read_thermo(data, 'Time')
    msd = read_thermo(data, 'c_msd[4]')

    return time, msd

