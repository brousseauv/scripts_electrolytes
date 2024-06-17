#!/usr/bin/env python

from ase.io import read
import pandas as pd
from abipy.core.structure import Structure

def abistruct_to_ase(struct):

    return struct.to_ase_atoms()

def ase_to_abistruct(atoms):

    return Structure.from_ase_atoms(atoms)


def create_mdlogger_dataframe(fname):
    ''' Read output data from ASe MDLogger output file and store in DataFrame '''

    cols = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, nrows=1).columns
    df = pd.read_csv(fname, delimiter=' ', skipinitialspace=True)

    return df

