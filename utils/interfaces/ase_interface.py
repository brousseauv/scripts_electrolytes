#!/usr/bin/env python

from ase.io import read

''' Convert Abipy Structure object to ASE Atoms object '''

def abistruct_to_ase(struct):

    return struct.to_ase_atoms

