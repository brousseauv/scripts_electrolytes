#!/usr/bin/env python

from ase.io import read
from abipy.core.structure import Structure

def abistruct_to_ase(struct):

    return struct.to_ase_atoms()

def ase_to_abistruct(atoms):

    return Structure.from_ase_atoms(atoms)
