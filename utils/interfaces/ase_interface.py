#!/usr/bin/env python

from ase.io import read


def abistruct_to_ase(struct):

    return struct.to_ase_atoms()
