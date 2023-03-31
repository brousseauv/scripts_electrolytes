#! /usr/bin/env/python

import argparse
import os
from abipy.core.structure import Structure
from abipy.dynamics.hist import HistFile
from ase.db import connect
import numpy as np
from utils.interfaces.ase_interface import abistruct_to_ase
from utils.constants import ha_to_ev, bohr_to_ang

# FIX ME: should this better be done as classes?

def db_from_hist(args):

    hist = HistFile(args.fname)

    structures = hist.structures
    energies = hist.etotals
    forces = hist.reader.read_variable('fcart')
    stresses = hist.reader.read_variable('strten')

    db = create_database(args)

    dblst = list(range(0, hist.num_steps, args.mdskip))

    for i in dblst:
        
        atoms = convert_structure(args, structures[i])
        energy = energies[i]  # in eV
        current_forces = forces[i, :, :] * ha_to_ev/bohr_to_ang  # in eV/ang
        current_stress = stresses[i, :] * ha_to_ev/(bohr_to_ang**3)  # in eV/ang^3
        # Je suis rendue à tester ça, mais ça devrait fonctionner correctement...
        add_to_database(args, atoms, energy, current_forces, current_stress)
        

# make one function for the ASE db format and one for the MTP-MLIP .cfg format

def db_from_gsr(args):

    # Find all GSR.nc files in PATH
    fname = [f for f in os.listdir(args.path) if f.endswith("GSR.nc")]

    # NOTE: ase.io.read() can read the format abinit-gsr
    # I could also simply use the to_ase_atoms class method

def create_database(args):

    if args.format == 'ase':
        return  connect(args.dbname)

    elif args.format == 'mtp':
        raise NotImplementedError('MTP database creation not yet implemented')

def convert_structure(args, struct):

    if args.format == 'ase':
        return abistruct_to_ase(struct)

    elif args.format == 'mtp':
        raise NotImplementedError('MTP database creation not yet implemented')

def add_to_database(args, atoms, energy, forces, stresses):

    if args.format == 'ase':
        db.write(atoms, data = {'energy': energy, 'forces': current_forces})
    elif args.format == 'mtp':
        raise NotImplementedError('MTP database creation not yet implemented')

def create_parser():

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--dbname", default=None, help="Database name")
    data = parser.add_mutually_exclusive_group(required=True)
    data.add_argument("--path", help="Path to the calculation folder")
    data.add_argument("--fname", help="HIST.nc file name")

    parser.add_argument("--source", choices=['hist', 'gsr'], help="Data source file type ('hist' for AIMD runs, 'gsr' for independent configuration)",
            required=True)
    parser.add_argument("--format", choices=['ase', 'mtp'], help="Output format for the database", required=True)
    parser.add_argument("--mdskip", type=int, default=10, help="Database will include every 'mdskip' configuration")

    return parser

def check_parser(args, parser):

    if args.source =='gsr' and args.path is None:
        parser.error("--source 'gsr' requires --path argument")
    if args.source =='hist' and args.fname is None:
        parser.error("--source 'hist' requires --fname argument")


def main(args):

    if args.source == 'hist':
        db_from_hist(args)

    elif args.source == 'gsr':
        raise NotImplementedError("Dataase creation from GSR.nc files is not yet implemented,")

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)
    main(args)
