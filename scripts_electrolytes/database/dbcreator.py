#! /usr/bin/env/python

import argparse
import os
from abipy.dynamics.hist import HistFile
from ase.db import connect
from scripts_electrolytes.interfaces.ase_interface import abistruct_to_ase
from scripts_electrolytes.interfaces.mtp_interface import abistruct_to_cfg
from scripts_electrolytes.utils.constants import ha_to_ev, bohr_to_ang


# FIX ME: would this better work as classes?
'''
    These functions create a databse of atomic configurations extracted from Abinit HIST.nc files
    (eventually: and a directory containing multiple GSR.nc files)
    and converts them either in ASE .db format (to be used with SchNetPack)
    or in .cfg format (to be used with MTP/mlip-2 code).

    Simply call python dbcreator.py --<OPTION1> <value1> -- <OPTION2> <value2> etc.

    Options: dbname: filename for the creater database

        path (required, mutually exclusive with fname): Path to the calculation folder containing the configurations

        fname (required, mutually exclusive with path): filename containing the configurations

        source (required): Data source file type.  Can be 'hist' for AIMD runs or 'gsr' for independent configuration.

        format(required): Output format for the database.  Can be either 'mtp' or 'ase'.

        mdskip: Integer; elect each 'mdskip' configuration in the AIMD trajectory (to prevent having too may correlated configurations).  
                Default = 10
        
        overwrite: Boolean; indicates if the database should be overwritten in case the filename already exists.
                   Default = False

        ex: the following command creates a database called mydatabase in .cfg format from a calc_HIST.nc file in subdirectory aimd/,
            selecting one every 50 configurations:

            python dbcreator.py --dbname mydatabase.cfg --fname aimd/calc_HIST.nc --source 'hist' --format='mtp' --mdskip 50

    For help about these options on the command line, type 
        python dbcreator.py --help

'''

def db_from_hist(args):

    hist = HistFile(args.fname)

    structures = hist.structures
    energies = hist.etotals
    forces = hist.reader.read_value('fcart')
    stresses = hist.reader.read_value('strten')

    db = create_database(args)
    dblst = list(range(0, hist.num_steps, args.mdskip))

    for i in dblst:

        atoms = convert_structure(args, structures[i])
        energy = energies[i]  # in eV
        current_forces = forces[i, :, :] * ha_to_ev/bohr_to_ang  # in eV/ang
        current_stress = stresses[i, :] * ha_to_ev/(bohr_to_ang**3)  # in eV/ang^3
        add_to_database(args, db, atoms, energy, current_forces, current_stress)


# make one function for the ASE db format and one for the MTP-MLIP .cfg format

def db_from_gsr(args):

    raise NotImplementedError('db creation from GSR files not yet implemented')
    # Find all GSR.nc files in PATH
    # fname = [f for f in os.listdir(args.path) if f.endswith("GSR.nc")]

    # NOTE: ase.io.read() can read the format abinit-gsr
    # I could also simply use the to_ase_atoms class method


def create_database(args):

    if args.format == 'ase':
        return connect(args.dbname)

    elif args.format == 'mtp':
        if args.dbname.endswith('.cfg'):
            return open(args.dbname, 'w')
        else:
            return open('{}.cfg'.format(args.dbname), 'w')


def convert_structure(args, struct):

    if args.format == 'ase':
        return abistruct_to_ase(struct)

    elif args.format == 'mtp':
        return struct


def add_to_database(args, db, atoms, energy, forces, stresses):

    if args.format == 'ase':
        db.write(atoms, data={'energy': energy, 'forces': forces, 'stresses': stresses})
    elif args.format == 'mtp':
        abistruct_to_cfg(db, atoms, energy=energy, forces=forces, stresses=stresses)


def create_parser():

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--dbname", default=None, help="Database name")
    data = parser.add_mutually_exclusive_group(required=True)
    data.add_argument("--path", help="Path to the calculation folder")
    data.add_argument("--fname", help="HIST.nc file name")

    parser.add_argument("--source", choices=['hist', 'gsr'], help="""Data source file type ('hist' for AIMD runs, 'gsr'
            for independent configuration)""", required=True)
    parser.add_argument("--format", choices=['ase', 'mtp'], help="Output format for the database", required=True)
    parser.add_argument("--mdskip", type=int, default=10, help="Database will include every 'mdskip' configuration")
    parser.add_argument("--overwrite", type=bool, default=False, help="Should an existing database be overwritten or not")

    return parser


def check_parser(args, parser):

    if args.source == 'gsr' and args.path is None:
        parser.error("--source 'gsr' requires --path argument")
    if args.source == 'hist' and args.fname is None:
        parser.error("--source 'hist' requires --fname argument")


def check_db_exists(fname, delete):

    if os.path.exists(os.path.join(os.getcwd(), fname)):
        if delete:
            os.remove(fname)
        else:
            raise FileExistsError("""{} file already exists. Either choose another name or use --overwrite True
                    keyword.""".format(os.path.join(os.getcwd(), fname)))

def main(args):

    check_db_exists(args.dbname, args.overwrite)

    if args.source == 'hist':
        db_from_hist(args)

    elif args.source == 'gsr':
        raise NotImplementedError("Dataase creation from GSR.nc files is not yet implemented,")


if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)
    main(args)
