#! /usr/bin/env/python

import argparse
import os
from abipy.dynamics.hist import HistFile
from abipy.abilab import abiopen
from ase.db import connect
from scripts_electrolytes.interfaces.ase_interface import abistruct_to_ase
from scripts_electrolytes.interfaces.mtp_interface import abistruct_to_cfg
from scripts_electrolytes.utils.constants import ha_to_ev, bohr_to_ang, gpa_to_evang3
import glob

# FIX ME: would this better work as classes?
'''
    These functions create a databse of atomic configurations extracted from Abinit HIST.nc files
    (eventually: and a directory containing multiple GSR.nc files)
    and converts them either in ASE .db format (to be used with SchNetPack)
    or in .cfg format (to be used with MTP/mlip-2 code).

    Simply call python dbcreator.py --<OPTION1> <value1> -- <OPTION2> <value2> etc.
    or load one of the classes, providing the required arguments

    Options: dbname: filename for the creater database

        path (required, mutually exclusive with fname): Path to the calculation folder containing the configurations

        fname (required, mutually exclusive with path): filename containing the configurations

        source (required): Data source file type.  Can be 'hist' for AIMD runs or 'gsr' for independent configuration.

        format(required): Output format for the database.  Can be either 'mtp' or 'ase'.

        mdskip: Integer; elect each 'mdskip' configuration in the AIMD trajectory (to prevent having too may correlated configurations).  
                Default = 10

        initstep: Integer, index of the first configuration selected in the database.
                  Defaut = 0
        
        overwrite: Boolean; indicates if the database should be overwritten in case the filename already exists.
                   Default = False

        ex: the following command creates a database called mydatabase in .cfg format from a calc_HIST.nc file in subdirectory aimd/,
            selecting one every 50 configurations:

            python dbcreator.py --dbname mydatabase.cfg --fname aimd/calc_HIST.nc --source 'hist' --format='mtp' --mdskip 50

    For help about these options on the command line, type 
        python dbcreator.py --help

'''
class DbCreator:

    def __init__(self, dbname, mdskip, initstep, overwrite, append):

        self.dbname = dbname
        self.mdskip = mdskip
        self.initstep = initstep
        self.overwrite = overwrite
        self.append = append


    def db_from_hist(self, fname):

        hist = HistFile(fname)

        structures = hist.structures
        energies = hist.etotals
        forces = hist.reader.read_value('fcart')
        stresses = hist.reader.read_value('strten')

        db = self.create_database()
        dblst = list(range(self.initstep, hist.num_steps, self.mdskip))

        for i in dblst:

            atoms = self.convert_structure(structures[i])
            energy = energies[i]  # in eV
            current_forces = forces[i, :, :] * ha_to_ev/bohr_to_ang  # in eV/ang
            current_stress = stresses[i, :] * ha_to_ev/(bohr_to_ang**3)  # in eV/ang^3
            self.add_to_database(db, atoms, energy, current_forces, current_stress)


    def db_from_gsr(self, path):

        # List of all GSR.nc files in path
        gsr_list = glob.glob(os.path.join(path, '**/*GSR.nc') ,recursive=True)
        # print(gsr_list)

        db = self.create_database()
        for fname in gsr_list:
            gsr = abiopen(fname)
            structure = gsr.structure
            atoms = self.convert_structure(structure)
            energy = gsr.energy
            current_forces = gsr.cart_forces  # already in eV/ang
            current_stress = gsr.reader.read_value('cartesian_stress_tensor') * ha_to_ev/(bohr_to_ang**3) # convert from ha_bohr3 to eV/ang^3 
            self.add_to_database(db, atoms, energy, current_forces, current_stress)


class MtpDbCreator(DbCreator):

    def __init__(self, dbname=None, mdskip=10, initstep=0, overwrite=False, append=False):

        super(MtpDbCreator, self).__init__(dbname, mdskip, initstep, overwrite, append)
        self.check_db_exists()


    def check_db_exists(self):

        if not self.dbname.endswith('.cfg'):
            self.dbname = '{}.cfg'.format(self.dbname)

        if os.path.exists(os.path.join(os.getcwd(), self.dbname)):
            if self.overwrite:
                os.remove(self.dbname)
            elif self.append:
                return
            else:
                raise FileExistsError("""{} file already exists. Either choose another name or use --overwrite True or --append True keywords.""".format(
                                       os.path.join(os.getcwd(), self.dbname)))


    def create_database(self):
        if self.append:
            return open(self.dbname, 'a')
        else:
            return open(self.dbname, 'w')


    def convert_structure(self, struct):
        return struct


    def add_to_database(self, db, atoms, energy, forces, stresses):
        abistruct_to_cfg(db, atoms, energy=energy, forces=forces, stresses=stresses)



class AseDbCreator(DbCreator):

    def __init__(self, dbname=None, mdskip=10, initstep=0, overwrite=False, append=False):

        super(AseDbCreator, self).__init__(dbname, mdskip, initstep, overwrite, append)
        self.check_db_exists()


    def check_db_exists(self):

        if not self.dbname.endswith('.db'):
            self.dbname = '{}.db'.format(self.dbname)

        if os.path.exists(os.path.join(os.getcwd(), self.dbname)):
            if self.overwrite:
                os.remove(self.dbname)
            elif self.append:
                return
            else:
                raise FileExistsError("""{} file already exists. Either choose another name or use --overwrite True keyword.""".format(
                                       os.path.join(os.getcwd(), self.dbname)))

    
    def create_database(self):
        # FIX ME: test if this appends to the db.
        return connect(self.dbname)


    def convert_structure(self, struct):
        return abistruct_to_ase(struct)


    def add_to_database(self, db, atoms, energy, forces, stresses):
        db.write(atoms, data={'energy': energy, 'forces': forces, 'stresses': stresses})



###################################

# make one function for the ASE db format and one for the MTP-MLIP .cfg format





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
    parser.add_argument("--initstep", type=int, default=0, help="Index of the first configuration selected")
    parser.add_argument("--overwrite", type=bool, default=False, help="Should an existing database be overwritten or not")
    parser.add_argument("--append", type=bool, default=False, help="Should data be appended to existing database or not")

    return parser


def check_parser(args, parser):

    if args.source == 'gsr' and args.path is None:
        parser.error("--source 'gsr' requires --path argument")
    if args.source == 'hist' and args.fname is None:
        parser.error("--source 'hist' requires --fname argument")



def main(args):

    if args.format == 'mtp':
        db = MtpDbCreator(dbname=args.dbname, mdskip=args.mdskip, initstep=args.initstep, overwrite=args.overwrite, append=args.append)
    elif args.format == 'ase':
        db = AseDbCreator(dbname=args.dbname, mdskip=args.mdskip, initstep=args.initstep, overwrite=args.overwrite, append=args.append)

    if args.source == 'hist':
        db.db_from_hist(args.fname)

    elif args.source == 'gsr':
        db.db_from_gsr(args.path)


if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)
    main(args)
