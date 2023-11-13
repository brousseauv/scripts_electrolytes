#! /usr/bin/env/python

import os
from ase.db import connect
from scripts_electrolytes.interfaces.ase_interface import abistruct_to_ase
from scripts_electrolytes.interfaces.mtp_interface import abistruct_to_cfg
from scripts_electrolytes.interfaces.lammps_interface import abistruct_to_xyz
from scripts_electrolytes.utils.constants import ha_to_ev, bohr_to_ang, gpa_to_evang3
from scripts_electrolytes.database.db_reader import MtpDbReader, AseDbReader, XyzDbReader

'''
    These classes merge databases of atomic configurations previously created with DbCreator.

    Simply call python dbmerger.py --<OPTION1> <value1> -- <OPTION2> <value2> etc.
    or load one of the classes, providing the required arguments

    Options: merged_dbname: filename for the merged database

        filenames: list of paths or filenames to be merged, i.e. [db1, db2, ...]

        format(required): Format of the databases.  Can be either 'mtp', 'ase' or 'xyz'.
                          The 'xyz' format is mostly for visualization purposes with Ovito.

        append: Boolean; indicates if the initial database should be appended in case the filename already exists.
                   Default = False

        atomic_numbers: List of integers specifying atomic numbers in the same order as MTP species (for MTP cfg format only)

        ex: the following command merges databases called mydatabase and anotherdatabase in .cfg format into a new file called merged_database:

            python dbmerger.py --merged_dbname merged_database.cfg --filenames [mydatabase.cfg, anotherdatabase.cfg] --format='mtp'

    For help about these options on the command line, type 
        python dbmerger.py --help

'''
class DbMerger:

    def __init__(self, dbname, filenames, append):

        self.dbname = dbname
        self.append = append
        self.filenames = filenames # check the correct list format


    def merge_db(self):

        if not self.append:
            os.system('cp {} {}'.format(self.filenames[0], self.dbname))

        newdb = self.open_database()

        for db in self.filenames[1:]:
            data = self.read_database(db)
            data.load_database()

            for idx, struct in enumerate(data.structures):
                atoms = self.convert_structure(struct)
                self.add_to_database(newdb, atoms, data.energy[idx], data.forces[idx], data.stresses[idx])


class MtpDbMerger(DbMerger):

    def __init__(self, merged_dbname, filenames, append=False, atomic_numbers=None):

        super(MtpDbMerger, self).__init__(merged_dbname, filenames, append)
        self.check_db_exists()

        if not atomic_numbers:
            raise ValueError('Must define a list for atomic_numbers')
        self.atomic_numbers = atomic_numbers


    def check_db_exists(self):

        if not self.dbname.endswith('.cfg'):
            self.dbname = '{}.cfg'.format(self.dbname)

        if os.path.exists(os.path.join(os.getcwd(), self.dbname)):
            if self.append:
                return
            else:
                raise FileExistsError("""{} file already exists. Either choose another name or use --append True keyword.""".format(
                                       os.path.join(os.getcwd(), self.dbname)))


    def open_database(self):
        return open(self.dbname, 'a')


    def read_database(self, fname):
        data = MtpDbReader(fname, atomic_numbers=self.atomic_numbers)
        return data


    def convert_structure(self, struct):
        return struct


    def add_to_database(self, db, atoms, energy, forces, stresses):
        abistruct_to_cfg(db, atoms, energy=energy, forces=forces, stresses=stresses)



class AseDbMerger(DbMerger):

    def __init__(self, merged_dbname, filenames, append=False):

        super(AseDbMerger, self).__init__(merged_dbname, filenames, append)
        self.check_db_exists()


    def check_db_exists(self):

        if not self.dbname.endswith('.db'):
            self.dbname = '{}.db'.format(self.dbname)

        if os.path.exists(os.path.join(os.getcwd(), self.dbname)):
            if self.append:
                return
            else:
                raise FileExistsError("""{} file already exists. Either choose another name or use --append True keyword.""".format(
                                       os.path.join(os.getcwd(), self.dbname)))

    
    def open_database(self):
        # FIX ME: test if this appends to the db.
        return connect(self.dbname)


    def read_database(self, fname):
        data = AseDbReader(fname)
        return data


    def convert_structure(self, struct):
        return abistruct_to_ase(struct)


    def add_to_database(self, db, atoms, energy, forces, stresses):
        db.write(atoms, data={'energy': energy, 'forces': forces, 'stresses': stresses})



class XyzDbMerger(DbMerger):

    def __init__(self, merged_dbname, filenames, append=False):

        super(XyzDbMerger, self).__init__(merged_dbname, filenames, append)
        self.check_db_exists()


    def check_db_exists(self):

        if not self.dbname.endswith('.xyz'):
            self.dbname = '{}.xyz'.format(self.dbname)

        if os.path.exists(os.path.join(os.getcwd(), self.dbname)):
            if self.append:
                return
            else:
                raise FileExistsError("""{} file already exists. Either choose another name or use --append True keyword.""".format(
                                       os.path.join(os.getcwd(), self.dbname)))


    def open_database(self):
        return open(self.dbname, 'a')


    def read_database(self, fname):
        data = XyzDbReader(fname)
        return data


    def convert_structure(self, struct):
        return struct


    def add_to_database(self, db, atoms, energy, forces, stresses):
        abistruct_to_xyz(db, atoms, energy=energy, forces=forces, stresses=stresses)

