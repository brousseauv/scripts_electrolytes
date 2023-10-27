import os
from .db_creator import MtpDbCreator, AseDbCreator, XyzDbCreator
from .db_reader import AseDbReader, MtpDbReader, XyzDbReader, DumpDbReader


''' 
    This class converts databases between different formats.
    Input format: "mtp", "ase", "xyz" or "dump"
    Output_format: "mtp", "ase" or "xyz"
'''


class DbConverter:

    def __init__(self, fname, input_format = None, output_format=None, dbname=None, overwrite=False, atomic_numbers=None, start=0, every=1):

        self.fname = fname
        self.overwrite = overwrite

        self.check_input_format(input_format)
        self.set_output_format(output_format)
        self.set_dbname(dbname)
        self.atomic_numbers = atomic_numbers
        self.every = every
        self.start = start

    def check_input_format(self, fmt):

        if fmt not in ['mtp', 'ase', 'xyz', 'dump']:
            raise ValueError('input_format must be either "mtp", "ase", "xyz" or "dump", but I got "{}"'.format(fmt))
        if fmt == 'mtp' and not self.fname.endswith('.cfg'):
            raise Exception('Input file {} does not have .cfg extension expected from MTP format')
        if fmt == 'ase' and not self.fname.endswith('.db'):
            raise Exception('Input file {} does not have .db extension expected from ASE format')
        if fmt == 'xyz' and not self.fname.endswith('.xyz'):
            raise Exception('Input file {} does not have .xyz extension expected from XYZ format')
        if fmt == 'dump' and not self.fname.endswith('.dump'):
            if not self.fname.endswith('dmp'):
                raise Exception('Input file {} does not have .dump or .dmp extension expected from LAMMPS-dump format')

        self.input_format = fmt


    def set_output_format(self, fmt):
       
        if not fmt:
            raise ValueError('output_format should be specified')
        if fmt not in ['mtp', 'ase', 'xyz']:
            raise ValueError('output_format shoud be either "mtp", "ase" or "xyz", but I got "{}"'.format(fmt))

        if fmt == self.input_format:
            raise Exception('Input file {} already has {} format. Why bother converting it?!?'.format(self.fname, fmt))

        self.output_format = fmt


    def set_dbname(self, dbname):

        if dbname is not None:
            if self.output_format == 'mtp' and os.path.splitext(dbname)[1] != '.cfg':
                self.dbname = '{}.cfg'.format(os.path.splitext(dbname)[0])
            elif self.output_format == 'ase' and os.path.splitext(dbname)[1] != '.db':
                self.dbname = '{}.db'.format(os.path.splitext(dbname)[0])
            elif self.output_format == 'xyz' and os.path.splitext(dbname)[1] != '.xyz':
                self.dbname = '{}.xyz'.format(os.path.splitext(dbname)[0])
            else:
                self.dbname = dbname

        else:
            root = os.path.splitext(self.fname)[0]
            if self.output_format == 'mtp':
                self.dbname = '{}.cfg'.format(root)
            elif self.output_format == 'ase':
                self.dbname = '{}.db'.format(root)
            elif self.output_format == 'xyz':
                self.dbname = '{}.xyz'.format(root)


    def convert_database(self):

        if self.output_format == 'mtp':
            newdb = MtpDbCreator(self.dbname, overwrite=self.overwrite)
        elif self.output_format == 'ase':
            newdb = AseDbCreator(self.dbname, overwrite=self.overwrite)
        elif self.output_format == 'xyz':
            newdb = XyzDbCreator(self.dbname, overwrite=self.overwrite)

        out = newdb.create_database()

        # read structures into abipy format, energy, forces, stresses
        if self.input_format == 'mtp':
            data = MtpDbReader(self.fname, self.atomic_numbers)
        elif self.input_format == 'ase':
            data = AseDbReader(self.fname)
        elif self.input_format == 'xyz':
            data = XyzDbReader(self.fname)
        elif self.input_format == 'dump':
            data = DumpDbReader(self.fname, self.atomic_numbers)

        data.load_database()
        # for now, we should have structures in abistruct format, energy, forces and stresses

            # ASE : will need to check. Get the number of configs and find how to loop on them
            # MTP : easier to break the text into single configs? (as a text object or a rewritten temp text file?)
        for idx, struct in enumerate(data.structures[self.start::self.every]):
            atoms = newdb.convert_structure(struct)
            newdb.add_to_database(out, atoms, data.energy[idx], data.forces[idx], data.stresses[idx])
