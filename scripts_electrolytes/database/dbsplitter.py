import argparse
import numpy as np
from ase.db import connect

'''
    These classes split a given database in ASE .db or MTP .cfg format in two distinct databases.
    The size of the first splitted database can be specified either from the number of configurations
    or the fraction of the total number of configurations.

    Syntax:
        python dbsplitter.py <dbname> --<OPTION> <value>

    dbname: Path to the database to be splitted.

    Options:
        
        nsplit (required, mutually exclusive with split_fraction): Integer, number of configurations in the 
                first splitted database.

        splitfraction (required, mutually exclusive with nsplit): Float between 0 and 1, fraction of the 
                       total number of configurations for the first database, rounded to the lowest integer.

        appendtxt: String containing the names that should be appended to the splitted databases.
                   Should be a single string, with the two entries separated by a space.
                   Default: "1 2"
        
    For help about these options on the command line, type 
        python dbsplitter.py --help
'''
class DbSplitter:

    def __init__(self, dbname, nsplit, split_fraction, appendtxt):

        self.dbname = dbname
        self.nsplit = int(nsplit) if nsplit is not None else None
        self.split_fraction = float(split_fraction) if split_fraction is not None else None
        self.appendtxt = appendtxt


class AseDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction, appendtxt):

        super(AseDbSplitter, self).__init__(dbname, nsplit, split_fraction, appendtxt)


    def split_data(self):

        split = self.appendtxt.split(' ')
        out_db1, out_db2 = (self.dbname.split('.db')[0] + '_{}.db'.format(split[0]), 
                            self.dbname.split('.db')[0] + '_{}.db'.format(split[1]))

        with connect(self.dbname) as db, connect(out_db1) as db1, connect(out_db2) as db2:

            if self.split_fraction:
                ndata = int(np.floor(self.split_fraction*db.count()))
            else:
                ndata = self.split()
            idx = np.arange(1, db.count() + 1)
            np.random.shuffle(idx)
            idx1, idx2 = idx[:ndata], idx[ndata:]

            for idx in idx1:
                row = db.get(id=idx.item())
                db1.write(row)
            for idx in idx2:
                row = db.get(id=idx.item())
                db2.write(row)


class MtpDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction, appendtxt):

        super(MtpDbSplitter, self).__init__(dbname, nsplit, split_fraction, appendtxt)


    def split_data(self):

        token = 'BEGIN_CFG'
        configs = []
        current_config = []

        for line in open(self.dbname).readlines():
            if line.startswith(token) and current_config:
                configs.append(current_config)
                current_config = []
            current_config.append(line)
        configs.append(current_config)

        split = self.appendtxt.split(' ')
        out_db1, out_db2 = (self.dbname.split('.cfg')[0] + '_{}.cfg'.format(split[0]), 
                            self.dbname.split('.cfg')[0] + '_{}.cfg'.format(split[1]))

        if self.split_fraction:
            ndata = int(np.floor(self.split_fraction*len(configs)))
        else:
            ndata = self.nsplit

        idx = np.arange(0, len(configs))
        np.random.shuffle(idx)
        idx1, idx2 = idx[:ndata], idx[ndata:]

        cfg1, cfg2 = (open(out_db1, 'w'), 
                      open(out_db2, 'w'))
        for idx in idx1:
            self.write_config(cfg1, configs[idx]) 
        for idx in idx2:
            self.write_config(cfg2, configs[idx]) 

        cfg1.close()
        cfg2.close()


    def write_config(self, g, config):

        for line in config:
            g.write(line)

def create_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument("dbname", help="path to the database file to be splitted")
    data = parser.add_mutually_exclusive_group(required=True)
    data.add_argument("--nsplit", type=int, help="number of entries in the 1st extracted database", default=None)
    data.add_argument("--fsplit", type=float, help="fraction of entries in the 1st extracted database, between 0 and 1",
                      default=None)
    parser.add_argument("--appendtxt", default="1 2", help="Text to be appended to the splitted databases. Single string separated by a space")
    return parser

def check_parser(args, parser):

    if args.nsplit is None and args.fsplit is None:
        parser.error("--nsplit or --fsplit should be passed as arguments")

    if args.fsplit is not None:
        if args.fsplit > 1.0 or args.fsplit < 0.0:
            parser.error("--fsplit should be a float between 0 and 1")

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)

    if args.dbname.endswith(".db"):
        dbspl = AseDbSplitter(args.dbname, args.nsplit, args.fsplit, args.appendtxt)    
    elif args.dbname.endswith(".cfg"):
        dbspl = MtpDbSplitter(args.dbname, args.nsplit, args.fsplit, args.appendtxt)    
    dbspl.split_data() 
