import argparse
import numpy as np
from ase.db import connect

class DbSplitter:

    def __init__(self, dbname, nsplit, split_fraction):

        self.dbname = dbname
        self.nsplit = int(nsplit) if nsplit is not None else None
        self.split_fraction = float(split_fraction) if split_fraction is not None else None


class AseDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction):

        super(AseDbSplitter, self).__init__(dbname, nsplit, split_fraction)


    def split_data(self):

        out_db1, out_db2 = (self.dbname.strip('.db') + '_1.db', 
                            self.dbname.strip('.db') + '_2.db')

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

    def __init__(self, dbname, nsplit, split_fraction):

        super(MtpDbSplitter, self).__init__(dbname, nsplit, split_fraction)


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

        out_db1, out_db2 = (self.dbname.strip('.cfg') + '_1.cfg', 
                            self.dbname.strip('.cfg') + '_2.cfg')
       
        if self.split_fraction:
            ndata = int(np.floor(self.split_fraction*len(configs)))
        else:
            ndata = self.nsplit
        idx = np.arange(0, len(configs))
        np.random.shuffle(idx)
        idx1, idx2 = idx[:ndata], idx[ndata:]

        cfg1, cfg2 = (open(out_db1, 'w'), 
                      open(out_db2, 'w'))
        print(idx, ndata, len(configs))
        print(idx1)
        print(idx2)
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
    data.add_argument("--split_fraction", type=float, help="fraction of entries in the 1st extracted database, between 0 and 1",
                      default=None)
    return parser


if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()

    if args.dbname.endswith(".db"):
        dbspl = AseDbSplitter(args.dbname, args.nsplit, args.split_fraction)    
    elif args.dbname.endswith(".cfg"):
        dbspl = MtpDbSplitter(args.dbname, args.nsplit, args.split_fraction)    
    dbspl.split_data() 
