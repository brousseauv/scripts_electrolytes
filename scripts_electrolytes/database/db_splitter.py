import numpy as np
import os
from ase.db import connect
from ..interfaces.mtp_interface import split_cfg_configs
from ..interfaces.partn_interface import split_xyz_configs

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

        seed: Value of the random seed, in case one wants to reproduce a given split/shuffle.
              Default: None
        
    For help about these options on the command line, type 
        python dbsplitter.py --help
'''
class DbSplitter:

    def __init__(self, dbname, nsplit, split_fraction, appendtxt, seed):

        self.dbname = dbname
        self.nsplit = int(nsplit) if nsplit is not None else None
        self.split_fraction = float(split_fraction) if split_fraction is not None else None
        self.appendtxt = appendtxt
        if seed is not None:
            self.seed = int(seed)
        else:
            self.seed = seed


class AseDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction, appendtxt, seed):

        super(AseDbSplitter, self).__init__(dbname, nsplit, split_fraction, appendtxt, seed)


    def split_data(self):

        split = self.appendtxt.split(' ')
        out_db1, out_db2 = (self.dbname.split('.db')[0] + '_{}.db'.format(split[0]), 
                            self.dbname.split('.db')[0] + '_{}.db'.format(split[1]))

        with connect(self.dbname) as db, connect(out_db1) as db1, connect(out_db2) as db2:

            if self.split_fraction:
                ndata = int(np.floor(self.split_fraction*db.count()))
            else:
                ndata = self.nsplit
            idx = np.arange(1, db.count() + 1)

            if self.seed:
                np.random.seed(self.seed)
            np.random.shuffle(idx)
            idx1, idx2 = idx[:ndata], idx[ndata:]

            for idx in idx1:
                row = db.get(id=idx.item())
                db1.write(row)
            for idx in idx2:
                row = db.get(id=idx.item())
                db2.write(row)


class MtpDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction, appendtxt, seed):

        super(MtpDbSplitter, self).__init__(dbname, nsplit, split_fraction, appendtxt, seed)


    def split_data(self):

        split = self.appendtxt.split(' ')
        out_db1, out_db2 = (os.path.basename(self.dbname).split('.cfg')[0] + '_{}.cfg'.format(split[0]), 
                            os.path.basename(self.dbname).split('.cfg')[0] + '_{}.cfg'.format(split[1]))

        configs = split_cfg_configs(self.dbname)

        if self.split_fraction:
            ndata = int(np.floor(self.split_fraction*len(configs)))
        else:
            ndata = self.nsplit

        idx = np.arange(0, len(configs))

        if self.seed:
            np.random.seed(self.seed)
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


class XyzDbSplitter(DbSplitter):

    def __init__(self, dbname, nsplit, split_fraction, appendtxt, seed):

        super(XyzDbSplitter, self).__init__(dbname, nsplit, split_fraction, appendtxt, seed)


    def split_data(self):

#        raise NotImplementedError('XyzDbSplitted not yet implemented.')
        split = self.appendtxt.split(' ')
        out_db1, out_db2 = (os.path.basename(self.dbname).split('.xyz')[0] + '_{}.xyz'.format(split[0]), 
                            os.path.basename(self.dbname).split('.xyz')[0] + '_{}.xyz'.format(split[1]))

        configs = split_xyz_configs(self.dbname)

        if self.split_fraction:
            ndata = int(np.floor(self.split_fraction*len(configs)))
        else:
            ndata = self.nsplit

        idx = np.arange(0, len(configs))

        if self.seed:
            np.random.seed(self.seed)
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

