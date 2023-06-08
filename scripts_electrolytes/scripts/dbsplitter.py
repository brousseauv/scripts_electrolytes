import argparse
from scripts_electrolytes.database.db_splitter import AseDbSplitter, MtpDbSplitter

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
