import argparse
from scripts_electrolytes.database.db_creator import MtpDbCreator, AseDbCreator


'''
    This scripts creates a databse of atomic configurations extracted from Abinit HIST.nc files
    or from a directory containing multiple GSR.nc files,
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
