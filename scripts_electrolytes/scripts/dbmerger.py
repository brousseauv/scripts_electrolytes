import argparse
from scripts_electrolytes.database.db_merger import MtpDbMerger, AseDbMerger, XyzDbMerger


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
def create_parser():

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--merged_dbname", default=None, required=True, help="Merged database name (new name or name of first database in filenames)")
    parser.add_argument("--filenames", default=None, type=list_of_strings, help="List of filenames to be merges, --filenames=db1,db2,db3,etc.", required=True)

    parser.add_argument("--format", choices=['ase', 'mtp', 'xyz'], help="Input/output format of the databases", required=True)
    parser.add_argument("--append", type=bool, default=False, help="Should data be appended to existing database or not")
    parser.add_argument("--atomic_numbers", default=None, type=list_of_integers, help="List of atomic numbers, ordered as in cfg file")
    return parser


def check_parser(args, parser):

    if args.format == 'mtp' and args.atomic_numbers is None:
        parser.error("--format 'mtp' requires --atomic_numbers argument")


def list_of_integers(arg):
    return list(map(int, arg.split(',')))


def list_of_strings(arg):
    return list(map(str, arg.split(',')))


def main(args):

    if args.format == 'mtp':
        db = MtpDbMerger(merged_dbname=args.merged_dbname, filenames=args.filenames, append=args.append, atomic_numbers=args.atomic_numbers)
    elif args.format == 'ase':
        db = AseDbMerger(merged_dbname=args.merged_dbname, filenames=args.filenames, append=args.append)
    elif args.format == 'xyz':
        db = XyzDbMerger(merged_dbname=args.merged_dbname, filenames=args.filenames, append=args.append)

    db.merge_db()


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)
    main(args)
