import argparse
from scripts_electrolytes.database.db_converter import DbConverter


'''
    This scripts converts a databse of atomic configurations from a given format to another, 
    keeping information about energy, forces and stresses.

    Simply call python dbconverter.py --<OPTION1> <value1> -- <OPTION2> <value2> etc.

    Options: 

        in_dbname (required): filename of the database to be converted

        input_format(required): Input format of the database.  Can be either 'mtp', 'ase', 'dump' or 'xyz'.

        output_format(required): Output format of the database.  Can be either 'mtp', 'ase' or 'xyz'.

        out_dbname: filename of the converted database 
                   Default: strips indbname and replaces the file extension

        overwrite: Boolean; indicates if thei converted database should be overwritten in case the filename already exists.
                   Default = False

        atomic_numbers: List of integer atomic numbers (ordered as the species type in the dump file), for dump file conversion.
                        Default: None

        start: Integer, index of the first configuration to convert
               Default: 0

        every: Integer, convert every "Every" configuration (i.e. data.structures[start::every])
               Default: 1 (all)

        ex: the following command converts a database called mydatabase in .cfg format to .xyz format

            python dbconverter.py --in_dbname mydatabase.cfg  --input_format='mtp' --output_format='xyz'

    For help about these options on the command line, type 
        python dbcreator.py --help

'''
def create_parser():

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--in_dbname", default=None, help="Input database name", required=True)
    parser.add_argument("--out_dbname", default=None, help="Output database name")

    parser.add_argument("--input_format", choices=['ase', 'mtp', 'xyz', 'dump', 'ncdump'], help="Input format for the database", required=True)
    parser.add_argument("--output_format", choices=['ase', 'mtp', 'xyz'], help="Output format for the database", required=True)
    parser.add_argument("--overwrite", type=bool, default=False, help="Should an existing database be overwritten or not")
    parser.add_argument("--atomic_numbers", default=None, type=list_of_integers, help="List of atomic numbers, ordered as in dump/cfg file")
    parser.add_argument("--start", type=int, default=0, help="Index of the first configuration to convert")
    parser.add_argument("--every", type=int, default=1, help="Convert every 'Every' configuration in the db")

    return parser


def list_of_integers(arg):
    return list(map(int, arg.split(',')))

def check_parser(args, parser):

    if args.input_format == 'dump' and not args.atomic_numbers:
        parser.error('Must provide a list for atomic_numbers')
    if args.input_format == 'ncdump' and not args.atomic_numbers:
        parser.error('Must provide a list for atomic_numbers')
    if args.input_format == 'mtp' and not args.atomic_numbers:
        parser.error('Must provide a list for atomic_numbers')


def main(args):

    if args.atomic_numbers is not None:
        db = DbConverter(fname=args.in_dbname, input_format=args.input_format, output_format=args.output_format, dbname=args.out_dbname, 
                         overwrite=args.overwrite, atomic_numbers=args.atomic_numbers, every=args.every, start=args.start)
    else:
        db = DbConverter(fname=args.in_dbname, input_format=args.input_format, output_format=args.output_format, dbname=args.out_dbname, 
                         overwrite=args.overwrite, every=args.every, start=args.start)

    db.convert_database()


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    check_parser(args, parser)
    main(args)
