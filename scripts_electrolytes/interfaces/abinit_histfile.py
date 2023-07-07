from abipy.dynamics.hist import HistFile
from abipy.core.structure import Structure
from ..utils.constants import ha_to_ev, bohr_to_ang
from ..database.db_creator import MtpDbCreator, AseDbCreator
from scripts_electrolytes.utils.constants import ha_to_ev, bohr_to_ang

''' Small utilities to work with Abinit output files'''


def extract_config_from_hist(fname, idx=0):

    hist = HistFile(fname)
    structs = hist.structures

    return structs[idx]


def convert_config(config, fmt, out_fname=None):

    new_config = config.to(fmt=fmt, filename=out_fname)

    
def extract_neb_trajectory(fname, fmt=None, output='neb'):

    ''' Extract the relaxed NEB trajectory, as well as 
        energy - forces - stresses, and dump into a database
    '''

    if not fmt:
        raise ValueError('Database output format unspecified. Choose fmt = "mtp" or fmt = "ase".')

    hist = HistFile(fname)

    structures = read_relaxed_neb_structures(hist)
    energy, forces, stresses = read_neb_efs(hist)

    if fmt == 'mtp':
        if not output.endswith('.cfg'):
            output = '{}.cfg'.format(output)
        nebtraj = MtpDbCreator(dbname=output)
    elif fmt == 'ase':
        if not output.endswith('.db'):
            output = '{}.db'.format(output)
        nebtraj = AseDbCreator(dbname=output)
    else:
        raise ValueError('Unknown output format. Choose fmt = "mtp" or fmt = "ase".')

    db = nebtraj.create_database()

    for s, struct in enumerate(structures):
        atoms = nebtraj.convert_structure(struct)
        nebtraj.add_to_database(db, atoms, energy[s], forces[s, :, :], stresses[s, :])

    
def read_relaxed_neb_structures(data):
   
    ''' This is adapted from Abipy dynamics.hist class, suited for
        multiple images present in NEB calculations'''

    rprimd = data.reader.read_value('rprimd')
    xred = data.reader.read_value('xred')
    nimage = data.reader.read_dimvalue('nimage')

    num_pseudos = data.reader.read_dimvalue("npsp")
    ntypat = data.reader.read_dimvalue("ntypat")
    if num_pseudos != ntypat:
        raise NotImplementedError("Alchemical mixing is not supported, num_pseudos != ntypat")

    znucl, typat = data.reader.read_value("znucl"), data.reader.read_value("typat").astype(int)

    structures = []

    for img in range(nimage):
        s = Structure.from_abivars(
                xred=xred[-1, img],
                rprim=rprimd[-1, img],
                acell=3 * [1.0],
                znucl=znucl,
                typat=typat,
            )
        structures.append(s)

    return structures


def read_neb_efs(data):

    e = data.reader.read_value('etotal')[-1, :] * ha_to_ev  # in eV
    f = data.reader.read_value('fcart')[-1, :, :, :] * ha_to_ev/bohr_to_ang  # in eV/ang
    s = data.reader.read_value('strten')[-1, :, :] * ha_to_ev/(bohr_to_ang**3)  # in eV/ang^3

    return e, f, s
