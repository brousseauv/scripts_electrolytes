from abipy.abilab import Structure
import json
from abipy.abio.inputs import AbinitInput
from abipy.abio.abivars import is_abivar

def poscar_to_abivars(vasp_fname):

    structure = Structure.from_file(vasp_fname)
    out = structure.to_abivars()
    # or maybe structure.to_abistring?
    # this does not work as MTP does not know the atomic species...
    return out

def load_abivars(fname):

    # extract the pseudo variables from the dictionnary, to work with AbinitInput class
    abivars = json.load(open(fname))

    check_abivars(abivars)

    abipseudos = {'pseudos': abivars['pseudos'],
                  'pp_dirpath': abivars['pp_dirpath']}
    abivars.pop('pseudos', None)
    abivars.pop('pp_dirpath', None)

    abivars['indata_prefix'] = None
    abivars['outdata_prefix'] = None
    abivars['tmpdata_prefix'] = None

    # deactivate default outputs, except for gsr
    prtkeys = ["prtden", "prteig", "prtwf"]
    for key in prtkeys:
        if not key in abivars:
            abivars[key] = 0
    abivars['prtgsr'] = 1

    return abivars, abipseudos


def check_abivars(myvars):
    fixvars = []
    for key in myvars.keys():
        if not is_abivar(key):
            fixvars.append(key)
            #print('{} is not an Abinit variable'.format(key))
    if len(fixvars)>0:
        raise ValueError('{} are not Abinit variables'.format(fixvars))


def input_from_dict(struct, myvars, mypseudos):

    abi_input = AbinitInput(struct, mypseudos['pseudos'], pseudo_dir=mypseudos['pp_dirpath'], abi_kwargs=myvars)
    abi_input.write()


def abivars_to_abistruct(myvars):

    return Structure.from_abivars(myvars)
