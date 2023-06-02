from abipy.abilab import Structure
import json
from abipy.abio.inputs import AbinitInput

def poscar_to_abivars(vasp_fname):

    structure = Structure.from_file(vasp_fname)
    out = structure.to_abivars()
    # or maybe structure.to_abistring?
    # this does not work as MTP does not know the atomic species...
    return out

def load_abivars(fname):

    return json.load(open(fname))

def input_from_dict(struct, myvars):
    # write an abinit input file from a dictionnaryo
    #with open('input.abi', 'w') as f:
    #    for key in data.keys():
    #        f.write('{} {}\n'.format(key, data[key]))
    pseudos = myvars['pseudos']
    pp_dirpath = myvars['pp_dirpath']
    myvars.pop('pseudos', None)
    myvars.pop('pp_dirpath', None)
    myvars['indata_prefix'] = None
    myvars['outdata_prefix'] = None
    myvars['tmpdata_prefix'] = None

    abi_input = AbinitInput(struct, pseudos, pseudo_dir=pp_dirpath, abi_kwargs=myvars)
    abi_input.write()

