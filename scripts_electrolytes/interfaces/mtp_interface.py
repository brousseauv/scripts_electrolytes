#!/usr/bin/env python
import numpy as np
import subprocess as subp
import re
import pandas as pd
from ..utils.constants import ang_to_bohr
from typing import Any, Dict, List, Optional, TextIO, Tuple
from collections import defaultdict

def abistruct_to_cfg(db, struct, energy=None, forces=None, stresses=None):

    # Write configuration in the .cfg format from MLIP-2 package
    # all the properties related to struct object ( an Abipy.core.Structure object)
    db.write("BEGIN_CFG\n")
    db.write(" Size\n    {}\n".format(struct.num_sites))
    db.write(" Supercell")
    for i in range(3):
        db.write("\n    {0[0]} {0[1]} {0[2]}".format(['{:13.6f}'.format(a) for a in struct.lattice.matrix[i, :]]))
        #  FIX ME: Check units in MLIP!!! Most likely angstrom as they work with VASP...
    db.write("\n")

    if forces is not None:
        db.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for idx, site in enumerate(struct):
            db.write("    {:10.0f} {:4.0f}".format(idx + 1, struct.types_of_species.index(site.specie)))
            db.write("  {0[0]} {0[1]} {0[2]}".format(['{:13.6f}'.format(x) for x in site.coords]))
            db.write("  {0[0]} {0[1]} {0[2]}\n".format(['{:11.6f}'.format(f) for f in forces[idx, :]]))
    else:
        db.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z\n")
        for idx, site in enumerate(struct):
            db.write("    {:10.0f} {:4.0f}".format(idx + 1, struct.types_of_species.index(site.specie)))
            db.write("  {0[0]} {0[1]} {0[2]}\n".format(['{:13.6f}'.format(x) for x in site.coords]))

    if energy is not None:
        db.write(" Energy\n    {:20.12f}\n".format(energy))
    if stresses is not None:
        db.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        db.write("     {0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]}\n".format(['{:11.5f}'.format(strs) for strs in stresses]))
    db.write(" Feature   ESF_by\tABINIT\n")
    # FIX ME: not sure this will be working during reading?!? if not, for now, change to vasp...
    db.write("END_CFG\n\n")


def read_errors(fname, runtype, version='mlip2'):

    check_runtype(runtype)

    if version not in ['mlip2', 'mlip3']:
        raise ValueError('"version" should either be "mlip2" or "mlip3" but I got {}'.format(version))

    try:
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
    except FileNotFoundError:
        print('File {} not found.'.format(fname))

    train_block = []
    valid_block = []

    if version == 'mlip2':

        train_block_on = False
        valid_block_on = False

        for line in lines:
            if 'TRAIN ERRORS' in line:
                train_block_on = True
                valid_block_on = False
            if 'VALIDATION ERRORS' in line:
                train_block_on = False
                valid_block_on = True

            if train_block_on:
                train_block.append(line)
            elif valid_block_on:
                valid_block.append(line)

    elif version == 'mlip3':
        if runtype == 'all':
            raise ValueError('With MLIP3, cannot read train errors and valid errors in the same file.')
        if runtype == 'train':
            train_block = []
            for line in lines:
                train_block.append(line)
        if runtype == 'valid':
            valid_block = []
            for line in lines:
                valid_block.append(line)
       
    if not train_block and runtype == 'train':
        raise Exception('Train errors block is empty')
    if not valid_block and runtype == 'valid':
        raise Exception('Validation errors block is empty')

    if runtype == 'train' or runtype == 'all':
        train_block = split_block(train_block, version)
        train_data = read_block(train_block)

    if runtype == 'valid' or runtype == 'all':
        valid_block = split_block(valid_block, version)
        valid_data = read_block(valid_block)

    if runtype == 'train':
        return train_data
    elif runtype == 'valid':
        return valid_data
    else:
        return train_data, valid_data


def check_runtype(runtype):
    ''' Checks the runtype keyword is correct'''
    keys = ['train', 'valid', 'all']
    if runtype not in keys:
        raise ValueError('Runtype should be {} but I received {}'.format(keys, runtype))


def split_block(block, version):
    ''' Splits an Errors report block into strings from individual variables'''

    block = "".join(block)
    split_block = []

    block = block.split('Energy:')[1]
    block = block.split('Energy per atom:')
    split_block.append(block[0])
    block = block[1].split('Forces:')
    split_block.append(block[0])

    if version == 'mlip2':
        block = block[1].split('Stresses (in eV):')
    elif version == 'mlip3':
        block = block[1].split('Stresses (in energy units):')
    split_block.append(block[0])
    if version == 'mlip2':
        block = block[1].split('Stresses (in GPa):')
    elif version == 'mlip3':
        block = block[1].split('Virial stresses (in pressure units):')
    block = block[1].split('_______________________________________________')
    split_block.append(block[0])

    return split_block


def read_block(block):
    ''' Extracts error data from splitted Errors report block'''

    if len(block) != 4:
        raise ValueError('Some error data is missing from the output file.')

    keys = ['energy', 'energy_per_atom', 'forces', 'stresses']
    data = {}
    for i, iblock in enumerate(block):
        data_block = []
        for line in iblock.split('\n'):
            if '=' in line:
                data_block.append(float(line.split('=')[-1]))
        data[keys[i]] = data_block

    return data


def read_mv_grade(fname, verbose=False):
    ''' Reads the MV_grade from a .cfg database and returns a dictionnary containing
        extrapolation grade values, and some statistics about the database
    '''

    data = {}
    status, output = subp.getstatusoutput('grep MV_grade {}'.format(fname))
    output = np.asarray([i.split('\t')[-1] for i in output.split('\n')], dtype=float)

    data['values'] = output
    data['max'] = max(output)
    data['min'] = min(output)
    data['argmax'] = np.argmax(output)
    data['argmin'] = np.argmin(output)
    data['mean'] = np.mean(output)
    data['median'] = np.median(output)
    data['stdev'] = np.std(output)

    if verbose:
        print('For file:{}'.format(fname))
        print('    gamma max = {:.2f} (step {})'.format(data['max'], data['argmax']))
        print('    gamma min = {:.2f} (step {})'.format(data['min'], data['argmin']))
        print('    gamma mean = {:.2f}'.format(data['mean']))
        print('    gamma stdev = {:.2f}'.format(data['stdev']))
        print('    gamma median = {:.2f}'.format(data['median']))

    return data


def read_nbh_mv_grade(fname, verbose=False):
    ''' Reads the atom and global (max) MV_grade from a .cfg database and returns a dictionnary containing
        extrapolation grade values, and some statistics about the database
    '''

    data = {}
    status, output = subp.getstatusoutput('grep MV_grade {}'.format(fname))
    output = np.asarray([i.split('\t')[-1] for i in output.split('\n')], dtype=float)

    data['values'] = output
    data['max'] = max(output)
    data['min'] = min(output)
    data['argmax'] = np.argmax(output)
    data['argmin'] = np.argmin(output)
    data['mean'] = np.mean(output)
    data['median'] = np.median(output)

    data['stdev'] = np.std(output)

    if verbose:
        print('For file:{}'.format(fname))
        print('    gamma max = {:.2f} (step {})'.format(data['max'], data['argmax']))
        print('    gamma min = {:.2f} (step {})'.format(data['min'], data['argmin']))
        print('    gamma mean = {:.2f}'.format(data['mean']))
        print('    gamma stdev = {:.2f}'.format(data['stdev']))
        print('    gamma median = {:.2f}'.format(data['median']))

    # Treat individual atom grade
    configs = split_cfg_configs(fname)
    nconfig = len(configs)

    for c, chunk in enumerate(configs):
        # strip atom lines
        # add to columns
        # concat df

        natom = read_natom(chunk) 
        energy, stresses, idx = read_config_properties(chunk)
        atomdata, header = read_atomdata(chunk, idx)
        split_header = re.sub('\s{2,}', ' ', re.sub('\n', '', header)).split(': ')[1].split(' ')
        split_atomdata = [re.sub('\s{2,}', ' ', re.sub('\n', '', row)) for row in atomdata]

        current_df = pd.DataFrame(data=[row.split() for row in split_atomdata], columns=split_header)
        current_df['cfg_index'] = c 
        try:
            atom_df = pd.concat([atom_df, current_df])
        except:
            atom_df = current_df.copy()

    atom_df.reset_index(drop=True, inplace=True)
    return data, atom_df

def split_cfg_configs(fname):

    token = 'BEGIN_CFG'
    configs = []
    current_config = []

    for line in open(fname).readlines():
        if line.startswith(token) and current_config:
            configs.append(current_config)
            current_config = []
        current_config.append(line)
    configs.append(current_config)
    return configs


def convert_chunk_to_abivars(data, atomic_numbers):

    ''' Converts a .cfg configuration into an abivars dict, and read EFS if applicable '''

    natom = read_natom(data) 

    lattice = read_lattice(data)

    energy, stresses, idx = read_config_properties(data)

    atomdata, header = read_atomdata(data, idx)

    typat, xcart, forces = split_atomic_data(atomdata, header, natom, lattice) 

    abivars = {'natom': natom,
               'ntypat': len(atomic_numbers),
               'znucl': atomic_numbers,
               'typat': typat,
               'acell': np.ones((3)),
               'rprim': lattice * ang_to_bohr,
               'xangst': xcart
              }

    return abivars, energy, forces, stresses


def read_natom(data):

    idx = data.index(' Size\n')
    return int(data[idx+1].split('\n')[0])


def read_lattice(data):

    idx = data.index(' Supercell\n')
    latt = np.zeros((3,3))
    
    for i in range(3):
        latt[i, :] = re.sub(r"\s+", " ", data[idx+i+1].strip()).split(' ')
    return latt


def read_config_properties(data):

    try:
        idx = data.index(' Energy\n')
        energy = data[idx+1].strip()
        my_idx = idx
    except:
        energy = None

    try:
        idx = [i for i, item in enumerate(data) if re.search('PlusStress:', item)][0]
        stress = np.asarray(re.sub(r"\s+", " ", data[idx+1].strip()).split(' '), dtype=np.double)
        try:
            my_idx
        except:
            my_idx = idx
    except:
        stress = None

    try:
        my_idx
    except:
        my_idx = [i for i, item in enumerate(data) if re.search('Feature', item)][0]

    return energy, stress, my_idx


def read_atomdata(data, idx):

    start_idx = [i for i, item in enumerate(data) if re.search('AtomData:', item)][0]
    header = data[start_idx]

    return data[start_idx+1:idx], header


def split_atomic_data(data, header, natom, latt):

    if len(data) != natom:
        raise ValueError('AtomData does not contain natom={} lines. Something went wrong, check your data.'.format(natom))

    # check header to see if forces are present
    if header.find('fx') != -1:
        forces = np.zeros((natom, 3), dtype=float)
    else:
        forces = None

    # check header to see which kind of coordinates are printed
    if header.find('cartes_x') != -1:
        pos_are_cart = True
    elif header.find('direct_x') != -1:
        pos_are_cart = False
        latt_inv = np.linalg.inv(latt.T)
    else:
        raise ValueError('Could not find either cartesian or reduced coordinates in AtomData header. Check your data.')


    header = re.sub(r'\s+', ' ', header.strip().split('AtomData:')[1].lstrip()).split(' ')

    # find index for atom type, coordinates (and forces)
    idx_typat = header.index('type')
    if pos_are_cart:
        idx_pos = header.index('cartes_x')
    else:
        idx_pos = header.index('direct_x')

    if forces is not None:
        idx_forces = header.index('fx')


    typat = np.zeros((natom), dtype=int)
    pos = np.zeros((natom, 3), dtype=float)
    # Retrieve atom type, positions and forces from AtomData block
    for a, atom in enumerate(data):
        atom = re.sub(r'\s+', ' ', atom.strip()).split(' ')
        typat[a] = int(atom[idx_typat])+1

        if pos_are_cart:
            pos[a, :] = atom[idx_pos:idx_pos+3]
        else:
            vec = np.asarray(atom[idx_pos:idx_pos+3])
            pos[a, :] = np.matmul(latt_inv, vec)

        if forces is not None:
            forces[a, :] = atom[idx_forces:idx_forces+3]

    return  typat, pos, forces

def read_cfgs_with_nbh_grade(filename: str, nbh_grade: bool=True, elements: list=None) -> pd.DataFrame:
    '''
    Read cfg file with atomic neighborhood MaxVol grade, 
    and convert to dataframe
    Adapted from diffusion_for_multi_scale_molecular_dynamics.crystal_diffusion.models.mtp.read_cfgs
    to account for cfg files containing numbered atomic index instead of species
    '''

    data_pool = []
    with open(filename, 'rt') as f:
        lines = f.read()

    block_pattern = re.compile("BEGIN_CFG\n(.*?)\nEND_CFG", re.S)
    size_pattern = re.compile("Size\n(.*?)\n SuperCell", re.S | re.I)
    energy_pattern = re.compile("Energy\n(.*?)\n (?=PlusStress|Stress)", re.S)

    for block in block_pattern.findall(lines):
        d = {"outputs": {}}
        size_str = size_pattern.findall(block)[0]
        size = int(size_str.lstrip())

        if nbh_grade:
            position_pattern = re.compile("nbh_grades\n(.*?)\n Energy", re.S)
            try:
                position_str = position_pattern.findall(block)[0]
            except:
                position_pattern = re.compile("nbh_grades\n(.*?)\n Feature", re.S)
                position_str = position_pattern.findall(block)[0]
        else:
            position_pattern = re.compile("fz\n(.*?)\n Energy", re.S)

        position = np.array(list(map(formatify, position_str.split("\n"))))
        species = np.array([elements[j] for j in position[:, 1].astype(np.int64)])

        if len(position[0])>6:
            forces = position[:, 5:8].tolist()
        else:
            forces = None

        try:
            energy_str = energy_pattern.findall(block)[0]
            energy = float(energy_str.lstrip())
        except:
            energy = None

        d["outputs"]["energy"] = energy
        d["num_atoms"] = size
        d["outputs"]["position"] = position[:, 2:5].tolist()
        d["outputs"]["forces"] = forces
        d["outputs"]["species"] = species

        if nbh_grade:
            nbh_grade_values = position[:, -1].tolist()
            d["outputs"]["nbh_grades"] = nbh_grade_values

        data_pool.append(d)
    df = convert_to_dataframe(docs=data_pool)
    return df

def convert_to_dataframe(docs: List[Dict[str, Any]]) -> pd.DataFrame:
    '''
    Convert a list of docs into DataFrame usable for computing metrics and analysis.
    Taken from from diffusion_for_multi_scale_molecular_dynamics.crystal_diffusion.models.mtp
    Written hereonly to limit the size of virtual environments...
    '''
    df = defaultdict(list)

    for s_idx, d in enumerate(docs):
        n_atom = d['num_atoms']
        outputs = d["outputs"]
        pos_arr = np.array(outputs["position"])
        assert n_atom == pos_arr.shape[0], "Number of positions do not match number of atoms"
        if outputs['forces'] is not None:
            force_arr = np.array(outputs["forces"])
            assert n_atom == force_arr.shape[0], "Number of forces do not match number of atoms"

        for i, x in enumerate(['x', 'y', 'z']):
            df[x] += pos_arr[:, i].tolist()
            if outputs['forces'] is not None:
                df[f'f{x}'] += force_arr[:, i].tolist()
            else:
                df[f'f{x}'] += [None] * n_atom

        if outputs['energy'] is not None:
            df['energy'] += [outputs['energy']] * n_atom  # copy the value to all atoms
        else:
            df['energy'] += [None] * n_atom

        if "nbh_grades" in outputs.keys():
            nbh_grades = outputs["nbh_grades"]
            assert n_atom == len(nbh_grades), "Number of gamma values do not match number of atoms"
            df['nbh_grades'] += nbh_grades

        df['atom_index'] += list(range(n_atom))
        df['structure_index'] += [s_idx] * n_atom

    df = pd.DataFrame(df)
    return df

def formatify(string: str) -> List[float]:
        """Convert string to a list of float."""
        return [float(s) for s in string.split()]
