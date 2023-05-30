#!/usr/bin/env python
import numpy as np
import subprocess as subp


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


def read_errors(fname, runtype):

    check_runtype(runtype)
    try:
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
    except FileNotFoundError:
        print('File {} not found.'.format(fname))

    train_block = []
    valid_block = []
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

    if not train_block and runtype == 'train':
        raise Exception('Train errors block is empty')
    if not valid_block and runtype == 'test':
        raise Exception('Validation errors block is empty')

    if runtype == 'train' or runtype == 'all':
        train_block = split_block(train_block)
        train_data = read_block(train_block)

    if runtype == 'valid' or runtype == 'all':
        valid_block = split_block(valid_block)
        valid_data = read_block(valid_block)

    if runtype == 'train':
        return train_data
    elif runtype == 'test':
        return valid_data
    else:
        return train_data, valid_data


def check_runtype(runtype):
    ''' Checks the runtype keyword is correct'''
    keys = ['train', 'test', 'all']
    if runtype not in keys:
        raise ValueError('Runtype should be {} but I received {}'.format(keys, runtype))


def split_block(block):
    ''' Splits an Errors report block into strings from individual variables'''

    block = "".join(block)
    split_block = []

    block = block.split('Energy:')[1]
    block = block.split('Energy per atom:')
    split_block.append(block[0])
    block = block[1].split('Forces:')
    split_block.append(block[0])
    block = block[1].split('Stresses (in eV):')
    split_block.append(block[0])
    block = block[1].split('Stresses (in GPa):')
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
    data['stdev'] = np.std(output)

    if verbose:
        print('For file:{}'.format(fname))
        print('    gamma max = {:.2f} (step {})'.format(data['max'], data['argmax']))
        print('    gamma min = {:.2f} (step {})'.format(data['min'], data['argmin']))
        print('    gamma mean = {:.2f}'.format(data['mean']))
        print('    gamma stdev = {:.2f}'.format(data['stdev']))

    return data
