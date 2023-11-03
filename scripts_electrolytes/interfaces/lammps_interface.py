import numpy as np
import pandas as pd
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.md.analysis import DiffusionCoefficient
import os

''' Some functions to treat the outputs from a LAMMPS run'''

def extract_thermo(fname, out='thermo.dat'):

    f = open(fname, 'r')
    g = open(out, 'w')
    start = False
    
    for line in f.readlines():
        if line.find('Loop time of') != -1:
            break
        if start:
            g.write(line)
        if line.find('Per MPI rank') != -1:
            start = True

    return out

def create_thermo_dataframe(fname):

    cols = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, nrows=1).columns
    df = pd.read_csv(fname, delimiter=' ', skipinitialspace=True, usecols=cols[:-1])

    return df


def read_thermo(data, key):

    data[key] = pd.to_numeric(data[key])
    return data[key].values


def read_msd_from_thermo(data):

    time = read_thermo(data, 'Time')
    step = read_thermo(data, 'Step')[:2]
    msd = read_thermo(data, 'c_msd[4]')
    temp = read_thermo(data, 'Temp')[0]
    timestep = (time[1]-time[0])/(step[1]-step[0])

    return time, timestep, msd, temp


def read_traj_from_dump(fname, atomic_numbers, which=':'):

    traj= ase_read(fname, format='lammps-dump-text', index=which)

    # as the lammps dump outputs only atom id (1,2,3...) and not type, ASE sees H, He, Li...
    # So, convert atom id to atomic masses
    for frame in traj:
        for a in range(len(frame.numbers)):
            frame.numbers[a] = atomic_numbers[frame.numbers[a]-1]
    return traj


def read_config_from_dump(fname, atomic_numbers, which=-1):

    atoms = ase_read(fname, format='lammps-dump-text', index=which)

    # as the lammps dump outputs only atom id (1,2,3...) and not type, ASE sees H, He, Li...
    # So, convert atom id to atomic masses
    for a in range(len(atoms.numbers)):
        atoms.numbers[a] = atomic_numbers[atoms.numbers[a]-1]
    return atoms

def read_neb_logfile(fname, rescale_energy):

    ''' Reads a log.lammps main log output file and extracts the converged results '''
    f = open(fname, 'r')
    lines = f.readlines()
    data = lines[-1]

    data = data.split()
    forward_barrier = float(data[6])
    backward_barrier  = float(data[7])
    reaction_coordinate_length = float(data[8])
    reaction_coordinate = np.array(data[9::2], dtype=float)
    energy = np.array(data[10::2], dtype=float)
    if rescale_energy:
        energy -= energy[0]  # Set the 0 of energy at the initial configuration

    return forward_barrier, backward_barrier, reaction_coordinate_length, reaction_coordinate, energy


def read_natoms(fname):
    # This reads the number of atoms in a lammps simulation from a log.lammps.X file
    # where the initial structure was read from a .lmp file
    f = open(fname, 'r')
    lines = f.readlines()

    found = False
    for line in lines:
        if found:
            natoms = int(line.split('atoms')[0])
            break
        else:
            if line.find('reading atoms') != -1:
                found = True
    return natoms


def slice_trajectory_from_dump(fname, out_rootname='traj', atomic_numbers=None, nskip=10):

    ''' Reads every nskip configuration of a dump trajectory file and writes it in .xyz format '''
    if not atomic_numbers:
        raise Exception('Must provide a list of atomic numbers')
    if not isinstance(nskip, int):
        raise TypeError('nskip should be an integer, but I got a {}'.format(type(nskip)))
    if not os.path.exists(fname):
        raise ValueError('File {} does not exist'.format(fname))

    which = '::{}'.format(nskip)

    trajectory = read_traj_from_dump(fname, atomic_numbers, which=which)

    out_fname = '{}.xyz'.format(out_rootname)
    ase_write(filename=out_fname, images=trajectory, append=True)


def abistruct_to_xyz(db, struct, energy=None, forces=None, stresses=None):

    '''add configuration in extended XYZ format'''
    db.write("{}\n".format(struct.num_sites))
    
    message = set_xyz_message(struct.lattice.matrix, energy, forces, stresses)
    db.write("{}\n".format(message))

    if forces is not None:
        for idx, site in enumerate(struct):
            db.write("{0} {1[0]} {1[1]} {1[2]} {2[0]} {2[1]} {2[2]}\n".format(
                '{:2s}'.format(site.specie),
                     ['{:16.8f}'.format(x) for x in site.coords],
                     ['{:16.8f}'.format(f) for f in forces[idx, :]]))
    else:
        for idx, site in enumerate(struct):
            db.write("{0} {1[0]} {1[1]} {1[2]}\n".format(
                '{:2s}'.format(site.specie),
                     ['{:16.8f}'.format(x) for x in site.coords]))



def set_xyz_message(latt, energy, forces, stresses):

    lattice = 'Lattice="'
    for i in range(3):
        lattice = ' '.join([lattice, '{0[0]} {0[1]} {0[2]}'.format(['{:.8f}'.format(a) for a in latt[i, :]])])
    lattice = lattice+'"'

    if forces is not None:
        properties = 'Properties=species:S:1:pos:R:3:forces:R:3'
    else:
        properties = 'Properties=species:S:1:pos:R:3'

    if energy is not None:
        ene = 'energy="{:.8f}"'.format(energy)
    if stresses is not None:
        mystress = 'stress="{0[0]} {0[5]} {0[4]} {0[5]} {0[1]} {0[3]} {0[4]} {0[3]} {0[2]}"'.format(['{:.8f}'.format(strs) for strs in stresses])

    if energy is not None and stresses is not None:
        data = ' '.join([ene, mystress])
    elif energy is not None:
        data = ene
    elif stresses is not None:
        data = mystress

    if energy is not None or stresses is not None:
        message = ' '.join([lattice, properties, data])
    else:
        message = ' '.join([lattice, properties])

    return message



