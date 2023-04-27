from abipy.dynamics.hist import HistFile

''' Small utilities to work with Abinit output files'''


def extract_config_from_hist(fname, idx=0):

    hist = HistFile(fname)
    structs = hist.structures

    return structs[idx]


def convert_config(config, fmt, out_fname=None):

    new_config = config.to(fmt=fmt, filename=out_fname)

    
