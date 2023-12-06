
def fix_species_in_xyz_mlip(fname, symbols):
   ''' When ran with LAMMPS and MLIP2 potentials, the .xyz files
       produced by pARTn does dot list atomic species by text, but 
       as a species index from 1 to nspecies.
       This function replaces these atomic indices by the correct 
       atomic string, provided as symbols (in list format).
   '''

   # for pARTn, there should be only one config in a sadX or minX file
   # so, no need to bother about multi-configurations for now.
   with open(fname, 'r+') as f:
       content = f.readlines()
       try:
           content[1] = content[1].replace('species:I', 'species:S')
       except:
           pass

       # Test if the first character is an integer or a letter
       try:
           isinstance(int(content[2][1]), int)

       except:
           pass
       else:
           for i, line in enumerate(content[2:]):
               content[2+i] = line.replace(line[1] , symbols[int(line[1])-1].capitalize(), 1)


       f.seek(0)
       f.writelines(content)


def split_xyz_configs(fname):
    
    token = 'Lattice='
    configs = []
    token_idx = []

    lines = open(fname).readlines()

    for iline, line in enumerate(lines):
        if line.startswith(token):
            token_idx.append(iline)

    for i, idx in enumerate(token_idx):
        try:
            # Read current block up to one line before next "Lattice" token
            current_config = lines[idx-1:token_idx[i+1]-1]
        except:
            # Last config has no next "Lattice" token
            current_config = lines[idx-1:]
        configs.append(current_config)
    return configs

