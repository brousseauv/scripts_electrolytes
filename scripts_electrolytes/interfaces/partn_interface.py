
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