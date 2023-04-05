#!/usr/bin/env python

def abistruct_to_cfg(db, struct, energy, forces, stresses):

    # Write configuration in the .cfg format from MLIP-2 package
    # all the properties related to struct object ( an Abipy.core.Structure object)
    # need to be checked! I am writing that without checking actually how to call these.
    db.write("BEGIN_CFG\n")
    db.write(" Size\n    {}\n".format(struct.num_sites))
    db.write(" Supercell")
    for i in range(3):
        db.write("\n    {0[0]} {0[1]} {0[2]}".format(['{:13.6f}'.format(a) for a in struct.lattice.matrix[i, :]]))
        #  FIX ME: Check units in MLIP!!!
    db.write("\n")

    db.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
    for idx, site in enumerate(struct):
        db.write("    {:10.0f} {:4.0f}".format(idx, struct.types_of_species.index(site.specie)))  # FIX ME: it has to be typat, not symbol
        db.write("  {0[0]} {0[1]} {0[2]}".format(['{:13.6f}'.format(x) for x in site.coords]))
        db.write("  {0[0]} {0[1]} {0[2]}\n".format(['{:11.6f}'.format(f) for f in forces[idx, :]]))

    db.write(" Energy\n    {:20.12f}\n".format(energy))
    db.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
    db.write("     {0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]}\n".format(['{:11.5f}'.format(strs) for strs in stresses]))
    db.write(" Feature   ESF_by\tABINIT\n")
    # FIX ME: not sure this will be working during reading?!? if not, for now, change to vasp...
    db.write("END_CFG\n\n")
