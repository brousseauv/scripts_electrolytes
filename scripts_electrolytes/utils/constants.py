import scipy.constants as cst

# Hartree to eV conversion
ha_to_ev = cst.physical_constants['Hartree energy in eV'][0]

# Bohr radius in Angstrom
bohr_to_ang = cst.physical_constants['Bohr radius'][0]/cst.angstrom
ang_to_bohr = 1./bohr_to_ang

# Atomic unit of time
atomic_timeunit = cst.physical_constants['atomic unit of time'][0]

# Hartree energy in eV
hartree_energy_ev = cst.physical_constants['Hartree energy in eV'][0]

# Boltzmann constant in eV/K
boltzmann_evK = cst.physical_constants['Boltzmann constant in eV/K'][0]

# Conversion factors for stress units
evang3_to_gpa = 160.21766208
gpa_to_evang3 = 1./evang3_to_gpa
