#!/usr/bin/bash
#SBATCH --job-name md-LSN
#SBATCH -t 0:15:00
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -A rrg-mousseau-ac

mlp='/home/broussev/software/mlip-3/bin/mlp'
lmp='/home/broussev/software/interface-lammps-mlip-3/lmp_mpi'
temp=800

module restore mlipmodule_gnu

MATNAME="Li9S3N"
LAUNCHDIR=$PWD
CONFIGDIR='/home/broussev/projects/rrg-cotemich-ac/broussev/calculs/Li9S3N/configs/singlevacancy'


lammpsinput='lammps.in'
fname='li9s3n_2x2x2_singlevac.lmp'

mlippot="/home/broussev/projects/rrg-cotemich-ac/broussev/calculs/Li9S3N/mtp/passive/800K/mlip3_nebparams_alpha16/103atom/weights_e1_f10/pot5/pot5.almtp"
lammpslog="lammps.log"

echo "Computing MD run for T="$temp"K"
struct=$CONFIGDIR/${fname[$d]}
seed=5486

# -v allows to pass variables to LAMMPS (see lammps.in file)
srun $lmp -v SEED $seed -v T $temp -v INITCONFIG $struct -v MLIPPOT $mlippot -log none -in $lammpsinput &>$lammpslog
