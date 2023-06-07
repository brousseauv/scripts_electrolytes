import os
import subprocess as subp
from ..interfaces.abinit_interface import poscar_to_abivars, load_abivars, input_from_dict
from ..interfaces.slurm_interface import SlurmWatcher
from ..database.dbcreator import MtpDbCreator

class OtfMtpTrainer:

    def __init__(self, mtp_path=None, init_mtp=None, init_train_db=None, abi_input=None,
                 dft_job_args=None, dft_job_script=None, username=None, train_job_args=None,
                 train_job_script=None, valid_db=None, submit=True):

        if not mtp_path:
            raise ValueError('Must provide path to MLP executable in mlp_path')
        else:
            self.mtp = mtp_path

        if not init_mtp:
            raise ValueError('Must provide initial mpt.pot in init_mtp')
        else:
            self.init_mtp = init_mtp

        if not init_train_db:
            raise ValueError('Must provide initial training set in init_train')
        else:
            self.init_train = init_train_db

        if not abi_input:
            raise ValueError('Must provide abinit variables set as abi_input')

        if submit:
            if not dft_job_script:
                msg = """DFT jobs will be submitted. Must provide path to a sample DFT submission 
                         script as dft_job_script, or use submit=False."""
                raise ValueError(msg)
            else:
                self.dft_jobscript = os.path.abspath(dft_job_script)

            if not train_job_script:
                msg = """Training jobs will be submitted. Must provide path to a sample training
                         submission script as train_job_script, or use submit=False."""
                raise ValueError(msg)
            else:
                self.train_jobscript = os.path.abspath(train_job_script)

            if not username:
                raise ValueError("When submit=True, must provide Slurm username as username for job monitoring.")
            else:
                self.username = username


        if dft_job_args:
            dft_job_args = self.set_job_args(dft_job_args)
        self.dft_job_args = dft_job_args 

        if train_job_args:
            train_job_args = self.set_job_args(train_job_args)
        self.train_job_args = train_job_args 

        self.valid_db = os.path.abspath(valid_db)
        self.submit = submit

        self.initiate_training()
        self.set_abivars(abi_input)


    def initiate_training(self):

        self.owd = os.getcwd()
        os.makedirs('0', exist_ok=True)
        os.system('cp {} 0/current.mtp'.format(self.init_mtp))
        os.system('cp {} 0/train.cfg'.format(self.init_train))
        self.iterstep = 1


    def set_abivars(self, fname):
        self.abivars, self.abipseudos = load_abivars(fname)


    def train_from_lammpsmd(self, lammps_path=None, md_nsteps=10000, lammps_input=None, lammps_struct=None, mlip_ini=None, temp=None,
                            atomic_species=None):

        self.set_variables(lammps_path, md_nsteps, lammps_input, lammps_struct, mlip_ini, temp, atomic_species)
        finished = False

        ##########
        ## THIS WILL BE THE WHILE LOOP
        while self.iterstep<5: # will be while fiished=False, with a condition on the number of preselected configurations
            print('Starting active learning step {}'.format(self.iterstep))
            # remove the while, do a single run before 
            iterdir = '{}'.format(self.iterstep)
            os.makedirs(iterdir, exist_ok=True) 
            self.iterdir = os.path.abspath(iterdir)
            os.chdir(self.iterdir)


            self.run(['echo "\nActive learning step {}">>iter_output.txt'.format(self.iterstep)])
            self.compute_als()
            self.mdrun_select()
            self.select_configs()
            nselect = self.check_configs()
            if nselect>0:
                self.launch_dft(nselect)
                self.watch_jobs('iter{}_config'.format(self.iterstep), msg='Watching DFT jobs...')
                # check jobs output, relaunch if necessary
                # self.check_dft_output -> is it finished, is it converged, does it have TIME LIMIT / oom handler etc. for relaunch
                # aggregate configs and add to training set
                self.collect_dft()
                # I can do this one then work on the check status and relaunch options.
                # retrain model
                self.train_model()
                self.watch_jobs('iter{}_train'.format(self.iterstep), msg='Watching training job...')

                self.iterstep += 1
                print('going to iterstep {}'.format(self.iterstep))
                os.chdir(self.owd)

            else:
                finished = True
                os.chdir(self.owd)

        self.run('echo "Active learning scheme completed after {} steps"'.format(self.iterstep))

    def compute_als(self):
        self.run('echo "  Active set construction...">>iter_output.txt')
        command = '{} calc-grade ../{}/current.mtp ../{}/train.cfg ../{}/train.cfg out.cfg --als-filename=state.als>>iter_output.txt'.format(
                  self.mtp, self.iterstep-1, self.iterstep-1, self.iterstep-1)
        self.run(command)
        self.run('rm out.cfg')


    def mdrun_select(self):
        self.run('echo "  Running MD trajectory...">>iter_output.txt')
        self.run(['cp ../{}/current.mtp prev.mtp'.format(self.iterstep-1)])
        command = '{} -v SEED 1 -v T {} -v NSTEP {} -v MLIP_INI {} -v STRUCT {} -log none -in {} &>lammps.log'.format(
                self.lammps, self.temperature, self.mdsteps, self.mlip_ini, self.lammps_struct, self.lammps_input)
        self.run(command)


    def select_configs(self):
        self.run('echo "  Selecting configurations...">>iter_output.txt')
        command = '{} select-add ../{}/current.mtp ../{}/train.cfg preselected.cfg add_to_train.cfg --als-filename=state.als>>iter_output.txt'.format(self.mtp, self.iterstep-1, self.iterstep-1)
        self.run(command)
        self.run('rm selected.cfg')
        self.run('rm state.als')


    def check_configs(self):
        nselect = int(subp.run('grep BEGIN_CFG add_to_train.cfg | wc -l', shell=True, capture_output=True).stdout)
        return nselect

    def launch_dft(self, njobs):
        self.run('echo "  Preparing and launching DFT calculations. This may take some time.">>iter_output.txt')
        os.makedirs('calc', exist_ok=True)
        self.calcdir = os.path.abspath('calc')
        os.chdir(self.calcdir)

        self.run('echo "  Processing POSCAR files...">>../iter_output.txt')
        command = '{} convert-cfg ../add_to_train.cfg POSCAR --output-format=vasp-poscar>>../iter_output.txt'.format(self.mtp)
        # UNCOMMENT WHEN READY
        self.run(command)
        
        if njobs == 1:
            self.launch_job(0, 'POSCAR')
        else:
            for j in range(njobs):
                self.launch_job(j, 'POSCAR{}'.format(j))

        os.chdir(self.iterdir) # returning to the iter directory

    def launch_job(self, j, config):

        # submit jobs. This would work better with a shell scripts, though...
        self.fix_poscar(config)
        struct = poscar_to_abivars(config)

        workdir = self.create_workdir(j)
        os.chdir(workdir)

        # do some stuff
        # create abi input
        self.write_abi_input(struct)
        # launch job
        if self.dft_job_args:
            self.dft_job_args = self.dft_job_args + ' --job-name=iter{}_config{}'.format(self.iterstep, j)
        else:
            self.dft_job_args = ' --job-name=iter{}_config{}'.format(self.iterstep, j)
        command = "sbatch {} {}".format(self.dft_job_args, self.dft_jobscript)
        self.run(command)
        os.chdir(self.calcdir)


    def watch_jobs(self, rootname):
        watcher = SlurmWatcher(rootname, self.username)
        watcher.watch()


    def collect_dft(self):
        # Here we are inside self.calcdir
        self.run('echo "  Collecting DFT results...">>iter_output.txt')
        self.run('cp {}/../{}/train.cfg .'.format(self.iterdir, self.iterstep-1))
        db = MtpDbCreator(dbname='train.cfg', append=True)
        db.db_from_gsr(self.calcdir)


    def fix_poscar(self, fname):
        # fix the POSCAR files produced by MTP as they do not include atomic species information
        # According to the VASP wiki, this line should always be the 6th line, as the 5 first are mandatory
        with open(fname, 'r+') as f:
            content = f.readlines()
            string = ' '+' '.join(self.atomic_species)+'\n'
            content.insert(5, string)
            f.seek(0)
            f.writelines(content)


    def train_model(self):
        # Start new training from the fitted potential from the previous iteration
        self.run('echo "  Training potential from updated training set. This may take some time.">>iter_output.txt')
        if self.train_job_args:
            self.train_job_args = self.train_job_args + ' --job-name=iter{}_train'.format(self.iterstep)
        else:
            self.train_job_args = ' --job-name=iter{}_train'.format(self.iterstep)

        # will need a train_job.sh example file and train_job_params, in opposition to dft job params (modify this after)
        # --update-mindist to be tested!!!
        #prev = os.path.abspath('prev.mtp')
        #train = os.path.abspath('train.cfg')
        if self.valid_db:
            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --valid-cfgs={} --update-mindist".format(self.mtp, self.valid_db)
        else:
            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --update-mindist".format(self.mtp)

        # copy the submission file to current dir
        # and replace a tagged keyword with the current string
        # FIX ME: this is a little sketchy but for now it works. 
        # Could be improved though
        # Perhaps keep only the header, i.e. the .sh script up to the mlp train command ?
        self.run('cp {} train.sh'.format(self.train_jobscript))
#        with open('train.sh', 'r+') as f:
#            content = f.read()
#            content = content.replace('mytraincommand', '{}'.format(traincommand))
#            f.seek(0)
#            f.write(content)

        # version 2: append command after header
        with open('train.sh', 'a') as f:
            f.write('srun {} >&train_$SLURM_JOB_ID.out'.format(traincommand))

        runcommand = "sbatch {} train.sh".format(self.train_job_args)
        self.run(runcommand)

    def write_abi_input(self, struct):
        input_from_dict(struct, self.abivars, self.abipseudos)


    def run(self, command):
        subp.run(command, shell=True)


    def create_workdir(self, j):
        string = 'config{}'.format(j)
        os.makedirs(string, exist_ok=True) # FIX ME: remove the exist_ok ?
        return os.path.abspath(string)


    def set_job_args(self, args):
        string = ' '.join(args)
        return string


    def set_variables(self,lammps_path, md_nsteps, lammps_input, lammps_struct, mlip_ini, temp, atomic_species):
        if not lammps_path:
            raise ValueError('Must provide path to LAMMPS executable as lammps_path')
        else:
            self.lammps = lammps_path

        if not lammps_input:
            raise ValueError('Must provide .in input file for LAMMPS as lammps_input')
        else:
            self.lammps_input = os.path.abspath(lammps_input)

        if not lammps_struct:
            raise ValueError('Must provide initial structure in .lmp format for LAMMPS as lammps_struct')
        else:
            self.lammps_struct = os.path.abspath(lammps_struct)

        if not temp:
            raise ValueError('Must provide MD temperature as temp')
        else:
            self.temperature = temp

        if not mlip_ini:
            raise ValueError('Must provide mlip.ini file as mlip_ini')
        else:
            self.mlip_ini = os.path.abspath(mlip_ini)

        if not atomic_species:
            raise ValueError("Must provide the list of atomic species using chemical symbols, in the same order as the MTP data as atomic_species")
        else:
            self.atomic_species = atomic_species

        self.mdsteps = md_nsteps


