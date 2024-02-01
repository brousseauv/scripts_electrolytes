import os
import logging
import json
import subprocess as subp
from ..interfaces.abinit_interface import poscar_to_abivars, load_abivars, input_from_dict
from ..interfaces.slurm_interface import SlurmWatcher
from ..database.db_creator import MtpDbCreator
from ..utils.time import when_is_now, increase_jobtime

class OtfMtpTrainer:

    def __init__(self, mtp_path=None, init_mtp=None, init_train_db=None, abi_input=None,
                 dft_job_args=None, dft_job_script=None, username=None, train_job_args=None,
                 train_job_script=None, valid_db=None, submit=True, abicommand=None,
                 restart_iterstep=None, stop_at_max_nsteps=False):

        '''
            Base class to train MTP models on-the-fly

            Input:
                mtp_path: path to the mlp executable

                init_mtp: initial MTP potential (can be empty if starting from scratch)

                init_train_db: path to the training database for the initial MTP model

                abi_input: json file containing Abinit variables

                dft_job_args: list of strings containing Slurm commands to be overriden compared to the dft_job_script
                              default=None

                dft_job_script: sample Slurm submission script (header) WITHOUT main command line for DFT jobs

                username: cluster username for job monitoring, when using submit=True option

                train_job_args: list of strings containing Slurm commands to be overriden compared to the train_job_script
                                default=None

                train_job_script: sample Slurm submission script (header) WITHOUT main command line for MTP training jobs

                valid_db: path to the validation database (optional)

                submit: whether so submit DFT and training jobs to the SLURM queue
                        Default: True

                abicommand: custom command to run Abinit jobs (if it differs from srun $abinit $inputfile>&$log)

                restart_iterstep: Continue OTF training at iterstep N (which was not completed on previous run)
                                  default: None (start from init_mtp and init_train_db)

                stop_at_max_nsteps: allows to stop the OTF procedure when the MD runs stops encouning configurations with gamma>gamma_break
                                    even if there are configurations with gamma>gamma_select
                                    intended for initial aggregation of configurations from potentials trained from a small number of configs
                                    Default:False (stop OTF procedure when no configurations are preselected)
        '''

        if not mtp_path:
            raise ValueError('Must provide path to MLP executable in mlp_path')
        else:
            if not os.path.exists(mtp_path):
                raise FileNotFoundError(f'mtp_path file {mtp_path} not found')
            self.mtp = mtp_path

        if not init_mtp:
            raise ValueError('Must provide initial mtp.pot or mtp.almtp in init_mtp')
        else:
            if not os.path.exists(init_mtp):
                raise FileNotFoundError(f'init_mtp file {init_mtp} not found')
            self.init_mtp = init_mtp

        if not init_train_db:
            raise ValueError('Must provide initial training set in init_train')
        else:
            if not os.path.exists(init_train_db):
                raise FileNotFoundError(f'init_train_db file {init_train_db} not found')
            self.init_train = init_train_db

        if not abi_input:
            raise ValueError('Must provide abinit variables set as abi_input')
        else:
            if not os.path.exists(abi_input):
               raise FileNotFoundError(f'abi_input file {abi_input} not found')
           
        if submit:
            if not dft_job_script:
                msg = """DFT jobs will be submitted. Must provide path to a sample DFT submission 
                         script as dft_job_script, or use submit=False."""
                raise ValueError(msg)
            else:
                if not os.path.exists(dft_job_script):
                    raise FileNotFoundError(f'dft_job_script file {dft_job_script} not found')
                self.dft_jobscript = os.path.abspath(dft_job_script)

            if not train_job_script:
                msg = """Training jobs will be submitted. Must provide path to a sample training
                         submission script as train_job_script, or use submit=False."""
                raise ValueError(msg)
            else:
                if not os.path.exists(train_job_script):
                    raise FileNotFoundError(f'train_job_script file {train_job_script} not found')
                self.train_jobscript = os.path.abspath(train_job_script)

            if not username:
                raise ValueError("When submit=True, must provide Slurm username as username for job monitoring.")
            else:
                self.username = username
        else:
            if not abicommand:
                raise ValueError("When submit=False, must provide abinit command as abicommand.")
            else:
                self.abicommand = abicommand


        if dft_job_args:
            dft_job_args = self.set_job_args(dft_job_args)
        self.dft_job_args = dft_job_args 

        if train_job_args:
            train_job_args = self.set_job_args(train_job_args)
        self.train_job_args = train_job_args 

        if not os.path.exists(valid_db):
            raise FileNotFoundError(f'valid_db file {valid_db} not found')
        self.valid_db = os.path.abspath(valid_db)
        self.submit = submit

        self.restart_iterstep = restart_iterstep
        self.stop_at_max_nsteps = stop_at_max_nsteps

        self.set_abivars(abi_input)
        
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


    def set_abivars(self, fname):
        self.abivars, self.abipseudos = load_abivars(fname)


    def set_mlip_flags(self, fname):
        self.mlip_flags = json.load(open(fname))


    def preselect_configs(self):
        ''' Calls all tasks to perform one learning-on-the-fly cycle '''

        # LAMMPS MD run to preselect extrapolative configurations
        self.mdrun_select()
        # Check if configurations were preselected during the MD run
        npreselect = self.check_configs(f'{self.preselect_fname}')

        return npreselect

    def run_otf_iteration(self):

        # Select configurations for DFT calculations
        self.select_configs()
        nselect = self.check_configs('add_to_train.cfg')
        # Evaluate selected configurations with DFT
        self.launch_dft(nselect)
        # Check job outputs and relaunch if necessary/requested
        dft_error = self.check_dft_output(nselect)
        
        if dft_error:
            self.relaunch_dft(nselect)

        # Collect energy/forces/stresses from DFT data
        self.collect_dft()
        # Retrain MTP potential
        self.train_model()


    def check_configs(self, fname):
        nselect = int(subp.run('grep BEGIN_CFG {} | wc -l'.format(fname), shell=True, capture_output=True).stdout)
        return nselect

    def launch_dft(self, njobs):
        logging.info('    {}: Running DFT calculcations...'.format(when_is_now()))
        self.run('echo "  Preparing and launching DFT calculations. This may take some time.">>iter_output.txt')
        os.makedirs('calc', exist_ok=True)
        self.calcdir = os.path.abspath('calc')
        os.chdir(self.calcdir)

        self.run('echo "  Processing POSCAR files...">>../iter_output.txt')
        command = self.set_convert_poscar_command()

        self.run(command)
        
        if njobs == 1:
            self.launch_job(0, 'POSCAR')
        else:
            for j in range(njobs):
                self.launch_job(j, 'POSCAR{}'.format(j))

        # Monitor job status if submitted to the queue
        if self.submit:
            self.watch_jobs('iter{}_config'.format(self.iterstep), msg='Watching DFT jobs...')


    def launch_job(self, j, config):
        # Must be run from calcdir
        self.fix_poscar(config)
        struct = poscar_to_abivars(config)

        workdir = self.create_workdir(j)
        os.chdir(workdir)

        self.write_abi_input(struct)
        
        if self.submit:
            if self.dft_job_args:
                dft_job_args = self.dft_job_args + ' --job-name=iter{}_config{}'.format(self.iterstep, j)
            else:
                dft_job_args = ' --job-name=iter{}_config{}'.format(self.iterstep, j)
            command = "sbatch {} {}".format(dft_job_args, self.dft_jobscript)
        else:
            command = "srun {} run.abi>& log".format(self.abicommand)
        self.run(command)
        os.chdir(self.calcdir)


    def watch_jobs(self, rootname, msg):
        os.chdir(self.iterdir)
        watcher = SlurmWatcher(rootname, self.username)
        watcher.watch(msg=msg)


    def collect_dft(self):
        os.chdir(self.iterdir)
        self.run('echo "  Collecting DFT results...">>iter_output.txt')
        self.run('cp {}/../{}/train.cfg .'.format(self.iterdir, self.iterstep-1))
        db = MtpDbCreator(dbname='train.cfg', append=True)
        db.db_from_gsr(self.calcdir)


    def train_model(self):
        # Start new training from the fitted potential from the previous iteration
        # Must run in iterdir
        os.chdir(self.iterdir)
        logging.info('    {}: Retraining MTP potential...'.format(when_is_now()))
        self.run('echo "  Training potential from updated training set. This may take some time.">>iter_output.txt')

        traincommand = self.set_traincommand()
        if self.mlip_flags:
            trainflags = ' '.join([f'--{key}={val}' for key, val in self.mlip_flags.items()])
            traincommand = ' '.join([traincommand, trainflags])
        # FIX ME: DELETE THIS ONCE TESTS ARE FINISHED
        logging.info(f' current training command:')
        logging.info(f'     {traincommand}')

        if self.submit:
            if self.train_job_args:
                train_job_args = self.train_job_args + ' --job-name=iter{}_train'.format(self.iterstep)
            else:
                train_job_args = ' --job-name=iter{}_train'.format(self.iterstep)

            self.run('cp {} train.sh'.format(self.train_jobscript))

            with open('train.sh', 'a') as f:
                f.write('srun {} >&train_$SLURM_JOB_ID.out'.format(traincommand))

            runcommand = "sbatch {} train.sh".format(train_job_args)
        else:
            runcommand = 'srun {} >&train.out'.format(traincommand)

        self.run(runcommand)
        # Monitor training job if submitted to the queue
        if self.submit:
            self.watch_jobs('iter{}_train'.format(self.iterstep), msg='Watching training job...')

        logging.info('    {}: ... training completed\n\n'.format(when_is_now()))
        self.run('echo "  Training completed\n\n">>iter_output.txt')

        if self.valid_db:
            self.compute_validation_errors()


    def check_dft_output(self, njobs):
        os.chdir(self.calcdir)
        self.failed_calc_index = []
        self.errormsg = []

        for j in range(njobs):
            # Check for some typical errors
            scf = self.check_process('grep -A2 ScfConvergenceWarning config{}/log*'.format(j))
            if scf.returncode == 0:
                self.failed_calc_index.append(j)
                self.errormsg.append('scfconv')
                # update nstep in abivars dictionnary
                nstep = int(str(scf.stdout).split('nstep ')[1].split(' ')[0])
                self.abivars['nstep'] = int(1.5*nstep)
                continue

            oom = self.check_process('grep  "out-of-memory handler" config{}/log*'.format(j))
            if oom.returncode == 0:
                self.failed_calc_index.append(j)
                self.errormsg.append('memory')
                continue

            time = self.check_process('grep  "TIME LIMIT" config{}/log*'.format(j))
            if time.returncode == 0:
                self.failed_calc_index.append(j)
                self.errormsg.append('timelimit')
                continue

            end_ok = self.check_process('grep  "Calculation completed" config{}/log*'.format(j))
            if end_ok.returncode != 0:
                self.failed_calc_index.append(j)
                self.errormsg.append('unknown')
                continue

        if len(self.failed_calc_index)>0:
            return True
        else:
            return False


    def check_process(self, command):
        output = subp.run(command, shell=True, capture_output=True)
        return output


    def relaunch_dft(self, njobs):

        os.chdir(self.calcdir) # launching must be done from calcdir
        if not self.relaunch:
            raise Exception('Some exceptions occured in iterstep {} for configs {}: {}.Stopping OTF procedure.'.format(
                            self.iterstep, self.failed_calc_index, self.errormsg))

        nrelaunch = 0
        while nrelaunch < self.max_relaunch:
            os.chdir(self.calcdir)
            for j, err in zip(self.failed_calc_index, self.errormsg):
                if err == 'scfconv':
                    try:
                        self.launch_job(j, 'POSCAR{}'.format(j))
                    except:
                        self.launch_job(j, 'POSCAR')

                if self.submit:  # "oom", "time" should only be relevant for submitted jobs. "unknown", well... check it out!
                    # relaunch options, relaunch job if submit, modifying parameters as needed
                    # errormsg : in "memory", "timelimit', 'unknown'
                    if err == 'timelimit':
                        # ADD A CHECK IF -t OR --time is in dft_job_args
                        try:
                            output = self.check_process("grep -- ' -t ' {}".format(self.dft_jobscript))
                            time = str(output.stdout).split(' -t ')[1].split('\\n')[0]
                        except:
                            output = self.check_process("grep -- ' --time ' {}".format(self.dft_jobscript))
                            time = str(output.stdout).split(' --time ')[1].split('\\n')[0]

                        jobtime = increase_jobtime(time)
                        arg = '--time {}'.format(jobtime)
                        self.relaunch_dft_job(j, arg)
                    
                    if err == 'memory':
                        raise NotImplementedError('Automatic memory increase not implemented. Do it or fix it by hand.')
                        try:
                            output = self.check_process("grep -- ' -m ' {}".format(self.dft_jobscript))
                            time = str(output.stdout).split(' -t ')[1].split('\\n')[0]
                        except:
                            output = self.check_process("grep -- ' --time ' {}".format(self.dft_jobscript))
                            time = str(output.stdout).split(' --time ')[1].split('\\n')[0]

                        jobtime = increase_jobtime(time, 1.25**nrelaunch)
                        arg = '--time {}'.format(jobtime)
                        self.relaunch_dft_job(j, arg)

            os.chdir(self.iterdir) # returning to the iter directory for job monitoring
            # Monitor job status if submitted to the queue
            if self.submit:
                self.watch_jobs('iter{}_config'.format(self.iterstep), msg='Watching DFT jobs, relaunch={}...'.format(nrelaunch+1))

            dft_error = self.check_dft_output(njobs) #-> break if not ok
            if dft_error:
                nrelaunch += 1
            else:
                return
        
        raise Exception('The number of DFT relaunch in iterstep {} reached max_relaunch = {}'.format(self.iterstep, self.max_relaunch))


    def relaunch_dft_job(self, idx, arg):
        # reset the dft_job_args with new arg and submit
        os.chdir('config{}'.format(idx))
        if self.dft_job_args:
            dft_job_args = self.dft_job_args + ' --job-name=iter{}_config{} {}'.format(self.iterstep, idx, arg)
        else:
            dft_job_args = ' --job-name=iter{}_config{} {}'.format(self.iterstep, idx, arg)
        command = "sbatch {} {}".format(dft_job_args, self.dft_jobscript)
        self.run(command)
        os.chdir(self.calcdir)


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


    def set_lammps_variables(self,lammps_path, md_nsteps, lammps_input, lammps_struct, temp, atomic_species, relaunch):
        if not lammps_path:
            raise ValueError('Must provide path to LAMMPS executable as lammps_path')
        else:
            if not os.path.exists(lammps_path):
                raise FileNotFoundError(f'lammps_path file {lammps_path} not found')
            self.lammps = lammps_path

        if not lammps_input:
            raise ValueError('Must provide .in input file for LAMMPS as lammps_input')
        else:
            if not os.path.exists(lammps_input):
                raise FileNotFoundError(f'lammps_input file {lammps_input} not found')
            self.lammps_input = os.path.abspath(lammps_input)

        if not lammps_struct:
            raise ValueError('Must provide initial structure in .lmp format for LAMMPS as lammps_struct')
        else:
            if not os.path.exists(lammps_struct):
                raise FileNotFoundError(f'lammps_struct file {lammps_struct} not found')
            self.lammps_struct = os.path.abspath(lammps_struct)

        if not temp:
            raise ValueError('Must provide MD temperature as temp')
        else:
            self.temperature = temp

        if not atomic_species:
            raise ValueError("Must provide the list of atomic species using chemical symbols, in the same order as the MTP data as atomic_species")
        else:
            self.atomic_species = atomic_species

        self.mdsteps = md_nsteps
        self.max_relaunch = 2
        self.relaunch = relaunch


