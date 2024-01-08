import os
import logging
import subprocess as subp
from ..utils.time import when_is_now
from .otf import OtfMtpTrainer

class OtfMtp3Trainer(OtfMtpTrainer):
    '''
        mlip_flags: json file containing custom training variables for MTP training

        training_mode: whether the MTP potential is trained on the full configuration ('cfg') on or atomic
                       neighborhoods ('nbh') 
                       default: 'cfg' (only mode supported)

        See OtfMtpTrainer for description of other input variables
    '''

    def __init__(self, mtp_path=None, init_mtp=None, init_train_db=None, abi_input=None,
                 dft_job_args=None, dft_job_script=None, username=None, train_job_args=None,
                 train_job_script=None, valid_db=None, submit=True, abicommand=None,
                 restart_iterstep=None, mlip_flags=None, training_mode='cfg'):

        super(OtfMtp3Trainer, self).__init__(mtp_path=mtp_path, init_mtp=init_mtp, init_train_db=init_train_db,
                                             abi_input=abi_input, dft_job_args=dft_job_args,
                                             dft_job_script=dft_job_script, username=username,
                                             train_job_args=train_job_args, train_job_script=train_job_script,
                                             valid_db=valid_db, submit=submit, abicommand=abicommand,
                                             restart_iterstep=restart_iterstep)

        self.preselect_fname = 'preselected.cfg'
        self.mlip_path = 'prev.almtp'

        self.set_mlip_variables(mlip_flags)
        self.check_training_mode(training_mode)

        if not self.restart_iterstep:
            self.initiate_training()

        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


    def initiate_training(self):

        self.owd = os.getcwd()
        os.makedirs('0', exist_ok=True)
        os.system('cp {} 0/current.almtp'.format(self.init_mtp))
        os.system('cp {} 0/train.cfg'.format(self.init_train))
        self.iterstep = 1


    def train_from_lammpsmd(self, lammps_path=None, md_nsteps=10000, lammps_input=None, lammps_struct=None, temp=None,
                            atomic_species=None, relaunch=True, preselect_fname=None):

        self.set_lammps_variables(lammps_path, md_nsteps, lammps_input, lammps_struct, temp, atomic_species, relaunch)

        if preselect_fname:
            self.preselect_fname = preselect_fname

        finished = False

        if self.restart_iterstep:
            self.iterstep = self.restart_iterstep
            self.owd = os.getcwd()

#        while self.iterstep<4: # will be while fiished=False, with a condition on the number of preselected configurations
        while finished == False:  

            logging.info('{}: Starting OTF learning step {}'.format(when_is_now(), self.iterstep))
            iterdir = '{}'.format(self.iterstep)
            os.makedirs(iterdir, exist_ok=True) 
            self.iterdir = os.path.abspath(iterdir)
            os.chdir(self.iterdir)

            self.run(['echo "\nOTF learning step {}">>iter_output.txt'.format(self.iterstep)])
            npreselect = self.preselect_configs()

            if npreselect>0:
                self.run_otf_iteration()
                self.iterstep += 1
                os.chdir(self.owd)

            else:
                finished = True
                self.run('cp prev.almtp {}/final.almtp'.format(self.owd))
                self.run('cp ../{}/train.cfg {}/final.cfg'.format(self.iterstep-1, self.owd))
                self.run('echo "  No extrapolative configuration were preselected.">>iter_output.txt')
                os.chdir(self.owd)

        logging.info("{}: OTF learning scheme completed after {} steps".format(when_is_now(), self.iterstep))

#    def compute_als(self):
#        logging.info('    {}: Creating active learning state...'.format(when_is_now()))
#        self.run('echo "  Active set construction...">>iter_output.txt')
#        command = '{} calc-grade ../{}/current.mtp ../{}/train.cfg ../{}/train.cfg out.cfg --als-filename=state.als>>iter_output.txt'.format(
#                  self.mtp, self.iterstep-1, self.iterstep-1, self.iterstep-1)
#        self.run(command)
#        self.run('rm out.cfg')


    def mdrun_select(self):
#        # Compute active learning state
#        self.compute_als()

        # Run MD trajectory in selection mode
        logging.info('    {}: Running MD trajectory...'.format(when_is_now()))
        self.run('echo "  Running MD trajectory...">>iter_output.txt')
        self.run(['cp ../{}/current.almtp prev.almtp'.format(self.iterstep-1)])
        self.run(f'touch {self.preselect_fname}')
        command = 'srun {} -v SEED 1 -v T {} -v NSTEP {} -v MLIPPATH {} -v STRUCT {} -log none -in {} &>lammps.log'.format(
                self.lammps, self.temperature, self.mdsteps, self.mlip_path, self.lammps_struct, self.lammps_input)
        self.run(command)


    def select_configs(self):
        self.run('echo "  Selecting configurations...">>iter_output.txt')
        command = (f'{self.mtp} select_add ../{self.iterstep-1}/current.almtp ../{self.iterstep-1}/train.cfg {self.preselect_fname} add_to_train.cfg>>iter_output.txt')
        self.run(command)
        #FIX ME: will this file still be created?!?
#        self.run('rm selected.cfg')
#        self.run('rm state.als')


    def set_traincommand(self):
        #if self.valid_db:
        #    traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --valid-cfgs={} --update-mindist".format(self.mtp, self.valid_db)
        #else:
        traincommand = "{} train prev.almtp train.cfg --save_to=current.almtp --update-mindist --al_mode={}".format(self.mtp, self.training_mode)
        return traincommand

    def set_convert_poscar_command(self):
        string = ','.join([str(a) for a in self.atomic_species])
        command = '{} convert ../add_to_train.cfg POSCAR --input_format=txt --output_format=poscar --type_order={}>>../iter_output.txt'.format(self.mtp, string)
        return command

    def fix_poscar(self, fname):
        return

    def compute_validation_errors(self):
        command = f'{self.mtp} check_errors current.almtp {self.valid_db} --report_to="valid_errors.log"'
        self.run(command)
        return
#
#    def check_configs(self, fname):
#        nselect = int(subp.run('grep BEGIN_CFG {} | wc -l'.format(fname), shell=True, capture_output=True).stdout)
#        return nselect

#    def launch_dft(self, njobs):
#        logging.info('    {}: Running DFT calculcations...'.format(when_is_now()))
#        self.run('echo "  Preparing and launching DFT calculations. This may take some time.">>iter_output.txt')
#        os.makedirs('calc', exist_ok=True)
#        self.calcdir = os.path.abspath('calc')
#        os.chdir(self.calcdir)
#
#        self.run('echo "  Processing POSCAR files...">>../iter_output.txt')
#        command = '{} convert-cfg ../add_to_train.cfg POSCAR --output-format=vasp-poscar>>../iter_output.txt'.format(self.mtp)
#        # UNCOMMENT WHEN READY
#        self.run(command)
#        
#        if njobs == 1:
#            self.launch_job(0, 'POSCAR')
#        else:
#            for j in range(njobs):
#                self.launch_job(j, 'POSCAR{}'.format(j))
#
#        # Monitor job status if submitted to the queue
#        if self.submit:
#            self.watch_jobs('iter{}_config'.format(self.iterstep), msg='Watching DFT jobs...')


#    def launch_job(self, j, config):
#        # Must be run from calcdir
#        self.fix_poscar(config)
#        struct = poscar_to_abivars(config)
#
#        workdir = self.create_workdir(j)
#        os.chdir(workdir)
#
#        self.write_abi_input(struct)
#        
#        if self.submit:
#            if self.dft_job_args:
#                dft_job_args = self.dft_job_args + ' --job-name=iter{}_config{}'.format(self.iterstep, j)
#            else:
#                dft_job_args = ' --job-name=iter{}_config{}'.format(self.iterstep, j)
#            command = "sbatch {} {}".format(dft_job_args, self.dft_jobscript)
#        else:
#            command = "srun {} run.abi>& log".format(self.abicommand)
#        self.run(command)
#        os.chdir(self.calcdir)


#    def watch_jobs(self, rootname, msg):
#        os.chdir(self.iterdir)
#        watcher = SlurmWatcher(rootname, self.username)
#        watcher.watch(msg=msg)


#    def collect_dft(self):
#        os.chdir(self.iterdir)
#        self.run('echo "  Collecting DFT results...">>iter_output.txt')
#        self.run('cp {}/../{}/train.cfg .'.format(self.iterdir, self.iterstep-1))
#        db = MtpDbCreator(dbname='train.cfg', append=True)
#        db.db_from_gsr(self.calcdir)


#    def train_model(self):
#        # Start new training from the fitted potential from the previous iteration
#        # Must run in iterdir
#        os.chdir(self.iterdir)
#        logging.info('    {}: Retraining MTP potential...'.format(when_is_now()))
#        self.run('echo "  Training potential from updated training set. This may take some time.">>iter_output.txt')
#
#        if self.valid_db:
#            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --valid-cfgs={} --update-mindist".format(self.mtp, self.valid_db)
#        else:
#            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --update-mindist".format(self.mtp)
#
#
#        if self.submit:
#            if self.train_job_args:
#                train_job_args = self.train_job_args + ' --job-name=iter{}_train'.format(self.iterstep)
#            else:
#                train_job_args = ' --job-name=iter{}_train'.format(self.iterstep)
#
#            # copy the submission file to current dir
#            # and replace a tagged keyword with the current string
#            # FIX ME: this is a little sketchy but for now it works. 
#            # Could be improved though
#            # Perhaps keep only the header, i.e. the .sh script up to the mlp train command ?
#            self.run('cp {} train.sh'.format(self.train_jobscript))
#    #        with open('train.sh', 'r+') as f:
#    #            content = f.read()
#    #            content = content.replace('mytraincommand', '{}'.format(traincommand))
#    #            f.seek(0)
#    #            f.write(content)
#
#            # version 2: append command after header
#            with open('train.sh', 'a') as f:
#                f.write('srun {} >&train_$SLURM_JOB_ID.out'.format(traincommand))
#
#            runcommand = "sbatch {} train.sh".format(train_job_args)
#        else:
#            runcommand = 'srun {} >&train.out'.format(traincommand)
#
#        self.run(runcommand)
#        # Monitor training job if submitted to the queue
#        if self.submit:
#            self.watch_jobs('iter{}_train'.format(self.iterstep), msg='Watching training job...')


    def fix_poscar(self, fname):
        # fix the POSCAR files produced by MTP as they do not include atomic species information
        # According to the VASP wiki, this line should always be the 6th line, as the 5 first are mandatory
        with open(fname, 'r+') as f:
            content = f.readlines()
            string = ' '+' '.join(self.atomic_species)+'\n'
            content.insert(5, string)
            f.seek(0)
            f.writelines(content)


#    def check_dft_output(self, njobs):
#        os.chdir(self.calcdir)
#        self.failed_calc_index = []
#        self.errormsg = []
#
#        for j in range(njobs):
#            # Check for some typical errors
#            scf = self.check_process('grep -A2 ScfConvergenceWarning config{}/log*'.format(j))
#            if scf.returncode == 0:
#                self.failed_calc_index.append(j)
#                self.errormsg.append('scfconv')
#                # update nstep in abivars dictionnary
#                nstep = int(str(scf.stdout).split('nstep ')[1].split(' ')[0])
#                self.abivars['nstep'] = int(1.5*nstep)
#                continue
#
#            oom = self.check_process('grep  "out-of-memory handler" config{}/log*'.format(j))
#            if oom.returncode == 0:
#                self.failed_calc_index.append(j)
#                self.errormsg.append('memory')
#                continue
#
#            time = self.check_process('grep  "TIME LIMIT" config{}/log*'.format(j))
#            if time.returncode == 0:
#                self.failed_calc_index.append(j)
#                self.errormsg.append('timelimit')
#                continue
#
#            end_ok = self.check_process('grep  "Calculation completed" config{}/log*'.format(j))
#            if end_ok.returncode != 0:
#                self.failed_calc_index.append(j)
#                self.errormsg.append('unknown')
#                continue
#
#        if len(self.failed_calc_index)>0:
#            return True
#        else:
#            return False


#    def check_process(self, command):
#        output = subp.run(command, shell=True, capture_output=True)
#        return output


#    def relaunch_dft(self, njobs):
#
#        os.chdir(self.calcdir) # launching must be done from calcdir
#        if not relaunch:
#            raise Exception('Some exceptions occured in iterstep {} for configs {}: {}.Stopping OTF procedure.'.format(
#                            self.iterstep, self.failed_calc_index, self.errormsg))
#
#        nrelaunch = 0
#        while nrelaunch < self.max_relaunch:
#            os.chdir(self.calcdir)
#            for j, err in zip(self.failed_calc_index, self.errormsg):
#                if err == 'scfconv':
#                    try:
#                        self.launch_job(j, 'POSCAR{}'.format(j))
#                    except:
#                        self.launch_job(j, 'POSCAR')
#
#                if self.submit:  # "oom", "time" should only be relevant for submitted jobs. "unknown", well... check it out!
#                    # relaunch options, relaunch job if submit, modifying parameters as needed
#                    # errormsg : in "memory", "timelimit', 'unknown'
#                    if err == 'timelimit':
#                        # ADD A CHECK IF -t OR --time is in dft_job_args
#                        try:
#                            output = self.check_process("grep -- ' -t ' {}".format(self.dft_jobscript))
#                            time = str(output.stdout).split(' -t ')[1].split('\\n')[0]
#                        except:
#                            output = self.check_process("grep -- ' --time ' {}".format(self.dft_jobscript))
#                            time = str(output.stdout).split(' --time ')[1].split('\\n')[0]
#
#                        jobtime = increase_jobtime(time)
#                        arg = '--time {}'.format(jobtime)
#                        self.relaunch_dft_job(j, arg)
#                    
#                    if err == 'memory':
#                        raise NotImplementedError('Automatic memory increase not implemented. Do it or fix it by hand.')
#                        try:
#                            output = self.check_process("grep -- ' -m ' {}".format(self.dft_jobscript))
#                            time = str(output.stdout).split(' -t ')[1].split('\\n')[0]
#                        except:
#                            output = self.check_process("grep -- ' --time ' {}".format(self.dft_jobscript))
#                            time = str(output.stdout).split(' --time ')[1].split('\\n')[0]
#
#                        jobtime = increase_jobtime(time, 1.25**nrelaunch)
#                        arg = '--time {}'.format(jobtime)
#                        self.relaunch_dft_job(j, arg)
#
#            os.chdir(self.iterdir) # returning to the iter directory for job monitoring
#            # Monitor job status if submitted to the queue
#            if self.submit:
#                self.watch_jobs('iter{}_config'.format(self.iterstep), msg='Watching DFT jobs, relaunch={}...'.format(nrelaunch+1))
#
#            dft_error = self.check_dft_output(njobs) #-> break if not ok
#            if dft_error:
#                nrelaunch += 1
#            else:
#                return
#        
#        raise Exception('The number of DFT relaunch in iterstep {} reached max_relaunch = {}'.format(self.iterstep, self.max_relaunch))


#    def relaunch_dft_job(self, idx, arg):
#        # reset the dft_job_args with new arg and submit
#        os.chdir('config{}'.format(idx))
#        if self.dft_job_args:
#            dft_job_args = self.dft_job_args + ' --job-name=iter{}_config{} {}'.format(self.iterstep, idx, arg)
#        else:
#            dft_job_args = ' --job-name=iter{}_config{} {}'.format(self.iterstep, idx, arg)
#        command = "sbatch {} {}".format(dft_job_args, self.dft_jobscript)
#        self.run(command)
#        os.chdir(self.calcdir)


#    def write_abi_input(self, struct):
#        input_from_dict(struct, self.abivars, self.abipseudos)


    def run(self, command):
        subp.run(command, shell=True)


#    def create_workdir(self, j):
#        string = 'config{}'.format(j)
#        os.makedirs(string, exist_ok=True) # FIX ME: remove the exist_ok ?
#        return os.path.abspath(string)


#    def set_job_args(self, args):
#        string = ' '.join(args)
#        return string


    def set_mlip_variables(self, training_flags):

        if training_flags:
            self.set_mlip_flags(training_flags)
        else:
            self.training_flags = None


    def check_training_mode(self, mode):

        if mode == 'cfg':
            self.training_mode = 'cfg'
        elif mode =='nbh':
            raise NotImplementedError('neighborhood OTF training mode not implemented yet.')

#    def set_lammps_variables(self,lammps_path, md_nsteps, lammps_input, lammps_struct, temp, atomic_species, relaunch):
#        if not lammps_path:
#            raise ValueError('Must provide path to LAMMPS executable as lammps_path')
#        else:
#            self.lammps = lammps_path
#
#        if not lammps_input:
#            raise ValueError('Must provide .in input file for LAMMPS as lammps_input')
#        else:
#            self.lammps_input = os.path.abspath(lammps_input)
#
#        if not lammps_struct:
#            raise ValueError('Must provide initial structure in .lmp format for LAMMPS as lammps_struct')
#        else:
#            self.lammps_struct = os.path.abspath(lammps_struct)
#
#        if not temp:
#            raise ValueError('Must provide MD temperature as temp')
#        else:
#            self.temperature = temp
#
#
#        if not atomic_species:
#            raise ValueError("Must provide the list of atomic species using chemical symbols, in the same order as the MTP data as atomic_species")
#        else:
#            self.atomic_species = atomic_species
#
#        self.mdsteps = md_nsteps
#        self.max_relaunch = 2
#        self.relaunch = relaunch


