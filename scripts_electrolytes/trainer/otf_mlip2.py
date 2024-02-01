import os
import logging
import subprocess as subp
from ..utils.time import when_is_now
from .otf import OtfMtpTrainer

class OtfMtp2Trainer(OtfMtpTrainer):
    '''
        mlip_flags: json file containing custom training variables for MTP training

        mlip_ini: path to the mlip.ini file containing info about the path to the initial potential and 
                  extrapolation selection flags

        See OtfMtpTrainer for description of other input variables
    '''
    
    def __init__(self, mtp_path=None, init_mtp=None, init_train_db=None, abi_input=None,
                 dft_job_args=None, dft_job_script=None, username=None, train_job_args=None,
                 train_job_script=None, valid_db=None, submit=True, abicommand=None,
                 restart_iterstep=None, mlip_flags=None, mlip_ini=None, stop_at_max_nsteps=False):

        super(OtfMtp2Trainer, self).__init__(mtp_path=mtp_path, init_mtp=init_mtp, init_train_db=init_train_db,
                                             abi_input=abi_input, dft_job_args=dft_job_args,
                                             dft_job_script=dft_job_script, username=username,
                                             train_job_args=train_job_args, train_job_script=train_job_script,
                                             valid_db=valid_db, submit=submit, abicommand=abicommand,
                                             restart_iterstep=restart_iterstep, stop_at_max_nsteps=stop_at_max_nsteps)

        self.preselect_fname = 'preselected.cfg'

        self.set_mlip_variables(mlip_ini, mlip_flags)

        if not self.restart_iterstep:
            self.initiate_training()

        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


    def initiate_training(self):

        self.owd = os.getcwd()
        os.makedirs('0', exist_ok=True)
        os.system('cp {} 0/current.mtp'.format(self.init_mtp))
        os.system('cp {} 0/train.cfg'.format(self.init_train))
        self.iterstep = 1


    def train_from_lammpsmd(self, lammps_path=None, md_nsteps=10000, lammps_input=None, lammps_struct=None, temp=None,
                            atomic_species=None, relaunch=True, dry_run=False):

        self.set_lammps_variables(lammps_path, md_nsteps, lammps_input, lammps_struct, temp, atomic_species, relaunch)
        if dry_run:
            print('Dry run finished, all paths were checked. Proceed with training.')
            quit()

        finished = False

        if self.restart_iterstep:
            self.iterstep = self.restart_iterstep
            self.owd = os.getcwd()

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

                if self.stop_at_max_nsteps and self.max_nsteps_reached:
                    # In stop_at_max_nsteps mode, use the last trained potential from the 1st completed MD run
                    finished = True
                    self.run('cp current.mtp {}/final.mtp'.format(self.owd))
                    self.run('cp train.cfg {}/final.cfg'.format(self.owd))
                    self.run('echo "  Extrapolative confugurations were preselected but max MD steps was reached.\n  Stopping per user request.">>iter_output.txt')
                    os.chdir(self.owd)
                else:
                    os.chdir(self.owd)

            else:
                finished = True
                self.run('cp prev.mtp {}/final.mtp'.format(self.owd))
                self.run('cp ../{}/train.cfg {}/final.cfg'.format(self.iterstep-1, self.owd))
                self.run('echo "  No extrapolative configuration were preselected.">>iter_output.txt')
                os.chdir(self.owd)

        logging.info("{}: OTF learning scheme completed after {} steps".format(when_is_now(), self.iterstep))

    def compute_als(self):
        logging.info('    {}: Creating active learning state...'.format(when_is_now()))
        self.run('echo "  Active set construction...">>iter_output.txt')
        command = '{} calc-grade ../{}/current.mtp ../{}/train.cfg ../{}/train.cfg out.cfg --als-filename=state.als>>iter_output.txt'.format(
                  self.mtp, self.iterstep-1, self.iterstep-1, self.iterstep-1)
        self.run(command)
        self.run('rm out.cfg')


    def mdrun_select(self):

        self.compute_als()

        logging.info('    {}: Running MD trajectory...'.format(when_is_now()))
        self.run('echo "  Running MD trajectory...">>iter_output.txt')
        self.run(['cp ../{}/current.mtp prev.mtp'.format(self.iterstep-1)])
        self.run('touch preselected.cfg')
        command = '{} -v SEED 1 -v T {} -v NSTEP {} -v MLIP_INI {} -v STRUCT {} -log none -in {} &>lammps.log'.format(
                self.lammps, self.temperature, self.mdsteps, self.mlip_ini, self.lammps_struct, self.lammps_input)
        self.run(command)

        # Check if nsteps was reached
        out = self.check_process(f'grep "Breaking threshold exceeded" lammps.log')
        if out.returncode == 0:
            self.max_nsteps_reached = False
            self.run('echo "        Breaking threshold exceeded">>iter_output.txt')
        else:
            self.max_nsteps_reached = True
            self.run('echo "        Max MD steps was reached without exceeding beaking threshold">>iter_output.txt')


    def select_configs(self):
        self.run('echo "  Selecting configurations...">>iter_output.txt')
        command = (f'{self.mtp} select-add ../{self.iterstep-1}/current.mtp ../{self.iterstep-1}/train.cfg {self.preselect_fname} add_to_train.cfg --als-filename=state.als>>iter_output.txt')
        self.run(command)
        self.run('rm selected.cfg')
        self.run('rm state.als')


    def set_traincommand(self):
        if self.valid_db:
            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --valid-cfgs={} --update-mindist".format(self.mtp, self.valid_db)
        else:
            traincommand = "{} train prev.mtp train.cfg --trained-pot-name=current.mtp --update-mindist".format(self.mtp)
        return traincommand

    def set_convert_poscar_command(self):
        command = '{} convert-cfg ../add_to_train.cfg POSCAR --output-format=vasp-poscar>>../iter_output.txt'.format(self.mtp)
        return command

    def fix_poscar(self, fname):
        # fix the POSCAR files produced by MTP as they do not include atomic species information
        # According to the VASP wiki, this line should always be the 6th line, as the 5 first are mandatory
        with open(fname, 'r+') as f:
            content = f.readlines()
            string = ' '+' '.join(self.atomic_species)+'\n'
            content.insert(5, string)
            f.seek(0)
            f.writelines(content)

    def compute_validation_errors(self):
        return


    def fix_poscar(self, fname):
        # fix the POSCAR files produced by MTP as they do not include atomic species information
        # According to the VASP wiki, this line should always be the 6th line, as the 5 first are mandatory
        with open(fname, 'r+') as f:
            content = f.readlines()
            string = ' '+' '.join(self.atomic_species)+'\n'
            content.insert(5, string)
            f.seek(0)
            f.writelines(content)


    def run(self, command):
        subp.run(command, shell=True)


    def set_mlip_variables(self, mlip_ini, training_flags):

        if not mlip_ini:
            raise ValueError('Must provide mlip.ini file as mlip_ini')
        else:
            if not os.path.exists(mlip_ini):
                raise FileNotFoundError(f'mlip_ini file {mlip_ini} not found')
            self.mlip_ini = os.path.abspath(mlip_ini)

        if training_flags:
            if not os.path.exists(mlip_flags):
                raise FileNotFoundError(f'mlip_flags file {training_flags} not found')
            self.set_mlip_flags(training_flags)
        else:
            self.training_flags = None
