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
                 restart_iterstep=None, mlip_flags=None, training_mode='cfg', stop_at_max_nsteps=False):

        super(OtfMtp3Trainer, self).__init__(mtp_path=mtp_path, init_mtp=init_mtp, init_train_db=init_train_db,
                                             abi_input=abi_input, dft_job_args=dft_job_args,
                                             dft_job_script=dft_job_script, username=username,
                                             train_job_args=train_job_args, train_job_script=train_job_script,
                                             valid_db=valid_db, submit=submit, abicommand=abicommand,
                                             restart_iterstep=restart_iterstep, stop_at_max_nsteps=stop_at_max_nsteps)

        self.preselect_fname = 'preselected.cfg'
        self.mlip_path = 'prev.almtp'

        self.set_mlip_variables(mlip_flags)
        self.check_training_mode(training_mode)

        if not self.restart_iterstep:
            self.initiate_training()

        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


    def check_paths(self):
        ''' Check that required filepaths actually exist. Do not start training. '''
        
        if not os.path.exists(self.init_mtp):
            raise FileNotFoundError(f'init_mtp file {self.init_mtp} not found')

    def initiate_training(self):

        self.owd = os.getcwd()
        os.makedirs('0', exist_ok=True)
        os.system('cp {} 0/current.almtp'.format(self.init_mtp))
        os.system('cp {} 0/train.cfg'.format(self.init_train))
        self.iterstep = 1


    def train_from_lammpsmd(self, lammps_path=None, md_nsteps=10000, lammps_input=None, lammps_struct=None, temp=None,
                            atomic_species=None, relaunch=True, preselect_fname=None, dry_run=False):

        self.set_lammps_variables(lammps_path, md_nsteps, lammps_input, lammps_struct, temp, atomic_species, relaunch)
        if dry_run:
            print('Dry run finished, all paths were checked. Proceed with training.')
            quit()

        if preselect_fname:
            self.preselect_fname = preselect_fname

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
                    self.run('cp current.almtp {}/final.almtp'.format(self.owd))
                    self.run('cp train.cfg {}/final.cfg'.format(self.owd))
                    self.run('echo "  Extrapolative configurations were preselected but max MD steps was reached.\n  Stopping per user request.">>iter_output.txt')
                    os.chdir(self.owd)
                else:
                    os.chdir(self.owd)

            else:
                finished = True
                self.run('cp prev.almtp {}/final.almtp'.format(self.owd))
                self.run('cp ../{}/train.cfg {}/final.cfg'.format(self.iterstep-1, self.owd))
                self.run('echo "  No extrapolative configuration were preselected.">>iter_output.txt')
                os.chdir(self.owd)

        logging.info("{}: OTF learning scheme completed after {} steps".format(when_is_now(), self.iterstep))


    def mdrun_select(self):
        # Run MD trajectory in selection mode
        logging.info('    {}: Running MD trajectory...'.format(when_is_now()))
        self.run('echo "  Running MD trajectory...">>iter_output.txt')
        self.run(['cp ../{}/current.almtp prev.almtp'.format(self.iterstep-1)])
        self.run(f'touch {self.preselect_fname}')
        command = 'srun {} -v SEED 1 -v T {} -v NSTEP {} -v MLIPPATH {} -v STRUCT {} -log none -in {} &>lammps.log'.format(
                self.lammps, self.temperature, self.mdsteps, self.mlip_path, self.lammps_struct, self.lammps_input)
        self.run(command)

        #Check if nsteps was reached
        out = self.check_process(f'grep "Breaking threshold exceeded" lammps.log')
        if out.returncode == 0:
            self.max_nsteps_reached = False
            self.run('echo "        Breaking threshold exceeded">>iter_output.txt')
        else:
            self.max_nsteps_reached = True
            self.run('echo "        Max MD steps was reached without exceeding breaking threshold">>iter_output.txt')


    def select_configs(self):
        self.run('echo "  Selecting configurations...">>iter_output.txt')
        command = (f'{self.mtp} select_add ../{self.iterstep-1}/current.almtp ../{self.iterstep-1}/train.cfg {self.preselect_fname} add_to_train.cfg>>iter_output.txt')
        self.run(command)

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


    def set_mlip_variables(self, training_flags):

        if training_flags:
            if not os.path.exists(training_flags):
                raise FileNotFoundError(f'mlip_flags file {training_flags} not found')
            self.set_mlip_flags(training_flags)
        else:
            self.training_flags = None


    def check_training_mode(self, mode):

        if mode == 'cfg':
            self.training_mode = 'cfg'
        elif mode =='nbh':
            raise NotImplementedError('neighborhood OTF training mode not implemented yet.')
