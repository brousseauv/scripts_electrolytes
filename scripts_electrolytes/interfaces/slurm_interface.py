from subprocess import run
from ..utils.time import when_is_now

class SlurmWatcher:

    def __init__(self, rootname, username):

        self.rootname = rootname
        self.username = username


    def watch(self, msg=None):
        if msg:
            run("echo {}>>status.txt".format(msg), shell=True)
        # this is the bash script
        # jobs_left=`squeue -u broussev| grep -v JOBID | grep $1 | wc -l`;  
        command = "squeue -u {}| grep -v JOBID | grep {} |wc -l".format(self.username, self.rootname)
        jobs_left = int(run(command, shell=True, capture_output=True).stdout)
        wait_time = 120

        #for j in range(20):
        while jobs_left>0:
            # Decide if it should be just overwrite (>) or append (>>)
            # I could also echo to the command line as it is intended to run from a terminal maybe using nohup
#            run("echo {}: {} jobs left>>status.txt".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), jobs_left), shell=True)
            run("echo {}: {} jobs left>>status.txt".format(when_is_now(), jobs_left), shell=True)

            run('sleep {}'.format(wait_time), shell=True)
            jobs_left = int(run(command, shell=True, capture_output=True).stdout)

        # so, get the number of jobs left
        # as long as it's more than 0, wait a bit and check again
        # to be tested


def write_slurm_submitfile_loop(args, precommands, command, nloop, calcdir, njobs=1):

    with open('job.sh', 'w') as f:

        f.write('#!/usr/bin/bash\n')
        # enumerate slurm args and write them
        for key, val in args.items():
            f.write(f'#SBATCH {key}={val}\n')
        if njobs>1:
            f.write(f'#SBATCH --array 0-{njobs-1}\n')
        f.write('\n')
        # write precommands
        for line in precommands:
            f.write(f'{line}\n')
        f.write('\n')

        f.write(f'cd {calcdir}\n\n')
        # define job index
        if njobs>1:
            # case n
            f.write(f'maxarray={njobs}\n')
            f.write(f'njobs={nloop}\n')
            f.write('array=$SLURM_ARRAY_TASK_ID\n')
            f.write(f'jobs_per_array=$((njobs/maxarray))\n\n')

            f.write(f'if [ $array -eq $((maxarray-1)) ] && [ $(( (array+1)*jobs_per_array-1 )) -ne $((njobs-1)) ]\n')
            f.write('then\n')
            f.write(f'  for i in $(seq $((0+array*jobs_per_array)) $((njobs-1 )) ); do\n')
            f.write(f'    cd $i\n')
            for line in command:
                f.write(f'    {line}\n')
            f.write('    cd ..\n')
            f.write('  done\n')

            f.write('else\n')
            f.write(f'  for i in $(seq $((0+array*jobs_per_array)) $(( (array+1)*jobs_per_array-1 )) ); do\n')
            f.write(f'    cd $i\n')
            for line in command:
                f.write(f'    {line}\n')
            f.write('    cd ..\n')
            f.write('  done\n')
            f.write('fi\n')

        else:
            # write the loop main command
            f.write(f'for i in $(seq $(()) {nloop-1}); do\n')
            f.write(f'  cd $i\n')
            for line in command:
                f.write(f'  {line}\n')
            f.write('  cd ..\n')
            f.write('done\n')

    f.close()
