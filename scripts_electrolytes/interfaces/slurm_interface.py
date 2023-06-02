from subprocess import run
from datetime import datetime

class SlurmWatcher:

    def __init__(self, rootname, username):

        self.rootname = rootname
        self.username = username


    def watch(self):
        # this is the bash script
        # jobs_left=`squeue -u broussev| grep -v JOBID | grep $1 | wc -l`;  
        command = "squeue -u {}| grep -v JOBID | grep {} |wc -l".format(self.username, self.rootname)
        jobs_left = int(run(command, shell=True, capture_output=True).stdout)
        wait_time = 120

        #for j in range(20):
        while jobs_left>0:
            run("echo {}: {} jobs left>>status.txt".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), jobs_left), shell=True)
            run('sleep {}'.format(wait_time), shell=True)
            jobs_left = int(run(command, shell=True, capture_output=True).stdout)

        # so, get the number of jobs left
        # as long as it's more than 0, wait a bit and check again
        # to be tested


