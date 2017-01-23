
import subprocess
import sys, os, shutil
import fileinput
import glob, time
from scipy.stats import poisson


class RunRAT(object):
  '''
  Packacke all the methods for running a set of Rat simulations into RunRat class. This
  package should be able to read in config files and alter the beamon 
  events parameter as well as the output files.
  It should be capable of spawning a certain number of child processes and feature a timeout and kill
  which should remove the failed process from the output directory.
  '''
  def __init__(self):
    self.rat = "rat"
    self.num_jobs = 1
    self.num_events = 35 # a week worth of muon events
    self.poisson = True # Draw a random possion distributed number for num_events
    self.master_config = False 
    self.out_dir=False
    self.out_file_name=False
    self.config_dir = False 	#/mac subdir in self.out_dir
    self.root_dir = False 	#/root subdir in self.out_dir
    self.proc_time_limit = 600  # proc time limit in seconds
    self.n_proc = 8		# number of allowed subprocesses
    self.killed_procs = 0

  def create_output_dirs(self):
     self.root_dir = self.out_dir+'/root/'
     self.config_dir = self.out_dir+'/mac/'
     self.log_dir = self.out_dir+'/log/'

     if self.out_dir and self.out_file_name and self.master_config:
       if os.path.exists(self.out_dir):
	  print 'Error - the output directory exists and is probably not empty! Aborting'
          a_val = raw_input( 'Hit a to abort anything else to continue: ')
          if 'a' is a_val :
            sys.exit(0)
       else:
         os.makedirs(self.out_dir)
         os.makedirs(self.root_dir)
         os.makedirs(self.config_dir)
         os.makedirs(self.log_dir)
     else:
       print ('Error, one of self.out_dir, self.out_file or self.master_config not specified.', 
	     'Specify these variables first')

  def alter_config(self, n_evts, out_file, mac_file):
    '''
    Create a copy of the master_config stored as mac_file. Mofify the beamon and outfile
    directives.
    '''
    with open(self.master_config, 'r') as f:
      print 'Alter config file:', mac_file 
      if os.path.exists(mac_file):
        print 'Error wanted to write', mac_file, 'but file exists already - Aborting'
        sys.exit(0)
      if os.path.exists(out_file):
        print 'Error ', out_file, 'exists already - Aborting'
        sys.exit(0)
      with open(mac_file, 'w') as outf:
        for line in f.readlines():
 	  if 'beamOn' in line:
	    line = '/run/beamOn '+str(n_evts)
          if 'rat/procset file' in line:
	    line = '/rat/procset file "'+ out_file + '"'
            print line
          outf.write(line)

  def get_last_file(self):
    print self.config_dir
    l_files = os.listdir(self.config_dir)
    if len(l_files) == 0:
      return 0
    l_files = [x.split('_')[-1] for x in l_files]
    l_files = [int(x.split('.')[0]) for x in l_files]
    l_files = sorted(l_files)
    print l_files
    last_number = l_files[-1]
    return last_number+1


  def do_sim(self):
    self.create_output_dirs()
    n = self.get_last_file()
    proc_list = []
    for i in range(n, n+self.num_jobs):
      n_start = self.num_events
      if (self.poisson):
        n_start = poisson.rvs(self.num_events, 1)
      out_file= self.root_dir+self.out_file_name.replace('.root','_'+str(i)+'.root')
      mac_file= self.config_dir+(os.path.basename(self.master_config)).replace('.mac', '_'+str(i)+'.mac')
      log_file = self.log_dir+self.out_file_name.replace('.root','_'+str(i)+'.log')
      print 'do_sim New mac file', mac_file, out_file, n_start
      self.alter_config(n_start, out_file, mac_file)
      #Now one is ready to tun the simulation as a subproccess and multiple can be run in parallel in principle.
      (p,logf) = self.submit_process_run(mac_file, log_file)
      proc_list.append( (p, time.time(), out_file, mac_file, log_file) )
      self.clean_proc_list(proc_list)
      while (len(proc_list) >= self.n_proc):
        time.sleep(20)
        proc_list = self.clean_proc_list(proc_list)
        print 'Waiting for processes to complete', proc_list
    while (len(proc_list) > 0 ):
      time.sleep(60)
      print "Waiting for final job completion"
      proc_list = self.clean_proc_list(proc_list)
    print 'All procs finished or terminated', self.num_jobs
    print 'Num procs killed', self.killed_procs  
      
  def clean_proc_list(self, proc_list):
      clean_list =[]
      for proc in proc_list:
         if (proc[0]).poll() is None: # proc has not finished
            dt = time.time()-(proc[1])
            if dt < self.proc_time_limit:
	      clean_list.append(proc)
            else:
              (proc[0]).kill()
              print "Needed to kill a process, since it exceeded its time_limit", self.proc_time_limit
              print "Removing mac,log and root file", proc[2]
              os.remove(proc[2])
              os.remove(proc[3])
              os.remove(proc[4])
              self.killed_procs +=1
      return clean_list

  def  submit_process_run(self, mac_file, log_file ):   
      command = self.rat+' '+ mac_file
#      command = 'printenv'
#      command = 'echo $DYLD_LIBRARY_PATH'
      my_env = os.environ.copy()
#      subprocess.Popen(my_command, env=my_env, stdout=f, stderr = subprocess.STDOUT)
      print 'Submit command:', command
#      subprocess.call(command, shell=True, stdout=subprocess.PIPE)
      f = open(log_file, 'w')
      p = subprocess.Popen(command, env=my_env, shell=True, stdout=f, stderr = subprocess.STDOUT )
      return (p,f)


def main():
#  ana_dir = sys.argv[1]
#Configure the directories and simulation parameters
  num_jobs= 530
  n_evts = 35
  poisson = True
  mac_file = os.path.expandvars('${RATROOT}/mac/TheiaRnD_mcprod_TeO2_cosmics_v1.mac')
  out_file = 'TheiaRnD_mcprod_TeO2_rough_cosmics_v1_.root'
#  out_dir = os.path.expandvars('${RATROOT}/results/TheiaRnD_TeO2_rough_comiscs_1') 
  out_dir = '/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_comiscs_1'
  a_rat = RunRAT()
  a_rat.num_jobs = num_jobs
  a_rat.num_events = 35
  a_rat.n_proc = 16
  a_rat.poisson = poisson
  a_rat.master_config = mac_file
  a_rat.out_dir = out_dir
  a_rat.out_file_name= out_file
# run the simulation
  print 'Determined macro and output file location from environment:'
  print '\t mac_file \t', mac_file
  print '\t out_dir \t', out_dir
  raw_input('Hit enter to start the simulation with the other parameters you specified in the main method.')
  a_rat.do_sim()


if __name__ == '__main__':
    main()
