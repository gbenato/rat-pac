import subprocess
import sys, os
import time

# Setup your stndard qsub directories - Note they currently needto be on nfs, scratch space
out_dir = '/nfs/cuore1/scratch/schmidtb/qsub_script/output/'
submit_dir = '/nfs/cuore1/scratch/schmidtb/qsub_script/submit/'
#e_mail =NONE

def write_pbs_file( command_str, template='/nfs/cuore1/scratch/schmidtb/qsub_script/example1.pbs', job_name=False):
  '''
  Load a template pbs file, and modify userspecific parameters (input, output, e-mail, command).
  Save the modified file in submit_dir with file_name generated from current time.
  '''
  file_templ = open(template, 'r')
  f_name = submit_dir
  if job_name:
    f_name += job_name
  else:
    t_name =str(time.strftime("%Y_%m_%d_%H_%M_%S"))
    f_name +=  t_name +'_queue_py.pbs'
  pbs_file = open(f_name, 'w')
  for i in range(13):
    a_string = file_templ.readline()
    pbs_file.write(a_string)
  out_name = out_dir+os.path.basename(f_name)
  PBS_err = '#PBS -e localhost:'+ out_name+'.e\n'
  PBS_out = '#PBS -o localhost:'+ out_name+'.o'
  pbs_file.write(PBS_err)
  pbs_file.write(PBS_out)
  string_end = '#end_config'
  while string_end not in a_string:
    a_string = file_templ.readline()
    pbs_file.write(a_string)
  pbs_file.write(command_str)
  pbs_file.close() 
  time.sleep(1) # make sure file_names are unique
  return pbs_file.name



if __name__ == '__main__':
 #read_write_access() #testing to enable afs home directory read (write access)
 template='/nfs/cuore1/scratch/schmidtb/qsub_script/example1.pbs'
 print '**************'
 print 'queue:', sys.argv
 job_name = False
 try:
   idx = sys.argv.index('-n')
   job_name = sys.argv[idx+1]
   argv = sys.argv
   del argv[idx+1]
   del argv[idx]
   print 'Your job name is ', job_name
 except:
   argv = sys.argv
 print 'Your command is', argv[1:]
 if argv[2] == os.path.basename(argv[2]) and argv[2][0] != '-':
   argv[2] = os.getcwd()+'/'+argv[2]
   print 'Expanded argument path to', argv[2]
 my_command = argv[1:]
 command_str = " ".join(my_command)
 print 'You entered command ', command_str
 pbs_file = write_pbs_file(command_str, template, job_name)
 print 'Your pbs file:', pbs_file
 #raw_input('Hit enter to continue')
 qsub_command = ['qsub', '-q', 'cuore', '-V', pbs_file]
 # -q select submission queue, -V export entire environement to batch job
 subprocess.call(qsub_command)

