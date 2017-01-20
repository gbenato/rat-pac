'''
Have a little qsubb module that can tell once all qsubb jobs are finished.
Necessery in order to submit jobs once all preivous jobs are finished.

'''

import subprocess
import time
import os, glob

def check_num_qsub():
  res =  subprocess.check_output(['qstat', '-u', 'benschmidt'])
  list_res = res.splitlines()
  count_jobs = 0
  for line in list_res:
    if 'benschmidt' in line and not 'STDIN' in line:
      count_jobs += 1
  print 'The number of currently running jobs is', count_jobs
  #print 'Print output', type(res), res
  return count_jobs

def check_list_qsub(job_ids):
  res =  subprocess.check_output(['qstat','-f', '-u', 'benschmidt'])
  list_res = res.splitlines()
  count_jobs = 0
  for idx,line in enumerate(list_res):
#    print idx, line
    if 'Job_Name = ' in line:
      a_job = line.split('Job_Name = ')[1]
      if a_job in job_ids:
        count_jobs += 1
  print 'The number of currently running jobs of the list is', count_jobs,'/', len(job_ids) 
#  print 'Print output', type(res), res
  return count_jobs

def check_completed_pbs_out(job_ids):
  out_path = '/nfs/cuore1/scratch/schmidtb/qsub_script/output/'
  # I would better get out_path from my qsub script!
  n_success = 0 
  l_failed = []
  n_total = len(job_ids)
  for el in job_ids:
    b_success = False
    if not '.o' in el:
      el += '.o'
    if os.path.basename(el) == el:
      el = out_path+el
    f = open(el)
    l = f.readlines()
    if len(l)> 10:
      l = l[-10:]
    for line in l:
      if 'Exit successfully' in line:
        n_success +=1
        b_success = True
        break
    if not b_success:
      l_failed.append(el)
  print '**********************************'
  print 'The following jobs failed:'
  print l_failed
  print n_success, ' out of ', n_total, ' jobs successfull.'
  print 'Length of the failed list',len(l_failed)
  #splitting by ds:
  res_by_ds = {}
  for el in l_failed:
    ds = el[el.find('ds'):el.find('ds')+6]
    if res_by_ds.has_key(ds):
      res_by_ds[ds].append(el)
    else:
      res_by_ds[ds] = [el]
  for lists in res_by_ds.itervalues():
    print lists, len(lists)

def wait_for_qsub():
  t0 = time.time()
  time.sleep(1)
  n_proc = check_num_qsub()
  while n_proc != 0:
    time.sleep(30)
    print 'Waiting for qsub job completion', (time.time()-t0)/60.0, 'minutes'
    n_proc = check_num_qsub()
  print 'All processes finshed at time', time.time()-t0

def wait_for_jobs(a_list):
  t0 = time.time()
  time.sleep(1)
  n_proc = check_list_qsub(a_list)
  while n_proc != 0:
    time.sleep(30)
    print 'Waiting for qsub job completion', int((time.time()-t0)/60.0), 'minutes', int((time.time()-t0)%60), 'sec' 
    n_proc = check_list_qsub(a_list)
  print 'All processes finshed at time', time.time()-t0

def check_entire_proc(pbs_out='/nfs/cuore1/scratch/schmidtb/qsub_script/output/', proc_name ='bi_standard*' ):
    proc_end = ['pre', 'avg_pulse', 'noise_ds', 'noise_run', 'amp', 'stab', 'cal', 'energize', \
               	'co_time', 'shape', 'jitter', 'coinc', 'reduce']
    for proc in proc_end:
      print pbs_out+proc_name+proc+'.pbs.o'
      list_out_files = glob.glob(pbs_out+proc_name+proc+'.pbs.o')
      check_completed_pbs_out(list_out_files)



if __name__ == '__main__':
#  check_num_qsub()    
#  wait_for_qsub()
  check_entire_proc()
  print '\n\n ************************'
  print 'Failed jobs in bi_sat_stab'
  check_entire_proc('/nfs/cuore1/scratch/schmidtb/qsub_script/output/', 'bi_sat_stab*')
#  wait_for_jobs(['2016_05_25_00_27_34_queue_py.pbs', '2016_05_25_00_27_12_queue_py.pbs'])
  
