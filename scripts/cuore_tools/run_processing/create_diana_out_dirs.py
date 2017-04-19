'''
Call 
python create_diana_out_dirs.py a_path 
to create a new diana output folder structure
in the supplied path and copy the cfg files from your
main $CUORE_INSTALL into it as well.
'''

import os, sys
import shutil
import glob
import fileinput
import subprocess
 
def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest, symlinks=False, ignore=None)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

def alter_all_bi_tables(dir, bi_name):
  print "alter_all_bi_tables"
  f_list = glob.glob(dir+"/*.cfg")
  for a_file in f_list:
    print "Change BITable in ", a_file, "to ", bi_name
    alter_bi_table_in_cfg(a_file, bi_name)

def alter_bi_table_in_cfg(cfg_file, bi_name):
  # In principle it would be good to execute this function as part of the copy process
  # but it is not clear which cfg files will be needed --> just run it over all copied
  # .cfg files ... 
  print 'Alter config file:', cfg_file
  print 'Using BiTableName = ', bi_name
  import fileinput
  b_done = False
  for line in fileinput.input(cfg_file, inplace=True):
    if 'BiTableName' in line:
      line = 'BiTableName = '+bi_name+'\n'
      b_done = True
#    print "%s"  % ( line)
    sys.stdout.write (line)
  if b_done:
    return
  for line in fileinput.input(cfg_file, inplace=True):
    if 'filter RejectBadIntervals' in line:
      sys.stdout.write (line)
#      print "%s"  % ( line)
      line = 'BiTableName = '+bi_name+'\n'
      b_done = True
    sys.stdout.write(line)
#    print "%s"  % ( line)
  return


def alter_file_path(list_file, a_path):
  '''
  Input: str list_file 	- full path to list e.g. file Calibration_Production...
	 str a_path	- new path (usually curent directory of the list_file)
  
  Update the list_file for diana data processing
  '''
  print 'Alter file:', list_file
  print 'New PATH ', a_path
  set_path = True
  for line in fileinput.input(list_file, inplace=True):
    if set_path:
      line = a_path+'\n'
      set_path = False
#    print "%s"  % ( line)
    line = line.replace('Blinded', 'Production')
    sys.stdout.write (line)
  return


def main(bi_name="bad_channels", overwrite = True):
  '''
  Create directory structure needed for diana processing.
  Execute: python create_diana_out_dirs 'abs_path'
 
  User needs to adjust ds_lists, old_cfg location,... in source file before use!
  '''
  main_out = sys.argv[1]
  diana_sub = [main_out+'/output', main_out+'/avg', main_out+'/stab']
  diana_calib = [main_out+'/calib', main_out+'/calib/SeedFiles', main_out+'/calib/SeedFiles/corc']
  ds_list = [2061]
#  ds_list = [2070]
#  ds_list = [2049, 2061, 2064, 2067, 2073, 2076, 2079, 2085, 2088, 2091, 2097, 2100, 2103, 2109, 2118, 2124, 2130, 2133, 2139, 2148]
#  ds_list = [2049, 2061, 2064, 2067, 2070, 2073, 2076, 2079, 2085, 2088, 2091, 2097, 2100, 2103, 2109, 2118, 2124, 2130, 2133, 2139, 2148]
  new_cfg = main_out+'/cfg/CUORE0'
#  old_cfg = '/home/benschmidt/scratch/diana_output/bi_sat_noise/cfg/CUORE0'
  old_cfg = os.path.expandvars('$CUORE_INSTALL/cfg/CUORE0')
  cal_seed = 'seed_cal_coeffs_ds'
#  old_seed_d = '/home/benschmidt/scratch/diana_output/bi_sat_noise/calib/SeedFiles/' is this needed?
  old_production_list = '/warehouse/data/CUORE0/OfficialProcessed_v02.30/'
  print "Ready to create dirs"
  print main_out
  for el in diana_sub:
    print el
  print 'Ready to create subdirs for datasets', ds_list
  print 'Ready to copy cfg files'
  print old_cfg, '-->', new_cfg
  raw_input('Hit enter to proceed: ')
  if not os.path.exists(main_out):
    os.mkdir(main_out)
  for el in diana_calib:
    if not os.path.exists(el):
       os.mkdir(el)
  for el in diana_sub:
    if not os.path.exists(el):
      os.mkdir(el)
    for ds in ds_list:
      ds_path = el+'/ds'+str(ds)
      if not os.path.exists(ds_path):
        os.mkdir(ds_path)
  for a_ds in ds_list:
    # Seed files can and probably should be created from the DB instead of such a copy
    print 'Execute script to create seed files:', a_ds
    script_name = os.path.expandvars('$CUORE_INSTALL/pat/calibration/CreateSeedFileFromDB.C')
    seed_path = ', "'+main_out+'/calib/SeedFiles"'
    print seed_path
    cmd = script_name+'('+str(a_ds)+ seed_path + ', "2.30"' + ', "DB")'
    print 'cmd', cmd  
    subprocess.call(['root', '-q','-l', cmd])
#    raw_input('Created seed file from DB '+ seed_path)
#    cal_seed = cal_seed + str(ds)+'.txt'
#    old_seed = old_seed_d + cal_seed
#    new_seed = main_out+'/calib/SeedFiles/'+cal_seed
#    shutil.copyfile(old_seed, new_seed)
  for a_ds in ds_list: # cal and bg Production lists
    old_list_dir = old_production_list +'ds'+str(a_ds)
    new_list_dir = main_out +'/output/ds'+str(a_ds)
    old_cal_list_file =  old_list_dir + '/calibration_Production_ds'+str(a_ds)+'.list.original'
    new_cal_list_file = new_list_dir + '/calibration_Production_ds'+str(a_ds)+'.list'
    old_bg_list_file = old_cal_list_file.replace('calibration', 'background')
    new_bg_list_file = new_cal_list_file.replace('calibration', 'background')
    if not os.path.exists(old_cal_list_file):
      old_cal_list_file =  old_list_dir + '/calibration_Blinded_ds'+str(a_ds)+'.list'
      old_bg_list_file = old_cal_list_file.replace('calibration', 'background')   
    if not os.path.exists(new_cal_list_file) or overwrite:
      shutil.copyfile(old_cal_list_file, new_cal_list_file)
      shutil.copyfile(old_bg_list_file, new_bg_list_file)
      alter_file_path(new_cal_list_file, 'DATAPATH ' + new_list_dir)
      alter_file_path(new_bg_list_file, 'DATAPATH ' + new_list_dir)
#    raw_input('Created Production.lists in '+ new_list_dir)
  if not os.path.exists(new_cfg):
    copyDirectory(old_cfg, new_cfg)
  alter_all_bi_tables(new_cfg, bi_name)
  if not os.path.exists(main_out+'/cuts'):
    os.mkdir(main_out+'/cuts')
  for a_ds in ds_list: # copy txt templatefes from cfg folder:
    # cuts/EnerggyRanges_ds????.txt files (coincidences_timing.cfg)
    new_e_ranges = main_out+'/cuts/EnergyRanges_ds'+str(a_ds)+'.txt'
    if not os.path.exists(new_e_ranges):
      e_ranges = main_out+'/cfg/CUORE0/EnergyRanges.txt'      
      shutil.copyfile(e_ranges, new_e_ranges)
    # calib/Intervals_ds????.txt files (shape_coefficients.cfg)
    new_intervals = main_out+'/calib/Intervals_ds'+str(a_ds)+'.txt'
    if not os.path.exists(new_intervals):
      intervals = main_out+'/cfg/CUORE0/Intervals.txt'
      shutil.copyfile(intervals, new_intervals)


def clean_diana( verbose = True, fake = True):
  ana_dir = '/warehouse/data/CUORE0/ModMongoProc_2017_01_09'
  base_dir = ana_dir 
  prod_files = glob.glob(base_dir+'/output/ds????/Production_2*')
  for a_file in prod_files:
    if verbose:
	print 'Remove', a_file
    if not fake:
      os.remove(a_file)
  stab_files =  glob.glob(base_dir+'/stab/ds????/stab_parameters_*.txt')
  for a_file in	stab_files:
    if verbose:
	print 'Remove', a_file
    if not fake:
      os.remove(a_file)
  cal_files =  glob.glob(base_dir+'/calib/cal_coeffs_ds????.*')
  cal_gui_files = glob.glob(base_dir+'/calib/SeedFiles/gui_ds????.gui')
  for a_file in	cal_files:
    if verbose:
	print 'Remove', a_file
    if not fake:
      os.remove(a_file)
  for a_file in cal_gui_files:
    if verbose:
	print 'Remove', a_file
    if not fake:
      os.remove(a_file)
  avg_files =  glob.glob(base_dir+'/avg/ds????/average_*.root')
  for a_file in	avg_files:
    if verbose:
	print 'Remove', a_file
    if not fake:
      os.remove(a_file)



if __name__ == "__main__":
   main("bad_channels")
#  clean_all()

