import glob, os, sys
import subprocess

# Give a directory as input and loop over all the root files inside to sum up
# what can be detected in each PMT.
analysis_dir = '/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_comiscs_1/'
out_file = '/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_comiscs_1/PMT_out_v2.root'
ana_exe = os.path.expandvars('${RATROOT}/scripts/DrawPlotsForTeO2.exe')

analysis_list = ['/warehouse/rat_optics_simulation/results/TheiaRnD_noTeO2_cosmics_1/',
  '/warehouse/rat_optics_simulation/results/TheiaRnD_TeO2_cosmics_1/',
  '/warehouse/rat_optics_simulation/results/TheiaRnD_TeO2_cosmics_2/',
  '/warehouse/rat_optics_simulation/results/TheiaRnD_TeO2_cosmics_scint_1/',
  '/warehouse/rat_optics_simulation/results/TheiaRnD_TeO2_cosmics_scint_2/']

def run_list(ana_list, out_file = 'PMT_out_v2.root' ):
  for a_dir in ana_list:
    print a_dir, out_file
    out_file = a_dir+out_file
    run_analysis(a_dir, out_file)

def run_analysis(analysis_dir, out_file):
  list_sim = glob.glob(analysis_dir+'root/*.root')
  print list_sim
  for a_sim in list_sim:
    command = ana_exe + ' -i '+a_sim + ' -o '+out_file + ' -b'
    print 'Executing command', command
#    subprocess.call(command, shell=True,  executable='/usr/local/bin/zsh')
    subprocess.call(command, shell=True)
#    raw_input('Hit enter to continue')

if __name__ == '__main__':
  run_list(analysis_list )
#  run_analysis(analysis_dir, out_file)
