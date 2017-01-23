import glob, os, sys
import subprocess

# Give a directory as input and loop over all the root files inside to sum up
# what can be detected in each PMT.
analysis_dir = '/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_comiscs_1/'
out_file = '/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_comiscs_1/PMT_out.root'
ana_exe = os.path.expandvars('${RATROOT}/scripts/DrawPlotsForTeO2.exe')

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
  run_analysis(analysis_dir, out_file)
