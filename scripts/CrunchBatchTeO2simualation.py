import glob, os, sys
import subprocess

# Give a directory as input and loop over all the root files inside to sum up what can be detected in each PMT.
#analysis_dir = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_1/'
#analysis_dir = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_2/'
#analysis_dir = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_scint_1/'
analysis_dir = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_scint_2/'
#out_file = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_1/PMT_out.root'
#out_file = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_2/PMT_out.root'
out_file = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_scint_2/PMT_out.root'
ana_exe = '/Users/benschmidt/CUORE/analysis/chess-teo2/scripts/DrawPlotsForTeo2.exe'

def run_analysis(analysis_dir, out_file):
  list_sim = glob.glob(analysis_dir+'root/*.root')
  print list_sim
  for a_sim in list_sim:
    command = ana_exe + ' -i '+a_sim + ' -o '+out_file
    print 'Executing command', command
    subprocess.call(command, shell=True,  executable='/usr/local/bin/zsh')
#    raw_input('Hit enter to continue')

if __name__ == '__main__':
  run_analysis(analysis_dir, out_file)
