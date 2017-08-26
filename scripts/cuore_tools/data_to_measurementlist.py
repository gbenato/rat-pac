import os,glob

def list_data(a_file):
  print sorted([os.path.basename( x) for x in (glob.glob(a_file))])

if __name__ == "__main__":
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-source-watertarget*.root')
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-source-uva*_cut_Source.root')
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-source-teo2-polish-face3top-teflonhut*cut_Source.root')
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-source-teo2-polish-face3top*cut_Source.root')
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-source-teo2-polish-face1top*cut_Source.root')
#  list_data('/warehouse/rat_optics_simulation/meas_data/cuore-muon-teo2-polish-face3top-teflonhut*cut_Through-going-muons.root')
#  list_data('/warehouse/rat_optics_simulation/TheiaRnD_uvtacrylic_cosmics/root/*.root')
#  list_data('/warehouse/rat_optics_simulation/TheiaRnD_TeO2_rough_cosmics_1/root/*.root')
  list_data('/warehouse/rat_optics_simulation/TheiaRnD_uvtacrylic_90Y/root/*.root')
