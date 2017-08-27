import ROOT
import math
import numpy as np

data_dir = "/warehouse/rat_optics_simulation/meas_data/"
#data_file = data_dir + "cuore-muon-uvt-0_28Feb2017-161923_merge_cut_Cosmics.root"
data_file = data_dir + "cuore-source-watertarget_02May2017-094645_11_cut_Source.root"
data_file = data_dir + "cuore-source-uvt-1_15Mar2017-145354_1_cut_Source.root"
data_file = data_dir + "cuore-source-teo2-polish-face3top-teflonhut-0_12Apr2017-052757_87_cut_Source.root"
data_file = data_dir + "cuore-source-watertarget_17Aug2017-140714_merge_cut_Source.root"
pmt_ids = range(25)
trigger = 23

def open_file(data_file):
  dsr = ROOT.RAT.DSReader(data_file)
  t = dsr.GetT()
  n_events = t.GetEntries()
  runT = dsr.GetRunT()
  return (t, dsr, n_events)

def show_pmt_tree(runT):
  for idx in range(runT.GetEntries()):
      runT.GetEntry(idx)
      arun  = runT.run
      pminfo = arun.GetPMTInfo()
      #return pminfo
      print pminfo.GetPMTCount(), pminfo.GetType(1), pminfo.GetName()
      for i in range(pminfo.GetPMTCount()):
         print i, pminfo.GetDirection(i)
      raw_input("Hit enter to continue")
      return pminfo


def loop_tree(dsr, n_events, pmt_ids):
  c = ROOT.TCanvas('pmt_canvas','pmt traces', 1200,800)
#  cols = 4
#  rows = int(math.ceil(len(pmt_ids)/float(cols)))
#  print rows, cols, len(pmt_ids)
#  c.Divide(cols, rows)
  for i_ev in range(n_events):
     rds = dsr.GetEvent(i_ev)
     ev = rds.GetEV(0); #FIXME: so far get only first event
     print 'Event ID', ev.GetID();
     print 'Event time', ev.GetClockTime();
#    for (int ipmt = 0; ipmt < charge_cut_pmts.size(); ipmt++) {
#      int pmtID = charge_cut_pmts[ipmt];
#      RAT::DS::PMT *pmt = ev->GetPMTWithID(pmtID);
#      if(!pmt==NULL){
#        double charge = pmt->GetCharge();
     traces = []
     i_pmt = trigger
#     for i_pmt in range(ev.GetPMTCount()):
     pmt = ev.GetPMT(i_pmt)
     print 'pmtid', pmt.GetID(), '/', ev.GetPMTCount();
#     c.cd(pmt.GetID()+1)
#        int pmtType = pmtInfo->GetType(pmtID);
     vPMTDigitizedWaveform = pmt.GetWaveform()
     trigger_arr = np.array(vPMTDigitizedWaveform)
     vWaveformTimes = pmt.GetWaveformTime()
     PMTDigitizedWaveform = ROOT.TGraph()
     for isample in range(vPMTDigitizedWaveform.size()):
        PMTDigitizedWaveform.SetPoint(isample,vWaveformTimes[isample],vPMTDigitizedWaveform[isample])
#     trace_arr = PMTDigitizedWaveform.GetY()
     tx,fifty_val = get_trigger_time(trigger_arr, vWaveformTimes)
     trig_line = ROOT.TLine(tx, min(trigger_arr), tx, max(trigger_arr))
     trig_line.SetLineColor(ROOT.kRed)
     PMTDigitizedWaveform.Draw("AL")
     trig_line.Draw("L")
     trig_line.DrawLine(0,fifty_val, 200, fifty_val)
     traces.append(PMTDigitizedWaveform)
     c.Update()
     raw_input("Hit enter to continue")


def get_trigger_time(trigger_arr, time_arr):
   '''
   Make a s imilar algorithm to the TES regime detection algorithm from Brad.
   E.g. walk along the trace from the beginning and end and identify the flat parts.
   Then set 2 brakepoints and finally find the 50% trigger time
   based on the upper and lower plateaus
   '''
   left_arr = trigger_arr[:50]
   right_arr =trigger_arr[-50:]
#   left_arr = np.array(trigger_arr, 10)
#   right_arr = np.array(treigger_arr, trigger_vec.size()-10, trigger_vec.size())
   left_val = np.mean(left_arr)
   right_val = np.mean(right_arr)
   fifty_val = left_val+(right_val - left_val)/2.0
   print left_val, right_val, fifty_val
   trigger_arr = trigger_arr - fifty_val
   abs_trigger = abs(trigger_arr)
   print abs_trigger
   import operator
   min_index, min_value = min(enumerate(abs_trigger), key=operator.itemgetter(1))
   print min_index, min_value, time_arr[min_index]
   print min_index-1, trigger_arr[min_index-1], time_arr[min_index-1]
   if (min_index-1) > 0 and (min_index+1)>0:
     dy = (trigger_arr[min_index+1]-trigger_arr[min_index-1])/2.0
     dt = time_arr[min_index+1]-time_arr[min_index]
     zero_frac = trigger_arr[min_index+1]/dy
   min_pos = time_arr[min_index-1]+zero_frac*dt
   print min_pos
   return min_pos, fifty_val


def main():
  (t, dsr, n_events) = open_file(data_file)
  loop_tree(dsr, n_events, pmt_ids)
#  loop_tree(dsr, n_events)

def get_event(i_ev):
  (t, dsr, n_events) = open_file(data_file)
  rds = dsr.GetEvent(i_ev)
  ev = rds.GetEV(0); #FIXME: so far get only first event
  return ev, rds, dsr

if __name__ == '__main__':
  main()
