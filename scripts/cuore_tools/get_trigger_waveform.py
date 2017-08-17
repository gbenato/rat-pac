import ROOT
import math

data_dir = "/warehouse/rat_optics_simulation/meas_data/"
#data_file = data_dir + "cuore-muon-uvt-0_28Feb2017-161923_merge_cut_Cosmics.root"
data_file = data_dir + "cuore-source-watertarget_02May2017-094645_11_cut_Source.root"
data_file = data_dir + "cuore-source-uvt-1_15Mar2017-145354_1_cut_Source.root"
data_file = data_dir + "cuore-source-teo2-polish-face3top-teflonhut-0_12Apr2017-052757_87_cut_Source.root"
pmt_ids = range(25)


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
  cols = 4
  rows = int(math.ceil(len(pmt_ids)/float(cols)))
  print rows, cols, len(pmt_ids)
  c.Divide(cols, rows)
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
     for i_pmt in range(ev.GetPMTCount()):
        pmt = ev.GetPMT(i_pmt)
        print 'pmtid', pmt.GetID(), '/', ev.GetPMTCount();
        c.cd(pmt.GetID()+1)
#        int pmtType = pmtInfo->GetType(pmtID);
        vPMTDigitizedWaveform = pmt.GetWaveform();
        vWaveformTimes = pmt.GetWaveformTime();
        PMTDigitizedWaveform = ROOT.TGraph()
        for isample in range(vPMTDigitizedWaveform.size()):
            PMTDigitizedWaveform.SetPoint(isample,vWaveformTimes[isample],vPMTDigitizedWaveform[isample])
        PMTDigitizedWaveform.Draw("AL")
        traces.append(PMTDigitizedWaveform)
     c.Update()
     raw_input("Hit enter to continue")



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
