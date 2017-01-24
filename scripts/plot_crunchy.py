import ROOT

#out_dir = /Users/benschmidt/CUORE/analysis/chess-teo2/
out_dir = '/warehouse/rat_optics_simulation/'
ch_bottom = out_dir + 'results/TheiaRnD_TeO2_comiscs_1/PMT_out_v2.root'
ch_sides = out_dir + 'results/TheiaRnD_TeO2_comiscs_2/PMT_out_v2.root'
scint_sides = out_dir + 'results/TheiaRnD_TeO2_comiscs_scint_2/PMT_out_v2.root'
scint_bottom = out_dir + 'results/TheiaRnD_TeO2_comiscs_scint_1/PMT_out_v2.root'
ch_rough_bottom = out_dir + 'TheiaRnD_TeO2_rough_comiscs_1/PMT_out_v2.root'
scint_rough_bottom = out_dir + 'TheiaRnD_TeO2_rough_comiscs_scint_1/PMT_out_v2.root'
ch_noTeO2 = out_dir+'results/TheiaRnD_noTeO2_comiscs_1/PMT_out_v2.root'

f1 = ROOT.TFile(ch_bottom)
t1=f1.Get("t")
f2 = ROOT.TFile(ch_sides)
t2=f2.Get("t")
f3 = ROOT.TFile(scint_bottom)
t3=f3.Get("t")
f4 = ROOT.TFile(scint_sides)
t4=f4.Get("t")
f5 = ROOT.TFile(ch_rough_bottom)
t5=f5.Get("t")
f6 = ROOT.TFile(ch_noTeO2)
t6=f6.Get("t")
f7 = ROOT.TFile(ch_noTeO2)
t7=f7.Get("t")


c1 = ROOT.TCanvas()
c1.Divide(2,2)
c1.cd(1)
t1.Draw("pmt3/pmt1>>h1(100,0,1)")
t4.Draw("pmt3/pmt1>>h4(100,0,1)","", "same")
t5.Draw("pmt3/pmt1>>h5(100,0,1)","", "same")
t6.Draw("pmt3/pmt1>>h6(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
(ROOT.gPad.GetListOfPrimitives().At(2)).SetLineColor(ROOT.kGreen)
(ROOT.gPad.GetListOfPrimitives().At(3)).SetLineColor(ROOT.kCyan)
c1.cd(2)
t2.Draw("pmt3/pmt1>>h2(100,0,1)")
t3.Draw("pmt3/pmt1>>h3(100,0,1)","", "same")
t7.Draw("pmt3/pmt1>>h7(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
(ROOT.gPad.GetListOfPrimitives().At(2)).SetLineColor(ROOT.kGreen)
c1.cd(3)
t1.Draw("pmt2/pmt1>>h11(100,0,1)")
t4.Draw("pmt2/pmt1>>h41(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
c1.cd(4)
t2.Draw("pmt2/pmt1>>h21(100,0,1)")
t3.Draw("pmt2/pmt1>>h31(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
c1.Update()
raw_input('Hit enter to quit')


h_pe_ch = ROOT.TH1F("h_pe_ratio_ch","Ratio of detected pes with acrylic on sides/ acrylic on bottom", 100, 1,2 )
h_pe_scint = ROOT.TH1F("h_pe_ratio_scint_ch","Ratio of detected pes with acrylic on sides/ acrylic on bottom", 100, 1,2 )
h_sum_ch = ROOT.TH1F("h_pe_sum_ch","Sum of detected PEs acrylic on bottom (Cherenkov only)", 100, 0,10000 )
h_sum_ch_sides = ROOT.TH1F("h_pe_sum_ch_sides","Sum of detected pes with acrylic on sides (Cherenkov only)", 100, 0,10000 )
h_sum_scint = ROOT.TH1F("h_pe_scint","Sum of detected pes with acrylic on bottom", 100, 0,10000 )
h_sum_scint_sides = ROOT.TH1F("h_pe_scint_sides","Sum of detected pes with acrylic on sides", 100, 0,10000 )
h_sum_rough = ROOT.TH1F("h_sum_rough", "Sum of detected PEs for a rough TeO2 crystal", 100, 0, 10000)
h_sum_rough_scint = ROOT.TH1F("h_sum_rough_scint", "Sum of detected PEs for a rough scintillating TeO2 crystal", 100, 0, 10000)
h_sum_NoTeO2 = ROOT.TH1F("h_sum_NoTeO2", "Sum of detected PEs for no TeO2 crystal", 100, 0, 10000)

my_range = min([t1.GetEntries(), t2.GetEntries(), t3.GetEntries(), t4.GetEntries(), t5.GetEntries(), t6.GetEntries()])
for idx in range(my_range):
  t1.GetEntry(idx)
  t2.GetEntry(idx)
  t3.GetEntry(idx)
  t4.GetEntry(idx)
  t5.GetEntry(idx)
  t6.GetEntry(idx)
  #At this point I also need to have the number of muon events simulated, 
  # to scale back to 35 muons, but 1st I need to add that to my PMT out trees. 
  h_pe_ch.Fill((t2.pmt3+t2.pmt2+t2.pmt1)*t1.evts_started/t2.evts_started/float(t1.pmt3+t1.pmt2+t1.pmt1))
  h_pe_scint.Fill((t4.pmt3+t4.pmt2+t4.pmt1)*t3.evts_started/t4.evts_started/float(t3.pmt3+t3.pmt2+t3.pmt1))
  h_sum_ch.Fill((t1.pmt3+t1.pmt2+t1.pmt1)*sim_evts/t1.evts_started)
  h_sum_ch_sides.Fill((t2.pmt3+t2.pmt2+t2.pmt1)*sim_evts/t2.evts_started)
  h_sum_scint.Fill((t3.pmt3+t3.pmt2+t3.pmt1)*sim_evts/t3.evts_started)
  h_sum_scint_sides.Fill((t4.pmt3+t4.pmt2+t4.pmt1)*sim_evts/t4.evts_started)
  h_sum_rough.Fill((t5.pmt3+t5.pmt2+t5.pmt1)*sim_evts/t5.evts_started)
  h_sum_rough_scint.Fill((t7.pmt3+t7.pmt2+t7.pmt1)*sim_evts/t7.evts_started)
  h_sum_NoTeO2.Fill((t6.pmt3+t6.pmt2+t6.pmt1)*sim_evts/t6.evts_started)

c2 = ROOT.TCanvas()
c2.Divide(2,2)
c2.cd(1)
h_sum_ch.Draw()
h_sum_ch_sides.SetLineColor(ROOT.kRed)
h_sum_ch_sides.Draw("same")
h_sum_rough.SetLineColor(ROOT.kGreen)
h_sum_rough.Draw("same")
h_sum_NoTeO2.SetLineColor(ROOT.kCyan)
h_sum_NoTeO2.Draw("same")
c2.cd(2)
h_sum_scint.Draw()
h_sum_scint_sides.SetLineColor(ROOT.kRed)
h_sum_scint_sides.Draw("same")
h_sum_rough_scint.SetLineColor(ROOT.kGreen)
h_sum_rough_scint.Draw("same")
c2.cd(3)
h_pe_ch.Draw()
h_pe_scint.SetLineColor(ROOT.kRed)
h_pe_scint.Draw("same")

c2.Update()
raw_input('Hit enter to quit')
