import ROOT

ch_bottom = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_1/PMT_out.root'
ch_sides = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_2/PMT_out.root'
scint_sides = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_scint_2/PMT_out.root'
scint_bottom = '/Users/benschmidt/CUORE/analysis/chess-teo2/results/TheiaRnD_TeO2_comiscs_scint_1/PMT_out.root'


f1 = ROOT.TFile(ch_bottom)
t1=f1.Get("t")
f2 = ROOT.TFile(ch_sides)
t2=f2.Get("t")
f3 = ROOT.TFile(scint_bottom)
t3=f3.Get("t")
f4 = ROOT.TFile(scint_sides)
t4=f4.Get("t")

c1 = ROOT.TCanvas()
c1.Divide(2,2)
c1.cd(1)
t1.Draw("pmt3/pmt1>>h1(100,0,1)")
t4.Draw("pmt3/pmt1>>h4(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
c1.cd(2)
t2.Draw("pmt3/pmt1>>h2(100,0,1)")
t3.Draw("pmt3/pmt1>>h3(100,0,1)","", "same")
(ROOT.gPad.GetListOfPrimitives().At(1)).SetLineColor(ROOT.kRed)
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
for idx in range(t1.GetEntries()):
  t1.GetEntry(idx)
  t2.GetEntry(idx)
  t3.GetEntry(idx)
  t4.GetEntry(idx)
  h_pe_ch.Fill((t2.pmt3+t2.pmt2+t2.pmt1)/float(t1.pmt3+t1.pmt2+t1.pmt1))
  h_pe_scint.Fill((t4.pmt3+t4.pmt2+t4.pmt1)/float(t3.pmt3+t3.pmt2+t3.pmt1))
  h_sum_ch.Fill((t1.pmt3+t1.pmt2+t1.pmt1))
  h_sum_ch_sides.Fill((t2.pmt3+t2.pmt2+t2.pmt1))
  h_sum_scint.Fill((t3.pmt3+t3.pmt2+t3.pmt1))
  h_sum_scint_sides.Fill((t4.pmt3+t4.pmt2+t4.pmt1))
  

c2 = ROOT.TCanvas()
c2.Divide(1,2)
c2.cd(1)
h_sum_ch.Draw()
h_sum_ch_sides.SetLineColor(ROOT.kRed)
h_sum_ch_sides.Draw("same")
c2.cd(2)
h_sum_scint.Draw()
h_sum_scint_sides.SetLineColor(ROOT.kRed)
h_sum_scint_sides.Draw("same")
c2.Update()
raw_input('Hit enter to quit')
