import ROOT
import math
import numpy as np


class TimingSimulator(object):
  def __init__(self, n_pmts, n_events, tau, t_min = None, t_max = None ):
     self.n_pmts = 12
     self.time_arr = np.ones(self.n_pmts, dtype='float')
     self.time_arr_var = []
     self.n_events = n_events
     self.tau = tau
     self.t_min = t_min
     if self.t_min is None:
        self.t_min = -5*tau
     self.t_max = t_max
     if self.t_max is None:
        self.t_max = 5*tau
     print self.t_min, self.t_max
     self.t_residuals = ROOT.TH1F("h_time_residuals", "time residuals", 1000, self.t_min, self.t_max )
     self.rand = ROOT.TRandom3()
     self.scint_pe = 0.05
     self.ch_pe = 0.1
     self.min_pe = 1.0
     self.sigma_pe = 1.0

  def gen_gauss(self):
     self.time_arr_var = []
     for i in range(self.n_pmts):
#        print i, self.time_arr, self.rand.Gaus(0, self.tau)
        self.time_arr[i] = self.rand.Gaus(0, self.tau)
        if self.rand.Gaus(self.ch_pe, self.sigma_pe) > self.min_pe:
          self.time_arr_var.append(self.rand.Gaus(0, self.tau))

  def gen_exp(self):
     self.time_arr_var = []
     for i in range(self.n_pmts):
        self.time_arr[i] = self.rand.Exp( self.tau)
        if self.rand.Gaus(self.scint_pe,self.sigma_pe) > self.min_pe:
          self.time_arr_var.append(self.rand.Exp( self.tau))
 
  def gen_gauss_exp(self):
     self.time_arr_var = []
     for i in range(self.n_pmts):
        if self.rand.Gaus(self.ch_pe, self.sigma_pe) > self.min_pe:
          self.time_arr_var.append(self.rand.Gaus(0, self.tau))
        elif self.rand.Gaus(self.scint_pe,self.sigma_pe) > self.min_pe:
          self.time_arr_var.append(self.rand.Exp( self.tau))
 

  def run_gauss_and_exp(self):
     del self.t_residuals
     self.t_residuals = ROOT.TH1F("h_time_residuals", "time residuals", 1000, self.t_min, self.t_max )
     for i in range(self.n_events):
        self.gen_gauss_exp()
        a_median = ROOT.TMath.Median(len(self.time_arr_var), np.array(self.time_arr_var))
        for j in range(0,len(self.time_arr_var)):
           if self.time_arr_var[j] == a_median:
             continue
           self.t_residuals.Fill(self.time_arr_var[j] - a_median)
    

  def run_gauss(self):
     del self.t_residuals
     self.t_residuals = ROOT.TH1F("h_time_residuals", "time residuals", 1000, self.t_min, self.t_max )
     for i in range(self.n_events):
        self.gen_gauss()
        a_median = ROOT.TMath.Median(len(self.time_arr_var), np.array(self.time_arr_var))
        for j in range(0,len(self.time_arr_var)):
           if self.time_arr_var[j] == a_median:
             continue
           self.t_residuals.Fill(self.time_arr_var[j] - a_median)
        #print a_median
#        a_median = ROOT.TMath.Median(self.n_pmts, self.time_arr)
#        a_median = self.time_arr[0]
#        for j in range(0,self.n_pmts):
#           self.t_residuals.Fill(self.time_arr[j]- a_median)

  def run_exp(self):
     del self.t_residuals
     self.t_residuals = ROOT.TH1F("h_time_residuals", "time residuals", 1000, self.t_min, self.t_max )
     for i in range(self.n_events):
        self.gen_exp()
        a_median = ROOT.TMath.Median(len(self.time_arr_var), np.array(self.time_arr_var))
        for j in range(0,len(self.time_arr_var)):
           if self.time_arr_var[j] == a_median:
             continue
           self.t_residuals.Fill(self.time_arr_var[j] - a_median)

#        a_median = self.time_arr[0]
#        a_median = ROOT.TMath.Median(self.n_pmts, self.time_arr)
#        for j in range(0,self.n_pmts):
#           self.t_residuals.Fill(self.time_arr[j]- a_median)

  

def create_reference_distributions():
  f = ROOT.TFile("median_time_ref_distributions.root", "recreate")
  ts = TimingSimulator(12, 100000, 1)  
  ts.run_gauss()
  ts.t_residuals.Draw()
  ts.t_residuals.SetName("median_gauss_1")
  ts.t_residuals.Write()
  raw_input('That is 12 PMTS, 1000evts, gauss with width of 1, Enter to continue')
  ts.tau = 1
  ts.run_exp()
  ts.t_residuals.Draw()
  ts.t_residuals.SetName("median_exp_1")
  ts.t_residuals.Write()
  raw_input('That is 12 PMTS, 1000evts, exp with width of 1, Enter to continue')
  ts.tau = 4
  ts.run_exp()
  ts.t_residuals.Draw()
  ts.t_residuals.SetName("median_exp_4")
  ts.t_residuals.Write()
  raw_input('That is 12 PMTS, 1000evts, exp with width of 4, Enter to continue')
  ts.tau = 11
  ts.run_exp()
  ts.t_residuals.Draw()
  ts.t_residuals.SetName("median_exp_11")
  ts.t_residuals.Write()
  raw_input('That is 12 PMTS, 1000evts, exp with width of 4, Enter to continue')
  f.Clsoe()

def gauss_plus_exp(n_gaus, n_exp, t_gaus, t_exp, t_min, t_max):
  ts = TimingSimulator(12, n_gaus, t_gaus, t_min, t_max)  
  ts.run_gauss()
  ts.t_residuals.Draw()
  h = ts.t_residuals.Clone()
  h.SetName("gauss_plus_exp")
  ts.tau = t_exp
  ts.n_events = n_exp
  ts.run_exp()
  h.Add(ts.t_residuals)
  return h

def main():
  h2 = gauss_plus_exp(10000, 2000, 2, 10, -40, 40)
  h2.Draw()
  raw_input("Hit enter to continue")


if __name__ == '__main__':
  main()
