'''
Author: Benjamin Schmidt
Date: 2017/04/12

Analysis tool for PMT SPE and few PE spectrum fits for CHESS calibration and pedestal correction
'''

import ROOT
import datetime
import os, glob
import time, math
import numpy as np
import json


chess_pmt_name = ["Control Channel", "Light PMT 1", "Light PMT 2", "Light PMT 3" ,"Light PMT 4", "Light PMT 5", "Light PMT 6", \
	"North Floor Panel", "South Floor Panel", "North Side Panel", "East Side Panel", "Ring PMT 0", "Ring PMT 1", \
	"Ring PMT 2", "Ring PMT 3", "Ring PMT 4", "Ring PMT 5", "Ring PMT 6", "Ring PMT 7", "Ring PMT 8", "Ring PMT 9", \
	"Ring PMT 10", "Ring PMT 11" , "Top Muon Tag", "Bottom Muon Tag", "Trigger PMT"]



def fit_pmt4(  h,  spe_mean=1.0, spe_sigma=1.01, x_min=0.315 ):
    '''
    Single method to fit a 4 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :param x_min:
    :return:
    '''
    x_max = spe_mean*6
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int(spe_mean/dV))
    hFit = ROOT.TF1("hFit",("[0]*exp(-((x-[1])^2)/(2*[2]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*2*[2]^2))"
                       "+ [5]*exp(-((x-3*[1])^2)/(3*2*[2]^2))+ [6]*exp(-((x-4*[1])^2)/(4*2*[2]^2))"),x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "tpeN");
    hFit.SetParName(6, "qpeN");

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.4 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 3.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, 0.1 * spe_evts);
#    hFit.SetParLimits(5, 0, 1000);
    hFit.SetParameter(6, 0.01* spe_evts);

    h.Fit(hFit, 'R')
#    raw_input('Hit enter to quit the fit histogram')
    return hFit

def fit_pmt2(  h,  spe_mean=1.03, spe_sigma=1.01, x_min=0.315 ):
    '''
    Single method to fit a 4 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :param x_min:
    :return:
    '''
    x_max = spe_mean*4
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int(spe_mean/dV))
    fit_res = {}
    fit_formula = "[0]*exp(-((x-[1])^2)/(2*[2]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*2*[2]^2))"
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = x_max 
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 3.5); # scaling of the position of the second and third...

    h.Fit(hFit, 'R')
    raw_input('Hit enter to quit the fit histogram')
    fit_res["par_array"] = hFit.GetParameters()
    fit_res["par_errors"] = hFit.GetParErrors()
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    return fit_res


def fit_gaus(h):
    mean = h.GetMean()
    rms = h.GetRMS()
    dV = h.GetBinWidth(1)
    amp = h.GetBinContent(int(mean/dV))
    fit = ROOT.TF1("gaus", "gaus", mean - 2* rms, mean + 2 * rms)
    fit.SetParameter(0, amp)
    fit.SetParameter(1, mean)
    fit.SetParameter(2, rms)
    h.Fit(fit, 'R')
#    raw_input('Hit enter to continue')
    return fit

def fit_pmt3_bg(  h,  outdir):
    '''
    Single method to fit an exp bg and 3 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
    spe_mean = h.GetMean()
    spe_sigma = h.GetRMS()
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
    x_max = 3*spe_mean
    x_min_bin = h.GetMaximumBin()
    x_min = h.GetXaxis().GetBinCenter(x_min_bin)+.8*f_i.GetParameter(2)
    par6 = (f_i.GetParameter(2)*fwhm_sigma/2.0)/math.log(2.0)/2.0
    par5 = h.GetMaximum()*math.exp(x_min/par6)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))/2.0
    print 'spe_evts', spe_evts, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-x/[6]) + [0]*exp(-((x-[1])^2)/(2*[2]^2)) +  [3]*exp(-((x-[1]*[4])^2)/(2*2*[2]^2)) + [7]*exp(-((x-[1]*3/2*[4])^2)/(3*2*[2]^2))"
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = x_max 
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    hFit.SetParName(7, "tpeN");

    print 'Initial parameters for the fit'
    print 'spe_pos', spe_mean, 'spe_width', spe_sigma, 'spe_amp', spe_evts
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts)
    hFit.SetParameter(1, 3*spe_mean)
    hFit.SetParLimits(1, .5, 2.0)
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    hFit.SetParameter(7, 0.01 * spe_evts);
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    new_xmax = hFit.GetParameter(1)*hFit.GetParameter(4)*3/2.0
    hFit.SetRange(x_min, new_xmax)
    h.Fit(hFit, 'R')
#   Assume the fit is good and constrain the range such that we cut at the center of the last gaussian

    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'tpe_fit_gaussbg'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res

def fit_pmt3_gaussbg(  h,  outdir):
    '''
    Single method to fit an exp bg and 3 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)-4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
    spe_mean = 1.2
    spe_sigma = (f_i.GetParameter(2))
#    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
    x_max = 10
    x_min_bin = h.GetMaximumBin()
    x_min = -.7
    par6 = f_i.GetParameter(2)
    par5 = f_i.GetParameter(0)
    par7 = f_i.GetParameter(1)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))*dV
    print 'spe_evts', spe_evts, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/(2*([2]^2+[6]^2))) +  [3]*exp(-((x-[4]*[1])^2)/(2*(2*[2]^2+[6]^2))) + [8]*exp(-((x-[1]*3/2*[4])^2)/(2*(3*[2]^2+[6]^2)))"
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    hFit.SetParName(7, "bgPos");
    hFit.SetParName(8, "tpeN");

    print 'Initial parameters for the fit'
    print 'spe_pos', 2*spe_sigma, 'spe_width', spe_sigma, 'spe_amp', spe_evts
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    hFit.SetParameter(7, par7)
    hFit.SetParameter(8, 0.01 * spe_evts);
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    new_xmax = hFit.GetParameter(1)*hFit.GetParameter(4)*4.0/2.0
    hFit.SetRange(x_min, new_xmax)
    h.Fit(hFit, 'R')
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"]= new_xmax 
    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'tpe_fit_'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res


def fit_pmt4_gaussbg(  h,  outdir):
    '''
    Single method to fit an exp bg and 3 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)-4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
#    spe_mean = h.GetMean()
    spe_mean = 1.2
    spe_sigma = f_i.GetParameter(2) # I guess assuming that spe_mean is already at more or less the correct place I could do sth. more elaborate here
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
#    x_max = spe_mean+10*spe_sigma
#    x_min_bin = h.GetMaximumBin()
#    x_min = h.GetXaxis().GetBinCenter(x_min_bin)-2*f_i.GetParameter(2)
    x_min= -.7
    x_max =10
    par6 = f_i.GetParameter(2)
    par5 = f_i.GetParameter(0)
    par7 = f_i.GetParameter(1)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))*dV
    print 'spe_evts', spe_evts/dV, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/2/([2]^2+[6]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*(2*[2]^2+[6]^2))) + [8]*exp(-((x-[1]*3/2*[4])^2)/(2*(3*[2]^2+[6]^2))) + [9]*exp(-((x-[1]*4/2*[4])^2)/(2*(4*[2]^2+[6]^2)))"
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    hFit.SetParName(7, "bgPos");
    hFit.SetParName(8, "tpeN");
    hFit.SetParName(9, "qpeN");

    print 'Initial parameters for the fit'
    print 'spe_pos', 2*spe_sigma, 'spe_width', spe_sigma, 'spe_amp', spe_evts/10.0
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParLimits(1, .5, 3)
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    hFit.SetParameter(7, par7)
    hFit.SetParameter(8, 0.01 * spe_evts);
    hFit.SetParameter(9, 0.001 * spe_evts);
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    new_xmax = hFit.GetParameter(1)*hFit.GetParameter(4)*5.0/2.0
    hFit.SetRange(x_min, new_xmax)
    h.Fit(hFit, 'R')
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = new_xmax 
    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'qpe_fit_'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res

def fit_pmt4_gaussbg_const_width(  h,  outdir):
    '''
    Single method to fit an exp bg and 3 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)-4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
#    spe_mean = h.GetMean()
    spe_mean = 1.2
    spe_sigma = h.GetRMS()
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
#    x_max = spe_mean+10*spe_sigma
#    x_min_bin = h.GetMaximumBin()
#    x_min = h.GetXaxis().GetBinCenter(x_min_bin)-2*f_i.GetParameter(2)
    x_min= -.7
    x_max =10
    par6 = f_i.GetParameter(2)
    par5 = f_i.GetParameter(0)
    par7 = f_i.GetParameter(1)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))*dV
    print 'spe_evts', spe_evts/dV, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/(2*[2]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*[2]^2)) + [8]*exp(-((x-[1]*3/2*[4])^2)/(2*[2]^2)) + [9]*exp(-((x-[1]*4/2*[4])^2)/(2*[2]^2))"
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    hFit.SetParName(7, "bgPos");
    hFit.SetParName(8, "tpeN");
    hFit.SetParName(9, "qpeN");

    print 'Initial parameters for the fit'
    print 'spe_pos', 2*spe_sigma, 'spe_width', spe_sigma, 'spe_amp', spe_evts/10.0
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParLimits(1, .5, 3)
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    hFit.SetParameter(7, par7)
    hFit.SetParameter(8, 0.01 * spe_evts);
    hFit.SetParameter(9, 0.001 * spe_evts);
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    new_xmax = hFit.GetParameter(1)*hFit.GetParameter(4)*5.0/2.0
    hFit.SetRange(x_min, new_xmax)
    h.Fit(hFit, 'R')
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = new_xmax 
    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'qpe_fit_const_width'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res


def fit_pmt2_bg(  h,  outdir):
    '''
    Single method to fit an exp bg and 2 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
    spe_mean = 1.2
    spe_sigma = f_i.GetParameter(2)
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
    x_max = 3*spe_mean
    x_min_bin = h.GetMaximumBin()
    x_min = h.GetXaxis().GetBinCenter(x_min_bin)+.8*f_i.GetParameter(2)
    par6 = (f_i.GetParameter(2)*fwhm_sigma/2.0)/math.log(2.0)/2.0
    par5 = h.GetMaximum()*math.exp(x_min/par6)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))/2.0
    print 'spe_evts', spe_evts, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-x/[6]) + [0]*exp(-((x-[1])^2)/(2*[2]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*2*[2]^2))"
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = x_max 
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    print 'Initial parameters for the fit'
    print 'spe_pos', spe_mean, 'spe_width', spe_sigma, 'spe_amp', spe_evts
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'dpe_fit_'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res


def fit_pmt2_gaussbg(  h,  outdir):
    '''
    Single method to fit an exp bg and 2 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    c = ROOT.TCanvas()
    fit_res = {}
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
    spe_mean = 1.2
    spe_sigma = f_i.GetParameter(2)
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
    x_max = 10
    x_min_bin = h.GetMaximumBin()
    x_min = -.7
    par6 = f_i.GetParameter(2)
    par5 = f_i.GetParameter(0)
    par7 = f_i.GetParameter(7)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))*dV
    print 'spe_evts', spe_evts, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/(2*([2]^2+[6]^2))) +  [3]*exp(-((x-[4]*[1])^2)/(2*(2*[2]^2+[6]^2)))"
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "dpeN");
    hFit.SetParName(4, "dpeScaling");
    hFit.SetParName(5, "bgAmp");
    hFit.SetParName(6, "bgWidth");
    print 'Initial parameters for the fit'
    print 'spe_pos', spe_mean, 'spe_width', spe_sigma, 'spe_amp', spe_evts
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, 0.1 * spe_evts);
    hFit.SetParameter(4, 2);
#   hFit.FixParameter(4, 2);
    hFit.SetParLimits(4, 1.5, 2.5); # scaling of the position of the second and third...
    hFit.SetParameter(5, par5)
    hFit.SetParameter(6, par6)
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)
#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    new_xmax = hFit.GetParameter(1)*hFit.GetParameter(4)*3.0/2.0
    hFit.SetRange(x_min, new_xmax)
    h.Fit(hFit, 'R')
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = new_xmax 

    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'dpe_fit_gaussbg_'+h.GetName()+'.pdf'
      c.Print(plotf)
    print len(hFit.GetParameters())
    print 'spinach', hFit.GetNumberFreeParameters(),  hFit.GetParameters()
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res

def fit_pmt1_bg(h, outdir):
    '''
    Single method to fit an exp bg and 2 gauss model to the PMT spectrum
    :param h:
    :param spe_mean:
    :param spe_sigma:
    :return:
    '''
    fit_res = {}
    c = ROOT.TCanvas()
    f_i = fit_gaus(h)
    fit_res["noise_gaus_par"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParameters(), dtype=float))
    fit_res["noise_gaus_err"] = list(np.ndarray(f_i.GetNumberFreeParameters(), buffer = f_i.GetParErrors(), dtype=float))
    fit_res["noise_chi2"] = f_i.GetChisquare()/f_i.GetNDF()
    fit_res["tot_events"] = h.Integral()
    fit_res["bin_width"] = h.GetXaxis().GetBinWidth(1)
    print h.GetMean()
    h_xmin = h.GetXaxis().GetXmin()
    h.GetXaxis().SetRangeUser(f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax() )
    print f_i.GetParameter(1)+4*f_i.GetParameter(2), h.GetXaxis().GetXmax(), 'Above noise', h.GetMean(), h.GetRMS()
    spe_mean = h.GetMean()
    spe_sigma = h.GetRMS()
    h.GetXaxis().SetRangeUser(h_xmin, h.GetXaxis().GetXmax())

    fwhm_sigma = 2.3548
    x_max = 3*spe_mean
    x_min_bin = h.GetMaximumBin()
    x_min = h.GetXaxis().GetBinCenter(x_min_bin)+.8*f_i.GetParameter(2)
    par6 = (f_i.GetParameter(2)*fwhm_sigma/2.0)/math.log(2.0)/2.0
    par5 = h.GetMaximum()*math.exp(x_min/par6)
    dV = h.GetBinWidth(1)
    spe_evts = h.GetBinContent(int((spe_mean-h_xmin)/dV))/2.0
    print 'spe_evts', spe_evts, 'in bin ', int((spe_mean-h_xmin)/dV)
    fit_formula = "[3]*exp(-x/[4]) + [0]*exp(-((x-[1])^2)/(2*[2]^2)) "
    fit_res["formula"] = fit_formula
    fit_res["x_min"] = x_min
    fit_res["x_max"] = x_max 
    hFit = ROOT.TF1("hFit",fit_formula,x_min,x_max)

    # Define the parameters and set initial values
    hFit.SetParName(0, "speN");
    hFit.SetParName(1, "speMean");
    hFit.SetParName(2, "speSigma");
    hFit.SetParName(3, "bgAmp");
    hFit.SetParName(4, "bgWidth");
    print 'Initial parameters for the fit'
    print 'spe_pos', spe_mean, 'spe_width', spe_sigma, 'spe_amp', spe_evts
    print 'bg_amp', par5, 'bg_decay', par6, 'xmin', x_min, 'xmax', x_max

    hFit.SetParameter(0, spe_evts);
    hFit.SetParameter(1, spe_mean);
    hFit.SetParameter(2, spe_sigma);
    hFit.SetParameter(3, par5);
    hFit.SetParameter(4, par6);
    h.Draw()
    hFit.Draw("SAME")
    ROOT.gPad.SetLogy(1)

#    raw_input("Initial parameters drawn")

    h.Fit(hFit, 'R')
    c.Update()
    raw_input('Hit enter to quit the fit histogram')
    if outdir:
      plotf = outdir + 'spe_fit_'+h.GetName()+'.pdf'
      c.Print(plotf)
    fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
    fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))

    if hFit.GetNDF():
      fit_res["chi2"] = hFit.GetChisquare()/hFit.GetNDF()
    else:
      fit_res["chi2"] = -99
    del c
    return fit_res


def fit_pol1(g, out_dir):
  hFit = ROOT.TF1("pol1", "[0]+[1]*x", -60, 60)
  fit_res = {}
  c = ROOT.TCanvas()
  g.Draw("AP")
  g.Fit(hFit, 'R')
  out_name = out_dir+"charge_correction_"+g.GetName()+".pdf"
  c.Print(out_name)
  raw_input("Hit enter to continue")
  fit_res["formula"] = "[0]+[1]*x"
  fit_res["par_array"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParameters(), dtype=float))
  fit_res["par_errors"] = list(np.ndarray(hFit.GetNumberFreeParameters(), buffer = hFit.GetParErrors(), dtype=float))
  return fit_res

def get_avg_number_pes(h, spe_pos, x_min, x_max ):
    '''
    Calculate the avg. number of photoelectron in the given PMT spectrum. In case you analyze cobodaq data and use a
    low avg photo number the value will be very biased.
    You cannot acount for the noise e.g. 0pe events as you are triggering above the noise...
    You could however get the real avge pe values from the relative intensities and match up with poisson statistics
    :param h: PMT spectrum
    :param spe_pos: Photo electron peak position in V
    :param x_min: integration boundaries
    :param x_max: integration boundaries
    :return: float avg number of pes
    '''
    n_evts = 0
    idx = 0
    dV = h.GetBinWidth(1)
    x = h.GetBinCenter(idx)
    integral = 0.0
    while x < x_min:
        idx +=1
        x = h.GetBinCenter(idx)
    while x < x_max:
        integral += h.GetBinCenter(idx)*h.GetBinContent(idx)
        n_evts += h.GetBinContent(idx)
        x = h.GetBinCenter(idx)
        idx +=1
    avg_pe = integral/float(n_evts)/spe_pos
    return avg_pe


def get_gauss_integral(mean, sigma, amp, dV= 1.0):
   f = ROOT.TF1('my_gauss','gaus',0, 1 )
   f.SetParameter(0, amp)
   f.SetParameter(1, mean)
   f.SetParameter(2, sigma)
   return f.Integral(0,1)*1.0/dV



def read_chess_pmt_histos(chess_file_name):
    '''
    Read in the root histos for the 9 chess PMTs of the ring 1-inch cross PMTs
    plus while we are at it we could in principle also do the light yield pmts
    :return a dictionary of thistograms
    '''
    file = ROOT.TFile(chess_file_name, 'r')
    file_keys = file.GetListOfKeys()
    a_dict = {}
    for key in file_keys:
      print key.GetName()
      if 'h_charge' in key.GetName() and 'Ring PMT' in key.GetTitle():
        a_dict[key.GetName()] = {}
        a_dict[key.GetName()]['title'] = key.GetTitle()
        a_dict[key.GetName()]['hist'] = file.Get(key.GetName())
        a_dict[key.GetName()]['hist'].SetDirectory(0)
    print a_dict
    return a_dict


def read_chess_charge_correction_graphs(chess_file_name):
    '''
    Read in the root histos for the 9 chess PMTs of the ring 1-inch cross PMTs
    plus while we are at it we could in principle also do the light yield pmts
    :return a dictionary of thistograms
    '''
    file = ROOT.TFile(chess_file_name, 'r')
    file_keys = file.GetListOfKeys()
    a_dict = {}
    for key in file_keys:
      print key.GetName()
      if 'g_charge' in key.GetName() and 'Ring PMT' in key.GetTitle():
        a_dict[key.GetName()] = {}
        a_dict[key.GetName()]['title'] = key.GetTitle()
        a_dict[key.GetName()]['hist'] = file.Get(key.GetName())
    print a_dict
    return a_dict


def run_chess_calibration(chess_file, plot_dir = False, outf = 'CHESS_PMT_correction_2017_04_11.json', fit_key = 'fit_result_dpe', l_npe_fits = [1,2,3,4]):
    hists = read_chess_pmt_histos(chess_file)
    write_json({'test':'success'}, outf)
    print hists
    spe_mean = np.ones(26)
    noise_offset = np.zeros(26)
    for key, adict in hists.iteritems():
       print adict
       h = adict['hist']
       if 4 in l_npe_fits:
         print 'Fitting 4PE model with const+varying width component',h
         adict['fit_result_qpe_gaussbg'] = fit_pmt4_gaussbg(h, plot_dir)
         print 'Fitting 4PE model with const width component',h
         adict['fit_result_qpe_gaussbg_const_res'] = fit_pmt4_gaussbg_const_width(h, plot_dir)
       if 3 in l_npe_fits:
         adict['fit_result_tpe_gaussbg'] = fit_pmt3_gaussbg(h, plot_dir)
         adict['fit_result_tpe'] = fit_pmt3_bg(h, plot_dir)
       if 2 in l_npe_fits:
         adict['fit_result_dpe_gaussbg'] = fit_pmt2_gaussbg(h, plot_dir)
         adict['fit_result_dpe'] = fit_pmt2_bg(h, plot_dir)
       if 1 in l_npe_fits:
         adict['fit_result_spe'] = fit_pmt1_bg(h, plot_dir)
       adict['hist'] = 'balabala'
       spe_mean[int((h.GetName().split('_'))[-1])] = adict[fit_key]['par_array'][1]
       noise_offset[int((h.GetName().split('_'))[-1])] = adict[fit_key]['noise_gaus_par'][1]        
    print hists
    print 'spe_mean:', list(spe_mean)
    print 'noise_offset:', list(noise_offset)
    write_json(hists, outf)

def display_calib_correction(json_file, fit_key):
  ajson = read_json(json_file)
  spe_mean = np.ones(26)
  noise_offset = np.zeros(26)
  for key, adict in ajson.iteritems():
     spe_mean[int((key.split('_'))[-1])] = adict[fit_key]['par_array'][1]
     noise_offset[int((key.split('_'))[-1])] = adict[fit_key]['noise_gaus_par'][1]        
  print "spe means"
  print list(spe_mean)
  print "noise_offset"
  print list(noise_offset)


def run_chess_charge_correction(chess_file, plot_dir = False, outf = 'CHESS_charge_correction_2017_04_12.json'):
    graphs = read_chess_charge_correction_graphs(chess_file)
    write_json({'test':'success'}, outf)
    print graphs
    slopes = np.ones(26)
    offsets = np.zeros(26)
    for key, adict in graphs.iteritems():
       print adict
       h = adict['hist']
       print h
       adict['fit_result_pol1'] = fit_pol1(h, plot_dir)
       adict['hist'] = 'balabala'
       slopes[int((h.GetName().split('_'))[-1])] = adict['fit_result_pol1']['par_array'][1]
       offsets[int((h.GetName().split('_'))[-1])] = adict['fit_result_pol1']['par_array'][0]
    print graphs
    print 'Slopes:', list(slopes)
    print 'Offsets:', list(offsets)
    write_json(graphs, outf)



def get_dpefit_integral(json_file):
  # Read in the json file, it is sorted already
  pmt_fits = read_json(json_file)
  for key,val in pmt_fits.iteritems():
    print "Extract results for", key
    dpe = val["fit_result_dpe"]
    par_s = dpe["par_array"]
#    print par_s
    f_spe = ROOT.TF1("f_spe","[0]*exp(-((x-[1])^2)/(2*[2]^2))", 0, 5)
    f_dpe = ROOT.TF1("f_dpe","[0]*exp(-((x-[1]*[3])^2)/(2*2*[2]^2))", 0, 5)
    f_spe.SetParameter(0, par_s[0])
    f_spe.SetParameter(1, par_s[1])
    f_spe.SetParameter(2, par_s[2])
    f_dpe.SetParameter(0, par_s[3])
    f_dpe.SetParameter(1, par_s[1])
    f_dpe.SetParameter(2, par_s[2])
    f_dpe.SetParameter(3, par_s[4])
    print 'f_spe integral:', f_spe.Integral(0,5), 'f_dpe integral', f_dpe.Integral(0,5)
    dpe['integral'] = (f_spe.Integral(0,5) + 2*f_dpe.Integral(0,5))*1.0/dpe["bin_width"]
    dpe['avg_pes'] = dpe['integral']/dpe['tot_events']
  write_json(pmt_fits, json_file)
  return pmt_fits


def get_fit_integral(json_file, npe=2, fit_key = "fit_result_dpe"):
  # Read in the json file, it is sorted already
  pmt_fits = read_json(json_file)
  integ_min = -3.0
  integ_max = 10.0
  for key,val in pmt_fits.iteritems():
    print "Extract results for", key
    dpe = val[fit_key]
    par_s = dpe["par_array"]
#    print par_s
    f_spe = ROOT.TF1("f_spe","[0]*exp(-((x-[1])^2)/(2*([2]^2+[3]^2)))", integ_min, integ_max)
    f_dpe = ROOT.TF1("f_dpe","[0]*exp(-((x-[1]*[4])^2)/(2*(2*[2]^2+[3]^2)))", integ_min, integ_max)
    f_tpe = ROOT.TF1("f_tpe","[0]*exp(-((x-[1]*[4]*3.0/2.0)^2)/(2*(3*[2]^2+[3]^2)))", integ_min, integ_max)
    f_qpe = ROOT.TF1("f_qpe","[0]*exp(-((x-[1]*[4]*2.0)^2)/(2*(4*[2]^2+[3]^2)))", integ_min, integ_max)
    
    if npe == 1:
      f_spe.SetParameter(0, par_s[0])
      f_spe.SetParameter(1, par_s[1])
      f_spe.SetParameter(2, par_s[2])
      f_spe.SetParameter(3, 0.0) # par 3 only being used in a few models
      print 'f_spe integral:', f_spe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"]
      dpe['integral'] = f_spe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"]
      dpe['avg_pes'] = dpe['integral']/dpe['tot_events']
    elif npe == 2 and fit_key == 'fit_result_dpe_gaussbg':
# fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/(2*([2]^2+[6]^2))) +  [3]*exp(-((x-[4]*[1])^2)/(2*(2*[2]....
      f_spe.SetParameter(0, par_s[0])
      f_spe.SetParameter(1, par_s[1])
      f_spe.SetParameter(2, par_s[2])
      f_spe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(0, par_s[3])
      f_dpe.SetParameter(1, par_s[1])
      f_dpe.SetParameter(2, par_s[2])
      f_dpe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(4, par_s[4])
      print ('f_spe integral:', f_spe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
	 'f_dpe integral', f_dpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"]) 
      dpe['integral'] = (f_spe.Integral(integ_min,integ_max) + 2*f_dpe.Integral(integ_min,integ_max))*1.0/dpe["bin_width"]
      dpe['avg_pes'] = dpe['integral']/dpe['tot_events']
    elif npe == 3 and fit_key == 'fit_result_tpe_gaussbg':
#  fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/(2*([2]^2+[6]^2))) +  [3]*exp(-((x-[4]*[1])^2)/...
      f_spe.SetParameter(0, par_s[0])
      f_spe.SetParameter(1, par_s[1])
      f_spe.SetParameter(2, par_s[2])
      f_spe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(0, par_s[3])
      f_dpe.SetParameter(1, par_s[1])
      f_dpe.SetParameter(2, par_s[2])
      f_dpe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(4, par_s[4])
      f_tpe.SetParameter(0, par_s[8])
      f_tpe.SetParameter(1, par_s[1])
      f_tpe.SetParameter(2, par_s[2])
      f_tpe.SetParameter(3, par_s[6])
      f_tpe.SetParameter(4, par_s[4])
      print ('f_spe integral:', f_spe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
	 'f_dpe integral',2* f_dpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
          'f_tpe gauss bg integral', 3*f_tpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"])
      dpe['integral'] = (f_spe.Integral(integ_min,integ_max) + 2*f_dpe.Integral(integ_min,integ_max)+3*f_tpe.Integral(integ_min,integ_max))*1.0/dpe["bin_width"]
      dpe['avg_pes'] = dpe['integral']/dpe['tot_events']
    elif npe == 4 and fit_key == 'fit_result_qpe_gaussbg':
#    fit_formula = "[5]*exp(-(x-[7])^2/(2*[6]^2)) + [0]*exp(-((x-[1])^2)/2/([2]^2+[6]^2)) +  [3]*exp(-((x-[4]*[1])^2)/(2*(2*[2]^2+[6]^2))) + [8]*exp(-((x-[1]*3/2*[4])^2)/(2*(3*[2]^2+[6]^2))) + [9]*exp(-((x-[1]*4/2*[4])^2)/(2*(4*[2]^2+[6]^2)))"
      f_spe.SetParameter(0, par_s[0])
      f_spe.SetParameter(1, par_s[1])
      f_spe.SetParameter(2, par_s[2])
      f_spe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(0, par_s[3])
      f_dpe.SetParameter(1, par_s[1])
      f_dpe.SetParameter(2, par_s[2])
      f_dpe.SetParameter(3, par_s[6])
      f_dpe.SetParameter(4, par_s[4])
      f_tpe.SetParameter(0, par_s[8])
      f_tpe.SetParameter(1, par_s[1])
      f_tpe.SetParameter(2, par_s[2])
      f_tpe.SetParameter(3, par_s[6])
      f_tpe.SetParameter(4, par_s[4])
      f_qpe.SetParameter(0, par_s[9])
      f_qpe.SetParameter(1, par_s[1])
      f_qpe.SetParameter(2, par_s[2])
      f_qpe.SetParameter(3, par_s[6])
      f_qpe.SetParameter(4, par_s[4])
      print ('f_spe integral:', f_spe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
	 'f_dpe integral',2* f_dpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
          'f_tpe gauss bg integral', 3*f_tpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"],
           'f_qpe gauss bg integral', 4*f_qpe.Integral(integ_min,integ_max)*1.0/dpe["bin_width"]) 
      dpe['integral'] = (f_spe.Integral(integ_min,integ_max) + 2*f_dpe.Integral(integ_min,integ_max)+3*f_tpe.Integral(integ_min,integ_max)+4*f_qpe.Integral(integ_min,integ_max))*1.0/dpe["bin_width"]
      dpe['avg_pes'] = dpe['integral']/dpe['tot_events']
      
  write_json(pmt_fits, json_file)
  return pmt_fits


def plot_2D_from_fit(json_file, rat_data_file, fit_key = 'fit_result_dpe'):
  ROOT.gSystem.Load("libRATEvent")
  h_ringcandidate_npevspos = ROOT.TH2F("h_ringcandidate_npevspos","h_ringcandidate_npevspos",7,-30.*3.5, 30.*3.5, 7, -30.*3.5, 30.*3.5);

# Get the pmt positions from the actual rat root file
  dsreader = ROOT.RAT.DSReader(rat_data_file);
  tree = dsreader.GetT();
  runT = dsreader.GetRunT();
  run = ROOT.RAT.DS.Run();
  runT.SetBranchAddress("run",run);
  runT.GetEntry(0);
  pmtInfo=run.GetPMTInfo();    
  pmtTypeCount= 0
  pmt_pos = []
  ipmt = 0
  centerpos = ROOT.TVector3(0,0,0)
  for pmt_idx in range(pmtInfo.GetPMTCount()):
      pmtpos = pmtInfo.GetPosition(pmt_idx)
      if(pmtInfo.GetType(pmt_idx) == 1):
        pmt_pos.append(pmtpos)
        centerpos = centerpos + pmtpos
        pmtTypeCount+=1     
  centerpos = centerpos*(1./pmtTypeCount)

  #Read the json file and do the actual plotting
  pmt_fits = read_json(json_file)
  pmt_index_list =[]
  for key in pmt_fits.iterkeys():
    pmt_index_list.append(int(key.split('_')[-1]))
  min_idx = min( pmt_index_list)
  for key,val in pmt_fits.iteritems():
    print "Extract results for", key
    idx = int(key.split('_')[-1])
    idx = idx - min_idx
    print  idx
    h_ringcandidate_npevspos.Fill(pmt_pos[idx].X()-centerpos.X(), pmt_pos[idx].Y()-centerpos.Y(), val[fit_key]['avg_pes']) 
    a_bin =  h_ringcandidate_npevspos.FindBin(pmt_pos[idx].X()-centerpos.X(), pmt_pos[idx].Y()-centerpos.Y())
    print val[fit_key]['integral'], val[fit_key]['tot_events']
    if val[fit_key]['integral'] > 0:
      h_ringcandidate_npevspos.SetBinError(a_bin, math.sqrt(val[fit_key]['integral'])/val[fit_key]['tot_events'])
  h_ringcandidate_npevspos.Draw('colz TEXTE')
  raw_input("Hit enter to continue")
  return
#    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
#    h_ringcandidate_timevspos->Fill(pmtpos.X(), pmtpos.Y(), 1010.);
#    h_ringcandidate_npevspos->Fill(pmtpos.X(), pmtpos.Y(), 1010.);


def read_json(a_file = './output/Nov21Analysis.json'):
    with open(a_file) as json_data:
        d = json.load(json_data)
        return d
    return

def write_json(a_dict, a_file = './output/CHESS_PMT_correction.json'):
    with open(a_file, "w") as outf:
      json.dump(a_dict, outf, indent=4, sort_keys=True)    
    return




if __name__ ==  "__main__":
    ROOT.gStyle.SetOptFit(1)
    chess_dir = '/warehouse/rat_optics_simulation/meas_data/'
#    chess_dir = '/Users/benschmidt/CUORE/data/CHESS_data/'

# config for water data 
#    json_file = 'CHESS_PMT_correction_water_2017_05_04.json'
#    json_file = 'CHESS_PMT_correction_water_2017_05_016.json'
#    rat_file = '/warehouse/rat_optics_simulation/meas_data/cuore-source-watertarget_02May2017-094645_11_cut_Source.root'
#    main_fit_key = 'fit_result_qpe_gaussbg'
#    list_npe_fits = [3,4]
#    main_npe_fit = 4 
#    chess_plot_file = chess_dir+'plots/cuore-source-watertarget_02May2017_cut_Source_no_correction.root'
#    chess_cal_plots = chess_dir + 'plots/plots/'

#config for uvt data 
#    json_file = 'CHESS_PMT_correction_uvt_after_water_correction_2017_05_016.json'
#    rat_file = '/warehouse/rat_optics_simulation/meas_data/cuore-source-uvt-0_15Mar2017-134503_0_cut_Source.root'
#    main_fit_key = 'fit_result_tpe_gaussbg'
#    list_npe_fits = [3,4]
#    main_npe_fit = 3
#    chess_plot_file = chess_dir+'/plots/cuore-source-uvt_all_data_water_correction.root'
#    chess_cal_plots = chess_dir + 'plots/plots/uvt_with_water_correction/'

#config for uva data 
    json_file = 'CHESS_PMT_correction_uva_after_water_correction_2017_05_016.json'
    rat_file = '/warehouse/rat_optics_simulation/meas_data/cuore-source-uva-0_19Apr2017-171022_0_cut_Source.root'
    main_fit_key = 'fit_result_dpe_gaussbg'
    main_npe_fit = 2
    list_npe_fits = [1,2]
    chess_plot_file = chess_dir+'plots/cuore-source-uva-0_19Apr2017_after_water_correction.root'
    chess_cal_plots = chess_dir + 'plots/plots/uva_with_water_correction/'




#run the code
    run_chess_calibration(chess_plot_file, chess_cal_plots , json_file, main_fit_key, list_npe_fits)
#    run_chess_charge_correction(chess_plot_file, './plots/')    
#    get_dpefit_integral(json_file)
#    display_calib_correction(json_file, main_fit_key) # print lists of correction coeffs.
    raw_input('Hit enter to continue')
    get_fit_integral(json_file, main_npe_fit, main_fit_key)
    plot_2D_from_fit(json_file, rat_file, main_fit_key)
