//////////////////////////////////////////////////////////////////////
///////////////////////// TheiaRnD Tools /////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// Fitter.cpp created by Javier Caravaca. 1st February 2015
//
// Summary: Fit to measure the collection efficiency of a PMT
//
// Details: This code throws MC simulations changing the collection
//          efficiency of the PMT at each step and finds the best
//          fit value
//
//////////////////////////////////////////////////////////////////////
#include<iostream>
#include <stdio.h>
#include<fstream>
#include <stdlib.h>
#include <ctime>

#include<TROOT.h>
#include<TH1F.h>
#include<TH2D.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<TBrowser.h>
#include<TGraph.h>
#include<TLegend.h>
#include<TMath.h>
#include<TStyle.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include <RAT/DB.hh>

#define DEBUG false
#define BATCH true
#define FITTER "SIMPLEX" //Migrad, SIMPLEX
#define ENABLESR false
#define NBINS 500 //(500) number of bins in the pdfs
#define INITBIN 10 //(10) first bin to consider in the fit
#define LASTBIN 200 //(200) last bin to consider in the fit

namespace TheiaRnD{


class Fitter{

public:

  Fitter(int, char**);
  ~Fitter(){};
  void ParseArgs(int argc, char **argv);
  void GetDataPDFs();
  void GetMCPDFs();
  void GetMCPDFsWithCollEff(double);
  double Likelihood(const double*); //likelihood function
  double ChiSquare(const double*); //chisquare function
  void DoFit();
  void DrawPlots();
  std::vector<TH1F*> hpdf_mc; //prefit MC pdfs
  std::vector<TH1F*> hpdf_dt; //prefit DT pdfs
  std::vector<TH1F*> hpdf_mc_fit; //posfit MC pdfs
  std::vector<TH1F*> hpdf_dt_fit; //posfit DT pdfs
  TH1F* hpdf_mc_sum;
  TH1F* hpdf_mc_sum_fit;


protected:

  //Input files
  std::map<int, char*> fMCFiles; //0->90Sr, 1->90Y
  std::map<int, char*> fDATAFiles; //0->Backgrounds, 1-> Cerenkov
  std::string f_log; //log filename
  bool exist_srfile;
  bool exist_yfile;
  bool exist_datafile;
  bool exist_bkgfile;
  time_t time_now;
  tm *time_local;

  //  TMinuit *fMinuit;

  ROOT::Math::Minimizer* min;
  double fMinLikelihood;
  std::vector<double> fParFitted;
  std::vector<double> fParErrFitted;
  int nlikestep;

};

//Constructor
Fitter::Fitter(int argc, char **argv){

  //Init
  if(BATCH) gROOT->SetBatch();
  nlikestep = 0;
  exist_srfile = false;
  exist_yfile = false;
  exist_datafile = false;
  exist_bkgfile = false;
  
  //Get time
  time_now = time(0);
  time_local = localtime(&time_now);

  //Set log filename
  //sprintf(f_log,"./logs/TheiaRnD_log_%s_%s_%s.txt",time_local->tm_mday,time_local->tm_mon,1900+time_local->tm_year);
  f_log = Form("./logs/TheiaRnD_log.txt");
  
  //Parse arguments and get pdfs
  ParseArgs(argc, argv);
  GetDataPDFs();
  GetMCPDFs();

}

//Parser for main() arguments
void Fitter::ParseArgs(int argc, char **argv){

  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-b") {fDATAFiles[0] = argv[++i]; exist_bkgfile=true;} //data - background
    if(std::string(argv[i]) == "-d") {fDATAFiles[1] = argv[++i]; exist_datafile=true;} //data - cerenkov
  }
  
  //  if(!exist_srfile || !exist_yfile || !exist_bkgfile || !exist_datafile){
  if(!exist_bkgfile){
    std::cerr<<" Specify background with option: '-b'"<<std::endl;
    std::cerr<<" Performing fit without background!"<<std::endl;
    exit(0);
  }
  if(!exist_datafile){
    std::cerr<<" Specify input files with options: '-d'"<<std::endl;
    std::cerr<<" Fit cannot be perform! Drawing prefit histograms..."<<std::endl;
    exit(0);
  }
}

//Extract PDFs from real DATA input files and populate the TH1F vectors
void Fitter::GetDataPDFs(){

  std::cout<<" GetDataPDFs "<<std::endl;

  TGraph* gpdf_dt;
  TH1F* hdata;
  TH1F* hscale;
  for(int ifile=0; ifile<fDATAFiles.size(); ifile++){
    hdata = new TH1F(Form("hdata%i",ifile),"hdata",NBINS,0,100);
    hscale = new TH1F(Form("hscale%i",ifile),"hscale",NBINS,0,100); //count how many times we fill the bin and scale it, root is anoying...
    gpdf_dt = new TGraph(fDATAFiles[ifile],"%lg %lg",",");
    double x=0.;
    double y=0.;
    for(int ip=0; ip<gpdf_dt->GetN(); ip++){
      gpdf_dt->GetPoint(ip,x,y);
      hdata->Fill(x,y);
      hscale->Fill(x,1.);
    }
    hdata->Divide(hscale);
    //    hdata->Scale(1./hdata->Integral()); //Nomalize for shape only analysis
    hpdf_dt.push_back(hdata);
  }
  //  hpdf_dt[0]->Scale(1.6/4.7);//scale background
  hpdf_dt[0]->Scale(0.);//NO BACKGROUND

  //Initialize post fit histograms
  for(int ih=0; ih<hpdf_dt.size(); ih++)
    hpdf_dt_fit.push_back((TH1F*)hpdf_dt[ih]->Clone());

}


//Get MC pdfs
void Fitter::GetMCPDFs(){

  std::cout<<" GetMCPDFs "<<std::endl;

  GetMCPDFsWithCollEff(1.0);

  //Fill prefit histos
  for(int ih=0; ih<hpdf_mc_fit.size(); ih++){
    hpdf_mc.push_back(hpdf_mc_fit[ih]);
    hpdf_mc[ih]->SetName(Form("hpdf_mc_%d",ih));
  }
  
  //Fill sum histos: noise+90Sr+90Y
  hpdf_mc_sum = (TH1F*)hpdf_mc[0]->Clone();
  hpdf_mc_sum->Reset();
  for(int ih=0; ih<hpdf_mc.size(); ih++)
    hpdf_mc_sum->Add(hpdf_mc[ih]);
  hpdf_mc_sum->Add(hpdf_dt[0]);

}


//For the MC we need to relaunch the simulation for a set value
//of the efficiency_correction in order to float it and fit it.
//Will fill hpdf_mc_fit vector
void Fitter::GetMCPDFsWithCollEff(double collection_eff){

  //Prune histo vector
  for(int ifile=0; ifile<hpdf_mc_fit.size(); ifile++){
    hpdf_mc_fit[ifile]->Clear();
  }
  hpdf_mc_fit.clear();

  //Set the efficiency_correction according to fit step: we use a geometry template
  //where we subtitute the efficiency_correction for the desired value
  ifstream fgeo_template("./data/onepmt_template.geo"); //set template
  ofstream fgeo_used("./data/onepmt_fitter.geo"); //set outputfile
  std::string line;
  std::string valname_template = "CE_VALUE";
  size_t len = valname_template.length();
  while (getline(fgeo_template, line)){
    if(line.find("efficiency_correction")){
      size_t pos = line.find(valname_template);
      if (pos != std::string::npos)
	line.replace(pos, len, std::to_string(collection_eff));
    }
    fgeo_used<<line<<"\n"; //copy line by line
  }
  fgeo_used.close();
  //Run simulation
  char command_90Sr[10000];
  sprintf(command_90Sr,"rat -l ./logs/rat_TheiaRnD_90Sr_%f.log ./mac/TheiaRnD_90Sr.mac >> %s 2>&1", collection_eff, f_log.c_str());
  char command_90Y[10000];
  sprintf(command_90Y,"rat -l ./logs/rat_TheiaRnD_90Y_%f.log ./mac/TheiaRnD_90Y.mac >> %s 2>&1", collection_eff, f_log.c_str());
  if(ENABLESR) system(command_90Sr); //Run 90Sr
  system(command_90Y); //Run 90Y

  //Read file and extract PDFs
  fMCFiles[0] = "./results/TheiaRnD_90Sr_fitter.root"; //signal - strontium
  fMCFiles[1] = "./results/TheiaRnD_90Y_fitter.root"; //signal - ytrium

  RAT::DSReader *dsreader;
  TH1F *htemp;
  for(int ifile=0; ifile<fMCFiles.size(); ifile++){
    htemp = new TH1F(Form("hpdf_mc_fit_%i",ifile),"Charge",NBINS,0,100);
    if(!ENABLESR && ifile == 0){ //ignore Sr
      hpdf_mc_fit.push_back(htemp);
      continue;
    }
    dsreader = new RAT::DSReader(fMCFiles[ifile]);
    dsreader->GetT()->Draw(Form("ds.ev.pmt.charge>>hpdf_mc_fit_%i",ifile), "ds.ev.pmt.id==0");
    hpdf_mc_fit.push_back(htemp);
  }

  //Scale MC accoring to the DATA considering 1to1 rate for 90Sr and 90Y
  double norm = (hpdf_dt[1]->Integral(INITBIN,LASTBIN)-hpdf_dt[0]->Integral(INITBIN,LASTBIN))/(hpdf_mc_fit[0]->Integral(INITBIN,LASTBIN)+hpdf_mc_fit[1]->Integral(INITBIN,LASTBIN));
  hpdf_mc_fit[0]->Scale(norm);
  hpdf_mc_fit[1]->Scale(norm);
  
}

//Likelihood function
double Fitter::Likelihood(const double *par){

  double p_sr = par[0];
  double p_y = par[1];
  double p_bkg = par[2];
  double p_coleff = par[3];
  
  double likelihood=0.;
  double n_mc=0.;
  double n_dt=0.;

  //Fill mc PDFs with the appropriate collection efficiency
  GetMCPDFsWithCollEff(p_coleff);

  //Calculate likelihood
  for(int ibin=INITBIN; ibin<LASTBIN; ibin++){
    n_mc = p_sr*hpdf_mc_fit[0]->GetBinContent(ibin) +  p_y*hpdf_mc_fit[1]->GetBinContent(ibin) + p_bkg*hpdf_dt[0]->GetBinContent(ibin);
    n_dt = hpdf_dt[1]->GetBinContent(ibin);
    //    likelihood += n_mc - n_dt + n_dt*( log( TMath::Max(n_dt,0.000001) ) - log( TMath::Max(n_mc,0.000001) ) );
    if(n_mc>0 && n_dt>0)
      likelihood += n_mc - n_dt + n_dt*TMath::Log(n_dt/n_mc);
  }
  
  std::cout<<" Likelihood-> step: "<<nlikestep<<" "<<likelihood<<" "<<p_sr<<" "<<p_y<<" "<<p_bkg<<" "<<p_coleff<<std::endl;
  nlikestep++;
  return 2*likelihood;

}

//ChiSquare function
double Fitter::ChiSquare(const double *par){

  double p_sr = par[0];
  double p_y = par[1];
  double p_bkg = par[2];
  
  double chisquare=0.;
  double n_mc=0.;
  double n_dt=0.;

  //Calculate chisquear
  for(int ibin=INITBIN; ibin<LASTBIN; ibin++){
    n_mc = p_sr*hpdf_mc[0]->GetBinContent(ibin) +  p_y*hpdf_mc[1]->GetBinContent(ibin) + p_bkg*hpdf_dt[0]->GetBinContent(ibin);
    n_dt = hpdf_dt[1]->GetBinContent(ibin);
    if(n_mc>0 && n_dt>0)
      chisquare += pow((n_dt - n_mc),2)/n_mc;
  }
  
  std::cout<<" ChiSquare-> step: "<<chisquare<<" "<<p_sr<<" "<<p_y<<" "<<p_bkg<<std::endl;
  return chisquare;

}


  
//Do fit
void Fitter::DoFit(){

  min = ROOT::Math::Factory::CreateMinimizer("Minuit2", FITTER);

  //Set Function
  ROOT::Math::Functor f(this, &Fitter::Likelihood,4);
  //  ROOT::Math::Functor f(this, &Fitter::ChiSquare,4);
  min->SetFunction(f);
  
  // set tolerance , etc...
  //  min->SetMaxIterations(5);  // for GSL
  //  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxFunctionCalls(1000); // for Minuit/Minuit2 
  min->SetErrorDef(0.5);
  min->SetTolerance(1.0); //Minimization will stop when the estimated vertical
                          //distance to the minimum (EDM) is less than 0.001*tolerance*SetErrDef
  min->SetPrintLevel(3);

  // Set the free variables to be minimized!
  min->SetVariable(0,"90Sr",1.,0.01);
  min->SetVariableLimits(0,0.,20.);
  min->SetVariable(1,"90Y",1.,0.01);
  min->SetVariableLimits(1,0.,10.);
  min->SetVariable(2,"Noise",1.,0.01);
  min->SetVariableLimits(2,0.,10.);
  min->SetVariable(3,"Collection Efficiency",.5,0.5);
  min->SetVariableLimits(3,0.,1.);
  min->FixVariable(0); //Fix 90Sr norm
  min->FixVariable(1); //Fix 90Y norm
  min->FixVariable(2); //Fix background
  //  min->FixVariable(3); //Fix collection eff

  //  do the minimization
  min->Minimize();
  //  GetMCPDFsWithCollEff(0.6);

  fMinLikelihood = min->MinValue(); //Minimum likelihood
  const double *xs = min->X();
  for(int ipar=0; ipar<3; ipar++)
    fParFitted.push_back(xs[ipar]);

  //***Prepare best fit histograms
  hpdf_mc_fit[0]->Scale(fParFitted[0]);
  hpdf_mc_fit[1]->Scale(fParFitted[1]);
  hpdf_dt_fit[0]->Scale(fParFitted[2]);

  hpdf_mc_sum_fit = (TH1F*)hpdf_mc_fit[0]->Clone();
  hpdf_mc_sum_fit->Add(hpdf_mc_fit[1]);
  hpdf_mc_sum_fit->Add(hpdf_dt_fit[0]);

  
}

  
void Fitter::DrawPlots(){

  TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
  leg->AddEntry(hpdf_dt[1],"Data","PL");
  leg->AddEntry(hpdf_mc_sum,"MC","L");
  leg->AddEntry(hpdf_dt[0],"DT-Noise","F");
  leg->AddEntry(hpdf_mc[0],"MC-90Sr","L");
  leg->AddEntry(hpdf_mc[1],"MC-90Y","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  
  TCanvas *c_fit = new TCanvas("c_fit","c_fit",600,600);
  c_fit->Divide(2,1);
  //*********Before fit
  c_fit->cd(1);
  c_fit->cd(1)->SetLogy();
  if(DEBUG) std::cout<<" Draw data before fit "<<std::endl;
  //Data
  hpdf_dt[1]->SetLineWidth(2);
  hpdf_dt[1]->SetLineColor(1);
  hpdf_dt[1]->GetXaxis()->SetTitle("PE");
  hpdf_dt[1]->GetYaxis()->SetRangeUser(1.,1.e6);
  hpdf_dt[1]->Draw("");
  //Background
  for(int ih=0; ih<hpdf_dt.size()-1; ih++){
    hpdf_dt[ih]->SetLineWidth(1);
    hpdf_dt[ih]->SetFillColor(kGreen);
    hpdf_dt[ih]->SetFillStyle(3001);
    hpdf_dt[ih]->SetLineColor(kGreen+1);
    //    hpdf_dt[ih]->SetLineStyle(ih+1);
    hpdf_dt[ih]->Draw("same");
  }
  //MC
  if(DEBUG) std::cout<<" Draw MC before fit "<<std::endl;
  for(int ih=0; ih<hpdf_mc.size(); ih++){
    hpdf_mc[ih]->SetLineWidth(1);
    //    hpdf_mc[ih]->SetLineStyle(2);
    hpdf_mc[ih]->SetLineColor(kOrange+ih);
    hpdf_mc[ih]->Draw("same");
  }
  hpdf_mc_sum->SetLineWidth(2);
  hpdf_mc_sum->SetLineColor(kRed+1);
  hpdf_mc_sum->Draw("same");

  //***********After fit
  c_fit->cd(2);
  c_fit->cd(2)->SetLogy();
  //Data
  if(DEBUG) std::cout<<" Draw data after fit "<<std::endl;
  hpdf_dt[1]->SetLineWidth(2);
  hpdf_dt[1]->SetLineColor(1);
  hpdf_dt[1]->GetXaxis()->SetTitle("PE");
  hpdf_dt[1]->GetYaxis()->SetRangeUser(1.,1.e6);
  hpdf_dt[1]->Draw("");
  for(int ih=0; ih<hpdf_dt_fit.size()-1; ih++){
    hpdf_dt_fit[ih]->SetLineWidth(1);
    hpdf_dt_fit[ih]->SetFillColor(kGreen);
    hpdf_dt_fit[ih]->SetFillStyle(3001);
    hpdf_dt_fit[ih]->SetLineColor(kGreen);
    //    hpdf_dt_fit[ih]->SetLineStyle(ih+1);
    hpdf_dt_fit[ih]->Draw("same");
  }
  //MC
  if(DEBUG) std::cout<<" Draw MC after fit "<<std::endl;
  for(int ih=0; ih<hpdf_mc_fit.size(); ih++){
    hpdf_mc_fit[ih]->SetLineWidth(1);
    //    hpdf_mc_fit[ih]->SetLineStyle(2);
    hpdf_mc_fit[ih]->SetLineColor(kOrange+ih);
    hpdf_mc_fit[ih]->Draw("same");
  }
  hpdf_mc_sum_fit->SetLineWidth(2);
  hpdf_mc_sum_fit->SetLineColor(kRed+1);
  hpdf_mc_sum_fit->Draw("same");
  leg->Draw("same");

  TFile *f_out = new TFile(Form("./plots/TheiaRnD_fitter_%d_%d_%d.root",time_local->tm_mday,time_local->tm_mon,1900+time_local->tm_year),"RECREATE");
  f_out->cd();
  c_fit->Write();
  f_out->Close();
  
}

}//end TheiaRnD namespace


int main(int argc, char **argv){

  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy("App", &appargc, appargv);
  gStyle->SetOptStat(0);
  
  TheiaRnD::Fitter myfitter(argc, argv);
  myfitter.DoFit();
  myfitter.DrawPlots();

  if(!BATCH) dummy.Run();
  return 0;

  
}
