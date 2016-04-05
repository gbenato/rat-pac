#include<iostream>
#include<fstream>

#include<TF1.h>
#include<TH1F.h>
#include<TH2F.h>
#include<THStack.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<TBrowser.h>
#include<TGraph.h>
#include<TColor.h>
#include<TMath.h>
#include<TROOT.h>
#include "Math/Minimizer.h"

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include <RAT/DB.hh>
#include <RAT/DS/RunStore.hh>

#define REFTUBE 14
#define MCPHOTONLOOP false
#define NLOGENTRIES 10
#define FIT_LIMIT_MIN 50. //-60.
#define FIT_LIMIT_MAX 800.

//Constants
double cspeed = 300/1.4; // (mm/ns)/rindex
TVector3* target_pos = new TVector3(-400.,-400.,-200.);
TVector3 centerpos(0,0,0);

//Globals
RAT::DSReader *dsreader;
TTree *tree;
TTree *runT;
char *gInputFile = NULL;
char *gOutFile = NULL;
char *gTargetMaterial = "LAB";
RAT::DS::PMTInfo *pmtInfo;
bool isCherenkovData = false, isSPEData = false;

//Methods
void ParseArgs(int argc, char **argv);
void GetPMTInfo(char*);
void GetDBTables();
void GetHistos();
void DrawHistos();
void PrintHistos(char*);
void NormalizeHistos();
void GetPMTCalibration(); //Extract SPE and time delays

//SPE fit function
Double_t fmultigaus(Double_t *x, Double_t *par) {
  Double_t gaus_noise = 1./(par[2]*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[1] )/par[2] , 2.) );
  Double_t gaus_spe = 1./(par[5]*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[4] )/par[5] , 2.) );
  Double_t gaus_2pe = 1./(par[5]*sqrt(2.)*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[4]*2. )/(par[5]*sqrt(2.)), 2.) );
  Double_t gaus_3pe = 1./(par[5]*sqrt(3.)*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[4]*3. )/(par[5]*sqrt(3.)), 2.) );
  //  gaus_noise = 0.;
  //  return par[0]*gaus_noise + par[3]*gaus_spe;
  //  return par[0]*gaus_noise + par[3]*gaus_spe + par[6]*gaus_2pe + par[7]*gaus_3pe;
  return par[3]*gaus_spe + par[6]*gaus_2pe + par[7]*gaus_3pe;
  //  return par[0]*gaus_noise + par[3]*gaus_spe;
}

//Log Scale
void BinLogX(TH1*h)
{

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = TMath::Log10(axis->GetXmin());
   Axis_t to = TMath::Log10(axis->GetXmax());
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   delete new_bins;
}

void BinLogY(TH1*h)
{

   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from = TMath::Log10(axis->GetXmin());
   Axis_t to = TMath::Log10(axis->GetXmax());
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   delete new_bins;
}

//Global variables
std::vector<Color_t> pmtidtocolor;
std::vector<int> pmtidtopos;
std::vector<double> time_delay;
vector<double> qScintCorr; //Scintillation correction
vector<double> qScintCorrErr; // Scintillation correction error
vector<double> spe; //SPE

//// Histograms
//Event level
ULong64_t lastbursttime = 0;
TH1F* h_event_time;
TH1F* h_event_deltat;
TH1F* h_charge_total; //Total charge in the event
TH2F* h_pmt_chargevspos; //PMT charge vs PMT position
TH2F* h_pmt_npevspos; //NPE vs PMT position
TH2F* h_pmt_timevspos; //PMT charge vs PMT position
TH2F* h_pmt_countvspos; //Hit counter for normalization
TH1F* h_chi2; //Chi2 cher vs scint
TH2F* h_mcpmt_npevspos; //MC NPE vs PMT position

//Hit level (PMT)
int t_nbins = 100;
double t_min = -2.0;
double t_max = 5.0;
std::vector<char*> pmtname;
std::vector<double> q_xmin;
std::vector<double> q_xmax;
std::vector<double> npes_xmax;
std::vector<TF1*> f_spe;
std::vector<double> f_spe_mean;
std::vector<double> f_spe_sigma;
std::vector<double> f_timecal_mean;
std::vector<double> f_timecal_sigma;
std::vector<double> f_timecal_sigma_mc;
std::vector<TH1F*> h_npes; //Measured charge
std::vector<TH1F*> h_npes_ring;
std::vector<TH1F*> h_charge; //Measured charge
std::vector<TH1F*> h_charge_ring;
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_ring;
std::vector<TH1F*> h_time_bottom; //Measured time
std::vector<TH1F*> h_time_bottom_ring;
std::vector<TH1F*> h_time_trigger; //Measured time
std::vector<TH2F*> h_charge_vs_trigq; //PMT charge vs trigger charge
TH2F* h_charge_muontrigs; //Muon trigger charges correlation
TH1F* h_time_muontrigs; //DeltaT muon tags

//MC Truth
std::vector<TH1F*> h_mcpmt_npe; //PEs by PMT ID
std::vector<TH1F*> h_mcpmt_charge; //MC charge
std::vector<TH1F*> h_mcpmt_time; //MC FE time
std::vector<TH1F*> h_mcpmt_fetime; //MC FE time
std::vector<TH1F*> h_mctime;
std::vector<TH1F*> h_mctime_res;

int main(int argc, char **argv){

  //Init********
  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy("App", &appargc, appargv);
  ParseArgs(argc, argv);
  //************

  GetPMTInfo(gInputFile);
  GetDBTables();
  GetHistos();
  NormalizeHistos();
//  GetPMTCalibration();

  DrawHistos();
  if(gOutFile){
    PrintHistos(gOutFile);
  }

  dummy.Run();
  return 1;

}

void GetDBTables(){

  std::cout<<" Get DB scint correction tables for "<<gTargetMaterial<<std::endl;

  RAT::DB* db = RAT::DB::Get();
  db->Load("/Users/snoplus/Work/Chess/rat-pac/data/TheiaRnD/SCINTCORR.ratdb");
  RAT::DBLinkPtr dbScintCorr = db->GetLink("SCINTCORR",gTargetMaterial);
  qScintCorr = dbScintCorr->GetDArray("corr");
  qScintCorrErr = dbScintCorr->GetDArray("corr_err");
  //Calculate total charge and normalize
  double qtotal_smallpmts = 0.;
  for(int pmtid=0; pmtid<qScintCorr.size(); pmtid++){
    if(pmtInfo->GetType(pmtid)==1) qtotal_smallpmts += qScintCorr[pmtid];
  }
  for(int pmtid=0; pmtid<qScintCorr.size(); pmtid++){
    if(pmtInfo->GetType(pmtid)==1) qScintCorr[pmtid] /= qtotal_smallpmts;
    if(pmtInfo->GetType(pmtid)==1) qScintCorrErr[pmtid] /= qtotal_smallpmts;
  }

  db->Load("../data/PMTGAUSCHARGE.ratdb");
  RAT::DBLinkPtr dbSPE = db->GetLink("PMTGAUSCHARGE");
  spe = dbSPE->GetDArray("gaus_mean");

  db->Load("../data/PMTGAUSTIME.ratdb");
  RAT::DBLinkPtr dbTime = db->GetLink("PMTGAUSTIME");
  time_delay = dbTime->GetDArray("cable_delay");
  f_timecal_sigma_mc = dbTime->GetDArray("jitter");

}

void GetPMTInfo(char* inputfile){

  std::cout<<" GetPMTInfo "<<inputfile<<std::endl;

  //Init pmt positions
  dsreader = new RAT::DSReader(inputfile);

  std::cout<<" GetPMTInfo "<<std::endl;

  tree = dsreader->GetT();
  runT = dsreader->GetRunT();
  RAT::DS::Run *run = 0;
  runT->SetBranchAddress("run",&run);
  runT->GetEntry(0);
  pmtInfo = run->GetPMTInfo();

  //Get center of the small PMT array and locate
  int pmtTypeCount = 0;

  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    if(pmtInfo->GetType(ipmt) == 1) {
      centerpos = centerpos + pmtpos;
      pmtTypeCount++;
    }
  }
  centerpos = centerpos*(1./pmtTypeCount);

  Color_t mycolors[] = {1, 1, 1, 1, 1, 1,  kBlue, kRed,  kRed, 1, 1, 1,  kOrange, kBlue, kRed, kRed, kBlue, kOrange, kOrange, kBlue, kRed, kRed, kBlue, kOrange, 1}; //By Position (WATER)
  //Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kRed, kOrange, kBlue, kBlue, kOrange, kRed, kRed, kOrange, kBlue, 1}; //By Position (LAB)
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue-2, kOrange-2, kRed-2, kRed-1, kOrange-1, kBlue-1, kBlue, kOrange, kRed, kRed+1, kOrange+1, kBlue+1, 1}; //By Position
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kBlue, kBlue, kBlue, kOrange, kOrange, kOrange, kOrange, kRed, kRed, kRed, kRed, 1}; //By Digitizer group
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kCyan, kBlack, kGray, kGreen, kTeal, kAzure, kViolet, kPink, kYellow, 1}; //Individual

  int mypmtpos[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 1, 0, 0, 1, 2, 2, 1, 0, 3}; //In space
  //int mypmtpos[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3}; //In digit group

  pmtidtocolor.insert(pmtidtocolor.begin(), mycolors, mycolors + pmtInfo->GetPMTCount() );
  pmtidtopos.insert(pmtidtopos.begin(), mypmtpos, mypmtpos + pmtInfo->GetPMTCount() );

  //Plot axis limits
  double myqxmins[] = {-100,-100,-100,-100,-100,-100, 1.,1., 1.,1.,1.,1., -100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100, -10};
  q_xmin.insert(q_xmin.begin(), myqxmins, myqxmins + pmtInfo->GetPMTCount() );
  double myqxmaxs[] = {1500,1500,1500,1500,1500,1500, 5e4,5e4, 1e6,1e6,1e6,1e6, 400,400,400,400,400,400,400,400,400,400,400,400, 40};
  q_xmax.insert(q_xmax.begin(), myqxmaxs, myqxmaxs + pmtInfo->GetPMTCount() );
  double mynpesxmaxs[] = {10,10,10,10,10,10, 100,100, 100,100,100,100, 20,20,20,20,20,20,20,20,20,20,20,20, 10};
  npes_xmax.insert(npes_xmax.begin(), mynpesxmaxs, mynpesxmaxs + pmtInfo->GetPMTCount() );
  char* mypmtname[] = {"Light PMT 0", "Light PMT 1", "Light PMT 2", "Light PMT 3" ,"Light PMT 4", "Light PMT 5", "Top Muon Tag", "Bottom Muon Tag", "North Floor Panel", "South Floor Panel", "North Side Panel", "East Side Panel", "Ring PMT 0", "Ring PMT 1", "Ring PMT 2", "Ring PMT 3", "Ring PMT 4", "Ring PMT 5", "Ring PMT 6", "Ring PMT 7", "Ring PMT 8", "Ring PMT 9", "Ring PMT 10", "Ring PMT 11" , "Trigger PMT" };
  pmtname.insert(pmtname.begin(), mypmtname, mypmtname + pmtInfo->GetPMTCount() );

}

void GetHistos(){

  std::cout<<" Get PDFs "<<std::endl;

  //Event level
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_countvspos = new TH2F("h_pmt_countvspos","h_pmt_countvspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_chargevspos = new TH2F("h_pmt_chargevspos","h_pmt_chargevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_npevspos = new TH2F("h_pmt_npevspos","h_pmt_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_timevspos = new TH2F("h_pmt_timevspos","h_pmt_timevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);

  //Ring PMTs radius level
  for(int ih=0; ih<4; ih++){
    h_npes_ring.push_back(new TH1F(Form("h_npes_ring_%i",ih),"Ring Tube NPEs",50,0,npes_xmax[13]*4));
    h_charge_ring.push_back(new TH1F(Form("h_charge_ring_%i",ih),"Ring Tube Charges",25,q_xmin[13],q_xmax[13]*4));
    h_time_ring.push_back(new TH1F(Form("h_time_ring_%i",ih),"Ring Tube Event Time",t_nbins,t_min,t_max));
    h_time_bottom_ring.push_back(new TH1F(Form("h_time_bottom_ring_%i",ih),"Ring Tube Bottom Time",t_nbins,t_min,t_max));
  }

  for(int ih=0; ih<pmtInfo->GetPMTCount(); ih++){
    //MCTRUTH
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",200,0,200));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_time.push_back(new TH1F(Form("h_mcpmt_time_%i",ih),"h_mcpmt_time",200,5,10));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",200,5,10));
    //DAQ
    h_npes.push_back(new TH1F(Form("h_npes_%i",ih),"h_npes",npes_xmax[ih]*3,0,npes_xmax[ih]));
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),pmtname[ih],100,q_xmin[ih],q_xmax[ih]));
    if(ih>=6 && ih<=11){
      BinLogX(h_charge.back());
    }
    h_charge_vs_trigq.push_back(new TH2F(Form("h_charge_vs_trigq_%i",ih),"h_charge_vs_trigq",200,0,100,200,0,100));
    h_time.push_back(new TH1F(Form("h_time_%i",ih),"h_time",t_nbins,t_min,t_max));
    h_time_bottom.push_back(new TH1F(Form("h_time_bottom_%i",ih),"h_time_bottom",t_nbins,t_min,t_max));
    h_time_trigger.push_back(new TH1F(Form("h_time_trigger%i",ih),"h_time_trigger",t_nbins,t_min,t_max));
    f_spe.push_back(new TF1(Form("f_spe_%i",ih),fmultigaus,FIT_LIMIT_MIN,FIT_LIMIT_MAX,8));
  }
  h_charge_muontrigs = new TH2F("h_charge_muontrigs","Top Tag vs Bottom Tag Charges",100,q_xmin[6],q_xmax[6],100,q_xmin[7],q_xmax[7]);
  BinLogX(h_charge_muontrigs);
  BinLogY(h_charge_muontrigs);
  h_time_muontrigs = new TH1F("h_time_muontrigs","Bottom Tag Time - Top Tag Time",100,-5.,5.);
  h_event_time = new TH1F("h_event_time","Event Time",300,150,250);
  h_event_deltat = new TH1F("h_event_deltat","Event DeltaT",300,0,1e11);
  h_charge_total = new TH1F("h_charge_total","Total charge",200,-20,50000);
  h_chi2 = new TH1F("h_chi2","h_chi2",50,0,200);

  RAT::DS::Root *rds;
  RAT::DS::MC *mc;
  int nentries = tree->GetEntries();
  std::cout<<" Number of entries: "<<nentries<<std::endl;
  for(int ientry=0; ientry<nentries;++ientry){

    if(nentries>NLOGENTRIES && ientry%(nentries/NLOGENTRIES) == 0) std::cout<<" Entry "<<ientry<<std::endl;

    tree->GetEntry(ientry);
    rds = dsreader->GetEvent(ientry);
    mc = rds->GetMC();

    //MC**************
    //MCPMT loop
    for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
      RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      //Make grid display only for the small PMTs
      int pmtid = mcpmt->GetID();
      //count PE
      h_mcpmt_npe[pmtid]->Fill(mcpmt->GetMCPhotonCount());
      h_mcpmt_npevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),mcpmt->GetMCPhotonCount()/(double)nentries);
      h_mcpmt_charge[pmtid]->Fill(mcpmt->GetCharge());

      if(MCPHOTONLOOP){
        for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++){
          h_mcpmt_charge[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetCharge());
          h_mcpmt_time[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetHitTime());
          h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
        }
      } else{
        h_mcpmt_time[pmtid]->Fill(mcpmt->GetTime());
        h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetFrontEndTime());
      }
    } //end MCPMT loop
    //****************

    //DAQ EVENTS******
    //Event loop
    int nevents = rds->GetEVCount();
    for(int ievt=0; ievt<nevents; ievt++){

      RAT::DS::EV *ev = rds->GetEV(ievt);
      double chi2 = 0.;
      double qtotal = 0.;
      double qtotal_smallpmts = 0.;
      double npes = 0.;
      double charge = 0.;
      double qshort = 0.;
      double charge_trig = 0.;

      // Get muon tags, source trigger and muon panels charge and time
      double trigger_time = -9999.;
      double ring_time = -9999.;
      double ring_timeres = -9999.;
      double ring_charge = 0.;
      double topmuon_charge = 0.;
      double topmuon_time = 0.;
      double bottommuon_charge = 0.;
      double bottommuon_time = -9999.;
      double bottommuon_timeres = -9999.;
      double panel_charge[] = {0., 0., 0., 0.};
      RAT::DS::PMT *pmt = ev->GetPMTWithID(REFTUBE); //ref ring tube
      if(pmt!=NULL) {
        ring_charge = pmt->GetCharge();
        ring_time = pmt->GetTime();
        double ring_dist = (pmtInfo->GetPosition(REFTUBE) - *target_pos).Mag();
        double ring_tof = ring_dist/cspeed;
        ring_timeres = ring_time - ring_tof - time_delay[REFTUBE];
        //std::cout<<" ring_timeres "<<ring_timeres<<" "<<ring_time<<" "<<ring_tof<<" "<<time_delay[REFTUBE]<<std::endl;
      }
      pmt = ev->GetPMTWithID(24); //trigger
      if(pmt!=NULL) {
        charge_trig = pmt->GetCharge();
      }
      pmt = ev->GetPMTWithID(6); //top tag
      if(pmt!=NULL) {
        topmuon_charge = pmt->GetCharge();
        topmuon_time = pmt->GetTime();
      }
      pmt = ev->GetPMTWithID(7); //bottom tag
      if(pmt!=NULL) {
        bottommuon_charge = pmt->GetCharge();
        bottommuon_time = pmt->GetTime();
        double bottom_dist = (pmtInfo->GetPosition(7) - *target_pos).Mag();
        double bottom_tof = bottom_dist/cspeed;
        bottommuon_timeres = bottommuon_time - bottom_tof;
      }
      pmt = ev->GetPMTWithID(8); //pannels
      if(pmt!=NULL) panel_charge[0] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(9);
      if(pmt!=NULL) panel_charge[1] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(10);
      if(pmt!=NULL) panel_charge[2] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(11);
      if(pmt!=NULL) panel_charge[3] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(20);
      if(pmt!=NULL) trigger_time = pmt->GetTime();

      double event_time = -99999;
      std::vector<double> ringPMTTimes;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        charge = ev->GetPMT(ipmt)->GetCharge();
        qtotal += charge;
        if(pmttype==1) {
          double dist = (pmtInfo->GetPosition(pmtid) - *target_pos).Mag();
          double tof = dist/cspeed;
          if(ev->GetPMT(ipmt)->GetTime()>-9000) ringPMTTimes.push_back(ev->GetPMT(ipmt)->GetTime() - tof - time_delay[pmtid]);
        }
      }
      //Sort in ascending order
      std::sort(ringPMTTimes.begin(), ringPMTTimes.end());

      //Cuts for cherenkov imaging
      if(panel_charge[1]>100. || panel_charge[2]>100. || panel_charge[3]>100.) { //Veto cut
        lastbursttime = ev->GetClockTime();
      }
      if(isCherenkovData){
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) { //Veto cut
          continue;
        }
        if(topmuon_charge<200.0) continue;
        //if(bottommuon_charge<500.0) continue;
        // if(bottommuon_time - topmuon_time < 0. || bottommuon_time - topmuon_time > 1. ) continue; //Time cut
        // if(panel_charge[0]<500.) continue; //Veto cut
        // if(panel_charge[0]>8000.) continue; //Veto cut

        //Calculate event time
        // if(ringPMTTimes.size() > 1){
        //   event_time = TMath::KOrdStat((int)ringPMTTimes.size(), &ringPMTTimes[0], 1);
        // }
        if(ringPMTTimes.size() > 2){
          event_time = (ringPMTTimes[0] + ringPMTTimes[1] + ringPMTTimes[2])/3.;
        }
      }
      else if(isSPEData){ //Cuts for SPE
        event_time = ring_timeres;
        if(panel_charge[0]>50 || panel_charge[1]>50 || panel_charge[2]>50 || panel_charge[3]>50) continue;
        if(ring_time < 150) continue;
      }

      h_event_time->Fill(event_time);
      double deltat = (ev->GetClockTime() - lastbursttime) *2;
      // if(deltat < 10e9) continue;
      h_event_deltat->Fill(deltat); //2ns resolution

      //Fill Histograms
      std::vector<double> charge_ring(4,0.);
      std::vector<double> npes_ring(4,0.);
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        double dist = (pmtInfo->GetPosition(pmtid) - *target_pos).Mag();
        double tof = dist/cspeed;
        double pmttime = ev->GetPMT(ipmt)->GetTime();
        double pmtfcn = ev->GetPMT(ipmt)->GetFCN();
        charge = ev->GetPMT(ipmt)->GetCharge();
        npes = charge/spe[pmtid];
        qshort = ev->GetPMT(ipmt)->GetQShort();
        if(isCherenkovData){
          if(charge < 50) continue;
        }else if(isSPEData){
          //if(pmtfcn > 1000) continue;
          if(pmttime < 150) continue;
          //if(charge < 30) continue;
        }
        charge_ring[pmtidtopos[pmtid]] += charge;
        npes_ring[pmtidtopos[pmtid]] += npes;
        double timeres = pmttime - tof - time_delay[pmtid];
        h_time[pmtid]->Fill(timeres - event_time);
        // std::cout<<" ToF "<<pmtid<<": "<<tof<<" "<<dist<<std::endl;
        //if(pmtid != 13 && pmtid != 22 && pmtid != 19 && pmtid != 16) continue;
        h_time_bottom[pmtid]->Fill(timeres - bottommuon_time + 20);
        h_time_trigger[pmtid]->Fill(timeres - trigger_time);
        //Event level averaged
        if(timeres>-900) {
          if (pmtInfo->GetType(pmtid)==1){
            h_pmt_chargevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),charge);
            h_pmt_npevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),charge/spe[pmtid]);
            h_pmt_timevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),timeres - event_time);
            h_pmt_countvspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y());
          }
        }
        if(pmtInfo->GetType(pmtid)==1){
          qtotal_smallpmts += charge;
        }
        h_npes[pmtid]->Fill(npes);
        h_charge[pmtid]->Fill(charge);
//        h_charge[pmtid]->Fill(qshort);
        h_charge_vs_trigq[pmtid]->Fill(charge,charge_trig);
        //Compute chi2 for cher/scint
        if(pmtInfo->GetType(pmtid)==1){
          chi2 += pow( (charge - qScintCorr[pmtid])/qScintCorrErr[pmtid], 2.);
          //            std::cout<<" chi2 "<<ipmt<<" "<<pmtid<<" "<<chi2<<" "<<charge<<" "<<qScintCorr[pmtid]<<" "<<qScintCorrErr[pmtid]<<std::endl;
        }
      }
      h_charge_muontrigs->Fill(bottommuon_charge,topmuon_charge);
      h_time_muontrigs->Fill(bottommuon_time - topmuon_time);
      h_charge_total->Fill(qtotal_smallpmts);
      h_chi2->Fill(chi2);
      for(int irad=0; irad<4; irad++){
        if(charge_ring[irad] !=0 ) {
          h_charge_ring[irad]->Fill(charge_ring[irad]);
          h_npes_ring[irad]->Fill(npes_ring[irad]);
        }
      }
    } //end daq event loop

  } //end ds entry loop

  //Create Stacks for Ring Tubes
  for(int ipmt=0; ipmt<pmtInfo->GetPMTCount();ipmt++){
    h_time_ring[pmtidtopos[ipmt]]->Add(h_time[ipmt]);
    h_time_bottom_ring[pmtidtopos[ipmt]]->Add(h_time_bottom[ipmt]);
    h_time[ipmt]->SetLineColor(pmtidtocolor[ipmt]);
    //    h_time[ipmt]->SetFillColor(pmtidtocolor[ipmt]);
    //    h_time[ipmt]->SetFillStyle(3003+pmtidtopos[ipmt]);
    h_time_bottom[ipmt]->SetLineColor(pmtidtocolor[ipmt]);
    //    h_time_bottom[ipmt]->SetFillColor(pmtidtocolor[ipmt]);
    //    h_time_bottom[ipmt]->SetFillStyle(3003+pmtidtopos[ipmt]);
  }
}


void GetPMTCalibration(){

  double *pspe;
  double *ptime;
  double *pspe_err;
  double *ptime_err;
  for(int ipmt=0; ipmt<pmtInfo->GetPMTCount(); ipmt++){

//    if(pmtInfo->GetType(ipmt)!=1 && pmtInfo->GetType(ipmt)!=2) continue;

    f_spe[ipmt]->SetParameters(1.5e7,0.,10.,1e6,50.,50.,1e4,1e3);
    f_spe[ipmt]->SetParLimits(0.,1e6,1e8); //Noise norm
    f_spe[ipmt]->SetParLimits(1.,-30.,30.); //Noise mean
    f_spe[ipmt]->SetParLimits(2.,0.,30.); //Noise sigma
    f_spe[ipmt]->SetParLimits(3.,1e5,1e7); //SPE norm
    f_spe[ipmt]->SetParLimits(4.,20.,400.); //SPE mean
    f_spe[ipmt]->SetParLimits(5.,0.,200.); //SPE sigma
    f_spe[ipmt]->SetParLimits(6.,1e3,1e7); //2PE norm
    f_spe[ipmt]->SetParLimits(7.,1e2,2e6); //3PE norm

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Simplex");
    h_charge[ipmt]->Fit(Form("f_spe_%i",ipmt),"Q","",FIT_LIMIT_MIN,FIT_LIMIT_MAX);
    h_charge[ipmt]->Fit(Form("f_spe_%i",ipmt),"Q M","",FIT_LIMIT_MIN,FIT_LIMIT_MAX);
    pspe = f_spe[ipmt]->GetParameters();
    f_spe_mean.push_back(pspe[4]);
    f_spe_sigma.push_back(sqrt(pspe[5]*pspe[5] - pspe[2]*pspe[2]));

    h_time[ipmt]->Fit("gaus","Q M","",-0.7,0.7);
    if(h_time[ipmt]->GetFunction("gaus") == NULL){
      ptime = new double[3];
      ptime[0] = 0.; ptime[1] = 0.; ptime[2] = 0.;
      ptime_err = new double[3];
      ptime_err[0] = 0.; ptime_err[1] = 0.; ptime_err[2] = 0.;
      f_timecal_mean.push_back(0.);
      f_timecal_sigma.push_back(0.);
    } else{
      ptime = h_time[ipmt]->GetFunction("gaus")->GetParameters();
      ptime_err = h_time[ipmt]->GetFunction("gaus")->GetParErrors();
      f_timecal_mean.push_back(ptime[1]);
      double jittercorr = 1./sqrt(2.)*sqrt(ptime[2]*ptime[2] - f_timecal_sigma_mc[ipmt]*f_timecal_sigma_mc[ipmt]);
      f_timecal_sigma.push_back(jittercorr);
    }

    std::cout<<" PMT SPE "<<ipmt<<":"<<std::endl;
    std::cout<<"  |-> Noise norm: "<<pspe[0]<<" mean: "<<pspe[1]<<" sigma: "<<pspe[2]<<std::endl;
    std::cout<<"  |-> SPE norm: "<<pspe[3]<<" mean: "<<pspe[4]<<" sigma: "<<pspe[5]<<std::endl;
    std::cout<<"  |-> 2PE norm: "<<pspe[6]<<std::endl;
    std::cout<<"  |-> 3PE norm: "<<pspe[7]<<std::endl;
    std::cout<<"  |-> Delay: "<<ptime[1]<<" +- "<<ptime_err[1]<<" jitter: "<<ptime[2]<<std::endl;

  }

  ofstream ratdb_charge("PMTGAUSCHARGE.ratdb"), ratdb_time("PMTGAUSTIME.ratdb");

  //Charge ratdb
  ratdb_charge<<"{\n";
  ratdb_charge<<"  name: \"PMTGAUSCHARGE\", \n  valid_begin: [0,0], \n  valid_end: [0,0],\n";
  ratdb_charge<<"  gaus_mean: ["<<f_spe_mean.at(0);
  for(int ipmt=1; ipmt<f_spe_mean.size(); ipmt++) ratdb_charge<<", "<<f_spe_mean.at(ipmt);
  ratdb_charge<<"],\n";
  ratdb_charge<<"  gaus_sigma: ["<<f_spe_sigma.at(0);
  for(int ipmt=1; ipmt<f_spe_sigma.size(); ipmt++) ratdb_charge<<", "<<f_spe_sigma.at(ipmt);
  ratdb_charge<<"],\n";
  ratdb_charge<<"}\n";

  //Time ratdb
  ratdb_time<<"{\n";
  ratdb_time<<"  name: \"PMTGAUSTIME\", \n  valid_begin: [0,0], \n  valid_end: [0,0],\n";
  ratdb_time<<"  cable_delay: ["<<f_timecal_mean.at(0);
  for(int ipmt=1; ipmt<f_timecal_mean.size(); ipmt++) ratdb_time<<", "<<f_timecal_mean.at(ipmt);
  ratdb_time<<"],\n";
  ratdb_time<<"  transit_time: [63.0,63.0,63.0,63.0,63.0,63.0,9.0,9.0,20.0,20.0,20.0,20.0,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8,5.8],\n";
  ratdb_time<<"  tts: [1.23,1.23,1.23,1.23,1.23,1.23,0.212,0.212,1.0,1.0,1.0,1.0,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127,0.127],\n";
  ratdb_time<<"  jitter: ["<<f_timecal_sigma.at(0);
  for(int ipmt=1; ipmt<f_timecal_sigma.size(); ipmt++) ratdb_time<<", "<<f_timecal_sigma.at(ipmt);
  ratdb_time<<"],\n";
  ratdb_time<<"}\n";

}




//Draw histrograms
void DrawHistos(){

  TCanvas *c_chi2 = new TCanvas("c_chi2","c_chi2",400,400);
  h_chi2->Draw();

  //General plots
  TCanvas *c_ring_event = new TCanvas("c_ring_event","Ring Event",1200,900);
  c_ring_event->Divide(3,2);
  c_ring_event->cd(1);
  //h_pmt_chargevspos->Draw("colz text"); //Average charge per pmt vs position
  h_pmt_npevspos->Draw("colz text"); //Average charge per pmt vs position
  c_ring_event->cd(2);
  h_charge_total->Draw(); //Total charge in event
  c_ring_event->cd(3);
  h_pmt_timevspos->Draw("colz text"); //Average charge per pmt vs position
  c_ring_event->cd(4);
  h_event_time->Draw();
  c_ring_event->cd(5);
  h_event_deltat->Draw();
  TCanvas *c_charge[5];
  c_charge[0] = new TCanvas("c_charge_0","Charge Ring Tubes",1200,900);
  c_charge[1] = new TCanvas("c_charge_1","Charge Light Tubes",900,600);
  c_charge[2] = new TCanvas("c_charge_2","Muon Tags",600,600);
  c_charge[3] = new TCanvas("c_charge_3","Charge Trigger",300,300);
  c_charge[4] = new TCanvas("c_charge_4","Charge Vetos",600,600);
  c_charge[0]->Divide(4,3);
  c_charge[1]->Divide(3,2);
  c_charge[2]->Divide(2,2);
  c_charge[3]->Divide(1,1);
  c_charge[4]->Divide(2,2);
  TCanvas *c_npes[2];
  c_npes[0] = new TCanvas("c_npes_0","NPEs Ring Tubes",1200,900);
  c_npes[1] = new TCanvas("c_npes_1","NPEs Light Tubes",900,600);
  c_npes[0]->Divide(4,3);
  c_npes[1]->Divide(3,2);
  TCanvas *c_time[5];
  c_time[0] = new TCanvas("c_time_0","Time Ring Tubes",1200,900);
  c_time[1] = new TCanvas("c_time_1","Time Light Tubes",900,600);
  c_time[3] = new TCanvas("c_time_3","Time Trigger",300,300);
  c_time[4] = new TCanvas("c_time_4","Time Vetos",600,600);
  c_time[0]->Divide(4,3);
  c_time[1]->Divide(3,2);
  c_time[3]->Divide(1,1);
  c_time[4]->Divide(2,2);

  int cc0 = 0, cc1 = 0, cc2 = 0, cc3 = 0, cc4 = 0;
  for(int pmtid = 0; pmtid<pmtInfo->GetPMTCount(); pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype==1){//Ring tubes
      c_charge[0]->cd(++cc0)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
      h_charge[pmtid]->Draw();
      f_spe[pmtid]->Draw("L same");
      c_npes[0]->cd(cc0)->SetLogy();
      h_npes[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_npes[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
      h_npes[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
      h_npes[pmtid]->Draw();
      c_time[0]->cd(cc0);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==2){ //Light tubes
      c_charge[1]->cd(++cc1)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      f_spe[pmtid]->Draw("L same");
      c_npes[1]->cd(cc1)->SetLogy();
      h_npes[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_npes[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
      h_npes[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
      h_npes[pmtid]->Draw();
      c_time[1]->cd(cc1);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==3){ //Muon tags
      c_charge[2]->cd(++cc2)->SetLogx();
      c_charge[2]->cd(cc2)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
    } else if(pmttype==0){ //Trigger tube
      c_charge[3]->cd(++cc3)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      c_time[3]->cd(cc3);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==4){ //Muon panels
      c_charge[4]->cd(++cc4)->SetLogy();
      c_charge[4]->cd(cc4)->SetLogx();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      c_time[4]->cd(cc4);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    }
  }

  c_charge[2]->cd(3)->SetLogx();
  c_charge[2]->cd(3)->SetLogy();
  h_charge_muontrigs->Draw("colz");
  c_charge[2]->cd(4);
  h_time_muontrigs->Draw();


  TCanvas *c_time_ring = new TCanvas("c_time_ring","Ring Time respect to Event Time",900,900);
  h_time_ring[0]->SetLineColor(pmtidtocolor[12]);
  h_time_ring[1]->SetLineColor(pmtidtocolor[13]);
  h_time_ring[2]->SetLineColor(pmtidtocolor[14]);
  h_time_ring[0]->Draw();
  h_time_ring[1]->Draw("sames");
  h_time_ring[2]->Draw("sames");

  TCanvas *c_charge_ring = new TCanvas("c_charge_ring","Ring Tube Charges",900,900);
  h_charge_ring[0]->SetLineColor(pmtidtocolor[12]);
  h_charge_ring[1]->SetLineColor(pmtidtocolor[13]);
  h_charge_ring[2]->SetLineColor(pmtidtocolor[14]);
  h_charge_ring[0]->Draw();
  h_charge_ring[1]->Draw("sames");
  h_charge_ring[2]->Draw("sames");

  TCanvas *c_npes_ring = new TCanvas("c_npes_ring","Ring Tube NPEs",900,900);
  h_npes_ring[0]->SetLineColor(pmtidtocolor[12]);
  h_npes_ring[1]->SetLineColor(pmtidtocolor[13]);
  h_npes_ring[2]->SetLineColor(pmtidtocolor[14]);
  h_npes_ring[0]->Draw();
  h_npes_ring[1]->Draw("sames");
  h_npes_ring[2]->Draw("sames");

  TCanvas *c_time_bottom = new TCanvas("c_time_bottom","Ring Time respect to Bottom",900,900);
  h_time_bottom_ring[0]->SetLineColor(pmtidtocolor[12]);
  h_time_bottom_ring[1]->SetLineColor(pmtidtocolor[13]);
  h_time_bottom_ring[2]->SetLineColor(pmtidtocolor[14]);
  h_time_bottom_ring[0]->Draw();
  h_time_bottom_ring[1]->Draw("same");
  h_time_bottom_ring[2]->Draw("same");

  TCanvas *c_time_trigger = new TCanvas("c_time_trigger","Ring Time respect to Trigger",900,900);
  for(int pmtid = 0; pmtid<pmtInfo->GetPMTCount(); pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype!=1) continue;
    h_time_trigger[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    h_time_trigger[pmtid]->Draw("sames");
  }

  //MCTruth
  bool firstdrawn0=false, firstdrawn1=false, firstdrawn2=false;
  TCanvas *c_mc = new TCanvas("c_mc","c_mc",900,1000);
  c_mc->Divide(3,3);
  for(int pmtid = 0; pmtid<pmtInfo->GetPMTCount(); pmtid++){
    for(int pmtid = 0; pmtid<pmtInfo->GetPMTCount(); pmtid++){
      int pmttype = pmtInfo->GetType(pmtid);
      if(pmttype==1){//Fast tubes
        const char *opt = firstdrawn0 ? "sames" : "";
        c_mc->cd(1);
        h_mcpmt_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_charge[pmtid]->Draw(opt);
        c_mc->cd(2);
        h_mcpmt_npe[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_npe[pmtid]->Draw(opt);
        c_mc->cd(3);
        h_mcpmt_fetime[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_fetime[pmtid]->Draw(opt);
        firstdrawn0 = true;
      } else if(pmttype==2){ //Large tubes
        const char *opt = firstdrawn1 ? "sames" : "";
        c_mc->cd(4);
        h_mcpmt_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_charge[pmtid]->Draw(opt);
        c_mc->cd(5);
        h_mcpmt_npe[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_npe[pmtid]->Draw(opt);
        c_mc->cd(6);
        h_mcpmt_fetime[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_fetime[pmtid]->Draw(opt);
        firstdrawn1 = true;
      } else if(pmttype==0){ //Trigger tube
        const char *opt = firstdrawn2 ? "sames" : "";
        // c_mc->cd(7);
        // h_mcpmt_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        // h_mcpmt_charge[pmtid]->Draw(opt);
        // c_mc->cd(8);
        // h_mcpmt_npe[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        // h_mcpmt_npe[pmtid]->Draw(opt);
        c_mc->cd(9);
        h_mcpmt_fetime[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_fetime[pmtid]->Draw(opt);
        firstdrawn2 = true;
      }
    }
  }

  TCanvas *c_pmtmaps = new TCanvas("c_pmtmaps","c_pmtmaps",300,300);
  // c_pmtmaps->Divide(2,1);
  // c_pmtmaps->cd(1);
  h_mcpmt_npevspos->Draw("colz text");

}



//Output histrograms to file
void PrintHistos(char *filename){

  std::cout<<" Output histrograms to "<<filename<<std::endl;

  TFile *fout = new TFile(filename,"RECREATE");
  fout->cd();
  h_time_muontrigs->Write(); //For event counting and norm
  for(int ipmt=0; ipmt<pmtInfo->GetPMTCount();ipmt++){
    h_charge[ipmt]->Write();
    h_npes[ipmt]->Write();
    h_time[ipmt]->Write();
  }
  for(int ih=0; ih<4;ih++){
    h_charge_ring[ih]->Write();
    h_npes_ring[ih]->Write();
    h_time_ring[ih]->Write();
  }
  fout->Close();

}

void NormalizeHistos(){

  double norm = h_event_time->GetEntries();
  // h_pmt_chargevspos->Scale(1./norm);
  // h_pmt_npevspos->Scale(1./norm);
  //h_pmt_timevspos->Scale(1./norm);

  h_pmt_chargevspos->Divide(h_pmt_countvspos);
  h_pmt_npevspos->Divide(h_pmt_countvspos);
  h_pmt_timevspos->Divide(h_pmt_countvspos);

  //Fill with zeroes
  for (int ibin = 1; ibin < h_pmt_timevspos->GetXaxis()->GetNbins()+1; ibin++) {
    for (int jbin = 1; jbin < h_pmt_timevspos->GetYaxis()->GetNbins()+1; jbin++) {
      h_pmt_chargevspos->SetBinContent(ibin, jbin, -1000. + h_pmt_chargevspos->GetBinContent(ibin, jbin));
      h_pmt_timevspos->SetBinContent(ibin, jbin, -1000. + h_pmt_timevspos->GetBinContent(ibin, jbin));
    }
  }
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    if(pmtInfo->GetType(ipmt)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    h_pmt_chargevspos->Fill(pmtpos.X(), pmtpos.Y(), 1000.);
    h_pmt_timevspos->Fill(pmtpos.X(), pmtpos.Y(), 1000.);
  }

  h_pmt_chargevspos->SetMaximum(1000.);
  h_pmt_chargevspos->SetMinimum(-100.0);
  h_pmt_npevspos->SetMaximum(10.);
  h_pmt_npevspos->SetMinimum(0.0);
  h_pmt_timevspos->SetMaximum(1.5);
  h_pmt_timevspos->SetMinimum(-0.5);

}


void ParseArgs(int argc, char **argv){
  bool exist_inputfile = false;
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-i") {gInputFile = argv[++i];}
    if(std::string(argv[i]) == "-o")  {gOutFile = argv[++i];}
    if(std::string(argv[i]) == "-m")  {gTargetMaterial = argv[++i];}
    if(std::string(argv[i]) == "-cher") {isCherenkovData = true;}
    if(std::string(argv[i]) == "-spe") {isSPEData = true;}
  }

  if(argc<=1 || (!isCherenkovData && !isSPEData) || (isCherenkovData && isSPEData)){
    std::cerr<<" Usage: ./DrawChargeAndTime.exe (-cher || -spe) -i INPUT_FILE [-o OUTPUT_FILE -m TARGET_MATERIAL]"<<std::endl;
    exit(0);
  }
}
