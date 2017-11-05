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
#include<TRandom3.h>
#include "Math/Minimizer.h"

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include<RAT/DB.hh>
#include<RAT/DS/RunStore.hh>

#include "CHESSTools.h"

#define SMEARING 0.000 //0.214
#define CORRCH 24 //24 for source data // 0 - for cosmics as empty channel for correction (if negative do not apply correction)
#define REFCH 18// A single example channels ( to study/visualize the charge correction signal crosstalks)
#define REFTUBE 23 //Reference tube 19 - PMT8 (inner ring more statistics,2nd best noise resolution) for time measurements 
#define TIMENPELIMIT 1.2
//hmmm there is an discrepencay with CHESSTOOLS.h naming convention of the pmts
// there chan 23, 24, 25 are the top muon tag, bottom muon tag and trigger PMT
// but chan 0 is named control channel 
#define MCPHOTONLOOP false // false
#define NLOGENTRIES 10
#define FIT_LIMIT_MIN -60. //50.
#define FIT_LIMIT_MAX 800.
#define RING_CANDIDATE 5 //99
#define TOTAL_NPE_CUT 160 // -1 to deactivate, Cut events with secondaries, bundles...

//PMT IDs
#define PMTID_INNER 11
#define PMTID_MID 12
#define PMTID_OUTER 13


//Binning
int t_nbins = 40;
double t_min[2] = {-1.0, -1.0}; //ring, light //What unit?
double t_max[2] = {1.5, 1.5}; //ring, light

//double lab_bins[21] = {-2., -1.75, -1.5, -1.25, -1., -.75, -.5, -.25, .0, .25, .5, .75, 1., 1.5, 2., 3., 5., 8., 12., 16., 20.};
double lab_bins[26] = {-2., -1.75, -1.5, -1.25, -1., -.75, -.5, -.25, .0, .25, .5, .75, 1., 1.25, 1.5, 2., 3., 5., 6.5, 8., 10., 12., 14., 16., 18., 20.};

//double q_xmin[25] = {-10,-10,-10,-10,-10,-10,  0.1,0.1,0.1,0.1,  -100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,  1,1,  100}; //In PEs
//double q_xmax_water[] = {10,10,10,10,10,10,  1000,1000,1000,1000,  100,100,100,100,100,100,100,100,100,100,100,100,  10000,10000,  100}; //In PEs
//double q_xmax_water[] = {10,10,10,10,10,10,  1000,1000,1000,1000,  5,5,5,5,5,5,5,5,5,5,5,5,  10000,10000,  100}; //In PEs

double q_xmin[] = {0, -10,-10,-10,-10,-10,-10,  0.1,0.1,0.1,0.1,  -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,  1,1,  -5}; //In PEs
double q_xmax_water[] = {0, 10,10,10,10,10,10,  1000,1000,1000,1000,  10,10,10,10,10,10,10,10,10,10,10,10,  10000,10000,  100}; //In PEs
double q_xmax_lab[] = {0, 20,20,20,20,20,20,  1000,1000,1000,1000,  20,20,20,20,20,20,20,20,20,20,20,20,  10000,10000,  100}; //In PEs
double q_xmax_labppo[] = {0, 100,100,100,100,100,100,  1000,1000,1000,1000,  100,100,100,100,100,100,100,100,100,100,100,100,  10000,10000,  200}; // In PEs

//Constants
double cspeed = 300/1.4; // (mm/ns)/rindex
TVector3* target_pos = new TVector3(-398.0, -367.0, -220.91);
TVector3 centerpos(0,0,0);

//Globals
TRandom3 *trand3 = new TRandom3();
RAT::DSReader *dsreader;
TTree *tree;
TTree *runT;
string data_dir = "";
int total_npe_cut = -1;
char *gInputFile = NULL;
char *gOutFile = NULL;
//std::string gTargetMaterial = "WATER";
std::string gTargetMaterial = "LABPPO";
RAT::DS::PMTInfo *pmtInfo;
bool isMC = false;
bool isCosmicData = false, isSourceData = false;

double tof_fixed[3] = {0.625727, 0.535623, 0.473399};

//Methods
void ParseArgs(int argc, char **argv);
void GetPMTInfo(const char*);
void GetDBTables();
void InitHistos();
void GetHistos();
void DrawHistos();
void PrintHistos(char*);
void NormalizeHistos();

//Log Scale
void BinLogX(TH1*h);
void BinLogY(TH1*h);

//Global variables
std::vector<double> time_delay;
vector<double> qScintCorr; //Scintillation correction
vector<double> qScintCorrErr; // Scintillation correction error
vector<double> spe; //SPE
vector<double> charge_slopes; // channel based charge correction
vector<double> charge_offsets; // channel based charge correction
vector<double> spe_corrections; // channel based spe scale correction
vector<double> noise_offsets; // channel based noise position correction
vector<string> data_files; // data files as read from MEASUREMENTSFILES
vector<double> noise_cutoffs; // channel based 3-sigma threshold for (signal selection)

//// Histograms
//Event level
ULong64_t lastbursttime = 0;
TH1F* h_nhits_ring;
TH1F* h_nhits_light;
TH1F* h_event_time;
TH1F* h_event_deltat;
TH1F* h_qtotal_light; //Total charge in the event
TH1F* h_qtotal_ring; //Total charge in the event
TH1F* h_charge_r_inout; // Charge ratio from the 3 well calibrated inner PMTs over the 3 well calibrated outer PMTs 
TH2F* h_pmt_qratiovspos; //PMT charge vs PMT position
TH2F* h_pmt_chargevspos; //NPE vs PMT position
TH2F* h_pmt_timevspos; //PMT charge vs PMT position
TH2F* h_pmt_countvspos; //Hit counter for normalization
TH2F* h_mcpmt_npevspos; //MC NPE vs PMT position
TH1F* h_mc_npe; // MC NPE in all cross PMTs
TH2F* h_ringcandidate_npevspos;
TH2F* h_ringcandidate_timevspos;
TH2F* h_timevsnpe;
TH2F* h_timevsqratio;
TGraph* g_corr_channel; // corr_channel over time


std::vector<TH1F*> h_qratio_ring;
std::vector<TH1F*> h_charge; //Measured charged
std::vector<TH1F*> h_qratio; //Measured charge
std::vector<TH1F*> h_charge_ring;
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_ring;
std::vector<TGraph*> g_charge_correction; // plot charge in chx against charge in Corrch (cutting on NPE)
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
  //

  //Style
//  gROOT->SetMacroPath("/home/beschmidt/rat-pac/scripts/cuore_tools/");
  SetCHESSStyle();

  if(SMEARING != 0) {
    std::cout<<" SMEARING IS NOT ZERO!!! "<<std::endl;
    std::cout<<" Are you analyzing real data? y/n ...."<<std::endl;
    char input[256]; std::cin>>input;
    if (input[0] == 'y') {
      std::cout<<" You can't do that! Exiting..."<<std::endl;
      exit(0);
    }
    else if (input[0] == 'n') std::cout<<" OK! Sorry about that. "<<std::endl;
    else exit(0);
  }
  GetDBTables();
  GetPMTInfo(data_files[0].c_str());
  InitHistos();
  GetHistos();
  NormalizeHistos();

  DrawHistos();
  if(gOutFile){
    PrintHistos(gOutFile);
  }

  dummy.Run();
  return 1;

}

double 	GetEventTime(RAT::DS::EV *ev, double corr_charge, double min_npe = 1.0){
  //
  //Calculate the median of all cross PMTs with a significant NPE number.
  //
  vector<double> ev_time;
  for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
      int pmtid = ev->GetPMT(ipmt)->GetID();
      int pmttype = pmtInfo->GetType(pmtid);
      if(pmttype != 1) continue;
        double charge = ev->GetPMT(ipmt)->GetCharge();
        charge -= ((corr_charge-charge_offsets[pmtid])*charge_slopes[pmtid]);
        double npes = charge/spe[pmtid]-noise_offsets[pmtid];
        npes = npes / (spe_corrections[pmtid] - noise_offsets[pmtid]);
        if (npes > min_npe){
          double refring_time = ev->GetPMT(ipmt)->GetTime();
//          double refring_dist = (pmtInfo->GetPosition(pmtid) - *target_pos).Mag(); // Make sure whether this is pmtid or ipmt
//         double refring_tof = refring_dist/cspeed; // this can't really work if my reftube is the movable trigger cube
          double refring_tof = tof_fixed[pmtidtopos[pmtid]];
          double refring_timeres = refring_time - refring_tof - time_delay[pmtid];
	  ev_time.push_back(refring_timeres);
       }
  }
  if (ev_time.size()){
//    return TMath::Median(ev_time.begin(), ev_time.end());
    return TMath::Median(ev_time.size(), &ev_time[0]);
  }
  return -9999.9;
}

double GetTriggerTime(RAT::DS::EV *ev){
  //
  //Calculate the median of all cross PMTs with a significant NPE number.
  //
  vector<double> ev_time;
  std::vector<double> trace = ev->GetPMT(REFTUBE)->GetWaveform();
  int trace_length = trace.size();
  int n_avg = 50;
  double* begin = &trace[0];
  double* end = &trace[trace_length-n_avg-1];  
  double begin_val = TMath::Mean(n_avg, begin);
  double end_val = TMath::Mean(n_avg, end);
  double fifty_val = begin_val + (end_val-begin_val)/2.0;
  double temp_trace[trace_length];
  for (int i = 0 ; i < trace_length; i++){
     temp_trace[i] = abs(trace[i] - fifty_val)
  }
  long min_index = TMath::LocMin(temp_trace);
  if (min_index-1)>0 && min_index+1 < trace_length(){
    double dy = (trace[min_index+1]-trace[min_index-1])/2.0;
    double dt = 
  }
 
  }
  return -9999.9;
}

void GetDBTables(){

  RAT::DB* db = RAT::DB::Get();

//  db->Load("../data/PMTGAUSCHARGE.ratdb");
  db->Load("PMTGAUSCHARGE.ratdb"); // Do not change name
  RAT::DBLinkPtr dbSPE = db->GetLink("PMTGAUSCHARGE");
  spe = dbSPE->GetDArray("gaus_mean");
  charge_slopes = dbSPE->GetDArray("charge_correction_slopes");
  charge_offsets = dbSPE->GetDArray("charge_correction_offsets");
  spe_corrections = dbSPE->GetDArray("gaus_correction");
  noise_offsets = dbSPE->GetDArray("noise_offset");
  noise_cutoffs = dbSPE->GetDArray("noise_3_sigma");
  std::cout << "Red charge correction " << charge_slopes[25] << std::endl;
//  db->Load("../data/PMTGAUSTIME.ratdb");
  db->Load("PMTGAUSTIME.ratdb");
  RAT::DBLinkPtr dbTime = db->GetLink("PMTGAUSTIME");
  time_delay = dbTime->GetDArray("cable_delay");
  cout << "Reading MEASUREMENTFILES.ratdb" << endl;
  db->Load("MEASUREMENTFILES.ratdb");
  RAT::DBLinkPtr dbFiles = db->GetLink("MEASUREMENTFILES") ;
  data_files = dbFiles->GetSArray("data_files");
  data_dir = dbFiles->GetS("data_dir");
  total_npe_cut = dbFiles->GetI("total_npe_cut");

  cout << "Input files loaded: " <<  endl;
  for (int i = 0; i < data_files.size(); i++){
    data_files[i] = data_dir+data_files[i];
    cout << data_files[i] << endl;
  }
}

void GetPMTInfo(const char* inputfile){

  std::cout<<" GetPMTInfo "<<inputfile<<std::endl;

  //Init pmt positions
  dsreader = new RAT::DSReader(inputfile);
  for(int i = 1; i < data_files.size(); i++){ 
    dsreader->Add(data_files[i].c_str());
  }

  std::cout<<" GetPMTInfo "<<std::endl;

  tree = dsreader->GetT();
  std::cout<< tree << std::endl;
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

}

void InitHistos(){

  std::cout<<" Init Histograms "<<std::endl;

  //Set Q limit and Time limits
  int n_pmts = 25; 
  double q_xmax_used[n_pmts];
  double qtotal_xmax;
  double qrat_xmax = 0.;
  if(isCosmicData){
    if(gTargetMaterial=="WATER"){
      qtotal_xmax = 50.;
      qrat_xmax = 1.0;
      memcpy(q_xmax_used, q_xmax_water, (n_pmts+1)*sizeof(double));
      t_min[0] = -1.5;
      t_max[0] = 1.5;
    } else if(gTargetMaterial=="LAB"){
      qtotal_xmax = 100.;
      qrat_xmax = 0.5;
      memcpy(q_xmax_used, q_xmax_lab, (n_pmts+1)*sizeof(double));
      t_nbins = 500;
      t_min[0] = -2.0;
      t_max[0] = 10.0;
    }else if(gTargetMaterial=="LABPPO"){
      qtotal_xmax = 1500.;
      qrat_xmax = 0.1;
      memcpy(q_xmax_used, q_xmax_labppo, (n_pmts+1)*sizeof(double));
      t_min[0] = -2.0;
      t_max[0] = 5.0;
    }else if(gTargetMaterial=="WBLS1PCT"){
      qtotal_xmax = 100.;
      qrat_xmax = 0.5;
      memcpy(q_xmax_used, q_xmax_labppo, (n_pmts+1)*sizeof(double));
      t_nbins = 50;
      t_min[0] = -2.0;
      t_max[0] = 10.0;
    }else if(gTargetMaterial=="WBLS5PCT"){
      qtotal_xmax = 100.;
      qrat_xmax = 0.3;
      memcpy(q_xmax_used, q_xmax_labppo, (n_pmts+1)*sizeof(double));
      t_nbins = 50;
      t_min[0] = -2.0;
      t_max[0] = 10.0;
    }else if(gTargetMaterial=="WBLS10PCT"){
      qtotal_xmax = 200.;
      qrat_xmax = 0.3;
      memcpy(q_xmax_used, q_xmax_labppo, (n_pmts+1)*sizeof(double));
      t_nbins = 50;
      t_min[0] = -2.0;
      t_max[0] = 10.0;
      //npes_xmax[11] = 200.0;
    }
  } else if(isSourceData){
    qtotal_xmax = 50.;
    qrat_xmax = 1.0;
    memcpy(q_xmax_used, q_xmax_water, (n_pmts+1)*sizeof(double));
    t_nbins = 2000;
    t_min[0] = -100.;
    t_max[0] = 100.;
    t_min[1] = 0.;
    t_max[1] = 1000.;
  }


  //Event level
  g_corr_channel = new TGraph();
  g_corr_channel->SetName("g_corr_channel");
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_timevsnpe = new TH2F("h_timevsnpe", "h_timevsnpe", 400, -100, 100, 100, -5, 10);
  h_timevsqratio = new TH2F("h_timevsqratio", "h_timevsqratio", 400, -100, 100, 100, -1, 1.5);
  h_mc_npe = new TH1F("h_mc_npe","h_mc_npe",2500,0,500 );
  h_pmt_countvspos = new TH2F("h_pmt_countvspos","h_pmt_countvspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_qratiovspos = new TH2F("h_pmt_qratiovspos","h_pmt_qratiovspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_chargevspos = new TH2F("h_pmt_chargevspos","h_pmt_chargevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_timevspos = new TH2F("h_pmt_timevspos","h_pmt_timevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_ringcandidate_npevspos = new TH2F("h_ringcandidate_npevspos","h_ringcandidate_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_ringcandidate_timevspos = new TH2F("h_ringcandidate_timevspos","h_ringcandidate_timevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_event_time = new TH1F("h_event_time","Event Time",500,0,200);
  h_nhits_ring = new TH1F("h_nhits_ring","NHits Ring Tubes",40,0,20);
  h_nhits_light = new TH1F("h_nhits_light","NHits Light Tubes",40,0,20);
  h_event_deltat = new TH1F("h_event_deltat","Event DeltaT",300,1e8,1e12);
  h_qtotal_light = new TH1F("h_qtotal_light","Total Charge Light Tubes",200,q_xmin[3],q_xmax_used[3]*3);
  h_qtotal_ring = new TH1F("h_qtotal_ring","Total Charge Ring Tubes",200,q_xmin[15],q_xmax_used[15]*12);
  BinLogX(h_event_deltat);

  for (int ibin = 1; ibin < h_ringcandidate_npevspos->GetXaxis()->GetNbins()+1; ibin++) {
    for (int jbin = 1; jbin < h_ringcandidate_npevspos->GetYaxis()->GetNbins()+1; jbin++) {
      h_ringcandidate_timevspos->SetBinContent(ibin, jbin, -1000.);//(timeVsPos->GetMaximum() + timeVsPos->GetMinimum())/2.);
      h_ringcandidate_npevspos->SetBinContent(ibin, jbin, -1000.);//(timeVsPos->GetMaximum() + timeVsPos->GetMinimum())/2.);
    }
  }
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    if(pmtInfo->GetType(ipmt)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    h_ringcandidate_timevspos->Fill(pmtpos.X(), pmtpos.Y(), 1010.);
    h_ringcandidate_npevspos->Fill(pmtpos.X(), pmtpos.Y(), 1010.);
    // if(pmtidtopos[ipmt]!=1) h_ringcandidate_timevspos->Fill(pmtInfo->GetPosition(ipmt).X(),pmtInfo->GetPosition(ipmt).Y(),2.);
  }

  //Ring PMTs radius level
  for(int ih=0; ih<4; ih++){
    h_charge_ring.push_back(new TH1F(Form("h_charge_ring_%i",ih),"Ring Tubes Charges",100,-10,qtotal_xmax));
    h_charge_ring.back()->SetXTitle("NPEs");
    h_qratio_ring.push_back(new TH1F(Form("h_qratio_ring_%i",ih),"Ring Tubes Charge Ratio",30,-0.1,qrat_xmax));
    h_qratio_ring.back()->GetXaxis()->SetTitle("QRatio");
    h_time_ring.push_back(new TH1F(Form("h_time_ring_%i",ih), "Ring Tubes Event Time", t_nbins, ih==3? t_min[1]:t_min[0], ih==3? t_max[1]:t_max[0]));

    h_time_ring.back()->SetXTitle("Time Residuals (ns)");
  }

  //Individual PMTs
  for(int ih=0; ih<pmtInfo->GetPMTCount(); ih++){
    //MCTRUTH
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",200,0,200));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_time.push_back(new TH1F(Form("h_mcpmt_time_%i",ih),"h_mcpmt_time",200,5,10));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",200,5,10));
    //DAQ
 //   if (ih == 11){
 //     std::cout << "Create ring PMT 0 histogram" << std::endl;
 //     h_charge.push_back(new TH1F(Form("h_charge_%i",ih),pmtname[ih], 800, q_xmin[ih], q_xmax_used[ih]));
 //   }
 //   else
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),pmtname[ih], 200, q_xmin[ih], q_xmax_used[ih]));
    h_charge.back()->SetXTitle("NPEs");
    h_charge.back()->SetLineColor(pmtidtocolor[ih]);
    g_charge_correction.push_back(new TGraph());
    g_charge_correction.back()->SetName(Form("g_charge_correction_%i",ih));
    g_charge_correction.back()->SetTitle(pmtname[ih]);
    h_qratio.push_back(new TH1F(Form("h_qratio_%i",ih),pmtname[ih],30,-0.1,qrat_xmax));
    h_qratio.back()->SetXTitle("QRatio");
    h_qratio.back()->SetLineColor(pmtidtocolor[ih]);
    if((ih>=7 && ih<=10) || ih==23 || ih==24){
      BinLogX(h_charge.back());
    }
    // if(gTargetMaterial=="LAB")
    //  h_time.push_back(new TH1F(Form("h_time_%i",ih),pmtname[ih],t_nbins,lab_bins));
    // else
    h_time.push_back(new TH1F(Form("h_time_%i",ih),pmtname[ih],t_nbins, pmtidtopos[ih]==3? t_min[1]:t_min[0], pmtidtopos[ih]==3? t_max[1]:t_max[0]));
    h_time.back()->StatOverflows(kTRUE);
    h_time.back()->SetXTitle("Time Residuals (ns)");
    h_time.back()->SetLineColor(pmtidtocolor[ih]);
  }
  h_charge_muontrigs = new TH2F("h_charge_muontrigs","Top Tag vs Bottom Tag Charges",100,q_xmin[23],q_xmax_used[23],100,q_xmin[24],q_xmax_used[24]);
  BinLogX(h_charge_muontrigs);
  BinLogY(h_charge_muontrigs);
  h_time_muontrigs = new TH1F("h_time_muontrigs","Bottom Tag Time - Top Tag Time",100,-5.,5.);
  h_pmt_qratiovspos->SetMaximum(qrat_xmax);

}

void GetHistos(){

  std::cout<<" Get Histograms "<<std::endl;

  RAT::DS::Root *rds;
  RAT::DS::MC *mc;
  int clip_type = 0;
  int nentries = tree->GetEntries();
  // nentries = 2;
  std::cout<<" Number of entries: "<< nentries <<std::endl;
  for(int ientry=0; ientry<nentries;++ientry){

    if(nentries>NLOGENTRIES && ientry%(nentries/NLOGENTRIES) == 0) std::cout<<" Entry "<<ientry<<std::endl;

    tree->GetEntry(ientry);
    rds = dsreader->GetEvent(ientry);
    if(rds->ExistMC() == 1) {
      isMC = true;
      mc = rds->GetMC();
      clip_type = mc->GetID();

      //MC**************
      //MCPMT loop
      double mc_charge_cut = 0;
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
        RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
        int pmtid = mcpmt->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        if(pmttype == 1) {
           mc_charge_cut += mcpmt->GetMCPhotonCount();           
        }
      }
      h_mc_npe->Fill(mc_charge_cut);
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
        RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
        //Make grid display only for the small PMTs
        int pmtid = mcpmt->GetID();
        //count PE
        h_mcpmt_npe[pmtid]->Fill(mcpmt->GetMCPhotonCount());
        int pmttype = pmtInfo->GetType(pmtid);

        if(pmttype == 1 && (mc_charge_cut < total_npe_cut || total_npe_cut < 0)) {
          h_mcpmt_npevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),mcpmt->GetMCPhotonCount());
        }
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
      // ****************
    } //! ExistMC()

    //DAQ EVENTS******
    //Event loop
    int nevents = rds->GetEVCount();
    for(int ievt=0; ievt<nevents; ievt++){ // not sure over what we are looping over here.

      //Load event
      RAT::DS::EV *ev = rds->GetEV(ievt);

      // Get muon tags, source trigger and muon panels charge and time
      //Empty channel for Q correction
      double corr_charge = 0.;
      double ref_charge = 0.;
      RAT::DS::PMT *pmt = ev->GetPMTWithID(CORRCH); //ref ring tube
//      RAT::DS::PMT *pmt2 = ev->GetPMTWithID(CORRCH+1); //ref ring tube
//      if(pmt!=NULL && pmt2!=NULL && !isMC){
      if(pmt!=NULL && !isMC){
        corr_charge = pmt->GetCharge();
//        std::cout << corr_charge << std::endl;
//        corr_charge2 = pKmt2->GetCharge();
      }
      pmt  = ev->GetPMTWithID(REFCH); //ref ring tube
      if (pmt)
        ref_charge = pmt->GetCharge();
      
      if(ievt){
        std::cout << "WARNING the 2nd level data event loop is being used" << std::endl;
      }
      g_corr_channel->SetPoint(g_corr_channel->GetN(), ientry, corr_charge);
      //Reference Ring Tube
      double refring_time = -9999.;
      double refring_timeres = -9999.;

      pmt = ev->GetPMTWithID(REFTUBE); //ref ring tube
      if(pmt!=NULL) {
        refring_time = pmt->GetTime();
//        refring_time = 	GetEventTime();
        double refring_dist = (pmtInfo->GetPosition(REFTUBE) - *target_pos).Mag();
        double refring_tof = refring_dist/cspeed; // this can't really work if my reftube is the movable trigger cube
        refring_tof = tof_fixed[pmtidtopos[REFTUBE]];
        refring_timeres = refring_time - refring_tof - time_delay[REFTUBE];
        refring_time = refring_timeres;
//        std::cout<<" ring_timeres "<<refring_timeres<<" "<<refring_time<<" "<<refring_tof<<" "<<time_delay[REFTUBE]<<std::endl;
      }

//      refring_time = GetEventTime(ev, corr_charge, TIMENPELIMIT );
      //Trigger PMT
      double trigger_q = 0.;
      double trigger_time = -9999.;
      pmt = ev->GetPMTWithID(25);
      if(pmt!=NULL) {
        trigger_q = pmt->GetCharge();
        trigger_time = pmt->GetTime();
      }
      //Top Muon Tag
      double topmuon_charge = 0.;
      double topmuon_time = 0.;
      double topmuon_fft = 0.;
      pmt = ev->GetPMTWithID(23); //top tag
      if(pmt!=NULL) {
        topmuon_charge = pmt->GetCharge();
        topmuon_time = pmt->GetTime();
        topmuon_fft = pmt->GetFCN();
      }
      //Bottom Muon Tag
      double bottommuon_charge = 0.;
      double bottommuon_time = -9999.;
      double bottommuon_fft = 0.;
      pmt = ev->GetPMTWithID(24); //bottom tag
      if(pmt!=NULL) {
        bottommuon_charge = pmt->GetCharge();
        bottommuon_time = pmt->GetTime();
        bottommuon_fft = pmt->GetFCN();
      }
//      if (topmuon_charge < 10 or bottommuon_charge < 10)
//          continue;
      //Veto pannels
      double panel_charge[] = {0., 0., 0., 0.};
      pmt = ev->GetPMTWithID(6);
      if(pmt!=NULL) panel_charge[0] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(7);
      if(pmt!=NULL) panel_charge[1] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(8);
      if(pmt!=NULL) panel_charge[2] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(9);
      if(pmt!=NULL) panel_charge[3] = pmt->GetCharge();

      //Ring PMTs
      //Loop over all the ring hits, if one has a bad pedestal, reject the
      //event for the ring PMTs
      bool skipEvent = false;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        if(pmttype != 1) continue;
        double charge = ev->GetPMT(ipmt)->GetCharge();
        if(charge < -9900) {
          skipEvent = true;
          break;
        }
      }

      //Time wrt last event
      double deltat = (ev->GetDeltaT());
      //double deltat = (ev->GetClockTime() - lastbursttime) *2;
      // std::cout<<" Clock time "<<ev->GetClockTime()<<" "<<lastbursttime<<" "<<deltat<<std::endl;

      //Get event time: for cosmic that we use the median of the first 4 tubes
      //while for source data we use the trigger time
      double event_time = -9999.;
      if(isCosmicData){
        event_time = -99.9;
        //event_time = ev->GetEventTime(); // non existent in our version!
        //if(event_time < -9900.) continue; //Cosmics: Less than 3 ring hits, Source: no defined trigger time
        //event_time = bottommuon_time - 18.3; //DT
        //event_time = bottommuon_time - 18.6; //MC
        h_ringcandidate_timevspos->SetMaximum(10.);
        if(gTargetMaterial == "LAB") h_ringcandidate_timevspos->SetMaximum(10.);
        h_ringcandidate_timevspos->SetMinimum(-2.0);
      }
      else if(isSourceData){
        if(isMC) event_time = 91.; //MC
        else event_time = refring_time; //DT
      }

      //Retrieve PMT information and fill histograms
      std::vector<double> charge_ring(4,0.), qshort_ring(4,0.);
      ref_charge -= ((corr_charge-charge_offsets[REFCH])*charge_slopes[REFCH]);
      
      double qtotal_lightpmts = 0.;
      double qtotal_ringpmts = 0.;
      double npe_ref_light_tube = 0.;
      double qratio_ref_light_tube = 0.;
      int nhits_ring = 0, nhits_light = 0;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        double charge = ev->GetPMT(ipmt)->GetCharge();
        double npes = charge/spe[pmtid];
        if(pmttype==1){
//          std::cout << corr_charge << std::endl;
//          charge -= corr_charge;
          charge -= ((corr_charge-charge_offsets[pmtid])*charge_slopes[pmtid]);
//          std::cout << ((corr_charge-charge_offsets[pmtid])*charge_slopes[pmtid]) << std::endl;
//          npes = charge/spe[pmtid];
          npes = charge/spe[pmtid]-noise_offsets[pmtid];
          npes = npes / (spe_corrections[pmtid] - noise_offsets[pmtid]);
          if (pmtid == REFTUBE) 
             qratio_ref_light_tube = ev->GetPMT(ipmt)->GetQShort()/charge;
     	     npe_ref_light_tube = npes;
          if((npes < -3.0 || npes > total_npe_cut) && total_npe_cut > 10) {
           std::cout << "Warning RING PMT with weird NPEs " << npes << " for PMT " << pmtid << std::endl;
           std::cout << "for event " << ientry << " " << ievt << std::endl;
//           npes = 999;
          }
          nhits_ring++;
          if ( pmtid_reliable[pmtid] && pmtid != REFTUBE)
          	qtotal_ringpmts += npes;
        } else if(pmttype==2){
          nhits_light++;
          qtotal_lightpmts += npes;
        }
        //Reject event if bad pedestals
        //if(pmttype == 1 && skipEvent) continue;
        //Calculate total charge
      }
//      std::cout << "Total NPEs in RING " << qtotal_ringpmts << std::endl;

      if( total_npe_cut > 0 && qtotal_ringpmts > total_npe_cut){
        std::cout << "Total NPE cut met - skipping event in RING " << qtotal_ringpmts << std::endl;
        continue;
      }
//      if(true){
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);

        //Hit selection
        double charge = ev->GetPMT(ipmt)->GetCharge();
        double charge_0 = charge;
        //Apply charge correction
//        if(pmttype==1) charge -= corr_charge;
        if(pmttype==1) charge -= ((corr_charge-charge_offsets[pmtid])*charge_slopes[pmtid]);
        double qshort = ev->GetPMT(ipmt)->GetQShort();
        if(charge < -9900) continue; //This is a tube with a bad pedestal
        double npes = charge/spe[pmtid];
        if(pmttype==1){
          npes = charge/spe[pmtid]-noise_offsets[pmtid];
          npes = npes / (spe_corrections[pmtid] - noise_offsets[pmtid]);
        }
        h_charge[pmtid]->Fill(npes);
//        std::cout << corr_charge << std::endl;
//        if((abs(corr_charge) + abs(charge_0)) < 200){
//           if( npes < 1.8)
//           g_charge_correction[pmtid]->SetPoint(g_charge_correction[pmtid]->GetN(), corr_charge, charge);
           g_charge_correction[pmtid]->SetPoint(g_charge_correction[pmtid]->GetN(), ref_charge, charge);
//        }
        //There are some outliers with crazy qshort and null charge probably
        //due to the baseline drifts. I remove those with the cut below.
        if(qshort/charge<2.0 && qshort/charge>-2.0 ) {
          h_qratio[pmtid]->Fill(qshort/charge);
          if(pmtInfo->GetType(pmtid)==1){
//            if(1){
//              std::cout << "NPES for PMT "<< pmtid << ": " << npes << std::endl;
//              std::cout << "total npes " << qtotal_ringpmts << endl;
//            }
            h_pmt_chargevspos->Fill(pmtInfo->GetPosition(pmtid).X(), pmtInfo->GetPosition(pmtid).Y(), npes);
            h_pmt_qratiovspos->Fill(pmtInfo->GetPosition(pmtid).X(), pmtInfo->GetPosition(pmtid).Y(), qshort/charge);
          }
        }

        double pmttime = ev->GetPMT(ipmt)->GetTime();
        if(pmttime < -9900) continue; //This tube didn't cross threshold

        //Count hits

        //Calculate charges by radius
        if (pmtid_reliable[pmtid]){
          charge_ring[pmtidtopos[pmtid]] += npes;
          qshort_ring[pmtidtopos[pmtid]] += qshort/spe[pmtid];
        }
        //Compute time residuals
        double dist = (pmtInfo->GetPosition(pmtid) - *target_pos).Mag();
        double tof = dist/cspeed;
        tof = tof_fixed[pmtidtopos[pmtid]];
        double timecorr = pmttime - tof - time_delay[pmtid];
        double timeres = timecorr - event_time;
        if (SMEARING > 0) timeres += trand3->Gaus(0., SMEARING);
        // double pmtfcn = ev->GetPMT(ipmt)->GetFCN();
        //std::cout<<" ToF "<<pmtid<<": "<<tof<<" "<<dist<<std::endl;
        //if(pmtid != 13 && pmtid != 22 && pmtid != 19 && pmtid != 16) continue;
        //Event level averaged
        if( npe_ref_light_tube > noise_cutoffs[REFTUBE]){
          if (npes > noise_cutoffs[pmtid]){
            if (qshort/charge > 0.2 && qratio_ref_light_tube > 0.2)
              h_time[pmtid]->Fill(timeres);
        
//        h_time[pmtid]->Fill(pmttime-time_delay[pmtid]);
           if(pmttype==1 && pmtid != REFTUBE){
            h_pmt_timevspos->Fill(pmtInfo->GetPosition(pmtid).X(), pmtInfo->GetPosition(pmtid).Y(), timeres);
            if (pmtid_reliable_time[pmtid]){
              if (qshort/charge > 0.2 && qratio_ref_light_tube > 0.2){
                h_timevsnpe->Fill(timeres, npes);
                h_time_ring[pmtidtopos[pmtid]]->Fill(timeres);
              }
              h_timevsqratio->Fill(timeres, qshort/charge);
            }
          }
         }
        }

        if(ientry==RING_CANDIDATE){
          h_ringcandidate_npevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),charge/spe[pmtid]);
          //if(pmtidtopos[pmtid]==1)
          h_ringcandidate_timevspos->Fill(pmtInfo->GetPosition(pmtid).X(),pmtInfo->GetPosition(pmtid).Y(),timeres - event_time - 10.); //TWEAK
          // std::cout<<" time "<<pmtid<<" "<<timeres - event_time<<std::endl;
        }       
      } //end PMT loop

      //Fill EV-level histograms
      h_nhits_ring->Fill(nhits_ring);
      h_nhits_light->Fill(nhits_light);
      h_qtotal_light->Fill(qtotal_lightpmts);
      h_qtotal_ring->Fill(qtotal_ringpmts);
      h_event_time->Fill(event_time);
      h_event_deltat->Fill(deltat);
      h_charge_muontrigs->Fill(bottommuon_charge, topmuon_charge);
      h_time_muontrigs->Fill(bottommuon_time - topmuon_time);
      for(int irad=0; irad<4; irad++){
        if(charge_ring[irad] !=0 ) {
//          std::cout <<"Filling event in h_charge_ring"<< endl;
          h_charge_ring[irad]->Fill(charge_ring[irad]);
          h_qratio_ring[irad]->Fill(qshort_ring[irad]/charge_ring[irad]);
        }
        else {
//          std::cout << "charge_ring value is 0" << endl;
        }
      }

    } //end DS::EV loop

  } //end ds entry loop

  //Create Stacks for Ring Tubes - do fill in event loop to be able to apply
//  for(int ipmt=0; ipmt<pmtInfo->GetPMTCount();ipmt++){
//      h_time_ring[pmtidtopos[ipmt]]->Add(h_time[ipmt]);
//  }


} //end GetHistos


//Draw histrograms
void DrawHistos(){
  int nentries = tree->GetEntries();

  //Ring candidates
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    if(pmtInfo->GetType(ipmt)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    h_ringcandidate_npevspos->Fill(pmtpos.X(), pmtpos.Y(), 1.5);
  }
/*  TString s = "Charge correction channel";
  s += string(CORRCH);
*/
  TCanvas* c_corr_channel = new TCanvas("c_corr_channel", "blub");
  std::cout << "Printing the channel to correct on with x entries: " <<std::endl;
  std::cout << g_corr_channel->GetN() << std::endl;
  g_corr_channel->SetMarkerStyle(1);
  g_corr_channel->Draw("AP");
  c_corr_channel->Update();

  TCanvas *c_ring_candidate = new TCanvas("c_ring_candidate","Ring Candidate",1200,600);
  c_ring_candidate->Divide(2,1);
  c_ring_candidate->cd(1);
  h_ringcandidate_npevspos->Draw("colz text");
  c_ring_candidate->cd(2);
  h_ringcandidate_timevspos->Draw("colz text");
  c_ring_candidate->Update();

  //General plots
  TCanvas *c_ring_event = new TCanvas("c_ring_event","Ring Event",1200,1200);
  c_ring_event->Divide(3,3);
  c_ring_event->cd(1);
  h_pmt_chargevspos->Draw("colz text");
  c_ring_event->cd(2);
  h_pmt_qratiovspos->Draw("colz text");
  c_ring_event->cd(3);
  h_pmt_timevspos->Draw("colz text");
  c_ring_event->cd(4);
  h_event_time->Draw();
  c_ring_event->cd(5);
  h_qtotal_light->Draw();
  c_ring_event->cd(6);
  h_qtotal_ring->Draw();
  c_ring_event->cd(7)->SetLogy();
  c_ring_event->cd(7)->SetLogx();
  h_event_deltat->Draw();
  c_ring_event->cd(8);
  h_nhits_ring->GetXaxis()->SetTitle("NHits");
  h_nhits_ring->Draw();
  c_ring_event->cd(9);
  h_nhits_light->GetXaxis()->SetTitle("NHits");
  h_nhits_light->Draw();
  c_ring_event->Update();



  //Light tubes
  TCanvas *c_charge_light;
  c_charge_light = new TCanvas("c_charge_light","Charge Light Tubes",900,600);
  c_charge_light->Divide(3,2);

  TCanvas *c_time_light;
  c_time_light = new TCanvas("c_time_light","Time Light Tubes",900,600);
  c_time_light->Divide(3,2);

  int padCount = 1;
  for(int pmtid = 1; pmtid<=6; pmtid++){
    c_charge_light->cd(padCount)->SetLogy();
    h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    //h_charge[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    //h_charge[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_charge[pmtid]->Draw();
    c_time_light->cd(padCount);
    h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    //h_time[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    //h_time[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_time[pmtid]->GetXaxis()->SetLabelSize(0.08);
    h_time[pmtid]->Draw();
    padCount++;
  }

  //Pannels
  TCanvas *c_charge_vetos;
  c_charge_vetos = new TCanvas("c_charge_vetos","Charge Vetos",600,600);
  c_charge_vetos->Divide(2,2);

  TCanvas *c_time_vetos;
  c_time_vetos = new TCanvas("c_time_vetos","Time Vetos",600,600);
  c_time_vetos->Divide(2,2);

  padCount = 1;
  for(int pmtid = 7; pmtid<=10; pmtid++){
    c_charge_vetos->cd(padCount)->SetLogy();
    c_charge_vetos->cd(padCount)->SetLogx();
    h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    //h_charge[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    //h_charge[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_charge[pmtid]->Draw();
    c_time_vetos->cd(padCount);
    h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    //h_time[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    //h_time[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_time[pmtid]->GetXaxis()->SetLabelSize(0.08);
    h_time[pmtid]->Draw();
    padCount++;
  }

  //Ring tubes
  TCanvas *c_charge_ring;
  c_charge_ring = new TCanvas("c_charge_ring","Charge Ring Tubes",1200,1200);
  Crucify(c_charge_ring);

  TCanvas *c_charge_corr_ring;
  c_charge_corr_ring = new TCanvas("c_charge_corr_ring","Charge correction Ring Tubes",1200,1200);
  Crucify(c_charge_corr_ring);

  TCanvas *c_qratio_ring;
  c_qratio_ring = new TCanvas("c_qratio_ring","QRatio Ring Tubes",1200,1200);
  Crucify(c_qratio_ring);

  TCanvas *c_time_ring;
  c_time_ring = new TCanvas("c_time_ring","Time Ring Tubes",1200,1200);
  Crucify(c_time_ring);

  padCount = 1;
  for(int pmtid = 11; pmtid<=22; pmtid++){
    c_charge_ring->cd(padCount)->SetLogy();
    h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    // h_charge[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    // h_charge[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_charge[pmtid]->Draw();

    c_charge_corr_ring->cd(padCount);
    g_charge_correction[pmtid]->SetMarkerColor(pmtidtocolor[pmtid]);
    g_charge_correction[pmtid]->SetMarkerStyle(1);
    g_charge_correction[pmtid]->Draw("AP");
    g_charge_correction[pmtid]->GetXaxis()->SetTitle(Form("Charge corr channel %i",CORRCH));
    g_charge_correction[pmtid]->GetYaxis()->SetTitle(Form("Charge channel %i",pmtid));
    c_charge_corr_ring->Update();

    c_qratio_ring->cd(padCount);
    h_qratio[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    // h_qratio[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    // h_qratio[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_qratio[pmtid]->Draw();
    c_time_ring->cd(padCount)->SetLogy();
    h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    // h_time[pmtid]->SetFillColor(pmtidtocolor[pmtid]);
    // h_time[pmtid]->SetFillStyle(3003+pmtidtopos[pmtid]);
    h_time[pmtid]->GetXaxis()->SetLabelSize(0.08);
    h_time[pmtid]->Draw();
    padCount++;
  }
  std::cout << "Drawing muon tags" << std::endl;
  //Muon tags
  TCanvas *c_tags;
  c_tags = new TCanvas("c_tags","Cosmic tags",600,600);
  c_tags->Divide(2,2);

  c_tags->cd(1)->SetLogy();
  h_charge[23]->SetLineColor(pmtidtocolor[22]);
  h_charge[23]->SetFillColor(pmtidtocolor[22]);
  h_charge[23]->SetFillStyle(3003+pmtidtopos[22]);
  h_charge[23]->Draw();
  c_tags->cd(2)->SetLogy();
  h_charge[24]->SetLineColor(pmtidtocolor[22]);
  h_charge[24]->SetFillColor(pmtidtocolor[22]);
  h_charge[24]->SetFillStyle(3003+pmtidtopos[22]);
  h_charge[24]->Draw();
  c_tags->cd(3)->SetLogy();
  c_tags->cd(3)->SetLogx();
  h_charge_muontrigs->Draw("colz");
  c_tags->cd(4)->SetLogy();
  h_time_muontrigs->Draw();

  std::cout << "Drawing trigger PMT" << std::endl;
  //Source Trigger PMT
  TCanvas *c_trigger;
  c_trigger = new TCanvas("c_trigger","Trigger PMT", 600,300);
  c_trigger->Divide(2,1);

  c_trigger->cd(1);
  h_charge.back()->Draw();
  c_trigger->cd(2);
  h_time.back()->Draw();

  std::cout << "Updating canvases" << std::endl;

  c_charge_ring->Update();
  c_charge_light->Update();
  c_tags->Update();
  c_trigger->Update();
  c_charge_vetos->Update();
  c_time_ring->Update();
  c_time_light->Update();
  c_trigger->Update();
  c_time_vetos->Update();
  c_qratio_ring->Update();

  std::cout << "Drawing radial PMT groups" << std::endl;

  TCanvas *c_ring_tubes_rad = new TCanvas("c_charge_ring_rad","Ring Tube Radial Groups",1800,600);
  c_ring_tubes_rad->Divide(3,1);
  c_ring_tubes_rad->cd(1);
  h_charge_ring[0]->SetLineColor(pmtidtocolor[PMTID_INNER]);
  h_charge_ring[1]->SetLineColor(pmtidtocolor[PMTID_MID]);
  h_charge_ring[2]->SetLineColor(pmtidtocolor[PMTID_OUTER]);
  h_charge_ring[0]->Draw();
  h_charge_ring[1]->Draw("sames");
  h_charge_ring[2]->Draw("sames");
  c_ring_tubes_rad->cd(2);
  h_qratio_ring[0]->SetLineColor(pmtidtocolor[PMTID_INNER]);
  h_qratio_ring[1]->SetLineColor(pmtidtocolor[PMTID_MID]);
  h_qratio_ring[2]->SetLineColor(pmtidtocolor[PMTID_OUTER]);
  h_qratio_ring[0]->Draw();
  h_qratio_ring[1]->Draw("sames");
  h_qratio_ring[2]->Draw("sames");
  c_ring_tubes_rad->cd(3)->SetLogy();
  h_time_ring[0]->SetLineColor(pmtidtocolor[PMTID_INNER]);
  h_time_ring[1]->SetLineColor(pmtidtocolor[PMTID_MID]);
  h_time_ring[2]->SetLineColor(pmtidtocolor[PMTID_OUTER]);
  if(gTargetMaterial=="LAB"){
    // h_time_ring[0]->Scale(1./10.,"width");
    // h_time_ring[1]->Scale(1./10.,"width");
    // h_time_ring[2]->Scale(1./10.,"width");
    h_time_ring[0]->Draw();
  } else if(gTargetMaterial=="LABPPO"){
    h_time_ring[0]->Draw();
  } else{
    h_time_ring[1]->Draw();
  }
  h_time_ring[1]->Draw("sames");
  h_time_ring[0]->Draw("sames");
  h_time_ring[2]->Draw("sames");
  int maxbin = h_time_ring[0]->GetXaxis()->FindBin(0.0);
  double total = h_time_ring[0]->Integral(0,maxbin) + h_time_ring[1]->Integral(0,maxbin) + h_time_ring[2]->Integral(0,maxbin);
  std::cout<<" time 0 " << h_time_ring[0]->Integral(0,maxbin)/total<<std::endl;
  std::cout<<" time 1 " << h_time_ring[1]->Integral(0,maxbin)/total<<std::endl;
  std::cout<<" time 2 " << h_time_ring[2]->Integral(0,maxbin)/total<<std::endl;
  c_ring_tubes_rad->Update();


  TCanvas *c_event_times = new TCanvas("c_event_times", "c_event_times", 900, 1000);
  c_event_times->Divide(1,2);
  c_event_times->cd(1);
  h_timevsnpe->Draw("colz");
  cout << "h_timevsnpe -> GetEntries() " << h_timevsnpe->GetEntries();
  c_event_times->cd(2);
  h_timevsqratio->Draw("colz");
  c_event_times->Update();

  //MCTruth
  TCanvas *c_mc_npevspos = new TCanvas("c_mc_npevspos","c_mc_npevspos",900,1000);
//  std::cout << "nentries for error calc " << nentries << "\t" << h_mcpmt_npevspos->GetSize()<<std::endl;
  for (int a_bin = 1;  a_bin < h_mcpmt_npevspos->GetSize(); a_bin++){   
//    std::cout << h_mcpmt_npevspos->GetBinContent(a_bin) << std::endl;
    if (h_mcpmt_npevspos->GetBinContent(a_bin) > 0 ){
       std::cout << nentries << std::endl;
       h_mcpmt_npevspos->SetBinError(a_bin, sqrt(h_mcpmt_npevspos->GetBinContent(a_bin))/(double)nentries);
       h_mcpmt_npevspos->SetBinContent(a_bin, h_mcpmt_npevspos->GetBinContent(a_bin)/(double)nentries);
    }
  }
  h_mcpmt_npevspos->Draw("colz texte");

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
        h_mc_npe->Draw();
//        h_mcpmt_fetime[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
//        h_mcpmt_fetime[pmtid]->Draw(opt);
        firstdrawn1 = true;
      } else if(pmttype==0){ //Trigger tube
        const char *opt = firstdrawn2 ? "sames" : "";
        c_mc->cd(7);
        h_mcpmt_npevspos->Draw("colz text");
        // h_mcpmt_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        // h_mcpmt_charge[pmtid]->Draw(opt);
        c_mc->cd(8);
        h_mcpmt_npe[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_npe[pmtid]->Draw(opt);
        c_mc->cd(9);
        h_mcpmt_fetime[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_fetime[pmtid]->Draw(opt);
        firstdrawn2 = true;
      }
    }
  }
  c_mc->Update();
}



//Output histrograms to file
void PrintHistos(char *filename){

  std::cout<<" Output histrograms to "<<filename<<std::endl;

  TFile *fout = new TFile(filename,"RECREATE");
  fout->cd();
  h_time_muontrigs->Write(); //For event counting and norm
  for(int ipmt=0; ipmt<pmtInfo->GetPMTCount();ipmt++){
    g_charge_correction[ipmt]->Write();
    h_charge[ipmt]->Write();
    h_time[ipmt]->Write();
  }
  for(int ih=0; ih<4;ih++){
    h_charge_ring[ih]->Write();
    h_time_ring[ih]->Write();
  }
  h_qtotal_light->Write(); //Total charge in the event
  h_qtotal_ring->Write(); //Total charge in the event

  h_pmt_timevspos->Write();
  h_pmt_chargevspos->Write();
  h_ringcandidate_npevspos->Write();
  h_ringcandidate_timevspos->Write();
  fout->Close();

}

void NormalizeHistos(){

//  double norm = h_event_time->GetEntries();
  double norm = h_qtotal_ring->GetEntries();
  h_pmt_qratiovspos->Scale(1./norm);
  h_pmt_chargevspos->Sumw2();
  h_pmt_chargevspos->Scale(1./norm);
  h_pmt_timevspos->Scale(1./norm);

  // h_pmt_chargevspos->Divide(h_pmt_countvspos);
  //h_pmt_timevspos->Divide(h_pmt_countvspos);

  //Fill with zeroes
  for (int ibin = 1; ibin < h_pmt_timevspos->GetXaxis()->GetNbins()+1; ibin++) {
    for (int jbin = 1; jbin < h_pmt_timevspos->GetYaxis()->GetNbins()+1; jbin++) {
      h_pmt_qratiovspos->SetBinContent(ibin, jbin, -1000. + h_pmt_qratiovspos->GetBinContent(ibin, jbin));
      h_pmt_timevspos->SetBinContent(ibin, jbin, -1000. + h_pmt_timevspos->GetBinContent(ibin, jbin));
    }
  }
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    if(pmtInfo->GetType(ipmt)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    h_pmt_qratiovspos->Fill(pmtpos.X(), pmtpos.Y(), 1000.);
    h_pmt_timevspos->Fill(pmtpos.X(), pmtpos.Y(), 1000.);
  }

  h_pmt_qratiovspos->SetMinimum(0.);
  h_pmt_chargevspos->SetMaximum(40.);
  h_pmt_chargevspos->SetMinimum(0.0);
  h_pmt_timevspos->SetMaximum(1.5);
  if(gTargetMaterial == "LAB") h_pmt_timevspos->SetMaximum(10.);
  h_pmt_timevspos->SetMinimum(-0.5);
  h_pmt_qratiovspos->SetMarkerColor(0);
  h_pmt_chargevspos->SetMarkerColor(0);
  h_pmt_timevspos->SetMarkerColor(0);

  h_ringcandidate_npevspos->SetMaximum(20.);
  h_ringcandidate_npevspos->SetMinimum(0.);

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


void ParseArgs(int argc, char **argv){
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-t") {
      if(std::string(argv[i+1]) == "c") isCosmicData = true;
      else if (std::string(argv[i+1]) == "s") isSourceData = true;
    }
    if(std::string(argv[i]) == "-i") {gInputFile = argv[++i];}
    if(std::string(argv[i]) == "-o") {gOutFile = argv[++i];}
    if(std::string(argv[i]) == "-m") {gTargetMaterial = argv[++i];}
  }

  if(argc<=1 || (!isCosmicData && !isSourceData) || (isCosmicData && isSourceData) ){
    std::cerr<<" Usage: ./DrawChargeAndTime.exe -t c(osmic) || s(ource) -i INPUT_FILE [-o OUTPUT_FILE -m TARGET_MATERIAL]"<<std::endl;
    exit(0);
  }
}
