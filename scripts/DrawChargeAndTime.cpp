#include<iostream>
#include<fstream>

#include<TF1.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<TBrowser.h>
#include<TGraph.h>
#include<TColor.h>
#include<TMath.h>

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include <RAT/DB.hh>
#include <RAT/DS/RunStore.hh>

#define MCPHOTONLOOP false
#define NLOGENTRIES 10

//Constants
double cspeed = 300/1.4; // (mm/ns)/rindex
TVector3* target_pos = new TVector3(-400.,-400.,-200.);
TVector3 centerpos(0,0,0);

//Globals
RAT::DSReader *dsreader;
TTree *tree;
TTree *runT;

//Methods
char *gInputFileMC = NULL;
std::vector<char*> gInputFileDT;
char *gOutFile = NULL;
char *gTargetMaterial = "LAB";
void ParseArgs(int argc, char **argv);
void GetPMTInfo();
void GetDBTables();
void GetHistos();
void DrawHistos();
void PrintHistos(char*);
void NormalizeHistos();
void ExtractSPE();

//Fit function
Double_t fmultigaus(Double_t *x, Double_t *par) {
  Double_t gaus_noise = 1./(par[2]*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[1] )/par[2] , 2.) );
  Double_t gaus_spe = 1./(par[5]*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[4] )/par[5] , 2.) );
  Double_t gaus_2pe = 1./(par[5]*sqrt(2.)*sqrt(2.*3.14159)) * exp( -0.5 * pow( ( x[0] - par[4]*2. )/(par[5]*sqrt(2.)), 2.) );
  gaus_noise = 0.;
  return par[0]*gaus_noise + par[3]*gaus_spe + par[6]*gaus_2pe;
}

//Global variables
std::vector<TVector3*> pos_pmts; //PMT positions by ID
int npmts; // # PMTs
std::map<int,int> npmts_type; // # PMTs per type
std::vector<Color_t> pmtidtocolor;
std::vector<int> pmtidtopos;
std::vector<double> pmttime_delay;
vector<double> qScintCorr; //Scintillation correction
vector<double> qScintCorrErr; // Scintillation correction error
vector<double> spe; //SPE


//// Histograms
//Event level
TH1F* h_event_time;
TH1F* h_charge_total; //Total charge in the event
TH2F* h_pmt_chargevspos; //PMT charge vs PMT position
TH2F* h_pmt_npevspos; //NPE vs PMT position
TH2F* h_pmt_timevspos; //PMT charge vs PMT position
TH2F* h_pmt_countvspos; //Hit counter for normalization
TH1F* h_chi2; //Chi2 cher vs scint
TH2F* h_mcpmt_npevspos; //MC NPE vs PMT position

//Hit level (PMT)
std::vector<double> q_xmin;
std::vector<double> q_xmax;
std::vector<TF1*> f_spe;
std::vector<double> f_spe_mean;
std::vector<double> f_spe_sigma;
std::vector<TH1F*> h_charge; //Measured charge
std::vector<TH1F*> h_charge_res; //Measured charge geometry corrected
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_bottom; //Measured time
std::vector<TH1F*> h_time_trigger; //Measured time
std::vector<TH1F*> h_time_event;
std::vector<TH1F*> h_time_diff; //Time diff between EV-MC
std::vector<TH2F*> h_charge_vs_trigq; //PMT charge vs trigger charge
TH2F* h_charge_muontrigs; //Muon trigger charges correlation
TH1F* h_time_muontrigs;
std::vector<TH1F*> h_mcpmt_npe; //PEs by PMT ID
std::vector<TH1F*> h_mcpmt_charge; //MC charge
std::vector<TH1F*> h_mcpmt_time; //MC FE time
std::vector<TH1F*> h_mcpmt_fetime; //MC FE time

//MC
std::vector<TH1F*> h_mctime; //Measured time
std::vector<TH1F*> h_mctime_res; //Measured time


//Real data
//std::vector<TH1F*> h_dt_charge; //Measured charge
RAT::DS::PMTInfo *pmtInfo;

int main(int argc, char **argv){

  //Init********
  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy("App", &appargc, appargv);
  ParseArgs(argc, argv);
  //************

  GetPMTInfo();
  GetDBTables();
  GetHistos();
  NormalizeHistos();
  ExtractSPE();

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
  db->Load("/Users/snoplus/Work/TheiaRnD/rat-pac/data/TheiaRnD/SCINTCORR.ratdb");
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

}

void GetPMTInfo(){

  //Init pmt positions
  dsreader = new RAT::DSReader(gInputFileMC);
  tree = dsreader->GetT();
  runT = dsreader->GetRunT();
  RAT::DS::Run *run = 0;
  runT->SetBranchAddress("run",&run);
  runT->GetEntry(0);
  pmtInfo = run->GetPMTInfo();
  npmts = pmtInfo->GetPMTCount();
  npmts_type = pmtInfo->GetPMTTypeCount();

  for (size_t ipmt = 0; ipmt < npmts; ipmt++) {
    TVector3* pmt_temp = new TVector3(pmtInfo->GetPosition(ipmt)[0],
                                      pmtInfo->GetPosition(ipmt)[1],
                                      pmtInfo->GetPosition(ipmt)[2]);
    pos_pmts.push_back(pmt_temp);
  }

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

  Color_t mycolors[] = {1, 1, 1, 1, 1, 1,  1, 1,  kRed, 1, 1, 1,  kOrange, kBlue, kRed, kRed, kBlue, kOrange, kOrange, kBlue, kRed, kRed, kBlue, kOrange, 1}; //By Position (WATER)
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kRed, kOrange, kBlue, kBlue, kOrange, kRed, kRed, kOrange, kBlue, 1}; //By Position (LAB)
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue-2, kOrange-2, kRed-2, kRed-1, kOrange-1, kBlue-1, kBlue, kOrange, kRed, kRed+1, kOrange+1, kBlue+1, 1}; //By Position
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kBlue, kBlue, kBlue, kOrange, kOrange, kOrange, kOrange, kRed, kRed, kRed, kRed, 1}; //By Digitizer group
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kCyan, kBlack, kGray, kGreen, kTeal, kAzure, kViolet, kPink, kYellow, 1}; //Individual

  //int mypmtpos[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 1, 0, 0, 1, 2, 2, 1, 0, 3}; //In space
  int mypmtpos[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3}; //In digit group

  pmtidtocolor.insert(pmtidtocolor.begin(), mycolors, mycolors + npmts );
  pmtidtopos.insert(pmtidtopos.begin(), mypmtpos, mypmtpos + npmts );
  double mydelays[] = {.0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0,
                    -20.57914449553667,
                    -20.659866192027934,
                    -20.827556184469994,
                    -20.83027277200025,
                    -20.41955157599053,
                    -20.42802016610127,
                    -20.565785808047142,
                    -20.607075743765186,
                    -20.628455171240507,
                    -20.65168715965919,
                    -20.70112435087266,
                    -20.450410008571303,
                    .0};

  pmttime_delay.insert(pmttime_delay.begin(), mydelays, mydelays + npmts );

  //Plot axis limits
  double myqxmins[] = {-100,-100,-100,-100,-100,-100, -100,-100, -100,-100,-100,-100, -100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100, -100};
  q_xmin.insert(q_xmin.begin(), myqxmins, myqxmins + npmts );
  double myqxmaxs[] = {3000,3000,3000,3000,3000,3000, 100000,100000, 100000,100000,100000,100000, 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000, 3000};
  q_xmax.insert(q_xmax.begin(), myqxmaxs, myqxmaxs + npmts );

}

void GetHistos(){

  //MC
  std::cout<<" Get PDFs "<<std::endl;

  //Init histos
  //  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",27,-27*15,27*15,27,-27*15,27*15);
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_countvspos = new TH2F("h_pmt_countvspos","h_pmt_countvspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_chargevspos = new TH2F("h_pmt_chargevspos","h_pmt_chargevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_npevspos = new TH2F("h_pmt_npevspos","h_pmt_npevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  h_pmt_timevspos = new TH2F("h_pmt_timevspos","h_pmt_timevspos",7, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 7, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);

  for(int ih=0; ih<npmts; ih++){
    //MCTRUTH
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",200,0,200));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_time.push_back(new TH1F(Form("h_mcpmt_time_%i",ih),"h_mcpmt_time",200,5,10));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",200,5,10));
    //DAQ
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),"h_charge",500,q_xmin[ih],q_xmax[ih]));
    h_charge_res.push_back(new TH1F(Form("h_charge_res_%i",ih),"h_charge_res",50,0,50));
    h_charge_vs_trigq.push_back(new TH2F(Form("h_charge_vs_trigq_%i",ih),"h_charge_vs_trigq",200,0,100,200,0,100));
    h_time.push_back(new TH1F(Form("h_time_%i",ih),"h_time",200,-4,4));
    h_time_bottom.push_back(new TH1F(Form("h_time_bottom_%i",ih),"h_time_bottom",50,0.0,5));
    h_time_trigger.push_back(new TH1F(Form("h_time_trigger%i",ih),"h_time_trigger",200,-10.0,10.0));
    h_time_event.push_back(new TH1F(Form("h_time_event%i",ih),"h_time_event",50,-2.0,2.0));
    h_time_diff.push_back(new TH1F(Form("h_time_diff_%i",ih),"h_time_diff",100,-100,100));
    f_spe.push_back(new TF1(Form("f_spe_%i",ih),fmultigaus,0,50,7));
  }
  h_charge_muontrigs = new TH2F("h_charge_muontrigs","h_charge_muontrigs",100,0,10000,100,0,10000);
  h_time_muontrigs = new TH1F("h_time_muontrigs","h_time_muontrigs",100,-5,5);
  h_event_time = new TH1F("h_event_time","h_event_time",300,150,250);
  h_charge_total = new TH1F("h_charge_total","h_charge_total",200,-20,50000);
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
      h_mcpmt_npevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),mcpmt->GetMCPhotonCount()/(double)nentries);
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
      double charge = 0.;
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
      RAT::DS::PMT *pmt = ev->GetPMTWithID(20); //ref ring tube
      if(pmt!=NULL) {
        ring_charge = pmt->GetCharge();
        ring_time = pmt->GetTime();
        double ring_dist = (*pos_pmts[20] - *target_pos).Mag();
        double ring_tof = ring_dist/cspeed;
        ring_timeres = ring_time - ring_tof;
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
        double bottom_dist = (*pos_pmts[7] - *target_pos).Mag();
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

      double event_time = -9999;
      std::vector<double> ringPMTTimes;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        charge = ev->GetPMT(ipmt)->GetCharge();
        qtotal += charge;
        if(pmttype==1) {
          double dist = (*pos_pmts[pmtid] - *target_pos).Mag();
          double tof = dist/cspeed;
          if(ev->GetPMT(ipmt)->GetTime()>-9000) ringPMTTimes.push_back(ev->GetPMT(ipmt)->GetTime() - tof - pmttime_delay[pmtid]);
        }
      }
      //Sort in ascending order
      std::sort(ringPMTTimes.begin(), ringPMTTimes.end());

      //Cuts for cherenkov imaging
      // if(ringPMTTimes.size()<3) continue; //More than 3 hits
      // if(bottommuon_charge<2000.0 || topmuon_charge<2000.0) continue; //Charge cut
      // if(bottommuon_time - topmuon_time < 0.4 || bottommuon_time - topmuon_time > 0.8 ) continue; //Time cut
      // if(panel_charge[0]<2000.) continue; //Veto cut
      // if(panel_charge[1]>400. || panel_charge[2]>400. || panel_charge[3]>400.) continue; //Veto cut

      //Cuts for SPE
      // if(panel_charge[0]>200 || panel_charge[1]>200 || panel_charge[2]>200 || panel_charge[3]>200) continue;
      // if(ring_time < 180) continue;

      //Calculate event time
      event_time = TMath::KOrdStat((int)ringPMTTimes.size(), &ringPMTTimes[0], 1);
      // event_time = (ringPMTTimes[0] + ringPMTTimes[1] + ringPMTTimes[2])/3.;
      // event_time = ring_time;
      h_event_time->Fill(event_time);


      //Fill Histograms
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        double dist = (*pos_pmts[pmtid] - *target_pos).Mag();
        double tof = dist/cspeed;
        double pmttime = ev->GetPMT(ipmt)->GetTime();
        charge = ev->GetPMT(ipmt)->GetCharge();
        //if(pmttime < 170) continue;
        //if(charge > 200) continue;
        double timeres = pmttime - tof;
        h_time[pmtid]->Fill(timeres - ring_timeres);
        // std::cout<<" ToF "<<pmtid<<": "<<tof<<" "<<dist<<std::endl;
        //h_time_bottom[pmtid]->Fill(timeres - bottommuon_time - pmttime_delay[pmtid]);
        //h_time_event[pmtid]->Fill(timeres - event_time - pmttime_delay[pmtid]);
        if(pmtid != 13 && pmtid != 22 && pmtid != 19 && pmtid != 16) continue;
        h_time_bottom[pmtidtopos[pmtid]]->Fill(timeres - bottommuon_time - pmttime_delay[pmtid]);
        h_time_event[pmtidtopos[pmtid]]->Fill(timeres - event_time - pmttime_delay[pmtid]);
        h_time_trigger[pmtid]->Fill(timeres - trigger_time - pmttime_delay[pmtid]);
        //Event level averaged
        if(timeres>-900) {
          if (pmtInfo->GetType(pmtid)==1){
            h_pmt_chargevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),charge);
            h_pmt_npevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),charge/spe[pmtid]);
            h_pmt_timevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),timeres - event_time - pmttime_delay[pmtid]);
            h_pmt_countvspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y());
          }
        }
        if(pmtInfo->GetType(pmtid)==1){
          qtotal_smallpmts += charge;
        }
        h_charge[pmtid]->Fill(charge);
        double charge_res = charge/qScintCorr[pmtid]/20; //20cm (?)
        h_charge_res[pmtid]->Fill(charge_res);
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
    } //end daq event loop

  } //end ds entry loop

}

void ExtractSPE(){

  double params[6];
  for(int ipmt=0; ipmt<npmts; ipmt++){

//    if(pmtInfo->GetType(ipmt)!=1 && pmtInfo->GetType(ipmt)!=2) continue;

    f_spe[ipmt]->SetParameters(1.,0.,10.,1.,50.,50.,1.);
    f_spe[ipmt]->SetParLimits(0.,0.,9999999.); //Noise norm
    f_spe[ipmt]->SetParLimits(1.,-60.,60.); //Noise mean
    f_spe[ipmt]->SetParLimits(2.,0.,100.); //Noise sigma
    f_spe[ipmt]->SetParLimits(3.,0.,9999999.); //SPE norm
    f_spe[ipmt]->SetParLimits(4.,0.,200.); //SPE mean
    f_spe[ipmt]->SetParLimits(5.,0.,100.); //SPE sigma
    f_spe[ipmt]->SetParLimits(6.,0.,9999999.); //2PE norm

    h_charge[ipmt]->Fit(Form("f_spe_%i",ipmt),"Q","Q",-60,500);
    f_spe[ipmt]->GetParameters(params);
    f_spe_mean.push_back(params[4]);
    f_spe_sigma.push_back(params[5]);

    std::cout<<" SPE "<<ipmt<<":"<<std::endl;
    std::cout<<"  |-> Noise norm: "<<params[0]<<std::endl;
    std::cout<<"  |-> Noise mean: "<<params[1]<<std::endl;
    std::cout<<"  |-> Noise sigma: "<<params[2]<<std::endl;
    std::cout<<"  |-> SPE norm: "<<params[3]<<std::endl;
    std::cout<<"  |-> SPE mean: "<<params[4]<<std::endl;
    std::cout<<"  |-> SPE sigma: "<<params[5]<<std::endl;
    std::cout<<"  |-> 2PE norm: "<<params[6]<<std::endl;

  }

  std::cout<<" RATDB TABLE "<<":"<<std::endl;
  std::cout<<"{"<<std::endl;
  std::cout<<"  name: \"PMTGAUSCHARGE\", \n  valid_begin: [0,0], \n  valid_end: [0,0],"<<std::endl;
  std::cout<<"  gaus_mean: ["<<f_spe_mean.at(0);
  for(int ipmt=1; ipmt<f_spe_mean.size(); ipmt++) std::cout<<", "<<f_spe_mean.at(ipmt);
  std::cout<<"],"<<std::endl;
  std::cout<<"  gaus_sigma: ["<<f_spe_sigma.at(0);
  for(int ipmt=1; ipmt<f_spe_sigma.size(); ipmt++) std::cout<<", "<<f_spe_sigma.at(ipmt);
  std::cout<<"],"<<std::endl;
  std::cout<<"}"<<std::endl;


}




//Draw histrograms
void DrawHistos(){

  TCanvas *c_chi2 = new TCanvas("c_chi2","c_chi2",400,400);
  h_chi2->Draw();

  //General plots
  TCanvas *c_ring_event = new TCanvas("c_ring_event","Ring Event",900,900);
  c_ring_event->Divide(2,2);
  c_ring_event->cd(1);
  //h_pmt_chargevspos->Draw("colz text"); //Average charge per pmt vs position
  h_pmt_npevspos->Draw("colz text"); //Average charge per pmt vs position
  c_ring_event->cd(2);
  h_charge_total->Draw(); //Total charge in event
  c_ring_event->cd(3);
  h_pmt_timevspos->Draw("colz text"); //Average charge per pmt vs position
  c_ring_event->cd(4);
  h_event_time->Draw();
  std::cout<<" Total Q: "<<h_charge_total->Integral(20,200)<<std::endl;
  TCanvas *c_charge[5];
  c_charge[0] = new TCanvas("c_charge_0","Charge Ring Tubes",1200,900);
  c_charge[1] = new TCanvas("c_charge_1","Charge Light Tubes",900,600);
  c_charge[2] = new TCanvas("c_charge_2","Charge Muon Tags",600,300);
  c_charge[3] = new TCanvas("c_charge_3","Charge Trigger",300,300);
  c_charge[4] = new TCanvas("c_charge_4","Charge Vetos",600,600);
  c_charge[0]->Divide(4,3);
  c_charge[1]->Divide(3,2);
  c_charge[2]->Divide(2,2);
  c_charge[3]->Divide(1,1);
  c_charge[4]->Divide(2,2);
  int cc0 = 0, cc1 = 0, cc2 = 0, cc3 = 0, cc4 = 0;
  TCanvas *c_time[5];
  c_time[0] = new TCanvas("c_time_0","Time Ring Tubes",1200,900);
  c_time[1] = new TCanvas("c_time_1","Time Light Tubes",900,600);
  c_time[2] = new TCanvas("c_time_2","Time Muon Tags",600,300);
  c_time[3] = new TCanvas("c_time_3","Time Trigger",300,300);
  c_time[4] = new TCanvas("c_time_4","Time Vetos",600,600);
  c_time[0]->Divide(4,3);
  c_time[1]->Divide(3,2);
  c_time[2]->Divide(2,2);
  c_time[3]->Divide(1,1);
  c_time[4]->Divide(2,2);

  for(int pmtid = 0; pmtid<npmts; pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype==1){//Ring tubes
      c_charge[0]->cd(++cc0)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      f_spe[pmtid]->Draw("L same");
      c_time[0]->cd(cc0);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==2){ //Light tubes
      c_charge[1]->cd(++cc1)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      f_spe[pmtid]->Draw("L same");
      c_time[1]->cd(cc1);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==3){ //Muon tags
      c_charge[2]->cd(++cc2)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      c_time[2]->cd(cc2);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==0){ //Trigger tube
      c_charge[3]->cd(++cc3)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      c_time[3]->cd(cc3);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    } else if(pmttype==4){ //Muon panels
      c_charge[4]->cd(++cc4)->SetLogy();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw();
      c_time[4]->cd(cc4);
      h_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time[pmtid]->Draw();
    }
  }

  c_charge[2]->cd(3);
  h_charge_muontrigs->Draw("colz");
  c_charge[2]->cd(4);
  h_time_muontrigs->Draw();


  TCanvas *c_time_bottom = new TCanvas("c_time_bottom","c_time_bottom",900,900);
  h_time_bottom[0]->GetYaxis()->SetRangeUser(0,100);
  h_time_bottom[0]->SetLineColor(pmtidtocolor[12]);
  h_time_bottom[0]->Draw();
  h_time_bottom[1]->SetLineColor(pmtidtocolor[13]);
  h_time_bottom[1]->Draw("sames");
  h_time_bottom[2]->SetLineColor(pmtidtocolor[14]);
  h_time_bottom[2]->Draw("sames");
  // for(int pmtid = 0; pmtid<npmts; pmtid++){
  //   h_time_bottom[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
  //   h_time_bottom[pmtid]->Draw("sames");
  // }

  TCanvas *c_time_event = new TCanvas("c_time_event","c_time_event",900,900);
  h_time_event[0]->GetYaxis()->SetRangeUser(0,100);
  h_time_event[0]->SetLineColor(pmtidtocolor[12]);
  h_time_event[0]->Draw();
  h_time_event[1]->SetLineColor(pmtidtocolor[13]);
  h_time_event[1]->Draw("sames");
  h_time_event[2]->SetLineColor(pmtidtocolor[14]);
  h_time_event[2]->Draw("sames");
  // h_time_event[0]->Draw("");
  // for(int pmtid = 0; pmtid<npmts; pmtid++){
  //   //if(pmtid !=13 && pmtid != 16 && pmtid !=19 && pmtid !=22) continue;
  //   h_time_event[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
  //   h_time_event[pmtid]->Draw("sames");
  // }

  TCanvas *c_time_trigger = new TCanvas("c_time_trigger","c_time_trigger",900,900);
  for(int pmtid = 0; pmtid<npmts; pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype!=1) continue;
    h_time_trigger[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    h_time_trigger[pmtid]->Draw("sames");
  }

  //MCTruth
  bool firstdrawn0=false, firstdrawn1=false, firstdrawn2=false;
  TCanvas *c_mc = new TCanvas("c_mc","c_mc",900,1000);
  c_mc->Divide(3,3);
  for(int pmtid = 0; pmtid<npmts; pmtid++){
    for(int pmtid = 0; pmtid<npmts; pmtid++){
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
  for(int ipmt=0; ipmt<npmts;ipmt++){
    h_mcpmt_npe[ipmt]->Write();
    h_mcpmt_charge[ipmt]->Write();
    h_mcpmt_fetime[ipmt]->Write();
    h_charge[ipmt]->Write();
    h_time[ipmt]->Write();
    h_time_diff[ipmt]->Write();
    h_charge_vs_trigq[ipmt]->Write();
  }
  h_chi2->Write();
  h_charge_total->Write();
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
    if(std::string(argv[i]) == "-mc") {gInputFileMC = argv[++i]; exist_inputfile=true;}
    if(std::string(argv[i]) == "-dt") {gInputFileDT.push_back(argv[++i]); exist_inputfile=true;}
    if(std::string(argv[i]) == "-o")  {gOutFile = argv[++i];}
    if(std::string(argv[i]) == "-m")  {gTargetMaterial = argv[++i];}
  }

  if(!exist_inputfile){
    std::cerr<<" Usage: ./DrawChargeAndTime.exe -mc MCFILE [Optional: -dt DATAFILE]"<<std::endl;
    exit(0);
  }
}
