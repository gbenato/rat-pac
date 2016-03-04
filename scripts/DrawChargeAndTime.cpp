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

#define USERROOTLOOP true
#define MCPHOTONLOOP false
#define NLOGENTRIES 10
#define NBINS 20
#define MIN -500.0
#define MAX -300.0

//Constants
double cspeed = 300/1.4; // (mm/ns)/rindex
TVector3* target_pos = new TVector3(-400.,-400.,-200.);

//Methods
char *gInputFileMC = NULL;
std::vector<char*> gInputFileDT;
char *gOutFile = NULL;
char *gTargetMaterial = "LAB";
void ParseArgs(int argc, char **argv);
void GetPMTInfo();
void GetDBPlots();
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


//// Histograms
//Event level
TH1F* h_event_time;
TH1F* h_charge_total; //Total charge in the event
TH2F* h_pmt_qresvspos; //PMT charge vs PMT position
TH2F* h_pmt_timevspos; //PMT charge vs PMT position
TH1F* h_chi2; //Chi2 cher vs scint
TH2F* h_mcpmt_npevspos; //MC NPE vs PMT position

//Hit level (PMT)
std::vector<double> q_xmin;
std::vector<double> q_xmax;
std::vector<TF1*> f_spe;
std::vector<TH1F*> h_charge; //Measured charge
std::vector<TH1F*> h_charge_res; //Measured charge geometry corrected
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_res; //Measured time
std::vector<TH1F*> h_time_bottom; //Measured time
std::vector<TH1F*> h_time_trigger; //Measured time
std::vector<TH1F*> h_time_ring; //Measured time
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
  GetDBPlots();
  GetHistos();
  NormalizeHistos();
  ExtractSPE();

  DrawHistos();
  if(gOutFile){
    PrintHistos(gOutFile);
  }

  new TBrowser;
  dummy.Run();
  return 1;

}

void GetDBPlots(){

  std::cout<<" Get DB scint correction plots for "<<gTargetMaterial<<std::endl;

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

}

void GetPMTInfo(){

  //Init pmt positions
  RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
  TTree *runT = dsreader->GetRunT();
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

  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kRed, kOrange, kBlue, kBlue, kOrange, kRed, kRed, kOrange, kBlue, 1}; //By Position
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue-2, kOrange-2, kRed-2, kRed-1, kOrange-1, kBlue-1, kBlue, kOrange, kRed, kRed+1, kOrange+1, kBlue+1, 1}; //By Position
  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kBlue, kBlue, kBlue, kOrange, kOrange, kOrange, kOrange, kRed, kRed, kRed, kRed, 1}; //By Digitizer group
  //  Color_t mycolors[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, kBlue, kOrange, kRed, kCyan, kBlack, kGray, kGreen, kTeal, kAzure, kViolet, kPink, kYellow, 1}; //Individual
  int mypmtpos[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 1, 0, 0, 1, 2, 2, 1, 0, 3};
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
  // pmttime_delay = std::vector<double>(25,-20.);
  // pmttime_delay = std::vector<double>(25,0.);

  //Plot axis limits
  double myqxmins[] = {-100,-100,-100,-100,-100,-100, -100,-100, -100,-100,-100,-100, -100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100, -100};
  q_xmin.insert(q_xmin.begin(), myqxmins, myqxmins + npmts );
  double myqxmaxs[] = {3000,3000,3000,3000,3000,3000, 3000,3000, 3000,3000,3000,3000, 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000, 3000};
  q_xmax.insert(q_xmax.begin(), myqxmaxs, myqxmaxs + npmts );

}

void GetHistos(){

  //MC
  std::cout<<" GetMCPDFs "<<std::endl;

  RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
  //Init histos
  //  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",27,-27*15,27*15,27,-27*15,27*15);
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",NBINS,MIN,MAX,NBINS,MIN,MAX);
  h_pmt_qresvspos = new TH2F("h_pmt_qresvspos","h_pmt_qresvspos",NBINS,MIN,MAX,NBINS,MIN,MAX);
  h_pmt_timevspos = new TH2F("h_pmt_timevspos","h_pmt_timevspos",NBINS,MIN,MAX,NBINS,MIN,MAX);
  for(int ih=0; ih<npmts; ih++){
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",200,0,200));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_time.push_back(new TH1F(Form("h_mcpmt_time_%i",ih),"h_mcpmt_time",500,0,100));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",500,0,100));
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),"h_charge",500,q_xmin[ih],q_xmax[ih]));
    h_charge_res.push_back(new TH1F(Form("h_charge_res_%i",ih),"h_charge_res",50,0,50));
    h_charge_vs_trigq.push_back(new TH2F(Form("h_charge_vs_trigq_%i",ih),"h_charge_vs_trigq",200,0,100,200,0,100));
    h_time.push_back(new TH1F(Form("h_time_%i",ih),"h_time",400,200,240));
    h_time_res.push_back(new TH1F(Form("h_time_res_%i",ih),"h_time_res",450,190,240));
    h_time_bottom.push_back(new TH1F(Form("h_time_bottom_%i",ih),"h_time_bottom",150,-1.0,14));
    h_time_trigger.push_back(new TH1F(Form("h_time_trigger%i",ih),"h_time_trigger",200,-10.0,10.0));
    h_time_ring.push_back(new TH1F(Form("h_time_ring%i",ih),"h_time_ring",200,-10.0,10.0));
    h_time_event.push_back(new TH1F(Form("h_time_event%i",ih),"h_time_event",150,-1.0,14.0));
    h_time_diff.push_back(new TH1F(Form("h_time_diff_%i",ih),"h_time_diff",100,-100,100));
    f_spe.push_back(new TF1(Form("f_spe_%i",ih),fmultigaus,0,50,7));
  }
  h_charge_muontrigs = new TH2F("h_charge_muontrigs","h_charge_muontrigs",200,0,8000,200,0,8000);
  h_time_muontrigs = new TH1F("h_time_muontrigs","h_time_muontrigs",100,-5,5);
  h_event_time = new TH1F("h_event_time","h_event_time",200,200,250);
  h_charge_total = new TH1F("h_charge_total","h_charge_total",200,-20,500);
  h_chi2 = new TH1F("h_chi2","h_chi2",50,0,200);

  //Fill histos with loop
  if(USERROOTLOOP){

    RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
    TTree *tree = dsreader->GetT();
    RAT::DS::Root *rds;
    RAT::DS::MC *mc;
    int nentries = tree->GetEntries();
    std::cout<<" Number of entries: "<<nentries<<std::endl;
    for(int ientry=0; ientry<nentries;++ientry){

      if(nentries>NLOGENTRIES && ientry%(nentries/NLOGENTRIES) == 0) std::cout<<" Entry "<<ientry<<std::endl;
      //    if(ientry>1000000) break;

      tree->GetEntry(ientry);

      rds = dsreader->GetEvent(ientry);
      mc = rds->GetMC();
      //If no PMT continue to save time (in theory)

      //      if(mc->GetMCPMTCount()==0) continue;

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
        h_mcpmt_time[pmtid]->Fill(mcpmt->GetTime());
        h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetFrontEndTime());

        if(MCPHOTONLOOP){
          for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++){
            h_mcpmt_charge[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetCharge());
            h_mcpmt_time[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetHitTime());
            h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
          }
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

        //First PMT loop:
        // - Calculate event time
        // - Look for trigger PMT and save charge
        // - Calculate total charge
        // - Get charge at top muon tag
        double event_time = -9999;
        std::vector<double> ringPMTTimes;
        for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
          int pmtid = ev->GetPMT(ipmt)->GetID();
          int pmttype = pmtInfo->GetType(pmtid);
          if(pmtid==1){
            charge_trig = ev->GetPMT(ipmt)->GetCharge();
          }
          // if(pmtid==6){
          //   topmuon_charge = ev->GetPMT(ipmt)->GetCharge();
          // }
          // else if(pmtid==7){
          //   bottommuon_charge = ev->GetPMT(ipmt)->GetCharge();
          // }
          charge = ev->GetPMT(ipmt)->GetCharge();
          qtotal += charge;
          if(pmttype==1) {
            double dist = (*pos_pmts[pmtid] - *target_pos).Mag();
            double tof = dist/cspeed;
            if(ev->GetPMT(ipmt)->GetTime()>-9000) ringPMTTimes.push_back(ev->GetPMT(ipmt)->GetTime() - tof - pmttime_delay[pmtid]);

            qtotal_smallpmts += charge;
          }
        }

        //Calculate event time
        event_time = TMath::KOrdStat((int)ringPMTTimes.size(), &ringPMTTimes[0], 1);
        h_event_time->Fill(event_time);

        //Get muon tags, source trigger and muon panels charge and time
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
        RAT::DS::PMT *pmt = ev->GetPMTWithID(20);
        if(pmt!=NULL) {
          ring_charge = pmt->GetCharge();
          ring_time = pmt->GetTime();
          double ring_dist = (*pos_pmts[20] - *target_pos).Mag();
          double ring_tof = ring_dist/cspeed;
          ring_timeres = ring_time - ring_tof;
        }
        pmt = ev->GetPMTWithID(6);
        if(pmt!=NULL) {
          topmuon_charge = pmt->GetCharge();
          topmuon_time = pmt->GetTime();
        }
        pmt = ev->GetPMTWithID(7);
        if(pmt!=NULL) {
          bottommuon_charge = pmt->GetCharge();
          bottommuon_time = pmt->GetTime();
          double bottom_dist = (*pos_pmts[7] - *target_pos).Mag();
          double bottom_tof = bottom_dist/cspeed;
          bottommuon_timeres = bottommuon_time - bottom_tof;
        }

//        std::cout<<" muon_charges "<<bottommuon_charge<<" "<<topmuon_charge<<std::endl;

        pmt = ev->GetPMTWithID(8);
        if(pmt!=NULL) panel_charge[0] = pmt->GetCharge();
        pmt = ev->GetPMTWithID(9);
        if(pmt!=NULL) panel_charge[1] = pmt->GetCharge();
        pmt = ev->GetPMTWithID(10);
        if(pmt!=NULL) panel_charge[2] = pmt->GetCharge();
        pmt = ev->GetPMTWithID(11);
        if(pmt!=NULL) panel_charge[3] = pmt->GetCharge();
        pmt = ev->GetPMTWithID(20);
        if(pmt!=NULL) trigger_time = pmt->GetTime();

        //Cuts for cherenkov imaging
        //if(bottommuon_charge<400.0 || topmuon_charge<400.0 || bottommuon_time - topmuon_time < 0.4 || bottommuon_time - topmuon_time > 1. ) continue;// || panel_charge[1]>200. || panel_charge[2]>200. || panel_charge[3]>200.) continue;
        //if(bottommuon_charge<400.0 || topmuon_charge<400.0 || bottommuon_time - topmuon_time < 0.2 || bottommuon_time - topmuon_time > 0.8 ) continue;// || panel_charge[1]>200. || panel_charge[2]>200. || panel_charge[3]>200.) continue;
        //Cuts for SPE
        // if(ring_time < -9000 || panel_charge[0]>200 || panel_charge[1]>200 || panel_charge[2]>200 || panel_charge[3]>200) continue;

        for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
          int pmtid = ev->GetPMT(ipmt)->GetID();
          double dist = (*pos_pmts[pmtid] - *target_pos).Mag();
          double tof = dist/cspeed;
          charge = ev->GetPMT(ipmt)->GetCharge();
          h_charge[pmtid]->Fill(charge);
          double charge_res = charge/qScintCorr[pmtid]/20; //20cm (?)
          h_charge_res[pmtid]->Fill(charge_res);
          h_pmt_qresvspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),charge_res);
          h_charge_vs_trigq[pmtid]->Fill(charge,charge_trig);
          double pmttime = ev->GetPMT(ipmt)->GetTime();
          h_time[pmtid]->Fill(pmttime);
          double timeres = pmttime - tof;
          // std::cout<<" ToF "<<pmtid<<": "<<tof<<" "<<dist<<std::endl;
          h_time_res[pmtid]->Fill(timeres);
          h_time_bottom[pmtidtopos[pmtid]]->Fill(timeres - bottommuon_time - pmttime_delay[pmtid]);
          h_time_trigger[pmtid]->Fill(timeres - trigger_time - pmttime_delay[pmtid]);
          h_time_ring[pmtid]->Fill(timeres - ring_timeres - pmttime_delay[pmtid]);
          h_time_event[pmtidtopos[pmtid]]->Fill(timeres - event_time - pmttime_delay[pmtid]);
          if(timeres>-900) h_pmt_timevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),timeres);
          //Compute chi2 for cher/scint
          if(ev->GetPMT(ipmt)->GetType()==1){
            chi2 += pow( (charge - qScintCorr[pmtid])/qScintCorrErr[pmtid], 2.);
//            std::cout<<" chi2 "<<ipmt<<" "<<pmtid<<" "<<chi2<<" "<<charge<<" "<<qScintCorr[pmtid]<<" "<<qScintCorrErr[pmtid]<<std::endl;
          }
        }
        h_charge_muontrigs->Fill(bottommuon_charge,topmuon_charge);
        h_time_muontrigs->Fill(bottommuon_time - topmuon_time);
        h_charge_total->Fill(qtotal_smallpmts);
        // if(qtotal!=0){
        //   h_charge_total->Fill(qtotal);
        // }
//        std::cout<<" TOTAL chi2 "<<chi2<<std::endl;
        h_chi2->Fill(chi2);
      } //end daq event loop

    } //end ds entry loop

  } else { //if !USERROOTLOOP

    //Fill histos with Tree::Draw method
    RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
    int nentries = dsreader->GetT()->GetEntries();
    std::cout<<" Number of entries: "<<nentries<<std::endl;
    TTree *T = dsreader->GetT();
    for(int ipmt=0; ipmt<npmts; ipmt++){
      //EV
      std::cout<<"   Loading charge... "<<std::endl;
      T->Draw(Form("ds.ev.pmt.charge>>h_charge_%i",ipmt),Form("ds.ev.pmt.id==%i",ipmt));
      if (MCPHOTONLOOP) {
        std::cout<<"   Loading time... "<<std::endl;
        T->Draw(Form("ds.ev.pmt.time>>h_time_%i",ipmt),Form("ds.ev.pmt.id==%i",ipmt));
        //MC
        std::cout<<"   Loading mc charge... "<<std::endl;
        T->Draw(Form("ds.mc.pmt.GetCharge()>>h_mcpmt_charge_%i",ipmt),Form("ds.mc.pmt.id==%i",ipmt));
        T->Draw(Form("ds.mc.pmt.GetMCPhotonCount()>>h_mcpmt_npe_%i",ipmt),Form("ds.ev.Nhits()>0 && ds.mc.pmt.id==%i",ipmt));
        T->Draw(Form("ds.mc.pmt.photon.at(0).frontEndTime>>h_mcpmt_fetime_%i",ipmt),Form("ds.ev.Nhits()>0 && ds.mc.pmt.id==%i",ipmt));
        T->Draw(Form("ds.mc.pmt.photon.at(0).frontEndTime - ds.ev.pmt.time >> h_time_diff_%i",ipmt),Form("ds.ev.Nhits()>0 && ds.mc.pmt.id==%i",ipmt));
      }
    }
  }

  // //REAL DATA
  // if(gInputFileDT.size()>0){
  //   std::cout<<" GetDataPDFs "<<std::endl;
  //
  //   TGraph* gpdf_dt;
  //   TH1F* hdata;
  //   TH1F* hscale;
  //   for(int ifile=0; ifile<gInputFileDT.size(); ifile++){
  //     std::cout<<" gInputFileDT "<<ifile<<" "<<gInputFileDT[ifile]<<std::endl;
  //     hdata = new TH1F(Form("hdata%i",ifile),"hdata",200,0,100);
  //     hscale = new TH1F(Form("hscale%i",ifile),"hscale",200,0,100); //count how many times we fill the bin and scale it, root is anoying...
  //     gpdf_dt = new TGraph(gInputFileDT[ifile],"%lg %lg",",");
  //     double x=0.;
  //     double y=0.;
  //     for(int ip=0; ip<gpdf_dt->GetN(); ip++){
  //       gpdf_dt->GetPoint(ip,x,y);
  //       hdata->Fill(x,y);
  //       hscale->Fill(x,1.);
  //     }
  //     hdata->Divide(hscale);
  //     //    hdata->Scale(1./hdata->Integral()); //Nomalize for shape only analysis
  //     h_dt_charge.push_back(hdata);
  //   }
  // }


}

void ExtractSPE(){

  double params[6];
  for(int ipmt=0; ipmt<npmts; ipmt++){

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

    std::cout<<" SPE "<<ipmt<<":"<<std::endl;
    std::cout<<"  |-> Noise norm: "<<params[0]<<std::endl;
    std::cout<<"  |-> Noise mean: "<<params[1]<<std::endl;
    std::cout<<"  |-> Noise sigma: "<<params[2]<<std::endl;
    std::cout<<"  |-> SPE norm: "<<params[3]<<std::endl;
    std::cout<<"  |-> SPE mean: "<<params[4]<<std::endl;
    std::cout<<"  |-> SPE sigma: "<<params[5]<<std::endl;
    std::cout<<"  |-> 2PE norm: "<<params[6]<<std::endl;

  }

}




//Draw histrograms
void DrawHistos(){

  TCanvas *c_chi2 = new TCanvas("c_chi2","c_chi2",400,400);
  h_chi2->Draw();

  //General plots
  bool firstdrawn0 = false;
  bool firstdrawn1 = false;
  bool firstdrawn2 = false;
  bool firstdrawn3 = false;
  bool firstdrawn4 = false;
  TCanvas *c_event = new TCanvas("c_event","c_event",900,900);
  c_event->Divide(2,2);
  c_event->cd(1);
  h_pmt_qresvspos->Draw("colz text"); //Average charge per pmt vs position
  c_event->cd(2);
  h_event_time->Draw();
  c_event->cd(3);
  h_charge_total->Draw(); //Total charge in event
  c_event->cd(4);
  h_pmt_timevspos->Draw("colz text"); //Average charge per pmt vs position
  std::cout<<" Total Q: "<<h_charge_total->Integral(20,200)<<std::endl;
  TCanvas *c_charge[5];
  c_charge[0] = new TCanvas("c_charge_0","c_charge_0",900,900);
  c_charge[1] = new TCanvas("c_charge_1","c_charge_1",900,900);
  c_charge[2] = new TCanvas("c_charge_2","c_charge_2",900,900);
  c_charge[3] = new TCanvas("c_charge_3","c_charge_3",900,900);
  c_charge[4] = new TCanvas("c_charge_4","c_charge_4",900,900);
  TCanvas *c_charge_ind[3];
  c_charge_ind[0] = new TCanvas("c_charge_ind_0","c_charge_ind_0",1500,900);
  c_charge_ind[0]->Divide(4,3);
  c_charge_ind[1] = new TCanvas("c_charge_ind_1","c_charge_ind_1",1500,900);
  c_charge_ind[1]->Divide(3,2);
  TCanvas *c_time[5];
  c_time[0] = new TCanvas("c_time_0","c_time_0",900,900);
  c_time[1] = new TCanvas("c_time_1","c_time_1",900,900);
  c_time[2] = new TCanvas("c_time_2","c_time_2",900,900);
  c_time[3] = new TCanvas("c_time_3","c_time_3",900,900);
  c_time[4] = new TCanvas("c_time_4","c_time_4",900,900);

  for(int pmtid = 0; pmtid<npmts; pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype==1){//Ring tubes
      const char *opt = firstdrawn0 ? "sames" : "";
      c_charge[0]->cd();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw(opt);
      c_charge_ind[0]->cd(pmtid-11);
      c_charge_ind[0]->cd(pmtid-11)->SetLogy();
      //      h_charge[pmtid]->GetXaxis()->SetRangeUser(-100,400);
      h_charge[pmtid]->Draw();
      f_spe[pmtid]->Draw("L same");
      c_time[0]->cd();
      h_time_res[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time_res[pmtid]->Draw(opt);
      firstdrawn0 = true;
    } else if(pmttype==2){ //Light tubes
      const char *opt = firstdrawn1 ? "sames" : "";
      c_charge[1]->cd();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw(opt);
      c_time[1]->cd();
      h_time_res[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time_res[pmtid]->Draw(opt);
      firstdrawn1 = true;
    } else if(pmttype==3){ //Muon tags
      const char *opt = firstdrawn2 ? "sames" : "";
      c_charge[2]->cd();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw(opt);
      c_time[2]->cd();
      h_time_res[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time_res[pmtid]->Draw(opt);
      firstdrawn2 = true;
    } else if(pmttype==0){ //Trigger tube
      const char *opt = firstdrawn3 ? "sames" : "";
      c_charge[3]->cd();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw(opt);
      c_time[3]->cd();
      h_time_res[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time_res[pmtid]->Draw(opt);
      firstdrawn3 = true;
    } else if(pmttype==4){ //Muon panels
      const char *opt = firstdrawn4 ? "sames" : "";
      c_charge[4]->cd();
      h_charge[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_charge[pmtid]->Draw(opt);
      c_time[4]->cd();
      h_time_res[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
      h_time_res[pmtid]->Draw(opt);
      firstdrawn4 = true;
    }
  }

  TCanvas *c_time_bottom = new TCanvas("c_time_bottom","c_time_bottom",900,900);
  h_time_bottom[0]->GetYaxis()->SetRangeUser(0,20);
  h_time_bottom[0]->SetLineColor(kRed);
  h_time_bottom[0]->Draw();
  h_time_bottom[1]->SetLineColor(kBlue);
  h_time_bottom[1]->Draw("sames");
  h_time_bottom[2]->SetLineColor(kOrange);
  h_time_bottom[2]->Draw("sames");
  // for(int pmtid = 0; pmtid<npmts; pmtid++){
  //   h_time_bottom[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
  //   h_time_bottom[pmtid]->Draw("sames");
  // }

  TCanvas *c_time_trigger = new TCanvas("c_time_trigger","c_time_trigger",900,900);
  for(int pmtid = 0; pmtid<npmts; pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype!=1) continue;
    h_time_trigger[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    h_time_trigger[pmtid]->Draw("sames");
  }

  TCanvas *c_time_ring = new TCanvas("c_time_ring","c_time_ring",900,900);
  for(int pmtid = 0; pmtid<npmts; pmtid++){
    int pmttype = pmtInfo->GetType(pmtid);
    if(pmttype!=1) continue;
    h_time_ring[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
    h_time_ring[pmtid]->Draw("sames");
  }

  TCanvas *c_time_event = new TCanvas("c_time_event","c_time_event",900,900);
  h_time_event[0]->GetYaxis()->SetRangeUser(0,20);
  h_time_event[0]->SetLineColor(kRed);
  h_time_event[0]->Draw();
  h_time_event[1]->SetLineColor(kBlue);
  h_time_event[1]->Draw("sames");
  h_time_event[2]->SetLineColor(kOrange);
  h_time_event[2]->Draw("sames");

  TCanvas *c_muontrigs = new TCanvas("c_muontrigs","c_muontrigs",1200,600);
  c_muontrigs->Divide(2,1);
  c_muontrigs->cd(1);
  h_charge_muontrigs->Draw("colz");
  c_muontrigs->cd(2);
  h_time_muontrigs->Draw();

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
        h_mcpmt_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_time[pmtid]->Draw(opt);
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
        h_mcpmt_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_time[pmtid]->Draw(opt);
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
        h_mcpmt_time[pmtid]->SetLineColor(pmtidtocolor[pmtid]);
        h_mcpmt_time[pmtid]->Draw(opt);
        firstdrawn2 = true;
      }
    }
  }

  TCanvas *c_pmtmaps = new TCanvas("c_pmtmaps","c_pmtmaps",300,300);
  // c_pmtmaps->Divide(2,1);
  // c_pmtmaps->cd(1);
  h_mcpmt_npevspos->Draw("colz text");
  // c_pmtmaps->cd(2);
  // h_pmt_qresvspos->Draw("colz text");

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
    h_time_res[ipmt]->Write();
    h_time_diff[ipmt]->Write();
    h_charge_vs_trigq[ipmt]->Write();
  }
  h_chi2->Write();
  h_charge_total->Write();
  fout->Close();


}

void NormalizeHistos(){

  // for(int ipmt=0; ipmt<npmts;ipmt++){
  //   double norm = h_charge[ipmt]->Integral(2,100);
  //   h_charge[ipmt]->Scale(1./norm);
  //   // std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
  //   norm = h_mcpmt_charge[ipmt]->Integral(2,100);
  //   h_mcpmt_charge[ipmt]->Scale(1./norm);
  //   // std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
  // }

  // if(gInputFileDT.size()>0){
  //   double norm = h_dt_charge[0]->Integral(2,100);
  //   h_dt_charge[0]->Scale(1./norm);
  // }

  // for(int ipmt=0;ipmt<npmts;ipmt++){
  //   h_pmt_qresvspos->Fill(pos_pmts[ipmt]->X(),pos_pmts[ipmt]->Y(),h_charge[ipmt]->GetMean());
  // }

  //Per event normalization
  double norm = (double)h_pmt_qresvspos->GetEntries();
  // std::cout<<" h_pmt_qresvspos "<<norm<<std::endl;
  h_pmt_qresvspos->Scale(npmts/norm);

  norm = (double)h_pmt_timevspos->GetEntries();
  // std::cout<<" h_pmt_qresvspos "<<norm<<std::endl;
  h_pmt_timevspos->Scale(npmts/norm);

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
