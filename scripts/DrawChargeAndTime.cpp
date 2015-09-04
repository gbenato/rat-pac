#include<iostream>
#include<fstream>

#include<TH1F.h>
#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<TBrowser.h>
#include<TGraph.h>

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include <RAT/DB.hh>
#include <RAT/DS/RunStore.hh>

#define NPMTs 81
#define USERROOTLOOP true
#define DRAWONLYCHARGE false
#define NLOGENTRIES 20

//Methods
char *gInputFileMC = NULL;
std::vector<char*> gInputFileDT;
char *gOutFile = NULL;
void ParseArgs(int argc, char **argv);
void GetHistos();
void DrawHistos();
void PrintHistos(char*);
void NormalizeHistos();

//Global variables
std::vector < std::vector < double > > pos_pmts;


//Histograms
TH2F* h_mcpmt_npevspos;
std::vector<TH1F*> h_mcpmt_npe; //Number of PE
std::vector<TH1F*> h_mcpmt_charge; //MC charge
std::vector<TH1F*> h_mcpmt_fetime; //MC FE time
std::vector<TH1F*> h_charge; //Measured charge
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_diff; //Time diff between EV-MC
std::vector<TH2F*> h_charge_scat; //PMT charge vs trigger charge
TH1F* h_charge_total;
TH1F* h_mcpmt_fetime_total;

//Real data
std::vector<TH1F*> h_dt_charge; //Measured charge

int main(int argc, char **argv){

  //Init********
  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy("App", &appargc, appargv);
  ParseArgs(argc, argv);
  //************

  GetHistos();
  if(gInputFileDT.size()>0)
    NormalizeHistos();

  DrawHistos();
  if(gOutFile)
    PrintHistos(gOutFile);

  new TBrowser;
  dummy.Run();
  return 1;

}


void GetHistos(){

  //MC
  std::cout<<" GetMCPDFs "<<std::endl;
  //Init pmt positions

  RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
  TTree *runT = dsreader->GetRunT();
  RAT::DS::Run *run = 0;
  runT->SetBranchAddress("run",&run);
  runT->GetEntry(0);
  RAT::DS::PMTInfo *pmtInfo = run->GetPMTInfo();

  for (size_t ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    std::vector <double> pmt_temp;
    pmt_temp.push_back(pmtInfo->GetPosition(ipmt)[0]);
    pmt_temp.push_back(pmtInfo->GetPosition(ipmt)[1]);
    pmt_temp.push_back(pmtInfo->GetPosition(ipmt)[2]);
    pos_pmts.push_back(pmt_temp);
  }

  // pos_pmts[0][0] = -30.; pos_pmts[0][1] = 0.; pos_pmts[0][2] = 115.;
  // pos_pmts[1][0] = -20.; pos_pmts[1][1] = 0.; pos_pmts[1][2] = 115.;
  // pos_pmts[2][0] = -10.; pos_pmts[2][1] = 0.; pos_pmts[2][2] = 115.;
  // pos_pmts[3][0] = 0.; pos_pmts[3][1] = 0.; pos_pmts[3][2] = 115.;
  // pos_pmts[4][0] = 10.; pos_pmts[4][1] = 0.; pos_pmts[4][2] = 115.;
  // pos_pmts[5][0] = 20.; pos_pmts[5][1] = 0.; pos_pmts[5][2] = 115.;
  // pos_pmts[6][0] = 30.; pos_pmts[6][1] = 0.; pos_pmts[6][2] = 115.;
  // pos_pmts[7][0] = 0.; pos_pmts[7][1] = -30.; pos_pmts[7][2] = 115.;
  // pos_pmts[8][0] = 0.; pos_pmts[8][1] = -20.; pos_pmts[8][2] = 115.;
  // pos_pmts[9][0] = 0.; pos_pmts[9][1] = -10.; pos_pmts[9][2] = 115.;
  // pos_pmts[10][0] = 0.; pos_pmts[10][1] = 10.; pos_pmts[10][2] = 115.;
  // pos_pmts[11][0] = 0.; pos_pmts[11][1] = 20.; pos_pmts[11][2] = 115.;
  // pos_pmts[12][0] = 0.; pos_pmts[12][1] = 30.; pos_pmts[12][2] = 115.;
  // pos_pmts[13][0] = 20.; pos_pmts[13][1] = 20.; pos_pmts[13][2] = 115.;
  // pos_pmts[14][0] = -20.; pos_pmts[14][1] = -20.; pos_pmts[14][2] = 115.;


  //Init histos
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",9,-9*15,9*15,9,-9*15,9*15);
  for(int ih=0; ih<NPMTs; ih++){
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",10,0,10));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",1000,0,200));
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),"h_charge",200,0,100));
    h_charge_scat.push_back(new TH2F(Form("h_charge_scat_%i",ih),"h_charge_scat",200,0,100,200,0,100));
    h_time.push_back(new TH1F(Form("h_time_%i",ih),"h_time",1000,0,200));
    h_time_diff.push_back(new TH1F(Form("h_time_diff_%i",ih),"h_time_diff",100,-100,100));
  }
  h_charge_total = new TH1F("h_charge_total","h_charge_total",100,0,150);
  h_mcpmt_fetime_total  = new TH1F("h_mcpmt_fetime_total","h_mcpmt_fetime_total",1000,0,150);

  //Fill histos with loop
  if(USERROOTLOOP){

    RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
    TTree *tree = dsreader->GetT();
    int nentries = tree->GetEntries();
    std::cout<<" Number of entries: "<<nentries<<std::endl;
    for(int ientry=0; ientry<nentries;++ientry){

      if(nentries>NLOGENTRIES && ientry%(nentries/NLOGENTRIES) == 0) std::cout<<" Entry "<<ientry<<std::endl;
      //    if(ientry>1000000) break;
      tree->GetEntry(ientry);
      RAT::DS::Root *rds = dsreader->GetEvent(ientry);
      RAT::DS::MC *mc = rds->GetMC();
      //If no PMT continue to save time (in theory)
      if(mc->GetMCPMTCount()==0) continue;
      if(mc->GetMCPMTCount() > NPMTs){
	std::cout<<" Number of PMTs larger than expected: change the parameter NPMTs! "<<std::endl;
	std::cout<<" Exiting now.... "<<std::endl;
	exit(0);
      }
      //MC**************
      //MCPMT loop
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
      	RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      	int pmtid = mcpmt->GetID();
	//count PE
	h_mcpmt_npe[pmtid]->Fill(mcpmt->GetMCPhotonCount());
	h_mcpmt_npevspos->Fill(pos_pmts[pmtid][0],pos_pmts[pmtid][1],mcpmt->GetMCPhotonCount()/(double)nentries);

      	for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++){
      	  h_mcpmt_charge[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetCharge());
      	  h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
	  h_mcpmt_fetime_total->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
      	}
      } //end MCPMT loop

      //DAQ EVENTS*****
      //Event loop
      int nevents = rds->GetEVCount();
      for(int ievt=0; ievt<nevents; ievt++){
	RAT::DS::EV *ev = rds->GetEV(ievt);
	double totalcharge = 0.;
	double charge = 0.;
	double charge_trig = 0.;
	//Look for trigger PMT and save charge
	for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
	  int trigid = ev->GetPMT(ipmt)->GetID();
	  if(trigid==1){
	    charge_trig = ev->GetPMT(ipmt)->GetCharge();
	    break;
	  }
	}
	for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
	  int pmtid = ev->GetPMT(ipmt)->GetID();
	  charge = ev->GetPMT(ipmt)->GetCharge();
	  h_charge[pmtid]->Fill(charge);
	  h_charge_scat[pmtid]->Fill(charge,charge_trig); //charge vs trigger charge
	  h_time[pmtid]->Fill(ev->GetPMT(ipmt)->GetTime());
	  totalcharge += charge;
	}
	if(totalcharge!=0){
	  h_charge_total->Fill(totalcharge);
	}
      } //end daq event loop
    } //end ds entry loop

  } else { //if !USERROOTLOOP

    //Fill histos with Tree::Draw method
    RAT::DSReader *dsreader = new RAT::DSReader(gInputFileMC);
    int nentries = dsreader->GetT()->GetEntries();
    std::cout<<" Number of entries: "<<nentries<<std::endl;
    TTree *T = dsreader->GetT();
    for(int ipmt=0; ipmt<NPMTs; ipmt++){
      //EV
      std::cout<<"   Loading charge... "<<std::endl;
      T->Draw(Form("ds.ev.pmt.charge>>h_charge_%i",ipmt),Form("ds.ev.pmt.id==%i",ipmt));
      if (!DRAWONLYCHARGE) {
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

  //REAL DATA
  if(gInputFileDT.size()>0){
    std::cout<<" GetDataPDFs "<<std::endl;

    TGraph* gpdf_dt;
    TH1F* hdata;
    TH1F* hscale;
    for(int ifile=0; ifile<gInputFileDT.size(); ifile++){
      std::cout<<" gInputFileDT "<<ifile<<" "<<gInputFileDT[ifile]<<std::endl;
      hdata = new TH1F(Form("hdata%i",ifile),"hdata",200,0,100);
      hscale = new TH1F(Form("hscale%i",ifile),"hscale",200,0,100); //count how many times we fill the bin and scale it, root is anoying...
      gpdf_dt = new TGraph(gInputFileDT[ifile],"%lg %lg",",");
      double x=0.;
      double y=0.;
      for(int ip=0; ip<gpdf_dt->GetN(); ip++){
	gpdf_dt->GetPoint(ip,x,y);
	hdata->Fill(x,y);
	hscale->Fill(x,1.);
      }
      hdata->Divide(hscale);
      //    hdata->Scale(1./hdata->Integral()); //Nomalize for shape only analysis
      h_dt_charge.push_back(hdata);
    }
  }


}


//Draw histrograms
void DrawHistos(){


  //DRAW PLOTS
  //Charge
  TCanvas *c_charge = new TCanvas("c_charge","c_charge",400*NPMTs,400);
  c_charge->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_charge->cd(ipmt+1);
    h_charge[ipmt]->Draw("");
  }
  if(gInputFileDT.size()>0){
    c_charge->cd(1);
    h_dt_charge[0]->SetLineColor(kRed);
    h_dt_charge[0]->Draw("same");
  }
  TCanvas *c_charge_scat = new TCanvas("c_charge_scat","c_charge_scat",400*NPMTs,400);
  c_charge_scat->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_charge_scat->cd(ipmt+1);
    h_charge_scat[ipmt]->Draw("colz");
  }
  TCanvas *c_mccharge = new TCanvas("c_mccharge","c_mccharge",400*NPMTs,400);
  c_mccharge->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_mccharge->cd(ipmt+1);
    h_mcpmt_charge[ipmt]->Draw("");
  }
  if(gInputFileDT.size()>0){
    c_mccharge->cd(1);
    h_dt_charge[0]->Draw("same");
  }
  TCanvas *c_time = new TCanvas("c_time","c_time",400*NPMTs,400);
  c_time->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_time->cd(ipmt+1);
    h_time[ipmt]->Draw("");
  }
  TCanvas *c_mcfetime = new TCanvas("c_mcfetime","c_mcfetime",400*NPMTs,400);
  c_mcfetime->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_mcfetime->cd(ipmt+1);
    h_mcpmt_fetime[ipmt]->Draw("");
  }
  TCanvas *c_time_diff = new TCanvas("c_time_diff","c_time_diff",400*NPMTs,400);
  c_time_diff->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_time_diff->cd(ipmt+1);
    h_time_diff[ipmt]->Draw("");
  }
  TCanvas *c_npe = new TCanvas("c_npe","c_npe",400*NPMTs,400);
  c_npe->Divide(NPMTs,1);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_npe->cd(ipmt+1);
    h_mcpmt_npe[ipmt]->Draw("");
  }
  TCanvas *c_npevspos = new TCanvas("c_npevspos","c_npevspos",400,400);
  h_mcpmt_npevspos->Draw("colz text");

}



//Output histrograms to file
void PrintHistos(char *filename){

  std::cout<<" Output histrograms to "<<filename<<std::endl;

  TFile *fout = new TFile(filename,"RECREATE");
  fout->cd();
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    h_mcpmt_npe[ipmt]->Write();
    h_mcpmt_charge[ipmt]->Write();
    h_mcpmt_fetime[ipmt]->Write();
    h_charge[ipmt]->Write();
    h_time[ipmt]->Write();
    h_time_diff[ipmt]->Write();
    h_charge_scat[ipmt]->Write();
  }
  h_charge_total->Write();
  h_mcpmt_fetime_total->Write();
  fout->Close();


}

void NormalizeHistos(){

  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    double norm = h_charge[ipmt]->Integral(2,100);
    h_charge[ipmt]->Scale(1./norm);
    std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
    norm = h_mcpmt_charge[ipmt]->Integral(2,100);
    h_mcpmt_charge[ipmt]->Scale(1./norm);
    std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
  }

  if(gInputFileDT.size()>0){
    double norm = h_dt_charge[0]->Integral(2,100);
    h_dt_charge[0]->Scale(1./norm);
  }

}


void ParseArgs(int argc, char **argv){
  bool exist_inputfile = false;
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-mc") {gInputFileMC = argv[++i]; exist_inputfile=true;}
    if(std::string(argv[i]) == "-dt") {gInputFileDT.push_back(argv[++i]); exist_inputfile=true;}
    if(std::string(argv[i]) == "-o")  {gOutFile = argv[++i];}
  }

  if(!exist_inputfile){
    std::cerr<<" Usage: ./DrawChargeAndTime.exe -mc MCFILE [Optional: -dt DATAFILE]"<<std::endl;
    exit(0);
  }
}
