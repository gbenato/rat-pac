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

#define USERROOTLOOP true
#define DRAWONLYCHARGE false
#define NLOGENTRIES 20
#define NBINS 20
#define MIN -500.0
#define MAX -300.0

//Methods
char *gInputFileMC = NULL;
std::vector<char*> gInputFileDT;
char *gOutFile = NULL;
char *gTargetMaterial = "WBLS10pct";
void ParseArgs(int argc, char **argv);
void GetDBPlots();
void GetHistos();
void DrawHistos();
void PrintHistos(char*);
void NormalizeHistos();

//Global variables
std::vector<TVector3*> pos_pmts; //PMT positions by ID
int npmts; // # PMTs
std::map<int,int> npmts_type; // # PMTs per type
std::vector<int> pmtidtotype; // PMT ID -> PMT Type
vector<double> qScintCorr; //Scintillation correction
vector<double> qScintCorrErr; // Scintillation correction error


//// Histograms
//Event level
TH1F* h_chi2;
TH1F* h_charge_total;
TH1F* h_mcpmt_fetime_total;
TH1F* h_mcpmt_pe_total;
TH2F* h_mcpmt_npevspos;
TH2F* h_pmt_qvspos;

//Hit level
std::vector<TH1F*> h_mcpmt_npe; //Number of PE by PMT ID
std::vector<TH1F*> h_mcpmt_charge; //MC charge
std::vector<TH1F*> h_mcpmt_fetime; //MC FE time
std::vector<TH1F*> h_charge; //Measured charge
std::vector<TH1F*> h_time; //Measured time
std::vector<TH1F*> h_time_diff; //Time diff between EV-MC
std::vector<TH2F*> h_charge_scat; //PMT charge vs trigger charge

//Real data
std::vector<TH1F*> h_dt_charge; //Measured charge

int main(int argc, char **argv){

  //Init********
  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy("App", &appargc, appargv);
  ParseArgs(argc, argv);
  //************

  GetDBPlots();
  GetHistos();
  NormalizeHistos();

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
  db->Load("/Users/snoplus/Work/snoing/install/rat-pac/data/TheiaRnD/SCINTCORR.ratdb");
  RAT::DBLinkPtr dbScintCorr = db->GetLink("SCINTCORR",gTargetMaterial);
  qScintCorr = dbScintCorr->GetDArray("corr");
  qScintCorrErr = dbScintCorr->GetDArray("corr_err");

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
  npmts = pmtInfo->GetPMTCount();
  npmts_type = pmtInfo->GetPMTTypeCount();

  for (size_t ipmt = 0; ipmt < npmts; ipmt++) {
    TVector3* pmt_temp = new TVector3(pmtInfo->GetPosition(ipmt)[0],
                                      pmtInfo->GetPosition(ipmt)[1],
                                      pmtInfo->GetPosition(ipmt)[2]);
    pos_pmts.push_back(pmt_temp);
    pmtidtotype.push_back(pmtInfo->GetType(ipmt));
  }

  //Init histos
  //  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",27,-27*15,27*15,27,-27*15,27*15);
  h_mcpmt_npevspos = new TH2F("h_mcpmt_npevspos","h_mcpmt_npevspos",NBINS,MIN,MAX,NBINS,MIN,MAX);
  h_pmt_qvspos = new TH2F("h_pmt_qvspos","h_pmt_qvspos",NBINS,MIN,MAX,NBINS,MIN,MAX);
  for(int ih=0; ih<npmts; ih++){
    h_mcpmt_npe.push_back(new TH1F(Form("h_mcpmt_npe_%i",ih),"h_mcpmt_npe",200,0,200));
    h_mcpmt_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",200,0,100));
    h_mcpmt_fetime.push_back(new TH1F(Form("h_mcpmt_fetime_%i",ih),"h_mcpmt_fetime",1000,0,200));
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),"h_charge",200,0,800));
    h_charge_scat.push_back(new TH2F(Form("h_charge_scat_%i",ih),"h_charge_scat",200,0,100,200,0,100));
    h_time.push_back(new TH1F(Form("h_time_%i",ih),"h_time",1000,0,200));
    h_time_diff.push_back(new TH1F(Form("h_time_diff_%i",ih),"h_time_diff",100,-100,100));
  }
  h_charge_total = new TH1F("h_charge_total","h_charge_total",100,0,150);
  h_mcpmt_fetime_total  = new TH1F("h_mcpmt_fetime_total","h_mcpmt_fetime_total",500,0,100);
  h_mcpmt_pe_total  = new TH1F("h_mcpmt_pe_total","h_mcpmt_pe_total",500,0,100);
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


      if(mc->GetMCPMTCount()==0) continue;

      //MC**************
      //MCPMT loop
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
        RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
        //Make grid display only for the small PMTs
        int pmtid = mcpmt->GetID();
        //count PE
        h_mcpmt_npe[pmtid]->Fill(mcpmt->GetMCPhotonCount());
        h_mcpmt_npevspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),mcpmt->GetMCPhotonCount()/(double)nentries);

        if(!DRAWONLYCHARGE){
          for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++){
            h_mcpmt_charge[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetCharge());
            h_mcpmt_fetime[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
            if(mcpmt->GetType()==1) h_mcpmt_fetime_total->Fill(mcpmt->GetMCPhoton(iph)->GetFrontEndTime());
            if(mcpmt->GetType()==1) h_mcpmt_pe_total->Fill(mcpmt->GetMCPhoton(iph)->GetHitTime());
          }
        }

      } //end MCPMT loop

      //DAQ EVENTS*****
      //Event loop
      int nevents = rds->GetEVCount();
      for(int ievt=0; ievt<nevents; ievt++){
        RAT::DS::EV *ev = rds->GetEV(ievt);
        double chi2 = 0.;
        double totalcharge = 0.;
        double charge = 0.;
        double charge_trig = 0.;
        //First PMT loop:
        // - Look for trigger PMT and save charge
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
          h_pmt_qvspos->Fill(pos_pmts[pmtid]->X(),pos_pmts[pmtid]->Y(),charge);
          h_charge_scat[pmtid]->Fill(charge,charge_trig); //charge vs trigger charge
          h_time[pmtid]->Fill(ev->GetPMT(ipmt)->GetTime());
          totalcharge += charge;
          //Compute chi2 for cher/scint
          if(ev->GetPMT(ipmt)->GetType()==1){
            chi2 += pow( (charge - qScintCorr[ipmt])/qScintCorrErr[ipmt], 2.);
            //            std::cout<<" chi2 "<<ipmt<<" "<<chi2<<" "<<charge<<" "<<qScintCorr[ipmt]<<" "<<qScintCorrErr[ipmt]<<std::endl;
          }
        }
        if(totalcharge!=0){
          h_charge_total->Fill(totalcharge);
        }
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

  //General plots
  TCanvas *c_total = new TCanvas("c_total","c_total",800,400);
  c_total->Divide(2,1);
  c_total->cd(1);
  h_mcpmt_fetime_total->Draw();
  c_total->cd(2);
  h_mcpmt_pe_total->Draw();

  int nvar = 7;
  std::map<int, TCanvas*> c_pmt;
  for(std::map<int,int>::iterator itype = npmts_type.begin(); itype != npmts_type.end(); itype++){

    //Init canvases
    int npmts_used = npmts_type[itype->first];
    c_pmt[itype->first] = new TCanvas(Form("c_pmt_type%d",itype->first),Form("c_pmt_type%d",itype->first),400*npmts_used,nvar*400);
    c_pmt[itype->first]->Divide(npmts_used,nvar);

    //Draw
    int icd = 1;
    for(int ipmt=0; ipmt<npmts; ipmt++){
      if(pmtidtotype[ipmt] == itype->first){
        c_pmt[itype->first]->cd(icd+npmts_used*0);
        h_charge[ipmt]->Draw("");
        c_pmt[itype->first]->cd(icd+npmts_used*1);
        h_mcpmt_charge[ipmt]->Draw("");
        c_pmt[itype->first]->cd(icd+npmts_used*2);
        h_charge_scat[ipmt]->Draw("colz");
        c_pmt[itype->first]->cd(icd+npmts_used*3);
        h_time[ipmt]->Draw("");
        c_pmt[itype->first]->cd(icd+npmts_used*4);
        h_mcpmt_fetime[ipmt]->Draw("");
        c_pmt[itype->first]->cd(icd+npmts_used*5);
        h_time_diff[ipmt]->Draw("");
        c_pmt[itype->first]->cd(icd+npmts_used*6);
        h_mcpmt_npe[ipmt]->Draw("");
        icd++;
      }
    }
    if(gInputFileDT.size()>0){
      c_pmt[itype->first]->cd(1*1);
      h_dt_charge[0]->SetLineColor(kRed);
      h_dt_charge[0]->Draw("same");
      c_pmt[itype->first]->cd(1*2);
      h_dt_charge[0]->Draw("same");
    }

  }// end PMT type loop

  TCanvas *c_npevspos = new TCanvas("c_npevspos","c_npevspos",400,400);
  c_npevspos->Divide(2,1);
  c_npevspos->cd(1);
  h_mcpmt_npevspos->Draw("colz text");
  c_npevspos->cd(2);
  h_pmt_qvspos->Draw("colz text");

  TCanvas *c_chi2 = new TCanvas("c_chi2","c_chi2",400,400);
  h_chi2->Draw();

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
    h_charge_scat[ipmt]->Write();
  }
  h_chi2->Write();
  h_charge_total->Write();
  h_mcpmt_fetime_total->Write();
  fout->Close();


}

void NormalizeHistos(){

  for(int ipmt=0; ipmt<npmts;ipmt++){
    double norm = h_charge[ipmt]->Integral(2,100);
    h_charge[ipmt]->Scale(1./norm);
    // std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
    norm = h_mcpmt_charge[ipmt]->Integral(2,100);
    h_mcpmt_charge[ipmt]->Scale(1./norm);
    // std::cout<<" norm "<<ipmt<<" "<<norm<<std::endl;
  }

  if(gInputFileDT.size()>0){
    double norm = h_dt_charge[0]->Integral(2,100);
    h_dt_charge[0]->Scale(1./norm);
  }

  double norm = (double)h_pmt_qvspos->GetEntries();
  // std::cout<<" h_pmt_qvspos "<<norm<<std::endl;
  h_pmt_qvspos->Scale(1./norm);

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
