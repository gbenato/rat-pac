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

#define NORM_ATT 1.e-3/1.4
#define NPMTs 2
#define LOOPTRACKS false

char * fInputFile = NULL;
void ParseArgs(int argc, char **argv);

int main(int argc, char **argv){

  int appargc = 0;
  char **appargv = NULL;  
  TApplication dummy("App", &appargc, appargv);
  ParseArgs(argc, argv);
  
  RAT::DB *db = RAT::DB::Get();
  db->Load("../data/OPTICS.ratdb");

  
  std::vector<TH1F*> h_MCPMT_charge; //MC charge
  std::vector<TH1F*> h_charge; //Measured charge
  for(int ih=0; ih<NPMTs; ih++){
    h_MCPMT_charge.push_back(new TH1F(Form("h_mcpmt_charge_%i",ih),"h_mcpmt_charge",100,0,50));
    h_charge.push_back(new TH1F(Form("h_charge_%i",ih),"h_charge",100,0,50));
  }
  TH1F* h_charge_total = new TH1F("h_charge_total","h_charge_total",100,0,50);
  
  TH1F* h_procinit = new TH1F("h_procinit","h_procinit",1,0,1);
  TH1F* h_proclast = new TH1F("h_proclast","h_proclast",1,0,1);
  TH1F* h_procinit_gm = new TH1F("h_procinit_gm","h_procinit_gm",1,0,1);
  TH1F* h_proclast_gm = new TH1F("h_proclast_gm","h_proclast_gm",1,0,1);
  TH1F* h_MCPMT_isdh = new TH1F("h_MCPMT_isdh","h_MCPMT_isdh",2,0,2);
  TH1F* h_MCPMT_time = new TH1F("h_MCPMT_time","h_MCPMT_time",100,1,5);
  TH1F* h_ntracks = new TH1F("h_ntracks","h_ntracks",100,0,5000);
  TH1F* h_ncp = new TH1F("h_ncp","h_ncp",500,0,5000);
  TH1F* h_npe = new TH1F("h_npe","h_npe",500,0,500);
  TH1F* h_e_length = new TH1F("h_e_length","h_e_length",300,0,50); //e- flight path
  TH2F* h_e_lengthvsnph = new TH2F("h_e_lengthvsnph","h_e_lengthvsnph",200,0,500,300,0,50);
  TH2F* h_e_lengthvsnpe = new TH2F("h_e_lengthvsnpe","h_e_lengthvsnpe",50,0,100,300,0,50);
  TH2F* h_e_lengthvsq = new TH2F("h_e_lengthvsq","h_e_lengthvsq",150,0,100,300,0,50);
  TH1F* h_ph_length = new TH1F("h_ph_length","h_ph_length",200,0,200);
  TH1F* h_cp_length = new TH1F("h_cp_length","h_cp_length",200,0,200);
  TH1F* h_cp_ke = new TH1F("h_cp_ke","h_cp_ke",300,0,3e-5);
  TH1F* h_cp_wl = new TH1F("h_cp_wl","h_cp_wl",300,0,1000);
  TH1F* h_ph_last = new TH1F("h_ph_last","h_ph_last",500,-4000,4000);
  TH2F* h_qvspe = new TH2F("h_qvspe","h_qvspe",100,0,100,150,0,100);

  RAT::DSReader *dsreader = new RAT::DSReader(fInputFile);
  int nentries = dsreader->GetT()->GetEntries();
  std::cout<<" Number of entries: "<<nentries<<std::endl;
  for(int ientry=0; ientry<nentries;++ientry){

    if(ientry%(10000) == 0) std::cout<<" Entry "<<ientry<<std::endl;
    RAT::DS::Root *rds = dsreader->GetEvent(ientry);
    RAT::DS::MC *mc = rds->GetMC();

    //***********MC TRUTH
    //track loop
    int ntracks = mc->GetMCTrackCount();
    double elength = 0; //generated e- length
    int ncerphotons = 0; //number of cerenkov photons
    if(LOOPTRACKS){
      h_ntracks->Fill(ntracks);
      //    std::cout<<" Number of tracks in Event "<<ientry<<": "<<ntracks<<std::endl;
      for (int itr = 0; itr < ntracks; itr++) {
	RAT::DS::MCTrack *track = mc->GetMCTrack(itr);
	RAT::DS::MCTrackStep *f_step = track->GetMCTrackStep(0);
	RAT::DS::MCTrackStep *l_step = track->GetLastMCTrackStep();
	
	//Fill histograms
	h_procinit->Fill(f_step->GetProcess().c_str(),1.);
	h_proclast->Fill(l_step->GetProcess().c_str(),1.);
	
	if(track->GetPDGCode() == 11) elength+=track->GetLength();
	else if(track->GetPDGCode() == 22 || track->GetPDGCode() == 0){ //gamma or optical photon
	  h_procinit_gm->Fill(f_step->GetProcess().c_str(),1.);
	  h_proclast_gm->Fill(l_step->GetProcess().c_str(),1.);
	  h_ph_last->Fill(l_step->GetEndpoint().x());
	  h_ph_length->Fill(track->GetLength());
	  if(f_step->GetProcess() == "Cerenkov"){
	    ncerphotons++;
	    h_cp_length->Fill(track->GetLength());
	    h_cp_ke->Fill(f_step->GetKE());
	    h_cp_wl->Fill(1.24/(f_step->GetKE()*1.e3));
	  }
	}
      }
      //end track loop
      h_e_length->Fill(elength);
      h_e_lengthvsnph->Fill(ncerphotons,elength);
      h_e_lengthvsnpe->Fill(mc->GetNumPE(),elength);
      h_ncp->Fill(ncerphotons);
    }
    
    //***********MC
    //MCPMT loop
    h_npe->Fill(mc->GetNumPE());
    for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
      RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      int pmtid = mcpmt->GetID();
      for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++){
	h_MCPMT_isdh->Fill(mcpmt->GetMCPhoton(iph)->IsDarkHit());
	h_MCPMT_charge[pmtid]->Fill(mcpmt->GetMCPhoton(iph)->GetCharge());
	h_MCPMT_time->Fill(mcpmt->GetMCPhoton(iph)->GetHitTime());
	// std::cout<<" IsDarkHit "<< i << mcpmt->GetMCPhoton(i)->isDarkHit <<std::endl;
      }
    }
    //end MCPMT loop


    //*********REAL EVENTS
    //Event loop
    int nevents = rds->GetEVCount();
    for(int ievt=0; ievt<nevents; ievt++){
      RAT::DS::EV *ev = rds->GetEV(ievt);

      double totalcharge = 0.;
      double charge = 0.;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
	int pmtid = ev->GetPMT(ipmt)->GetID();
	charge = ev->GetPMT(ipmt)->GetCharge();
	h_charge[pmtid]->Fill(charge);
	totalcharge += charge;
      }
      if(totalcharge!=0){
	h_e_lengthvsq->Fill(totalcharge,elength);
	h_charge_total->Fill(totalcharge);
	//if(elength<15) h_charge->Fill(totalcharge);
	h_qvspe->Fill(mc->GetNumPE(),totalcharge);
      }
    }
    
       
  }//end entry loop
  

  RAT::DBLinkPtr lAcrylic = db->GetLink("OPTICS", "acrylic_berkeley");
  std::vector<double> abs_x = lAcrylic->GetDArray("ABSLENGTH_value1");
  std::vector<double> abs_y = lAcrylic->GetDArray("ABSLENGTH_value2");
  TGraph *gAtt = new TGraph(abs_x.size(),&abs_x[0],&abs_y[0]);
  for (int i=0;i<gAtt->GetN();i++) gAtt->GetY()[i] *= NORM_ATT;
 
  RAT::DBLinkPtr lPC = db->GetLink("OPTICS", "photocathode_R7081_hqe");
  std::vector<double> qeff_x = lPC->GetDArray("EFFICIENCY_value1");
  std::vector<double> qeff_y = lPC->GetDArray("EFFICIENCY_value2");
  TGraph *gQEff = new TGraph(qeff_x.size(),&qeff_x[0],&qeff_y[0]);
  gQEff->SetLineColor(kOrange);
  gQEff->SetLineStyle(2);
  gQEff->SetLineWidth(2);

 
  TGraph *gB1 = new TGraph("../data/TheiaRnD/Acrylic_b1.csv","%lg %*lg %lg",",");
  TGraph *gB2 = new TGraph("../data/TheiaRnD/Acrylic_b2.csv","%lg %*lg %lg",",");
  TGraph *gB3 = new TGraph("../data/TheiaRnD/Acrylic_b3.csv","%lg %*lg %lg",",");
  TGraph *gB4 = new TGraph("../data/TheiaRnD/Acrylic_b4.csv","%lg %*lg %lg",",");
  TGraph *gB5 = new TGraph("../data/TheiaRnD/Acrylic_b5.csv","%lg %*lg %lg",",");
  TGraph *gB6 = new TGraph("../data/TheiaRnD/Acrylic_b6.csv","%lg %*lg %lg",",");
  TGraph *gB7 = new TGraph("../data/TheiaRnD/Acrylic_b7.csv","%lg %*lg %lg",",");
  TGraph *gB8 = new TGraph("../data/TheiaRnD/Acrylic_b8.csv","%lg %*lg %lg",",");
  gB1->SetLineColor(2);
  gB2->SetLineColor(3);
  gB3->SetLineColor(4);
  gB4->SetLineColor(5);
  gB5->SetLineColor(6);
  gB6->SetLineColor(7);
  gB7->SetLineColor(8);
  gB8->SetLineColor(9);






  
  //DRAW PLOTS
  //Charge
  TCanvas *c_charge = new TCanvas("c_charge","c_charge",1400,1400);
  c_charge->Divide(3,3);
  for(int ipmt=0; ipmt<NPMTs;ipmt++){
    c_charge->cd(ipmt+1);
    h_charge[ipmt]->Draw("");
  }
  // c_charge->cd(1);
  // dsreader->GetT()->Draw("ds.ev.pmt.charge>>test(150,0,100)");
  // h_charge->SetLineColor(kRed);
  // h_charge->Draw("");
  // c_charge->cd(2);
  // h_MCPMT_charge->Draw("same");
  //Process
  TCanvas *c_process = new TCanvas("c_process","c_process",1400,600);
  c_process->Divide(2,1);
  c_process->cd(1);
  h_procinit->Draw();
  c_process->cd(2);
  h_proclast->Draw();
  TCanvas *c_process_gm = new TCanvas("c_process_gm","c_process_gm",1400,600);
  c_process_gm->Divide(2,1);
  c_process_gm->cd(1);
  h_procinit_gm->Draw();
  c_process_gm->cd(2);
  h_proclast_gm->Draw();
  //IsDarkHit
  TCanvas *c_isdh = new TCanvas("c_isdh","c_isdh",600,600);
  c_isdh->cd();
  h_MCPMT_isdh->Draw();
  //Time
  TCanvas *c_time = new TCanvas("c_time","c_time",600,600);
  c_time->cd();
  h_MCPMT_time->Draw();
  //NTracks
  TCanvas *c_ntracks = new TCanvas("c_ntracks","c_ntracks",1400,600);
  c_ntracks->Divide(2,1);
  c_ntracks->cd(1);
  h_ntracks->Draw();
  c_ntracks->cd(2);
  //Generated electrons
  TCanvas *c_e = new TCanvas("c_e","c_e",1400,1400);
  c_e->Divide(3,2);
  c_e->cd(1);
  h_e_length->Draw();
  c_e->cd(2);
  h_npe->Draw();
  c_e->cd(3);
  h_e_lengthvsnph->Draw("colz");
  c_e->cd(4);
  h_e_lengthvsnpe->Draw("colz");
  c_e->cd(5);
  h_e_lengthvsq->Draw("colz");
  c_e->cd(5);
  h_qvspe->Draw();
  //  h_cp_ke->Draw();
  //Cerenkov photons
  TCanvas *c_cp = new TCanvas("c_cp","c_cp",1400,1400);
  c_cp->Divide(2,2);
  c_cp->cd(1);
  h_ncp->Draw();
  c_cp->cd(2);
  h_cp_length->GetYaxis()->SetRangeUser(0.,600.);
  h_cp_length->Draw();
  h_ph_length->SetLineColor(kRed);
  h_ph_length->Draw("same");
  c_cp->cd(3);
  h_cp_wl->Scale(1e-4);
  h_cp_wl->Draw();
  gAtt->Draw("same");
  gB1->Draw("same");
  gB2->Draw("same");
  gB3->Draw("same");
  gB4->Draw("same");
  gB5->Draw("same");
  gB6->Draw("same");
  gB7->Draw("same");
  gB8->Draw("same");
  gQEff->Draw("same");
  c_cp->cd(4);
  h_ph_last->Draw();


  //DRAW TABLE
  double ph_total = h_cp_length->GetEntries();
  double ph_att = h_cp_length->Integral(0,(int)100.*300./6000.);
  //  double ph_hitpmt = h_cp_length->Integral((int)600.*300./6000.,(int)850.*300./6000.); //OLDDB
  //  double ph_hitpmt_total = h_ph_last->Integral((int)(100.+4000)*500./8000.,(int)(220.+4000)*500./8000.); //OLDDB
  double ph_hitpmt_total = h_cp_length->Integral((int)580.*300./6000.,(int)1000.*300./6000.); //NEWDB
  double ph_elec = h_MCPMT_isdh->GetEntries();
  std::cout<<" # e- sttopped at acrylic: "<<h_e_length->Integral(0,20)<<std::endl;
  std::cout<<" # e- sttopped at PMT: "<<h_e_length->Integral(20,50)<<std::endl;
  std::cout<<" # Cherenkov photons: "<<ph_total<<" (in "<<ph_total/nentries<<" per event)"<<std::endl;
  std::cout<<" # photons attenuated: "<<ph_att<<" ("<<ph_att/ph_total*100<<"\%)"<<std::endl;
  std::cout<<" # photons hitting PC: "<<ph_hitpmt_total<<" ("<<ph_hitpmt_total/ph_total*100<<"\%)"<<std::endl;
  std::cout<<" # photo-electrons: "<<ph_elec<<" ("<<ph_elec/ph_hitpmt_total*100<<"\%)"<<std::endl;
  std::cout<<" 0 photo-electrons: "<<h_npe->GetBinContent(1)<<" ("<<h_npe->GetBinContent(1)/nentries*100<<"\%)"<<std::endl;
  std::cout<<" 1 photo-electrons: "<<h_npe->GetBinContent(2)<<" ("<<h_npe->GetBinContent(2)/nentries*100<<"\%)"<<std::endl;
  std::cout<<" Multi photo-electrons: "<<h_npe->Integral(3,10)<<" ("<<h_npe->Integral(3,10)/nentries*100<<"\%)"<<std::endl;
  std::cout<<" QEff: "<<ph_elec/ph_hitpmt_total*100<<"\%"<<std::endl;
  
  new TBrowser;
  dummy.Run();
  return 0;

}


void ParseArgs(int argc, char **argv){
  bool exist_inputfile = false;
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-i") {fInputFile = argv[++i]; exist_inputfile=true;}
  }
  if(!exist_inputfile){
    std::cerr<<" Specify input file with option: '-i'"<<std::endl;
    exit(0);
  }
}
