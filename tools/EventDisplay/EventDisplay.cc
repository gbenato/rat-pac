#include<TH2F.h>
#include<TCanvas.h>
#include<TPolyLine3D.h>
#include<TVector3.h>
#include<TMath.h>
#include<TGeoBBox.h>
#include<TGeoManager.h>
#include<TGeoMaterial.h>
#include<TGeoMedium.h>
#include<TGeoVolume.h>
#include<TLine.h>
#include<TPaveText.h>
#include<TStyle.h>

#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include<RAT/DB.hh>

#include "EventDisplay.hh"

bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

EventDisplay::EventDisplay(){

  //Init
  SetParameters();
  OpenFile(inputFileName);
  int appargc = 0;
  char **appargv = NULL;
  dummyApp = new TApplication("EventDisplay", &appargc, appargv);

  //Set canvas
  gStyle->SetGridWidth(1);
  canvas_event = new TCanvas("canvas_event", "Event", 1200, 800);
  canvas_event->Divide(2,2);
  canvas_event->cd(1)->SetPad(0., .33, .5, 1.);
  canvas_event->cd(2)->SetPad(.5, .33, 1., 1.);
  canvas_event->cd(3)->SetPad(.0, 0., .5, .33);
  canvas_event->cd(4)->SetPad(.5, 0., 1., .33);

  //Particle maps
  ParticleColor[11]=kGreen;   ParticleWidth[11]=1;   ParticleName[11]="Electron";
  ParticleColor[22]=kYellow;     ParticleWidth[22]=1;   ParticleName[22] = "Standard photon";
  ParticleColor[13]=kOrange;  ParticleWidth[13]=2;   ParticleName[13] = "Muon";
  ParticleColor[211]=kOrange; ParticleWidth[211]=2;  ParticleName[211]= "Pi+";
  ParticleColor[0]=kCyan+1;     ParticleWidth[0]=1;    ParticleName[0] = "Cherenkov photon"; //Indeed this is an optical photon, but I changed the definition
  ParticleColor[9999]=kRed-7;     ParticleWidth[9999]=1;    ParticleName[9999] = "Scintillation photon"; //Created by me, PDG number doesn't actually exist

  //Representation plane
  hxyplane["start"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);
  hxyplane["Cerenkov"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);
  hxyplane["Scintillation"] = new TH2F("hxyplane_scint","Track intersections with XY plane: Scintillation",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);

  SetGeometry();

  if(debugLevel > 0) std::cout<<" EventDisplay::EventDisplay - DONE "<<std::endl;

};


void EventDisplay::OpenFile(std::string inputfile){

  std::cout<<" EventDisplay >>> Opening file "<<inputfile<<" ..... "<<std::endl;
  dsreader = new RAT::DSReader(inputfile.c_str());
  nevents = dsreader->GetT()->GetEntries();

  if(debugLevel > 0) std::cout<<" DONE! "<<nevents<<" in file."<<std::endl;

};


void EventDisplay::CustomizeTrack(TPolyLine3D *track, RAT::DS::MCTrack *mctrack){

  RAT::DS::MCTrackStep *firststep = mctrack->GetMCTrackStep(0);
  if(firststep->GetProcess()=="Scintillation")
  mctrack->SetPDGCode(9999);

  //Set color
  //track->SetLineColor(ParticleColor[mctrack->GetPDGCode()]);
  track->SetLineColorAlpha(ParticleColor[mctrack->GetPDGCode()],0.2);
  //Set width
  track->SetLineWidth(ParticleWidth[mctrack->GetPDGCode()]);
  //Count particles
  ParticleCounter[ParticleName[mctrack->GetPDGCode()]] += 1;

}

void EventDisplay::LoadEvent(int ievt){

  if(debugLevel > 0) std::cout<<"EventDisplay::LoadEvent -- Loading event "<<ievt<<"......."<<std::endl;

  //Initialize
  SetParameters();
  //Objects
  rds = dsreader->GetEvent(ievt);
  mc = rds->GetMC();
  if(rds->ExistEV()) ev = rds->GetEV(0); //FIXME: so far get only first event
  else if(debugLevel > 0) std::cout<<"EventDisplay::LoadEvent -- EV do not exist! "<<std::endl;



  //Event features
  elength=0.; //e- lenght
  //Particles
  FirstProcessCounter.clear();
  LastProcessCounter.clear();
  ParticleCounter.clear();
  pl_tracks.clear();
  MCPMTWaveforms.clear();
  MCPMTDigitizedWaveforms.clear();
  PMTDigitizedWaveforms.clear();
  for (std::map<std::string,TH2F*>::iterator it=hxyplane.begin();it!=hxyplane.end();it++)
  it->second->Reset();
  npe.clear();
  if (finalTrack<=0 || finalTrack > mc->GetMCTrackCount()) finalTrack = mc->GetMCTrackCount();
  //Load tracks
  for (int itr = initialTrack; itr < finalTrack; itr++) {

    if(debugLevel > 0) std::cout<<"  Track "<<itr<<"/"<<finalTrack-1<<std::endl;

    RAT::DS::MCTrack *mctrack = mc->GetMCTrack(itr);
    //Create new track
    pl_tracks.resize(pl_tracks.size()+1);
    //Set PDGcode color code
    CustomizeTrack(&pl_tracks.back(), mctrack);
    //Measure electron length
    if(mctrack->GetPDGCode()==11) elength += mctrack->GetLength();
    //Count processes
    RAT::DS::MCTrackStep *firststep = mctrack->GetMCTrackStep(0);
    RAT::DS::MCTrackStep *laststep = mctrack->GetLastMCTrackStep();
    double last_pos[3];
    laststep->GetEndpoint().GetXYZ(last_pos);
    //      if (itr>=10 && itr<11) std::cout<<" Last pos: "<<last_pos[0]<<" "<<last_pos[1]<<" "<<last_pos[2]<<" "<<std::endl;
    FirstProcessCounter[firststep->GetProcess()] += 1;
    LastProcessCounter[laststep->GetProcess()] += 1;

    //Make XY-plane for ring display
    //Loop over all the steps
    int nsteps = mctrack->GetMCTrackStepCount();
    TVector3 top_pos(-9999.,-9999.,-9999.); //first interpolation point
    TVector3 bottom_pos(9999.,9999.,9999.); //second interpolation point
    TVector3 int_pos(9999.,9999.,9999.);//intersection point
    for (int istep = 0; istep < nsteps; istep++) {

      if(debugLevel > 0) std::cout<<"  |->Step "<<istep<<"/"<<nsteps-1<<" ";

      RAT::DS::MCTrackStep *step = mctrack->GetMCTrackStep(istep);
      const TVector3 endpointstep = step->GetEndpoint();
      pl_tracks.back().SetPoint(istep,endpointstep.X(),endpointstep.Y(),endpointstep.Z());

      if(debugLevel > 0) std::cout<<" -- DONE "<<std::endl;

      //Calculate intersection with XY plane
      // std::cout<<"step "<<istep<<" "<<top_pos.X()<<" "<<top_pos.Y()<<" "<<top_pos.Z()<<std::endl;
      // std::cout<<"step "<<istep<<" "<<bottom_pos.X()<<" "<<bottom_pos.Y()<<" "<<bottom_pos.Z()<<std::endl;
      if(mctrack->GetPDGCode()!=0 && mctrack->GetPDGCode()!=9999) continue; //only for OPTICAL photons

      if(debugLevel > 0) std::cout<<"    Intersectng Optical Photon "<<std::endl;

      if(bottom_pos.Z()!=-9999.){ //we haven't found the point yet
        if(endpointstep.Z()>intersection_zplane[2]){
          if(debugLevel > 0) std::cout<<"      Case 1: "<<std::endl;
          top_pos = endpointstep;
        }
        else if(top_pos.Z()!=-9999.){ //this is our guy

          if(debugLevel > 0) std::cout<<"      Case 2: "<<firststep->GetProcess()<<std::endl;

          bottom_pos = endpointstep;
          //Intersect!
          double lambda = (intersection_zplane[2] - top_pos.Z())/(bottom_pos.Z() - top_pos.Z());
          int_pos = top_pos + (bottom_pos - top_pos)*lambda;
          //	  std::cout<<"FILL IT! "<<int_pos.X()<<" "<<int_pos.Y()<<" "<<int_pos.Z()<<std::endl;
          if(firststep->GetProcess()=="Reemission") hxyplane["Scintillation"]->Fill(int_pos.X(),int_pos.Y());
          else hxyplane[firststep->GetProcess()]->Fill(int_pos.X(),int_pos.Y());
          bottom_pos.SetZ(-9999.);

        }
      }

      if(debugLevel > 0) std::cout<<"   EventDisplay::LoadEvent (Passed intersection) "<<std::endl;

    } //end step loop
  } //end track loop

  //Load photoelectrons
  //  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

  for (int ipmt = 0; ipmt < EDGeo->GetPMTCount(); ipmt++){
    npe[ipmt]=0;
  }
  for(int ipmt=0; ipmt<mc->GetMCPMTCount(); ipmt++){
    int pmtID = mc->GetMCPMT(ipmt)->GetID();
    npe[pmtID] = mc->GetMCPMT(ipmt)->GetMCPhotonCount();
  }

  if(drawPMTs){
    //Highlight PMT if was hit
    for(std::map<int, int>::iterator itpmt=npe.begin(); itpmt!=npe.end(); itpmt++){
      EDGeo->HitPMT(itpmt->first,itpmt->second);
    }
  }


#ifdef __WAVEFORMS_IN_DS__
  ////////////
  // PMT waveforms
  ////////////
  MCPMTWaveforms.resize(mc->GetMCPMTCount());
  MCPMTDigitizedWaveforms.resize(mc->GetMCPMTCount());
  double ymin=9999999.; //yaxis min limit analogue
  UShort_t ymax_d=0.; //yaxis max limit digital
  UShort_t ymin_d=9999999.; //yaxis min limit digital
  double ymax_temp=0.;
  double ymin_temp=0.;
  double xmax_temp=0.;//dummy
  double xmin_temp=0.;//dummy

  if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - GetMCPMTCount "<<mc->GetMCPMTCount()<<std::endl;

  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {

    RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(ipmt);

    //Set analogue graphs
    vMCPMTWaveforms[ipmt] = mcpmt->GetWaveform();
    for(int isample=0; isample<vMCPMTWaveforms[ipmt].size(); isample++){
      //      std::cout<<"waveform "<<isample<<" "<<vPMTWaveforms[ipmt][isample]<<std::endl;
      MCPMTWaveforms[ipmt].SetPoint(isample,isample,vMCPMTWaveforms[ipmt][isample]);
      ymin = TMath::Min(ymin,vMCPMTWaveforms[ipmt][isample]);
    }

    //Set digitized graphs
    vMCPMTDigitizedWaveforms[ipmt] = mcpmt->GetDigitizedWaveform();
    for(int isample=0; isample<vMCPMTDigitizedWaveforms[ipmt].size(); isample++){
      MCPMTDigitizedWaveforms[ipmt].SetPoint(isample,isample,vMCPMTDigitizedWaveforms[ipmt][isample]);
      ymax_d = TMath::Max(ymax_d,vMCPMTDigitizedWaveforms[ipmt][isample]);
      ymin_d = TMath::Min(ymin_d,vMCPMTDigitizedWaveforms[ipmt][isample]);
    }

  }

  //Set correct limits for drawing purposes
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    MCPMTWaveforms[ipmt].GetYaxis()->SetRangeUser(1.2*ymin,.5);
    MCPMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.99*ymin_d,1.01*ymax_d);
  }

  if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - GetPMTCount "<<ev->GetPMTCount()<<std::endl;

  if(!rds->ExistEV()) {
    if(debugLevel > 0) std::cout<<" EventDisplay::LoadEvent - DONE (EV do not exist) "<<std::endl;
    return;
  }

  PMTDigitizedWaveforms.resize(ev->GetPMTCount());
  for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {

    RAT::DS::PMT *pmt = ev->GetPMT(ipmt);
    vPMTDigitizedWaveforms[ipmt] = pmt->GetWaveform();
    if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - DigitWF: nsamples "<<vPMTDigitizedWaveforms[ipmt].size()<<std::endl;

    for(int isample=0; isample<vPMTDigitizedWaveforms[ipmt].size(); isample++){
      PMTDigitizedWaveforms[ipmt].SetPoint(isample,isample,vPMTDigitizedWaveforms[ipmt][isample]);

      if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - Digit WF: "
      <<"sample "<<isample<<" "<<vPMTDigitizedWaveforms[ipmt][isample]<<std::endl;

      ymax_d = TMath::Max(ymax_d,vPMTDigitizedWaveforms[ipmt][isample]);
      ymin_d = TMath::Min(ymin_d,vPMTDigitizedWaveforms[ipmt][isample]);
    }
  }

  //Set correct limits for drawing purposes
  for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {
    PMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.99*ymin_d,1.01*ymax_d);
  }

#endif

  if(debugLevel > 0) std::cout<<" EventDisplay::LoadEvent - DONE "<<std::endl;

}

//Define parameters
void EventDisplay::SetParameters(){

  //Get link to ratdb file
  RAT::DB* db = RAT::DB::Get();
  db->Load("ED.json");
  dbED = db->GetLink("EVENTDISPLAY");

  //Flags
  debugLevel = dbED->GetI("debug_level");
  drawGeometry = dbED->GetI("draw_geo");
  drawPMTs = dbED->GetI("draw_pmts");
  initialTrack = dbED->GetI("initial_track");
  finalTrack = dbED->GetI("final_track");
  event_option = dbED->GetS("event_option");
  event_number = dbED->GetI("event_number");
  intersection_zplane = dbED->GetDArray("intersection_zplane");

  //Analysis file
  inputFileName = dbED->GetS("input_file");

  //Geometry files
  geoFileName = dbED->GetS("geo_file");

  //Validate parameters
  if(event_number<-1) std::cout<<" EventDisplay >>> Event by event mode (Event navigation disabled) "<<std::endl;
  if(!fexists(geoFileName.c_str())) {std::cout<<" EventDisplay >>> "<<geoFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(!fexists(inputFileName.c_str())) {std::cout<<" EventDisplay >>> "<<inputFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(drawGeometry) std::cout<<" EventDisplay >>> Draw geometry in "<<geoFileName<<std::endl;
  else std::cout<<" EventDisplay >>> Draw geometry disabled "<<std::endl;
  if(drawPMTs) std::cout<<" EventDisplay >>> Draw PMTs from the PMTInfo header "<<std::endl;
  else std::cout<<" EventDisplay >>> Draw PMTs disabled "<<std::endl;



}


//Define experiment geometry
void EventDisplay::SetGeometry(){

  if(debugLevel > 0) std::cout<<" EventDisplay::SetGeometry "<<std::endl;

  if(drawPMTs) EDGeo = new EventGeometry(geoFileName, inputFileName);
  else EDGeo = new EventGeometry(geoFileName);

  if(debugLevel > 0) std::cout<<" EventDisplay::SetGeometry - DONE "<<std::endl;

}

void EventDisplay::DumpEventInfo(int ievt){

  std::cout<<"********EVENT "<<ievt<<"/"<<nevents<<"********"<<std::endl;
  for (std::map<int,int>::iterator it=npe.begin();it!=npe.end();it++){
    if(it==npe.begin()) std::cout<<"Number of PE"<<std::endl;
    if(it->second!=0) std::cout<<"ID: "<<it->first<<" -> "<<it->second<<std::endl;
  }

  std::cout<<"Electron lenght: "<<elength<<" mm"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"    INITIAL PROCESSES "<<std::endl;
  for (std::map<std::string,int>::iterator it=FirstProcessCounter.begin();it!=FirstProcessCounter.end();it++){
    if(it->second!=0) std::cout<<it->first<<" == "<<it->second<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<"    END PROCESSES   "<<std::endl;
  for (std::map<std::string,int>::iterator it=LastProcessCounter.begin();it!=LastProcessCounter.end();it++){
    if(it->second!=0) std::cout<<it->first<<" == "<<it->second<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<"    TRACKS        "<<std::endl;
  for (std::map<std::string,int>::iterator it=ParticleCounter.begin();it!=ParticleCounter.end();it++){
    if(it->second!=0) std::cout<<it->first<<" == "<<it->second<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<"    WAVEFORMS      "<<std::endl;
  std::cout<<" Number of Waveforms "<<mc->GetMCPMTCount()<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<std::endl;
  std::cout<<" Press any key to go to next event "<<std::endl;

}


void EventDisplay::DisplayEvent(int ievt){

  this->LoadEvent(ievt);
  if(event_option == "cherenkov" && !this->IsCerenkov()) return;
  if(event_option == "pe" && !this->IsPE()) return;
  this->DumpEventInfo(ievt);
  if(event_number>=0) this->DumpDisplayInfo();


  if(debugLevel > 0) std::cout<<"Display canvas 1 "<<std::endl;

  canvas_event->SetTitle(Form("Event %d/%d",ievt,nevents));

  //3D display
  canvas_event->cd(1);
  if(drawGeometry) EDGeo->DrawGeometry();
  for (int itr = 0; itr < finalTrack - initialTrack; itr++) {
    if(itr==0) pl_tracks[itr].Draw("LINE");
    pl_tracks[itr].Draw("LINE same");
  }

  if(debugLevel > 0) std::cout<<"Display canvas 2 "<<std::endl;

  //2D display
  canvas_event->cd(2);
  hxyplane["start"]->SetLineColor(ParticleColor[0]);
  hxyplane["Cerenkov"]->SetLineColor(ParticleColor[0]);
  hxyplane["Scintillation"]->SetLineColor(ParticleColor[9999]);
  if(hxyplane["Cerenkov"]->GetEntries()>0)
  hxyplane["Cerenkov"]->Draw("box");
  else if(hxyplane["Scintillation"]->GetEntries()>0)
  hxyplane["Scintillation"]->Draw("box");
  else
  hxyplane["start"]->Draw("box");
  for (std::map<std::string,TH2F*>::iterator it=hxyplane.begin();it!=hxyplane.end();it++){
    it->second->Draw("box same");
  }
  if(drawPMTs) EDGeo->DrawPMTMap();

#ifdef __WAVEFORMS_IN_DS__
  //Waveforms
  if(debugLevel > 0) std::cout<<"Display canvas 3 "<<std::endl;

  canvas_event->cd(3);
  //Analogue waveforms
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    if(ipmt==0){
      MCPMTWaveforms[ipmt].Draw("AP");
      MCPMTWaveforms[ipmt].GetXaxis()->SetTitle("t(ns)");
      MCPMTWaveforms[ipmt].GetYaxis()->SetTitle("V");
    }
    MCPMTWaveforms[ipmt].SetLineColor(ipmt+1);
    MCPMTWaveforms[ipmt].Draw("LINE same");
    //      PMTDigitizedWaveforms[ipmt].SetLineColor(kRed);
    //      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
  }

  if(debugLevel > 0) std::cout<<"Display canvas 4 "<<std::endl;
  //Digitized waveforms
  canvas_event->cd(4);
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    if(ipmt==0){
      MCPMTDigitizedWaveforms[ipmt].SetTitle("Non-triggered event");
      MCPMTDigitizedWaveforms[ipmt].Draw("AP");
      MCPMTDigitizedWaveforms[ipmt].GetXaxis()->SetTitle("sample");
      MCPMTDigitizedWaveforms[ipmt].GetYaxis()->SetTitle("ADC counts");
    }
    MCPMTDigitizedWaveforms[ipmt].SetLineColor(ipmt+1);
    MCPMTDigitizedWaveforms[ipmt].Draw("LINE same");
    //      MCPMTDigitizedWaveforms[ipmt].ComputeRange(xmin_temp,xmax_temp,ymin_temp,ymax_temp);
  }
  if(rds->ExistEV()){
    for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {
      if(ipmt==0){
        PMTDigitizedWaveforms[ipmt].SetTitle("Triggered event");
        PMTDigitizedWaveforms[ipmt].Draw("AP");
        PMTDigitizedWaveforms[ipmt].GetXaxis()->SetTitle("sample");
        PMTDigitizedWaveforms[ipmt].GetYaxis()->SetTitle("ADC counts");
      }
      PMTDigitizedWaveforms[ipmt].SetLineColor(ipmt+1);
      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
      //      PMTDigitizedWaveforms[ipmt].ComputeRange(xmin_temp,xmax_temp,ymin_temp,ymax_temp);
    }
  }
#endif

  //Wait for user action
  canvas_event->Modified();
  canvas_event->Update();
  canvas_event->WaitPrimitive();

  if(event_number>=0) dummyApp->Run();

  if(debugLevel > 0) std::cout<<" EventDisplay::DisplayEvent - DONE "<<std::endl;

}

bool EventDisplay::IsCerenkov(){

  return FirstProcessCounter["Cerenkov"]>0;

}

bool EventDisplay::IsPE(){

  int npe_total = 0;
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++)
  npe_total += npe[mc->GetMCPMT(ipmt)->GetID()];

  return npe_total>0;

}

void EventDisplay::DumpDisplayInfo(){

  std::cout<<"***NAVIGATION CONTROL***"<<std::endl;
  std::cout<<"Left: H"<<std::endl;
  std::cout<<"Right: L"<<std::endl;
  std::cout<<"Up: U"<<std::endl;
  std::cout<<"Down: I"<<std::endl;
  std::cout<<"Zoom in: J"<<std::endl;
  std::cout<<"Zoom out: K"<<std::endl;

}

void EventDisplay::Open(){

  //Display events
  if(event_number>=0) this->DisplayEvent(event_number);
  else{
    for(int ievt=0; ievt<nevents ; ievt++){
      this->DisplayEvent(ievt);
    }
  }
}
