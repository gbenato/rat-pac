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
  hxyplane["start"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,(-1)*intersection_zplane[0],intersection_zplane[0],1000,(-1)*intersection_zplane[1],intersection_zplane[1]);
  hxyplane["Cerenkov"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,(-1)*intersection_zplane[0],intersection_zplane[0],1000,(-1)*intersection_zplane[1],intersection_zplane[1]);
  hxyplane["Scintillation"] = new TH2F("hxyplane_scint","Track intersections with XY plane: Scintillation",1000,(-1)*intersection_zplane[0],intersection_zplane[0],1000,(-1)*intersection_zplane[1],intersection_zplane[1]);

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

  if(debugLevel > 0) std::cout<<"Loading event "<<ievt<<"......."<<std::endl;

  //Initialize
  //Objects
  rds = dsreader->GetEvent(ievt);
  mc = rds->GetMC();
  //Event features
  elength=0.; //e- lenght
  //Particles
  FirstProcessCounter.clear();
  LastProcessCounter.clear();
  ParticleCounter.clear();
  pl_tracks.clear();
  PMTWaveforms.clear();
  PMTDigitizedWaveforms.clear();
  for (std::map<std::string,TH2F*>::iterator it=hxyplane.begin();it!=hxyplane.end();it++)
    it->second->Reset();
  npe.clear();
  if (finalTrack<=0 || finalTrack > mc->GetMCTrackCount()) finalTrack = mc->GetMCTrackCount();
  //Load tracks
  for (int itr = initialTrack; itr < finalTrack; itr++) {
    
    if(debugLevel > 0) std::cout<<"  Track "<<itr<<std::endl;
    
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
      
      if(debugLevel > 0) std::cout<<"  |->Step "<<istep<<std::endl;

      RAT::DS::MCTrackStep *step = mctrack->GetMCTrackStep(istep);
      const TVector3 endpointstep = step->GetEndpoint();
      pl_tracks.back().SetPoint(istep,endpointstep.X(),endpointstep.Y(),endpointstep.Z());

      //Calculate intersection with XY plane
      // std::cout<<"step "<<istep<<" "<<top_pos.X()<<" "<<top_pos.Y()<<" "<<top_pos.Z()<<std::endl;
      // std::cout<<"step "<<istep<<" "<<bottom_pos.X()<<" "<<bottom_pos.Y()<<" "<<bottom_pos.Z()<<std::endl;
      if(mctrack->GetPDGCode()!=0 && mctrack->GetPDGCode()!=9999) continue; //only for OPTICAL photons

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
  for (int ipmt = 0; ipmt < 16; ipmt++){
    npe[ipmt]=0;
  }
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){
    int pmtID = mc->GetMCPMT(ipmt)->GetID();
    npe[pmtID] = mc->GetMCPMT(ipmt)->GetMCPhotonCount();
  }

  if(drawPMTs){
    //Highlight PMT if was hit
    for(int ipmt=0; ipmt<npe.size(); ipmt++){
      if(npe[ipmt]>0) EDGeo->HitPMT(ipmt,npe[ipmt]);
    }
  }
  
  ////////////
  // PMT waveforms
  ////////////
  PMTWaveforms.resize(mc->GetMCPMTCount());
  PMTDigitizedWaveforms.resize(mc->GetMCPMTCount());
  double ymin=9999999.; //yaxis min limit analogue
  int ymax_d=0.; //yaxis max limit digital
  int ymin_d=9999999.; //yaxis min limit digital
  double ymax_temp=0.;
  double ymin_temp=0.;
  double xmax_temp=0.;//dummy
  double xmin_temp=0.;//dummy


  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    
    RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(ipmt);

    //Set analogue graphs
#ifdef __WAVEFORMS_IN_DS__
    vPMTWaveforms[ipmt] = mcpmt->GetWaveform();
#endif
    for(int isample=0; isample<vPMTWaveforms[ipmt].size(); isample++){
      //      std::cout<<"waveform "<<isample<<" "<<vPMTWaveforms[ipmt][isample]<<std::endl;
      PMTWaveforms[ipmt].SetPoint(isample,isample,vPMTWaveforms[ipmt][isample]);
      ymin = TMath::Min(ymin,vPMTWaveforms[ipmt][isample]);
    }

    //Set digitized graphs
#ifdef __WAVEFORMS_IN_DS__
    vPMTDigitizedWaveforms[ipmt] = mcpmt->GetDigitizedWaveform();
#endif
    for(int isample=0; isample<vPMTDigitizedWaveforms[ipmt].size(); isample++){
      PMTDigitizedWaveforms[ipmt].SetPoint(isample,isample,vPMTDigitizedWaveforms[ipmt][isample]);
      ymax_d = TMath::Max(ymax_d,vPMTDigitizedWaveforms[ipmt][isample]);
      ymin_d = TMath::Min(ymin_d,vPMTDigitizedWaveforms[ipmt][isample]);
    }

  }
 
  //Set correct limits for drawing purposes
  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
    PMTWaveforms[ipmt].GetYaxis()->SetRangeUser(1.2*ymin,.5);
    PMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.99*ymin_d,1.01*ymax_d);
  }



  
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
  pmtInfoFileName = dbED->GetS("pmtinfo_file");

  //Validate parameters
  if(event_number<-1) std::cout<<" EventDisplay >>> Event by event mode (Event navigation disabled) "<<std::endl;
  if(!fexists(geoFileName.c_str())) {std::cout<<" EventDisplay >>> "<<geoFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(!fexists(pmtInfoFileName.c_str())) {std::cout<<" EventDisplay >>> "<<pmtInfoFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(drawGeometry) std::cout<<" EventDisplay >>> Draw geometry in "<<geoFileName<<std::endl;
  else std::cout<<" EventDisplay >>> Draw geometry disabled "<<std::endl;
  if(drawPMTs) std::cout<<" EventDisplay >>> Draw PMTs in "<<pmtInfoFileName<<std::endl;
  else std::cout<<" EventDisplay >>> Draw PMTs disabled "<<std::endl;


  
}


//Define experiment geometry
void EventDisplay::SetGeometry(){

  if(debugLevel > 0) std::cout<<" EventDisplay::SetGeometry "<<std::endl;
  
  if(drawPMTs) EDGeo = new EventGeometry(geoFileName, pmtInfoFileName);
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
     pl_tracks[0].Draw("LINE");
  for (int itr = 0; itr < finalTrack - initialTrack; itr++) {
    pl_tracks[itr].Draw("LINE same");
  }

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
  if(mc->GetMCPMTCount()>0){
    
    if(debugLevel > 0) std::cout<<"Display canvas 5 "<<std::endl;
    
    canvas_event->cd(3);
    PMTWaveforms[0].Draw("AP");
    PMTWaveforms[0].GetXaxis()->SetTitle("t(ns)");
    PMTWaveforms[0].GetYaxis()->SetTitle("V");
    for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
      PMTWaveforms[ipmt].SetLineColor(ipmt+1);
      PMTWaveforms[ipmt].Draw("LINE same");
      //      PMTDigitizedWaveforms[ipmt].SetLineColor(kRed);
      //      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
    }
    
    if(debugLevel > 0) std::cout<<"Display canvas 6 "<<std::endl;
    
    canvas_event->cd(4);
    PMTDigitizedWaveforms[0].Draw("AP");
    PMTDigitizedWaveforms[0].GetXaxis()->SetTitle("sample");
    PMTDigitizedWaveforms[0].GetYaxis()->SetTitle("ADC counts");
    for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
      PMTDigitizedWaveforms[ipmt].SetLineColor(ipmt+1);
      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
      //      PMTDigitizedWaveforms[ipmt].ComputeRange(xmin_temp,xmax_temp,ymin_temp,ymax_temp);
    }
  }
#endif
  
  /*
  //Draw grid
  int nlines = 2*XP_XSIDE/pmtwidth;
  TLine* xline;
  TLine* yline;
  for(int iline=0; iline<nlines; iline++){
    xline = new TLine(-XP_XSIDE,iline*pmtwidth-nlines/2.*pmtwidth,XP_XSIDE,iline*pmtwidth-nlines/2.*pmtwidth);
    xline->SetLineWidth(1.);
    xline->SetLineStyle(3);
    xline->SetLineColor(kGray);
    xline->Draw("same");
    yline = new TLine(iline*pmtwidth-nlines/2.*pmtwidth,-XP_YSIDE,iline*pmtwidth-nlines/2.*pmtwidth,XP_YSIDE);
    yline->SetLineWidth(1.);
    yline->SetLineStyle(3);
    yline->SetLineColor(kGray);
    yline->Draw("same");
  }

  //Draw pmts
  for(int ipmt=0; ipmt<vpmtbox.size();ipmt++)
    vpmtbox[ipmt].Draw("LINE same");
  */
  
  //Wait for user action
  canvas_event->Modified();
  canvas_event->Update();
  canvas_event->WaitPrimitive();

  //  if(event_number>=0) exit(0);
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
  for(int ievt=0; ievt<nevents ; ievt++){
    this->DisplayEvent(ievt);
  }

}
