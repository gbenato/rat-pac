//ROOT stuff
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
#include<TColor.h>
#include<TMarker.h>

//RAT Stuff
#include<RAT/DS/MC.hh>
#include<RAT/DS/MCTrack.hh>
#include<RAT/DS/MCTrackStep.hh>
#include<RAT/DS/MCPMT.hh>
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include<RAT/DB.hh>
#include<RAT/DS/Run.hh>

#include "EventDisplay.hh"

bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

EventDisplay::EventDisplay(std::string _inputFileName){

  //Init
  inputFileName = _inputFileName;
  SetParameters();
  OpenFile(inputFileName);
  int appargc = 0;
  char **appargv = NULL;
  dummyApp = new TApplication("EventDisplay", &appargc, appargv);

  //Set canvas
  gStyle->SetGridWidth(1);
  canvas_event = new TCanvas("canvas_event", "Event", 1600, 1000);
  canvas_event->Divide(2,5);
  //MC tracks and geometry
  canvas_event->cd(1)->SetPad(.0, .5, .33, 1.);
  //2D Ring plane
  canvas_event->cd(2)->SetPad(.33, .5, .66, 1.);
  //Charge vs Pos
  canvas_event->cd(3)->SetPad(.0, .0, .33, .5);
  canvas_event->cd(9)->SetPad(.05, .27, .15, .45);
  //Time vs Pos
  canvas_event->cd(4)->SetPad(.33, .0, .66, .5);
  //Ring PMT traces
  canvas_event->cd(5)->SetPad(.66, .75, .99, 1.);
  //Light PMT traces
  canvas_event->cd(6)->SetPad(.66, .5, .99, .75);
  //Muon PMT traces
  canvas_event->cd(7)->SetPad(.66, .25, .99, .5);
  //Spare
  canvas_event->cd(8)->SetPad(.66, .0, .99, .25);
  //Hidden
  canvas_event->cd(10)->SetPad(.99, .99, 1., 1.);

  //Particle maps
  ParticleColor[11]=kGreen;   ParticleWidth[11]=1;   ParticleName[11]="Electron";
  ParticleColor[22]=kYellow;     ParticleWidth[22]=1;   ParticleName[22] = "Standard photon";
  ParticleColor[13]=kOrange;  ParticleWidth[13]=2;   ParticleName[13] = "Muon-";
  ParticleColor[-13]=kOrange+1;  ParticleWidth[-13]=2;   ParticleName[-13] = "Muon+";
  ParticleColor[211]=kOrange-1; ParticleWidth[211]=2;  ParticleName[211]= "Pi+";
  ParticleColor[0]=kCyan+1;     ParticleWidth[0]=1;    ParticleName[0] = "Cherenkov photon"; //Indeed this is an optical photon, but I changed the definition
  ParticleColor[9999]=kRed-7;     ParticleWidth[9999]=1;    ParticleName[9999] = "Scintillation photon"; //Created by me, PDG number doesn't actually exist

  //Representation plane
  hxyplane["start"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);
  hxyplane["Cerenkov"] = new TH2F("hxyplane_cher","Track intersections with XY plane: Cherenkov",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);
  hxyplane["Scintillation"] = new TH2F("hxyplane_scint","Track intersections with XY plane: Scintillation",1000,intersection_zplane[0],intersection_zplane[1],1000,intersection_zplane[0],intersection_zplane[1]);

  //DAQ event
  hPmtTime = new TH1F("hPmtTime", "hPmtTime", 500, -500., 500.);
  chargeVsR = new TH1F("chargeVsR", "chargeVsR", 3, 0., 100.);
  chargeVsRScint = new TH1F("chargeVsRScint", "chargeVsRScint", 3, 0., 100.);
  chargeVsRCorr = new TH1F("chargeVsRCorr", "chargeVsRCorr", 3, 0., 100.);
  chargeVsR->GetYaxis()->SetLabelSize(.06);
  chargeVsRScint->GetYaxis()->SetLabelSize(.06);
  chargeVsRCorr->GetYaxis()->SetLabelSize(.06);



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


  /////////////////////////////////////
  // MC TRUTH
  /////////////////////////////////////

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
  for (std::map<std::string,TH2F*>::iterator it=hxyplane.begin();it!=hxyplane.end();it++) it->second->Reset();
  npe.clear();
  if (finalTrack<=0 || finalTrack > mc->GetMCTrackCount()) finalTrack = mc->GetMCTrackCount();
  //Load tracks
  for (int itr = initialTrack; itr < finalTrack; itr++) {

    if(debugLevel > 1) std::cout<<"  Track "<<itr<<"/"<<finalTrack-1<<std::endl;

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

      if(debugLevel > 1) std::cout<<"  |->Step "<<istep<<"/"<<nsteps-1<<" ";

      RAT::DS::MCTrackStep *step = mctrack->GetMCTrackStep(istep);
      const TVector3 endpointstep = step->GetEndpoint();
      pl_tracks.back().SetPoint(istep,endpointstep.X(),endpointstep.Y(),endpointstep.Z());

      if(debugLevel > 1) std::cout<<" -- DONE "<<std::endl;

      //Calculate intersection with XY plane
      // std::cout<<"step "<<istep<<" "<<top_pos.X()<<" "<<top_pos.Y()<<" "<<top_pos.Z()<<std::endl;
      // std::cout<<"step "<<istep<<" "<<bottom_pos.X()<<" "<<bottom_pos.Y()<<" "<<bottom_pos.Z()<<std::endl;
      if(mctrack->GetPDGCode()!=0 && mctrack->GetPDGCode()!=9999) continue; //only for OPTICAL photons

      if(debugLevel > 1) std::cout<<"    Intersectng Optical Photon "<<std::endl;

      if(bottom_pos.Z()!=-9999.){ //we haven't found the point yet
        if(endpointstep.Z()>intersection_zplane[2]){
          if(debugLevel > 1) std::cout<<"      Case 1: "<<std::endl;
          top_pos = endpointstep;
        }
        else if(top_pos.Z()!=-9999.){ //this is our guy

          if(debugLevel > 1) std::cout<<"      Case 2: "<<firststep->GetProcess()<<std::endl;

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

      if(debugLevel > 1) std::cout<<"   EventDisplay::LoadEvent (Passed intersection) "<<std::endl;

    } //end step loop
  } //end track loop

  //Load photoelectrons
  for (int ipmt = 0; ipmt < EDGeo->GetPMTCount(); ipmt++){
    npe[ipmt]=0;
  }
  for(int ipmt=0; ipmt<mc->GetMCPMTCount(); ipmt++){
    RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(ipmt);
    int pmtID = mcpmt->GetID();
    npe[pmtID] = mcpmt->GetMCPhotonCount();
  }

  if(drawPMTs){
    EDGeo->ResetHitPMTs();
    //Highlight PMT if was hit
    for(std::map<int, int>::iterator itpmt=npe.begin(); itpmt!=npe.end(); itpmt++){
      EDGeo->HitPMT(itpmt->first,itpmt->second);
    }
  }


#ifdef __WAVEFORMS_IN_DS__

  MCPMTWaveforms.resize(mc->GetMCPMTCount());
  MCPMTDigitizedWaveforms.resize(mc->GetMCPMTCount());
  double ymin=99999.; //yaxis min limit analogue
  UShort_t ymax_d=0.; //yaxis max limit digital
  UShort_t ymin_d=99999.; //yaxis min limit digital
  UShort_t ymax_d1=0.; //yaxis max limit digital
  UShort_t ymin_d1=99999.; //yaxis min limit digital
  UShort_t ymax_d2=0.; //yaxis max limit digital
  UShort_t ymin_d2=99999.; //yaxis min limit digital
  UShort_t ymax_d3=0.; //yaxis max limit digital
  UShort_t ymin_d3=99999.; //yaxis min limit digital
  double ymax_temp=0.;
  double ymin_temp=0.;
  double xmax_temp=0.;//dummy
  double xmin_temp=0.;//dummy

  for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {

    RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(ipmt);

    //Set analogue graphs
    vMCPMTWaveforms[ipmt] = mcpmt->GetWaveform();
    for(int ipoint=0; ipoint<vMCPMTWaveforms[ipmt].size(); ipoint++){
      //      std::cout<<" EventDisplay::LoadEvent - Waveform "<<itime<<" "<<vPMTWaveforms[ipmt][itime]<<std::endl;
      double itime = ipoint*0.1;
      MCPMTWaveforms[ipmt].SetPoint(ipoint,itime,vMCPMTWaveforms[ipmt][ipoint]);
      ymin = TMath::Min(ymin,vMCPMTWaveforms[ipmt][ipoint]);
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

#endif

  /////////////////////////////////////
  // DAQ EVENT
  /////////////////////////////////////

  if(!rds->ExistEV()) {
    if(debugLevel > 0) std::cout<<" EventDisplay::LoadEvent - DONE (EV do not exist) "<<std::endl;
    return;
  }

  TTree *runT = dsreader->GetRunT();
  RAT::DS::Run *run = 0;
  runT->SetBranchAddress("run",&run);
  runT->GetEntry(0);
  pmtInfo = run->GetPMTInfo();
  RAT::DS::DAQHeader *daqHeaderV1730 = run->GetDAQHeader("V1730");
  RAT::DS::DAQHeader *daqHeaderV1742 = run->GetDAQHeader("V1742");

  timeVsPos = new TH2F("timeVsPos", "timeVsPos", 1, 0., 1., 1, 0., 1.);
  timeVsPos->SetStats(0);
  chargeVsPos = new TH2F("chargeVsPos", "chargeVsPos", 1, 0., 1., 1, 0., 1.);
  chargeVsPos->SetStats(0);
  chargeVsPosScint = new TH2F("chargeVsPosScint", "chargeVsPosScint", 1, 0., 1., 1, 0., 1.);
  chargeVsPosScint->SetStats(0);
  chargeVsPosCorr = new TH2F("chargeVsPosCorr", "chargeVsPosCorr", 1, 0., 1., 1, 0., 1.);
  chargeVsPosCorr->SetStats(0);
  chargeVsR->Reset();
  chargeVsRScint->Reset();
  chargeVsRCorr->Reset();

  //Get center of the samll PMT array and locate
  TVector3 centerpos(0,0,0);
  int pmtTypeCount = 0;

  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    if(pmtInfo->GetType(ipmt) == 1) {
      centerpos = centerpos + pmtpos;
      pmtTypeCount++;
    }
  }
  centerpos = centerpos*(1./pmtTypeCount);
  timeVsPos->SetBins(21, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 21, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  chargeVsPos->SetBins(21, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 21, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  chargeVsPosScint->SetBins(21, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 21, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);
  chargeVsPosCorr->SetBins(21, centerpos.X()-30.*3.5, centerpos.X()+30.*3.5, 21, centerpos.Y()-30.*3.5, centerpos.Y()+30.*3.5);

  //Fill with zeroes
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    int pmtID = ev->GetPMT(ipmt)->GetID();
    //    if(pmtInfo->GetType(pmtID)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    timeVsPos->Fill(pmtpos.X(), pmtpos.Y(), 0.);
    chargeVsPos->Fill(pmtpos.X(), pmtpos.Y(), 0.);
    chargeVsPosCorr->Fill(pmtpos.X(), pmtpos.Y(), 0.);
    pmtCharge[ipmt] = 0.;
    pmtTime[ipmt] = 0.;
  }

  //Collect charges and times
  EDGeo->ResetHitPMTs();
  for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {
    int pmtID = ev->GetPMT(ipmt)->GetID();
    pmtCharge[pmtID] = ev->GetPMT(ipmt)->GetCharge();
    pmtTime[pmtID] = ev->GetPMT(ipmt)->GetTime();
    if(pmtTime[pmtID]<0){
      pmtTime[pmtID] = -400.;
    } else{
      EDGeo->HitPMT(pmtID,1); //If time>0 means that the WF crossed threshold
    }
    hPmtTime->Fill(pmtTime[pmtID]);
  }

  //Centroid
  centroid = ev->GetCentroid()->GetPosition();

  //Fill 2D plot
  for (int ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
    if(pmtInfo->GetType(ipmt)!=1) continue;
    TVector3 pmtpos = pmtInfo->GetPosition(ipmt);
    // PMT distance
    double XYdist = pow( (pmtpos.X()-centerpos.X()),2 ) + pow( (pmtpos.Y()-centerpos.Y()),2 );
    XYdist = sqrt(XYdist);
    // double dist = pow(70.,2) + pow(XYdist,2);
    // dist = sqrt(dist);

    // Geometry PMT charge correction
    // chargeVsPos->Fill(pmtpos.X(), pmtpos.Y(), pmtCharge[pmtID]*pow(dist,2) );
    timeVsPos->Fill(pmtpos.X(), pmtpos.Y(), pmtTime[ipmt]);
    chargeVsPos->Fill(pmtpos.X(), pmtpos.Y(), pmtCharge[ipmt]);
    chargeVsPosScint->Fill(pmtpos.X(), pmtpos.Y(), pmtGeoCorr[ipmt]);

    chargeVsPosCorr->Fill(pmtpos.X(), pmtpos.Y(), (pmtCharge[ipmt] - pmtGeoCorr[ipmt])/pmtGeoCorrErr[ipmt] );
    chargeVsR->Fill(XYdist, pmtCharge[ipmt]);
    chargeVsRScint->Fill(XYdist, pmtGeoCorr[ipmt]);
    chargeVsRCorr->Fill(XYdist, (pmtCharge[ipmt] - pmtGeoCorr[ipmt])/pmtGeoCorrErr[ipmt] );

  }

  // chargeVsR->Scale(1./4.); //PMT average
  // chargeVsRScint->Scale(1./4.); //PMT average
  // chargeVsRCorr->Scale(1./4.); //PMT average

  // chargeVsR->Scale(1./chargeVsR->Integral()); //Area norm
  // chargeVsRScint->Scale(1./chargeVsRScint->Integral()); //Area norm
  // chargeVsRCorr->Scale(1./chargeVsRCorr->Integral()); //Area norm

  //Perform KS test against scintillation light only
  double ks = chargeVsPos->KolmogorovTest(chargeVsPosScint);
  std::cout<<" KS 1 "<<ks<<std::endl;
  ks = chargeVsR->KolmogorovTest(chargeVsRScint);
  std::cout<<" KS 2 "<<ks<<std::endl;
  double chi2 = 0;
  for(int ibin=1; ibin<4; ibin++){
    //    chi2 += pow(chargeVsR->GetBinContent(ibin) - chargeVsRScint->GetBinContent(ibin),2)/chargeVsRScint->GetBinContent(ibin);
    chi2 += pow(chargeVsRCorr->GetBinContent(ibin),2);
  }
  std::cout<<" CHI2 "<<chi2<<std::endl;


#ifdef __WAVEFORMS_IN_DS__

  PMTDigitizedWaveforms.resize(ev->GetPMTCount());
  for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {

    RAT::DS::PMT *pmt = ev->GetPMT(ipmt);
    int pmtType = pmtInfo->GetType(ipmt);

    double timeStep = 0;
    double timeDelay = 0;

    if(pmtType==2){
      timeStep = daqHeaderV1730->GetDoubleAttribute("TIME_RES");
      timeDelay = daqHeaderV1730->GetDoubleAttribute("TIME_DELAY");
    } else if(pmtType==1 || pmtType==3){
      timeStep = daqHeaderV1742->GetDoubleAttribute("TIME_RES");
      timeDelay = daqHeaderV1742->GetDoubleAttribute("TIME_DELAY");
    }

    vPMTDigitizedWaveforms[ipmt] = pmt->GetWaveform();
    if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - DigitWF: PMT " << ipmt<<" nsamples "<<vPMTDigitizedWaveforms[ipmt].size()<<std::endl;

    for(int isample=0; isample<vPMTDigitizedWaveforms[ipmt].size(); isample++){
      PMTDigitizedWaveforms[ipmt].SetPoint(isample,isample*timeStep-timeDelay,vPMTDigitizedWaveforms[ipmt][isample]);

      if(debugLevel > 1) std::cout<<" EventDisplay::LoadEvent - Digit WF: sample "<<isample<<" "<<vPMTDigitizedWaveforms[ipmt][isample]<<std::endl;

      if(pmtType==1){
        ymax_d1 = TMath::Max(ymax_d1,vPMTDigitizedWaveforms[ipmt][isample]);
        ymin_d1 = TMath::Min(ymin_d1,vPMTDigitizedWaveforms[ipmt][isample]);
      } else if(pmtType==2){
        ymax_d2 = TMath::Max(ymax_d2,vPMTDigitizedWaveforms[ipmt][isample]);
        ymin_d2 = TMath::Min(ymin_d2,vPMTDigitizedWaveforms[ipmt][isample]);
      } else if(pmtType==3){
        ymax_d3 = TMath::Max(ymax_d3,vPMTDigitizedWaveforms[ipmt][isample]);
        ymin_d3 = TMath::Min(ymin_d3,vPMTDigitizedWaveforms[ipmt][isample]);
      }
      //Temporary
      ymin_d2 = 14400.;
    }
  }

  //Set correct limits for drawing purposes
  for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {
    int pmtType = pmtInfo->GetType(ipmt);
    if(pmtType==1){
      PMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.99*ymin_d1,1.01*ymax_d1);
    } else if(pmtType==2){
      PMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.999*ymin_d2,1.001*ymax_d2);
    } else if(pmtType==3){
      PMTDigitizedWaveforms[ipmt].GetYaxis()->SetRangeUser(0.99*ymin_d3,1.01*ymax_d3);
    }
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
  if(inputFileName == "") inputFileName = dbED->GetS("input_file");

  //Geometry files
  geoFileName = dbED->GetS("geo_file");

  //Geometry correction
  corrFileName = dbED->GetS("corr_file");
  targetMaterial = dbED->GetS("material");
  db->Load(corrFileName);

  //Validate parameters
  if(event_number<-1) std::cout<<" EventDisplay >>> Event by event mode (Event navigation disabled) "<<std::endl;
  if(!fexists(geoFileName.c_str())) {std::cout<<" EventDisplay >>> "<<geoFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(!fexists(inputFileName.c_str())) {std::cout<<" EventDisplay >>> "<<inputFileName<<" doesn't exist. Exit now!"<<std::endl; exit(0);}
  if(drawGeometry) std::cout<<" EventDisplay >>> Draw geometry in "<<geoFileName<<std::endl;
  else std::cout<<" EventDisplay >>> Draw geometry disabled "<<std::endl;
  if(drawPMTs) std::cout<<" EventDisplay >>> Draw PMTs from the PMTInfo header "<<std::endl;
  else std::cout<<" EventDisplay >>> Draw PMTs disabled "<<std::endl;

  //Geometry PMT correction
  dbCorr = db->GetLink("SCINTCORR",targetMaterial);
  pmtGeoCorr = dbCorr->GetDArray("corr");
  pmtGeoCorrErr = dbCorr->GetDArray("corr_err");

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

  //Custom palette
  Double_t r[]    = {1.0, 1.0, 0.0};
  Double_t g[]    = {0.0, 1.0, 0.0};
  Double_t b[]    = {0.0, 1.0, 1.0};
  Double_t stop[] = {0.0, .55, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(3, stop, r, g, b, 100);

  this->LoadEvent(ievt);
  if(debugLevel > 0) std::cout<<"Starting canvases "<<std::endl;
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
  if(drawPMTs) {
    EDGeo->DrawPMTMap(npe);
//    EDGeo->DrawPMTMap(pmtCharge);
  }

// #ifdef __WAVEFORMS_IN_DS__
//   //MC Analogue Waveforms
//   if(debugLevel > 0) std::cout<<"Display canvas 3 "<<std::endl;
//
//   canvas_event->cd(3);
//   //Analogue waveforms
//   for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
//     if(ipmt==0){
//       MCPMTWaveforms[ipmt].Draw("AP");
//       MCPMTWaveforms[ipmt].GetXaxis()->SetTitle("t(ns)");
//       MCPMTWaveforms[ipmt].GetYaxis()->SetTitle("V");
//     }
//     MCPMTWaveforms[ipmt].SetLineColor(ipmt+1);
//     MCPMTWaveforms[ipmt].Draw("LINE same");
//     //      PMTDigitizedWaveforms[ipmt].SetLineColor(kRed);
//     //      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
//   }
//
// #endif

  //Charge vs position
  if(rds->ExistEV()){
    if(debugLevel > 0) std::cout<<"Display canvas 3 and 9 "<<std::endl;
    // canvas_event->cd(3);
    // chargeVsPosCorr->GetZaxis()->SetRangeUser(-chargeVsPos->GetMaximum(), chargeVsPos->GetMaximum());
    // chargeVsPosCorr->Draw("colz");
    // TMarker marker;
    // marker.SetMarkerSize(3.);
    // marker.SetMarkerColor(kRed);
    // marker.SetMarkerStyle(30);
    // marker.DrawMarker(centroid.X(),centroid.Y());
    //
    // canvas_event->cd(6);
    // chargeVsR->SetLineWidth(3);
    // chargeVsRScint->SetLineColor(kRed);
    // chargeVsRScint->SetLineWidth(2);
    // chargeVsRScint->SetMinimum(0);
    // //    chargeVsRScint->SetMaximum(chargeVsRScint->GetMaximum()*1.5);
    // chargeVsRScint->Draw("");
    // chargeVsR->Draw("same");
    // // chargeVsRCorr->SetLineWidth(3);
    // // chargeVsRCorr->Draw("");

    //  gStyle->SetPalette(55);
    if(debugLevel > 0) std::cout<<"Display canvas 3 and 9"<<std::endl;
    //Charge and time
    canvas_event->cd(3);
    double rightmax = 1.1*timeVsPos->GetMaximum();
    double scale = gPad->GetUymax()/rightmax;
    if(scale > 0.) {
      timeVsPos->Scale(scale);
    }
    chargeVsPos->SetLineColor(kBlack);
    chargeVsPos->SetFillColor(kBlack);
    // chargeVsPos->SetFillStyle(3002);
    chargeVsPos->Draw("colz");
    chargeVsPos->SetMaximum(TMath::Max( chargeVsPos->GetMaximum(),chargeVsPos->GetMinimum() ) );
    chargeVsPos->SetMinimum(-chargeVsPos->GetMaximum());
    // timeVsPos->Draw("lego same");

    canvas_event->cd(9);
    chargeVsR->SetLineWidth(3);
    chargeVsR->Draw("");
  }

  if(debugLevel > 0) std::cout<<"Display canvas 5, 6 & 7"<<std::endl;
  //Digitized waveforms
  // for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
  //   if(ipmt==0){
  //     MCPMTDigitizedWaveforms[ipmt].SetTitle("MC event");
  //     MCPMTDigitizedWaveforms[ipmt].Draw("AP");
  //     MCPMTDigitizedWaveforms[ipmt].GetXaxis()->SetTitle("sample");
  //     MCPMTDigitizedWaveforms[ipmt].GetYaxis()->SetTitle("ADC counts");
  //   }
  //   MCPMTDigitizedWaveforms[ipmt].SetLineColor(ipmt+1);
  //   MCPMTDigitizedWaveforms[ipmt].Draw("LINE same");
  //   //      MCPMTDigitizedWaveforms[ipmt].ComputeRange(xmin_temp,xmax_temp,ymin_temp,ymax_temp);
  // }
  if(rds->ExistEV()){
    bool drawRingPMTs=false, drawLightPMTs=false, drawMuonPMTs=false;
    for (int ipmt = 0; ipmt < ev->GetPMTCount(); ipmt++) {
      int pmtID = ev->GetPMT(ipmt)->GetID();
      int pmtType = pmtInfo->GetType(ipmt);
      // std::cout<<" PMT "<<pmtID<<" "<<pmtType<<std::endl;
      if(pmtType == 1) canvas_event->cd(5);
      if(pmtType == 2) canvas_event->cd(6);
      if(pmtType == 3) canvas_event->cd(7);
      if(pmtType == 1 && drawRingPMTs==false || pmtType == 2 && drawLightPMTs==false || pmtType == 3 && drawMuonPMTs==false){
        PMTDigitizedWaveforms[ipmt].Draw("A LINE");
        PMTDigitizedWaveforms[ipmt].GetXaxis()->SetTitle("t(ns)");
        PMTDigitizedWaveforms[ipmt].GetYaxis()->SetTitle("ADC counts");
        if(pmtType == 1) {
          PMTDigitizedWaveforms[ipmt].SetTitle("Triggered event - Ring Tubes");
          drawRingPMTs=true;
        }
        if(pmtType == 2) {
          PMTDigitizedWaveforms[ipmt].SetTitle("Triggered event - Light Tubes");
          drawLightPMTs=true;
        }
        if(pmtType == 3) {
          PMTDigitizedWaveforms[ipmt].SetTitle("Triggered event - Muon Tags");
          drawMuonPMTs=true;
        }
      }
      PMTDigitizedWaveforms[ipmt].SetLineColor(ipmt+20);
      PMTDigitizedWaveforms[ipmt].Draw("LINE same");
      TMarker timeMarker;
      timeMarker.SetMarkerSize(1.);
      timeMarker.SetMarkerColor(ipmt+20);
      timeMarker.SetMarkerStyle(23);
      timeMarker.DrawMarker(pmtTime[pmtID], PMTDigitizedWaveforms[ipmt].GetYaxis()->GetXmax());

      //      PMTDigitizedWaveforms[ipmt].ComputeRange(xmin_temp,xmax_temp,ymin_temp,ymax_temp);
    }
  }

  // if(debugLevel > 0) std::cout<<"Display canvas 10 "<<std::endl;
  // canvas_event->cd(10);
  // hPmtTime->Draw();
  ////////////////////

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
