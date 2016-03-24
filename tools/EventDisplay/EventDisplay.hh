#ifndef __EventDisplay__
#define __EventDisplay__

//c-libs
#include<iostream>
#include<fstream>
#include<string>
//ROOT-libs
#include<TGraph.h>
#include<TH2F.h>
#include<TCanvas.h>
#include<TGeoVolume.h>
#include<TGeoBBox.h>
#include<TPaveText.h>
#include<TPolyLine3D.h>
#include<TApplication.h>
//RATPAC-libs
#include<RAT/DSReader.hh>
#include<RAT/DS/Root.hh>
#include<RAT/DS/MC.hh>
#include<RAT/DS/PMTInfo.hh>

#include<RAT/DB.hh>

#include"EventGeometry.hh"

class EventDisplay {
public:
  EventDisplay(std::string _inputFileName = "");
  ~EventDisplay(){};
  void Open();
  void OpenFile(std::string);
  int GetNEvents(){return nevents;};
  void DisplayEvent(int);
  bool LoadEvent(int);
  void DumpEventInfo(int);
  void DumpDisplayInfo();
  void SetGeometry();
  bool IsCut();
  bool IsCerenkov();
  bool IsTriggered();
  bool IsPE();
  void CustomizeTrack(TPolyLine3D*,RAT::DS::MCTrack*);

  void SetParameters();
  void SetIsCut(bool value){event_cut = value;};

protected:

  //Parameters
  RAT::DBLinkPtr dbED;
  RAT::DBLinkPtr dbCorr;
  int debugLevel;
  bool drawGeometry;
  bool drawPMTs;
  bool event_cut;
  std::string geoFileName;
  std::string pmtInfoFileName;
  std::string inputFileName;
  std::string corrFileName;
  std::string targetMaterial;
  int initialTrack;
  int finalTrack;
  std::string event_option;
  int event_number;
  std::vector<int> charge_cut_pmts;
  std::vector<double> charge_cut_lower;
  std::vector<double> charge_cut_higher;
  std::vector<double> intersection_zplane;
  RAT::DS::PMTInfo *pmtInfo;
  std::vector<double> time_delay;

  TApplication *dummyApp;

  //Reader stuff
  RAT::DSReader *dsreader;
  RAT::DS::Root *rds;
  RAT::DS::MC *mc;
  RAT::DS::EV *ev;

  //Maps
  std::map<std::string,int> FirstProcessCounter; //Process name - Counter
  std::map<std::string,int> LastProcessCounter; //Process name - Counter
  std::map<int,std::string> ParticleName; //PDG code - Particle name
  std::map<std::string,int> ParticleCounter; //PDG code - Counter
  std::map<int,Color_t> ParticleColor;
  std::map<int,int> ParticleWidth;
  int nevents;
  std::map<int,bool> hitpmts;
  std::vector<TPolyLine3D> pl_tracks;
  std::vector<TGraph> MCPMTWaveforms;
  std::vector<TGraph> MCPMTDigitizedWaveforms;
  std::vector<TGraph> PMTDigitizedWaveforms;
  std::map< int, std::vector<double> > vMCPMTWaveforms;
  std::map< int, std::vector<UShort_t> > vMCPMTDigitizedWaveforms;
  std::map< int, std::vector<UShort_t> > vPMTDigitizedWaveforms;
  std::map< int, std::vector<double> > vWaveformTimes;
  std::map<std::string,TH2F*> hxyplane;
  TCanvas *canvas_event;
  double elength;
  std::vector<Color_t> pmtidtocolor;
  std::map<int, int> npe; //number of photoelectrons per PMT
  std::map<int, double> pmtCharge; //measured PMT charge
  std::map<int, double> pmtQShort;
  std::map<int, double> pmtTime; //measured PMT time
  std::map<int, double> pmtTimeCorr; //measured PMT time
  std::vector<double> spe;

  //Geometry
  EventGeometry *EDGeo;
  TGeoVolume *vworld;
  std::map<int, TGeoBBox* > bpmt;
  std::map<int, TGeoVolume* > vpmt;
  std::vector<TPaveText> vpmtbox;

  //MC event
  TH1F * hMCPeTime;

  //DAQ event
  double event_id;
  double event_time;
  std::vector<double> pmtGeoCorr; // Geometry correction
  std::vector<double> pmtGeoCorrErr; //Error
  TH1F *hTime;
  TH2F *timeVsPos;
  TH2F *chargeVsPos;
  TH2F *chargeVsPosScint;
  TH2F *chargeVsPosCorr;
  TH1F *chargeVsR;
  TH1F *chargeVsRScint;
  TH1F *chargeVsRCorr;
  TH2F *npeVsPos;

  //Plot Limits
  UShort_t ymax_d;
  UShort_t ymin_d;
  UShort_t ymax_d1;
  UShort_t ymin_d1;
  UShort_t ymax_d2;
  UShort_t ymin_d2;
  UShort_t ymax_d3;
  UShort_t ymin_d3;
  UShort_t ymax_d4;
  UShort_t ymin_d4;
  UShort_t ymax_d5;
  UShort_t ymin_d5;

};

#endif
