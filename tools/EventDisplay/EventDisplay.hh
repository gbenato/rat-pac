//#define __WAVEFORMS_IN_DS__

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

#include<RAT/DB.hh>

#include"EventGeometry.hh"

class EventDisplay{
public:
  EventDisplay();
  ~EventDisplay(){};
  void Open();
  void OpenFile(std::string);
  int GetNEvents(){return nevents;};
  void DisplayEvent(int);
  void LoadEvent(int);
  void DumpEventInfo(int);
  void DumpDisplayInfo();
  void SetGeometry();
  bool IsCerenkov();
  bool IsPE();
  void CustomizeTrack(TPolyLine3D*,RAT::DS::MCTrack*);

  void SetParameters();

protected:
  
  //Parameters
  RAT::DBLinkPtr dbED;
  int debugLevel;
  bool drawGeometry;
  bool drawPMTs;
  std::string geoFileName;
  std::string pmtInfoFileName;
  std::string inputFileName;
  int initialTrack;
  int finalTrack;
  std::string event_option;
  int event_number;
  std::vector<double> intersection_zplane;

  TApplication *dummyApp;
  
  //Maps
  std::map<std::string,int> FirstProcessCounter; //Process name - Counter
  std::map<std::string,int> LastProcessCounter; //Process name - Counter
  std::map<int,std::string> ParticleName; //PDG code - Particle name
  std::map<std::string,int> ParticleCounter; //PDG code - Counter
  std::map<int,Color_t> ParticleColor;
  std::map<int,int> ParticleWidth;
  RAT::DSReader *dsreader;
  RAT::DS::Root *rds;
  RAT::DS::MC *mc;
  int nevents;
  std::map<int,bool> hitpmts;
  std::vector<TPolyLine3D> pl_tracks;
  std::vector<TGraph> PMTWaveforms;
  std::vector<TGraph> PMTDigitizedWaveforms;
  std::map< int, std::vector<double> > vPMTWaveforms;
  std::map< int, std::vector<int> > vPMTDigitizedWaveforms;
  std::map<std::string,TH2F*> hxyplane;
  TCanvas *canvas_event;

  double elength;
  std::map<int, int> npe; //number of photoelectrons per PMT

  //Geometry
  EventGeometry *EDGeo;
  TGeoVolume *vworld;
  std::map<int, TGeoBBox* > bpmt;
  std::map<int, TGeoVolume* > vpmt;
  std::vector<TPaveText> vpmtbox;

};

#endif
