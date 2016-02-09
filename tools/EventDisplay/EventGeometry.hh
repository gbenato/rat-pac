#ifndef __EventGeometry__
#define __EventGeometry__

#include<iostream>
#include<string>
#include<vector>
#include<TGeoVolume.h>
#include<TGeoMaterial.h>
#include<TGeoMedium.h>
#include<TGeoBBox.h>
#include<TGeoTube.h>
#include<TGeoSphere.h>
#include<TGeoMatrix.h>
#include<TPaveText.h>

class EDGeoBox{
public:
  EDGeoBox(std::string _name){SetName(_name); volume = NULL;};
  ~EDGeoBox();

  void SetName(std::string _name){ name = _name; };
  void SetMother(std::string _mother){ mother = _mother; };
  void SetSize(std::vector<double> _size){ size = _size; };
  void SetPos(std::vector<double> _pos){ pos = _pos; };

  std::string GetMother(){ return mother;};
  std::string GetName(){ return name;};
  std::vector<double> GetPos(){ return pos;};

  void AddVolume(TGeoVolume* vworld, std::vector<double> absPos){
    TGeoTranslation *trans_local = new TGeoTranslation(pos[0] + absPos[0], pos[1] + absPos[1], pos[2] + absPos[2]);
    vworld->AddNode(GetVolume(),1,trans_local);
  };

  TGeoVolume* GetVolume(){
    if(volume==NULL){
      TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
      TGeoMedium *med = new TGeoMedium("MED",1,mat);
      TGeoBBox *b = new TGeoBBox(name.c_str(), size[0], size[1], size[2]);
      volume = new TGeoVolume(name.c_str(),b,med);
      volume->SetLineWidth(1);
      volume->SetLineColor(1);
    }
    return volume;
  };

protected:

  std::string name;
  std::string mother;
  std::vector<double> size;
  std::vector<double> pos;
  TGeoVolume *volume;

};

class EDGeoTube{
public:
  EDGeoTube(std::string _name){SetName(_name); volume = NULL;};
  ~EDGeoTube();

  void SetName(std::string _name){ name = _name; };
  void SetMother(std::string _mother){ mother = _mother; };
  void SetRadius(double _r_min, double _r_max){ r_min = _r_min; r_max = _r_max; };
  void SetHeight(double _height){ height = _height; };
  void SetPos(std::vector<double> _pos){ pos = _pos; };

  std::string GetMother(){ return mother;};
  std::string GetName(){ return name;};
  std::vector<double> GetPos(){ return pos;};

  void AddVolume(TGeoVolume* vworld, std::vector<double> absPos){
    TGeoTranslation *trans_local = new TGeoTranslation(pos[0] + absPos[0], pos[1] + absPos[1], pos[2] + absPos[2]);
    vworld->AddNode(GetVolume(),1,trans_local);
  };

  TGeoVolume* GetVolume(){
    if(volume==NULL){
      TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
      TGeoMedium *med = new TGeoMedium("MED",1,mat);
      TGeoTube *b = new TGeoTube(name.c_str(), r_min, r_max, height);
      volume = new TGeoVolume(name.c_str(),b,med);
      volume->SetLineWidth(1);
      volume->SetLineColor(1);
    }
    return volume;
  };

protected:

  std::string name;
  std::string mother;
  double r_min;
  double r_max;
  double height;
  std::vector<double> pos;
  TGeoVolume *volume;

};


class EDGeoSphere{
public:
  EDGeoSphere(std::string _name){SetName(_name); volume = NULL;};
  ~EDGeoSphere();

  void SetName(std::string _name){ name = _name; };
  void SetMother(std::string _mother){ mother = _mother; };
  void SetRadius(double _r_min, double _r_max){ r_min = _r_min; r_max = _r_max; };
  void SetPos(std::vector<double> _pos){ pos = _pos; };

  std::string GetMother(){ return mother;};
  std::string GetName(){ return name;};
  std::vector<double> GetPos(){ return pos;};

  void AddVolume(TGeoVolume* vworld, std::vector<double> absPos){
    TGeoTranslation *trans_local = new TGeoTranslation(pos[0] + absPos[0], pos[1] + absPos[1], pos[2] + absPos[2]);
    vworld->AddNode(GetVolume(),1,trans_local);
  };

  TGeoVolume* GetVolume(){
    if(volume==NULL){
      TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
      TGeoMedium *med = new TGeoMedium("MED",1,mat);
      TGeoSphere *b = new TGeoSphere(name.c_str(),r_min,r_max);
      volume = new TGeoVolume(name.c_str(),b,med);
      volume->SetLineWidth(1);
      volume->SetLineColor(1);
    }
    return volume;
  };

protected:

  std::string name;
  std::string mother;
  double r_min;
  double r_max;
  std::vector<double> pos;
  TGeoVolume *volume;

};


class EDGeoPMT{
public:
  EDGeoPMT(std::string _name){SetName(_name); volume = NULL;};
  ~EDGeoPMT();

  void SetName(std::string _name){ name = _name; };
  void SetMother(std::string _mother){ mother = _mother; };
  void SetSize(std::vector<double> _size){ size = _size; };
  void SetPos(std::vector<double> _pos){ pos = _pos; };
  void SetNPE(int _npe){ npe = _npe; };
  void SetType(int _type){ type = _type; };

  std::string GetMother(){ return mother;};
  std::string GetName(){ return name;};
  std::vector<double> GetPos(){ return pos;};
  std::vector<double> GetSize(){ return size;};
  int GetNPE(){ return npe; };
  int GetType(){ return type; };

  void AddVolume(TGeoVolume* vworld, std::vector<double> absPos){
    TGeoTranslation *trans_local = new TGeoTranslation(pos[0] + absPos[0], pos[1] + absPos[1], pos[2] + absPos[2]);
    vworld->AddNode(GetVolume(),1,trans_local);
  };

  TGeoVolume* GetVolume(){
    if(volume==NULL){
      TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
      TGeoMedium *med = new TGeoMedium("MED",1,mat);
      TGeoBBox *b;
      if(type==0 || type==1) {
        b = new TGeoBBox(name.c_str(), size[0], size[1], size[2]);
      }
      else if(type == 2){
        b = new TGeoSphere(name.c_str(), 0.0, size[0]);
      }
      else if(type == 3){
        b = new TGeoTube(name.c_str(), 0, size[0]/2, size[1]);
      }
      else if(type == 4){
        b = new TGeoBBox(name.c_str(), size[0], size[1], size[2]);
      }
      else{
        std::cout<<" EDGeoPMT::GetVolume: type "<<type<<" not defined!"<<std::endl;
        exit(0);
      }
      volume = new TGeoVolume(name.c_str(),b,med);
      volume->SetLineWidth(1);
      volume->SetLineColor(1);
    }
    return volume;
  };

protected:

  std::string name;
  std::string mother;
  std::vector<double> size;
  std::vector<double> pos;
  int type;
  int npe;
  TGeoVolume *volume;

};



class EventGeometry{
public:
  EventGeometry(std::string,std::string analysisFile = "");
  ~EventGeometry(){};

  void AddNewBox(std::string,std::string,std::vector<double>,std::vector<double>);
  void AddNewTube(std::string,std::string,std::vector<double>,double,double,double);
  void AddNewSphere(std::string,std::string,std::vector<double>,double,double);
  void AddNewPMT(std::string,std::string,double,double,double,int);
  int GetPMTCount(){return pmts.size();};

  void BuildGeometry();
  void BuildPMTMap();
  void ResetHitPMTs();
  void HitPMT(int,int);
  void DrawGeometry();
  void DrawPMTMap(std::map<int,int>);
  void DrawPMTMap(std::map<int,double>);

  std::vector<double> GetAbsolutePositionFor(std::string);

protected:

  EDGeoBox *world;
  std::vector<EDGeoBox*> boxes;
  std::vector<EDGeoTube*> tubes;
  std::vector<EDGeoSphere*> spheres;
  std::vector<EDGeoPMT*> pmts;

  TGeoManager *tgeoman;
  TGeoMaterial *mat;
  TGeoMedium *med;

  std::map<std::string,TGeoVolume*> geoVolumes;
  std::map<std::string,std::string> geoHierarchy; //mother volume
  std::vector<std::string> geoOrder;

  //PMT map
  std::vector<TPaveText> vpmtbox;

};

#endif
