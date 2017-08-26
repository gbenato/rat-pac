#include"EventGeometry.hh"
#include<RAT/DSReader.hh>
#include<RAT/DB.hh>
#include<RAT/DS/Run.hh>

#include<TGeoManager.h>
#include<TGeoMaterial.h>
#include<TGeoMedium.h>
#include<TGeoBBox.h>
#include<TGeoArb8.h>


EventGeometry::EventGeometry(std::string dbGeoFile, std::string analysisFile){

   std::cout<<" EventGeometry::EventGeometry "<<std::endl;

  bool drawPMTs = false;
//  if(analysisFile!="")
//    drawPMTs = true;

  world = new EDGeoBox("world");
  RAT::DB* db = RAT::DB::Get();
  db->Load(dbGeoFile);
  std::cout<<" EventGeometry::EventGeometry DB loaded"<<std::endl;
  std::cout<<" EventGeometry::EventGeometry "<< analysisFile.c_str() <<std::endl;

  RAT::DSReader *dsreader = new RAT::DSReader(analysisFile.c_str());
  TTree *runT = dsreader->GetRunT();
  RAT::DS::Run *run = 0;
  runT->SetBranchAddress("run",&run);
  runT->GetEntry(0);
  RAT::DS::PMTInfo *pmtInfo = run->GetPMTInfo();

  RAT::DBLinkGroup mydblink = db->GetLinkGroup("GEO");
  std::cout<<" EventGeometry::EventGeometry for loop volumes"<<std::endl;

  for(RAT::DBLinkGroup::iterator it=mydblink.begin(); it!=mydblink.end();it++){

    std::string volName = it->first;
    std::cout<<" EventGeometry::EventGeometry "<<volName<<std::endl;
    RAT::DBLinkPtr dbGeo = db->GetLink("GEO",volName.c_str());
    std::string mother = dbGeo->GetS("mother");
    std::vector<double> pos = std::vector<double>(3,0);
    try{
      pos = dbGeo->GetDArray("position");
    }
    catch(RAT::DBNotFoundError){}

    //Define world
    if( volName == "world"){
      world->SetSize(dbGeo->GetDArray("size"));
      world->SetPos(pos);
      if(volName!="") geoHierarchy[volName] = mother;
    }
    //Define volumes
    else if( dbGeo->GetS("type") == "box" ){
      std::vector<double> size = dbGeo->GetDArray("size");
      AddNewBox(volName.c_str(), mother, pos, size);
      if(volName!="") geoHierarchy[volName] = mother;
    }
    else if( dbGeo->GetS("type") == "tube" ){
      double r_max = dbGeo->GetD("r_max");
      double r_min = 0;
      double height = dbGeo->GetD("size_z");
      try{
        r_min = dbGeo->GetD("r_min");
      }
      catch(RAT::DBNotFoundError){}
      AddNewTube(volName.c_str(), mother, pos, r_min, r_max, height);
      if(volName!="") geoHierarchy[volName] = mother;
    }
    else if( dbGeo->GetS("type") == "sphere" ){
      double r_max = dbGeo->GetD("r_max");
      double r_min = 0;
      try{
        r_min = dbGeo->GetD("r_min");
      }
      catch(RAT::DBNotFoundError){}
      AddNewSphere(volName.c_str(), mother, pos, r_min, r_max);
      if(volName!="") geoHierarchy[volName] = mother;
    }
  }

  //Build PMTs from PMTInfo
  if(drawPMTs){
    for (size_t ipmt = 0; ipmt < pmtInfo->GetPMTCount(); ipmt++) {
      std::string volName = "pmts";
      std::ostringstream volNameUsed;
      volNameUsed << volName <<"_"<<ipmt;
      TVector3 pmt_pos = pmtInfo->GetPosition(ipmt);
      AddNewPMT(volNameUsed.str().c_str(), "world", pmt_pos[0], pmt_pos[1], pmt_pos[2], pmtInfo->GetType(ipmt));
      if(volNameUsed.str()!="") geoHierarchy[volNameUsed.str()] = "world";
    }
  }
  std::cout<<" EventGeometry::EventGeometry call BuildGeometry()"<<std::endl;

  BuildGeometry();
  if(drawPMTs) BuildPMTMap();

  //  std::cout<<" EventGeometry::EventGeometry - DONE "<<std::endl;

}


void EventGeometry::BuildGeometry(){


  //Order hierarchy
  geoOrder.push_back("world");
  int counter = 0;
  while(counter<=geoOrder.size()){
    for(std::map<std::string,std::string>::iterator it=geoHierarchy.begin(); it != geoHierarchy.end(); it++){
      if(it->second == geoOrder[counter] && it->second!="") geoOrder.push_back(it->first);
    }
    counter++;
  }
   std::cout<<" EventGeometry::BuildGeometry "<<std::endl;
  //Print geo order
  for(int ivol=0; ivol<geoOrder.size();ivol++)
     std::cout<<" Volume "<<ivol<<": "<<geoOrder[ivol]<<std::endl;


  //Geometry
  tgeoman = new TGeoManager("box", "poza1");
  mat = new TGeoMaterial("Al", 26.98,13,2.7);
  med = new TGeoMedium("MED",1,mat);

  geoVolumes["world"] = world->GetVolume();
  gGeoManager->SetTopVolume(geoVolumes["world"]);


  for(int ivol=1; ivol<geoOrder.size();ivol++){
    for(int ibox=0; ibox<boxes.size(); ibox++){
      if(boxes[ibox]->GetName() == geoOrder[ivol]){
        std::vector<double> trans = this->GetAbsolutePositionFor(boxes[ibox]->GetName());
        boxes[ibox]->AddVolume(geoVolumes["world"],trans);
        geoVolumes[boxes[ibox]->GetName()] = boxes[ibox]->GetVolume();
        break;
      }
    }

    for(int itube=0; itube<tubes.size(); itube++){
      if(tubes[itube]->GetName() == geoOrder[ivol]){
        std::vector<double> trans = this->GetAbsolutePositionFor(tubes[itube]->GetName());
        tubes[itube]->AddVolume(geoVolumes["world"],trans);
        geoVolumes[tubes[itube]->GetName()] = tubes[itube]->GetVolume();
        break;
      }
    }

    for(int isphere=0; isphere<spheres.size(); isphere++){
      if(spheres[isphere]->GetName() == geoOrder[ivol]){
        std::vector<double> trans = this->GetAbsolutePositionFor(spheres[isphere]->GetName());
        spheres[isphere]->AddVolume(geoVolumes["world"],trans);
        geoVolumes[spheres[isphere]->GetName()] = spheres[isphere]->GetVolume();
        break;
      }
    }

    for(int ipmt=0; ipmt<pmts.size(); ipmt++){
      if(pmts[ipmt]->GetName() == geoOrder[ivol]){
        std::vector<double> trans(3,0);
        pmts[ipmt]->AddVolume(geoVolumes["world"],trans);
        geoVolumes[pmts[ipmt]->GetName()] = pmts[ipmt]->GetVolume();
        break;
      }
    }

  }

  tgeoman->CloseGeometry();

   std::cout<<" EventGeometry::BuildGeometry - DONE "<<std::endl;

}


void EventGeometry::BuildPMTMap(){

  //PMTs in XY plane
  for(int ipmt=0;ipmt<pmts.size();ipmt++){
    if(pmts[ipmt]->GetType() == 1){
      vpmtbox.push_back(TPaveText(pmts[ipmt]->GetPos()[0] + pmts[ipmt]->GetSize()[0] , pmts[ipmt]->GetPos()[1] + pmts[ipmt]->GetSize()[1] , pmts[ipmt]->GetPos()[0]-pmts[ipmt]->GetSize()[0] , pmts[ipmt]->GetPos()[1]-pmts[ipmt]->GetSize()[1] ));
    }
    else{
      vpmtbox.push_back(TPaveText(pmts[ipmt]->GetPos()[0] + pmts[ipmt]->GetSize()[0] , pmts[ipmt]->GetPos()[1] + pmts[ipmt]->GetSize()[1] , pmts[ipmt]->GetPos()[0]-pmts[ipmt]->GetSize()[0] , pmts[ipmt]->GetPos()[1]-pmts[ipmt]->GetSize()[1] ));
    }
    //customize
    vpmtbox[ipmt].SetFillColor(kGray);
    vpmtbox[ipmt].SetFillStyle(3002);
    vpmtbox[ipmt].SetLineColor(2);
    vpmtbox[ipmt].SetLineWidth(1);
    vpmtbox[ipmt].SetTextColor(kBlack);
    vpmtbox[ipmt].SetTextSize(0.04);
  }

}


void EventGeometry::DrawPMTMap(std::map<int,int> labels){

  for(int ipmt=0;ipmt<vpmtbox.size();ipmt++){
    vpmtbox[ipmt].Clear();
    vpmtbox[ipmt].AddText(Form("%d",labels[ipmt]));
  }

  for(int ipmt=0; ipmt<vpmtbox.size();ipmt++)
    vpmtbox[ipmt].Draw("LINE same");

}

void EventGeometry::DrawPMTMap(std::map<int,double> labels){

  //Draw grid
  // int nlines = 2*XP_XSIDE/pmtwidth;
  // TLine* xline;
  // TLine* yline;
  // for(int iline=0; iline<nlines; iline++){
  //   xline = new TLine(-XP_XSIDE,iline*pmtwidth-nlines/2.*pmtwidth,XP_XSIDE,iline*pmtwidth-nlines/2.*pmtwidth);
  //   xline->SetLineWidth(1.);
  //   xline->SetLineStyle(3);
  //   xline->SetLineColor(kGray);
  //   xline->Draw("same");
  //   yline = new TLine(iline*pmtwidth-nlines/2.*pmtwidth,-XP_YSIDE,iline*pmtwidth-nlines/2.*pmtwidth,XP_YSIDE);
  //   yline->SetLineWidth(1.);
  //   yline->SetLineStyle(3);
  //   yline->SetLineColor(kGray);
  //   yline->Draw("same");
  // }

  for(int ipmt=0;ipmt<vpmtbox.size();ipmt++){
    vpmtbox[ipmt].Clear();
    vpmtbox[ipmt].AddText(Form("%.1f",labels[ipmt]));
  }

  for(int ipmt=0; ipmt<vpmtbox.size();ipmt++)
    vpmtbox[ipmt].Draw("LINE same");

}


void EventGeometry::DrawGeometry(){

  std::cout<<" EventGeometry::DrawGeometry "<<pmts.size()<<std::endl;

  tgeoman->GetMasterVolume()->Draw();
  //tgeoman->GetMasterVolume()->Draw("ogle");

}

void EventGeometry::ResetHitPMTs(){

  //Reset PMT colors
  for(int ipmt = 0; ipmt<pmts.size(); ipmt++){
      pmts[ipmt]->GetVolume()->SetLineColor(1);
  }
}


void EventGeometry::HitPMT(int id, int npe){

  //Find PMT
  TGeoVolume* hitpmt = pmts[id]->GetVolume();
  if(npe>0) hitpmt->SetLineColor(kRed);
  pmts[id]->SetNPE(npe);

}

void EventGeometry::AddNewPMT(std::string name, std::string mother, double x_pos, double y_pos, double z_pos, int type){

  std::vector<double> pos(3,0);
  pos[0] = x_pos; pos[1] = y_pos; pos[2] = z_pos;
  pmts.push_back(new EDGeoPMT(name));
  EDGeoPMT* newpmt = pmts.back();
  newpmt->SetName(name);
  newpmt->SetMother(mother);
  newpmt->SetPos(pos);
  newpmt->SetType(type);
  if(type == 0 || type == 1){
    newpmt->SetSize(std::vector<double>(3,14.5));
  }
  else if(type == 2){
    newpmt->SetSize(std::vector<double>(3,125.0));
  }
  else if(type == 3){
    newpmt->SetSize(std::vector<double>(3,14.5));
  }
  else if(type == 4){
    newpmt->SetSize(std::vector<double>(3,14.5));
  }
  else{
    std::cout<<" EDGeoPMT::AddNewPMT: type "<<type<<" not defined!"<<std::endl;
    exit(0);
  }

  newpmt->SetNPE(0);

}


void EventGeometry::AddNewTube(std::string name, std::string mother, std::vector<double> pos, double r_min, double r_max, double height){

  tubes.push_back(new EDGeoTube(name));
  EDGeoTube* newtube = tubes.back();
  newtube->SetName(name);
  newtube->SetMother(mother);
  newtube->SetPos(pos);
  newtube->SetRadius(r_min,r_max);
  newtube->SetHeight(height);

}

void EventGeometry::AddNewSphere(std::string name, std::string mother, std::vector<double> pos, double r_min, double r_max){

  spheres.push_back(new EDGeoSphere(name));
  EDGeoSphere* newsphere = spheres.back();
  newsphere->SetName(name);
  newsphere->SetMother(mother);
  newsphere->SetPos(pos);
  newsphere->SetRadius(r_min,r_max);

}


void EventGeometry::AddNewBox(std::string name, std::string mother, std::vector<double> pos, std::vector<double> size){

  boxes.push_back(new EDGeoBox(name));
  EDGeoBox* newbox = boxes.back();
  newbox->SetName(name);
  newbox->SetMother(mother);
  newbox->SetPos(pos);
  newbox->SetSize(size);

}


std::vector<double> EventGeometry::GetAbsolutePositionFor(std::string volumeName){

  std::vector<double> pos_total(3,0);

  std::string currentVolume = volumeName;

  while(geoHierarchy[currentVolume] != "world"){

    for(int ibox=0; ibox<boxes.size(); ibox++){
      if(boxes[ibox]->GetName() == geoHierarchy[currentVolume]){
        std::vector<double> pos = boxes[ibox]->GetPos();
        pos_total[0] += pos[0];
        pos_total[1] += pos[1];
        pos_total[2] += pos[2];
        break;
      }
    }

    for(int itube=0; itube<tubes.size(); itube++){
      if(tubes[itube]->GetName() == geoHierarchy[currentVolume]){
        std::vector<double> pos = tubes[itube]->GetPos();
        pos_total[0] += pos[0];
        pos_total[1] += pos[1];
        pos_total[2] += pos[2];
        break;
      }
    }

    for(int isphere=0; isphere<spheres.size(); isphere++){
      if(spheres[isphere]->GetName() == geoHierarchy[currentVolume]){
        std::vector<double> pos = spheres[isphere]->GetPos();
        pos_total[0] += pos[0];
        pos_total[1] += pos[1];
        pos_total[2] += pos[2];
        break;
      }
    }

    for(int ipmt=0; ipmt<pmts.size(); ipmt++){
      if(pmts[ipmt]->GetName() == geoHierarchy[currentVolume]){
        std::vector<double> pos = pmts[ipmt]->GetPos();
        pos_total[0] += pos[0];
        pos_total[1] += pos[1];
        pos_total[2] += pos[2];
        break;
      }
    }

    currentVolume = geoHierarchy[currentVolume];
  }

  return pos_total;

}
