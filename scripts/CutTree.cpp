///////////////////////////////////////////
// Reduce the tree weight by skimming the branches that the user
// consider not important. Set the cut out level by CUTLEVEL:
// 0-> Store events with at least 1 MCPMT
// 1-> Store events with at least 1 PMT (triggered event)
///////////////////////////////////////////
#include<iostream>
#include<fstream>

#include<TTree.h>
#include<TFile.h>

#include<RAT/DS/Root.hh>
#include <RAT/DS/RunStore.hh>

//Parse args
char* filename = NULL;
void ParseArgs(int argc, char **argv);

int main(int argc, char **argv){

  //  gSystem->Load("$ROOTSYS/test/libEvent");

  ParseArgs(argc, argv);

  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(filename,"OPEN");
  TTree *oldtree = (TTree*)oldfile->Get("T");
  TTree *oldruntree = (TTree*)oldfile->Get("runT");

  Long64_t nentries = oldtree->GetEntries();
  std::cout<<" Cutting out tree with "<<nentries<<" entries "<<std::endl;
  RAT::DS::Root *ds = new RAT::DS::Root; //needs to be initialized
  oldtree->SetBranchAddress("ds",&ds);
  RAT::DS::Run *run = 0;
  oldruntree->SetBranchAddress("run",&run);
  oldruntree->GetEntry(0);
  RAT::DS::PMTInfo *pmtInfo = run->GetPMTInfo();

  ////Enable or disable branches
  // oldtree->SetBranchStatus("*",0);
  // oldtree->SetBranchStatus("calib",1);
  // oldtree->SetBranchStatus("ratVersion",1);
  // oldtree->SetBranchStatus("procResult",1);
  // oldtree->SetBranchStatus("ev.pmt",1);
  // oldtree->SetBranchStatus("mc.pmt",1);

  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile(Form("%s_cut.root",filename),"recreate");
  TTree *newtree = oldtree->CloneTree(0); //clone an empty tree
  TTree *newruntree = oldruntree->CloneTree(); //clone the whole tree
  //  TTree *newtree = oldtree->CloneTree(); //clone the whole tree

  for (Long64_t ient=0; ient<nentries; ient++) {

    if(ient%10000 == 0) std::cout<<"   Entry "<<ient<<std::endl;
    oldtree->GetEntry(ient);

    for(int iev=0; iev<ds->GetEVCount(); iev++){
      RAT::DS::EV *ev = ds->GetEV(iev);

      double topmuon_charge = 0.;
      double topmuon_time = 0.;
      double bottommuon_charge = 0.;
      double bottommuon_time = -9999.;
      double panel_charge[] = {0., 0., 0., 0.};
      RAT::DS::PMT *pmt = ev->GetPMTWithID(6);
      if(pmt!=NULL) {
        topmuon_charge = pmt->GetCharge();
        topmuon_time = pmt->GetTime();
      }
      pmt = ev->GetPMTWithID(7);
      if(pmt!=NULL) {
        bottommuon_charge = pmt->GetCharge();
        bottommuon_time = pmt->GetTime();
      }
      pmt = ev->GetPMTWithID(8);
      if(pmt!=NULL) panel_charge[0] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(9);
      if(pmt!=NULL) panel_charge[1] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(10);
      if(pmt!=NULL) panel_charge[2] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(11);
      if(pmt!=NULL) panel_charge[3] = pmt->GetCharge();

      int ringPMTs = 0, lightPMTs = 0;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        if(pmttype==1) {
          ringPMTs++;
        }
        if(pmttype==2) {
          lightPMTs++;
	}
      }

      //Cuts
      //Cuts for cherenkov imaging
      if(lightPMTs==0) continue; //More than 3 hits
      //if(ringPMTs<3) continue; //More than 3 hits
      //if(bottommuon_charge<200.0) continue;
      //if(topmuon_charge<200.0) continue;
      // if(bottommuon_time - topmuon_time < 0.4) continue;
      // if(bottommuon_time - topmuon_time > 0.8 ) continue;
      // if(panel_charge[0]<2000.) continue; //Veto cut
      // if(panel_charge[1]>200. || panel_charge[2]>200. || panel_charge[3]>200.) continue;

      //Cuts for SPE
      // if(ring_time < -9000 || panel_charge[0]>200 || panel_charge[1]>200 || panel_charge[2]>200 || panel_charge[3]>200) continue;

      newtree->Fill();
    }

  }

  newruntree->Print();
  newtree->Print();
  newfile->Write();
  delete oldfile;
  delete newfile;

}

void ParseArgs(int argc, char **argv){

  if(argc < 2){
    std::cout<<" Usage: ./CutTree.exe INPUTFILE "<<std::endl;
    exit(0);
  }
  filename = argv[1];

  std::cout<<" Compressing "<<filename<<" .......... "<<std::endl;

}
