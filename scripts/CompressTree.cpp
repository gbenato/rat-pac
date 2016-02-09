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

//Parse args
char* filename = NULL;
int cutlevel = 0;
void ParseArgs(int argc, char **argv);

int main(int argc, char **argv){

  //  gSystem->Load("$ROOTSYS/test/libEvent");

  ParseArgs(argc, argv);

  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(filename,"OPEN");
  TTree *oldtree = (TTree*)oldfile->Get("T");

  Long64_t nentries = oldtree->GetEntries();
  std::cout<<" Cutting out tree with "<<nentries<<" entries "<<std::endl;
  RAT::DS::Root *ds = new RAT::DS::Root; //needs to be initialized
  oldtree->SetBranchAddress("ds",&ds);

  ////Enable or disable branches
  // oldtree->SetBranchStatus("*",0);
  // oldtree->SetBranchStatus("calib",1);
  // oldtree->SetBranchStatus("ratVersion",1);
  // oldtree->SetBranchStatus("procResult",1);
  // oldtree->SetBranchStatus("ev.pmt",1);
  // oldtree->SetBranchStatus("mc.pmt",1);

  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile(Form("%s_cutlevel%d.root",filename,cutlevel),"recreate");
  TTree *newtree = oldtree->CloneTree(0); //clone an empty tree
  //  TTree *newtree = oldtree->CloneTree(); //clone the whole tree

  for (Long64_t ient=0; ient<nentries; ient++) {

    if(ient%10000 == 0) std::cout<<"   Entry "<<ient<<std::endl;
    oldtree->GetEntry(ient);

    if(cutlevel==0){

      RAT::DS::MC *mc = ds->GetMC();
      if (mc->GetMCPMTCount()==0) continue;
      newtree->Fill();

    } else if(cutlevel==1){

      for(int iev=0; iev<ds->GetEVCount(); iev++){
        RAT::DS::EV *ev = ds->GetEV(iev);
        if (ev->GetPMTCount()==0) continue;
        newtree->Fill();
      }

    }

  }

  newtree->Print();
  newfile->Write();
  delete oldfile;
  delete newfile;

}

void ParseArgs(int argc, char **argv){

  if(argc < 2){
    std::cout<<" Usage: ./CompressTree.exe INPUTFILE [CUTLEVEL (optional)] "<<std::endl;
    exit(0);
  }
  filename = argv[1];
  if(argc>2)
    cutlevel = std::stoi(argv[2]);

  std::cout<<" Skimming "<<filename<<" .......... Cutlevel = "<<cutlevel<<std::endl;

}
