///////////////////////////////////////////
// Select events or remove branches
///////////////////////////////////////////
#include<iostream>
#include<fstream>

#include<TTree.h>
#include<TFile.h>

#include<RAT/DS/Root.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DB.hh>

#define NLOGENTRIES 10

unsigned int nfiles;
char **filenames;
unsigned int fSelection;
char* selectionTypes[] = {
  "Cosmics",
  "Cosmics-Water",
  "Cosmics-LAB",
  "Cosmics-LABPPO",
  "Cosmics-WBLS1%",
  "Cosmics-WBLS5%",
  "Cosmics-WBLS10%",
  "Through-going-muons",
  "Source"
};
vector<double> spe; //SPE

void ParseArgs(int argc, char **argv);
void GetSPETable();

int main(int argc, char **argv){

  //  gSystem->Load("$ROOTSYS/test/libEvent");

  ParseArgs(argc, argv);
  GetSPETable();

  std::cout<<"Applying "<<selectionTypes[fSelection]<<" selection to "<<std::endl;
  std::cout<<nfiles<<" files selected "<<std::endl;
  
  for(int ifile=0; ifile<nfiles; ifile++){

  std::string filename = std::string(filenames[ifile]);
  
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile(filename.c_str(),"OPEN");
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
  std::size_t pos1 = filename.find_last_of(".");
  std::string newfilename = filename.substr(0,pos1) + "_cut_" + selectionTypes[fSelection] +".root";

  std::cout<<" Cutting "<<filename<<" to \n"
	   <<"    "<<newfilename<<" ... "<<std::endl;

  TFile *newfile = new TFile(newfilename.c_str(),"recreate");
  TTree *newtree = oldtree->CloneTree(0); //clone an empty tree
  TTree *newruntree = oldruntree->CloneTree(); //clone the whole tree
  //  TTree *newtree = oldtree->CloneTree(); //clone the whole tree

  for (Long64_t ient=0; ient<nentries; ient++) {

    if(nentries>NLOGENTRIES && ient%(nentries/NLOGENTRIES) == 0) std::cout<<" Entry "<<ient<<std::endl;
    oldtree->GetEntry(ient);

    for(int iev=0; iev<ds->GetEVCount(); iev++){
      RAT::DS::EV *ev = ds->GetEV(iev);

      double bottommuon_fft = 0., topmuon_fft = 0.;
      double topmuon_charge = 0.;
      double topmuon_time = 0.;
      double bottommuon_charge = 0.;
      double bottommuon_time = -9999.;
      double panel_charge[] = {0., 0., 0., 0.};
      RAT::DS::PMT *pmt = ev->GetPMTWithID(23);
      if(pmt!=NULL) {
        topmuon_charge = pmt->GetCharge();
        topmuon_time = pmt->GetTime();
        topmuon_fft = pmt->GetFCN();
      }
      pmt = ev->GetPMTWithID(24);
      if(pmt!=NULL) {
        bottommuon_charge = pmt->GetCharge();
        bottommuon_time = pmt->GetTime();
        bottommuon_fft = pmt->GetFCN();
      }
      pmt = ev->GetPMTWithID(7);
      if(pmt!=NULL) panel_charge[0] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(8);
      if(pmt!=NULL) panel_charge[1] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(9);
      if(pmt!=NULL) panel_charge[2] = pmt->GetCharge();
      pmt = ev->GetPMTWithID(10);
      if(pmt!=NULL) panel_charge[3] = pmt->GetCharge();

      int ringPMTs = 0, lightPMTs = 0;
      double qinner = 0., npeinner = 0.;
      double qinner_sh = 0., npeinner_sh = 0.;
      for(int ipmt=0; ipmt<ev->GetPMTCount(); ipmt++){
        int pmtid = ev->GetPMT(ipmt)->GetID();
        int pmttype = pmtInfo->GetType(pmtid);
        if(pmttype==1) {
          ringPMTs++;
          if(pmtid == 13 || pmtid == 14 || pmtid == 19 || pmtid == 20 ) {
            qinner += ev->GetPMT(ipmt)->GetCharge();
            qinner_sh += ev->GetPMT(ipmt)->GetQShort();
            npeinner += ev->GetPMT(ipmt)->GetCharge()/spe[pmtid];
            npeinner_sh += ev->GetPMT(ipmt)->GetQShort()/spe[pmtid];
          }
        }
        if(pmttype==2) {
          lightPMTs++;
        }
      }

      //Cuts
      //Basics
      //if(lightPMTs==0) continue; //More than 3 hits
      //if(ringPMTs==0) continue; //More than 3 hits

      //Cherenkov imaging
      switch(fSelection){
      case 0: //Cosmics Basic Cuts
        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	break;	

      case 1: //Cosmics Water

        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	//Cosmic water specific
        if(npeinner > 40) continue; //60
	break;
	
      case 2: //Cosmics LAB

        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	//Cosmic LAB specific
        if(npeinner > 50) continue; //60
	break;

      case 3: //Cosmics LABPPO
	
        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

        //Cosmic LABPPO Specific
        if(npeinner > 500) continue; //ring tubes q cut
        if(bottommuon_fft<18e3 || topmuon_fft<18e3) continue;
	break;

      case 4: //Cosmic WBLS 1%

	//Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	//Cosmic WBLS1% specific
        if(npeinner > 50) continue;
	break;

      case 5: //Cosmic WBLS 5%

        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	//Cosmic WBLS5% specific
        if(npeinner > 80) continue;
	break;

      case 6: //Cosmic WBLS 10%

        //Top & bottom: data cleaning
        if(bottommuon_charge<50.0) continue; //200
        if(topmuon_charge<50.0) continue; //200
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue; //500, 1000, 1100
        if(panel_charge[0]<2000.) continue;

	//Cosmic WBLS10% specific
        if(npeinner > 100) continue;
	break;

      case 7: //Through-going muons
        if(bottommuon_charge<200.0) continue;
        if(panel_charge[0]<4e3) continue;
        if(panel_charge[0]>2e4) continue;
        if(panel_charge[1]>500. || panel_charge[2]>1000. || panel_charge[3]>1100.) continue;
	break;

      case 8: //Source data
	if(panel_charge[0]>50 || panel_charge[1]>50 || panel_charge[2]>50 || panel_charge[3]>50) continue;
	//if(panel_charge[0]<500) continue;
	break;

      default:

	std::cout<<" Specify a valid option"<<std::endl;
	exit(0);

        //Rest
        // if(panel_charge[0]>10000.) continue; //Veto cut
        // if(bottommuon_charge<8000.0 || topmuon_charge<8000.0) continue; //High charge cut
        // if(bottommuon_time - topmuon_time < 0.) continue;
        // if(bottommuon_time - topmuon_time > 3.) continue;
        // if(npeinner > 600) continue; //qinner cut LABPPO
        // if(npeinner_sh > 50) continue; //qshort cut LAB
        // if(qinner > 1000) continue; //inner ring tubes q cut
        // if(ev->GetClockTime() < 140e12) continue;

      }

      //Fill tree
      newtree->Fill();

    }

  }

  
  //newruntree->Print();
  //newtree->Print();
  std::cout<<"Total number of selected events: "<<newtree->GetEntries()<<std::endl;
  newfile->Write();
  delete oldfile;
  delete newfile;

  }
  
}

void GetSPETable(){

  RAT::DB* db = RAT::DB::Get();
  db->Load("../data/PMTGAUSCHARGE.ratdb");
  RAT::DBLinkPtr dbSPE = db->GetLink("PMTGAUSCHARGE");
  spe = dbSPE->GetDArray("gaus_mean");

}

void ParseArgs(int argc, char **argv){

  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-i") {filenames = argv + ++i;}
    if(std::string(argv[i]) == "-o") {fSelection = std::atoi(argv[++i]);}
  }

  nfiles = 0;
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-i") {
      for(int ifile = i+1; ifile < argc; ifile++){
	nfiles++;
      }
      break;
    }
  }
  
  if(argc<=1){
    std::cout<<" Usage: ./Selection.exe -o SELECTION -i INPUTFILES "<<std::endl;
    std::cout<<" Posible Selections: "<<std::endl;
    for(int isel=0; isel<9; isel++){
      std::cout<<isel<<": "<<selectionTypes[isel]<<std::endl;
    }
    std::cout<<std::endl;
    exit(0);
  }

}
