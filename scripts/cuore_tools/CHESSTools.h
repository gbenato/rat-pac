#ifndef __CHESS_TOOLS__
#define __CHESS_TOOLS__

#include<TROOT.h>
#include<TCanvas.h>

//Globals
char* pmtname[26] = {"Control Channel", "Light PMT 1", "Light PMT 2", "Light PMT 3" ,"Light PMT 4", "Light PMT 5", "Light PMT 6", "North Floor Panel", "South Floor Panel", "North Side Panel", "East Side Panel", "Ring PMT 0", "Ring PMT 1", "Ring PMT 2", "Ring PMT 3", "Ring PMT 4", "Ring PMT 5", "Ring PMT 6", "Ring PMT 7", "Ring PMT 8", "Ring PMT 9", "Ring PMT 10", "Ring PMT 11" , "Top Muon Tag", "Bottom Muon Tag", "Trigger PMT" };
int pmtidtopos[] = {3,  3, 3, 3, 3, 3, 3,   3, 3, 3, 3,   0, 1, 2, 2, 1, 0, 0, 1, 2, 2, 1, 0,   3, 3,   3}; //In space
Color_t pmtidtocolor[] = {1,  1, 1, 1, 1, 1, 1, kRed, 1, 1, 1,  kOrange, kBlue, kRed, kRed, kBlue, kOrange, kOrange, kBlue, kRed, kRed, kBlue, kOrange, kBlue, kRed, 1}; //By Position (WATER)
//  Color_t pmtidtocolor[] =   {1, 1, 1, 1, 1, 1, kRed, 1, 1, 1,  kBlue, kOrange, kRed, kRed, kOrange, kBlue, kBlue, kOrange, kRed, kRed, kOrange, kBlue, kBlue, kRed, 1}; //By Position (LAB)

//Functions
void Crucify(TCanvas *canvas){

  canvas->Divide(3,4);
  canvas->cd(1)->SetPad(0./7., 3./7., 1./7., 4./7.);
  canvas->cd(2)->SetPad(1./7., 3./7., 2./7., 4./7.);
  canvas->cd(3)->SetPad(2./7., 3./7., 3./7., 4./7.);

  canvas->cd(4)->SetPad(4./7., 3./7., 5./7., 4./7.);
  canvas->cd(5)->SetPad(5./7., 3./7., 6./7., 4./7.);
  canvas->cd(6)->SetPad(6./7., 3./7., 7./7., 4./7.);

  canvas->cd(7)->SetPad(3./7., 6./7., 4./7., 7./7.);
  canvas->cd(8)->SetPad(3./7., 5./7., 4./7., 6./7.);
  canvas->cd(9)->SetPad(3./7., 4./7., 4./7., 5./7.);

  canvas->cd(10)->SetPad(3./7., 2./7., 4./7., 3./7.);
  canvas->cd(11)->SetPad(3./7., 1./7., 4./7., 2./7.);
  canvas->cd(12)->SetPad(3./7., 0./7., 4./7., 1./7.);


}

void SetCHESSStyle(){
  //Style
  gROOT->ProcessLine(".x CheSSStyle.C");
  gROOT->SetStyle("CHESS");
  gROOT->ForceStyle();
}

#endif
