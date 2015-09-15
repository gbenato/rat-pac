#include<iostream>
#include<string>
#include "EventDisplay.hh"

int main(int argc, char **argv){

  //Validate args
  if(argc > 2) {
    std::cout<<" Usage: ./EventDisplay [INPUTFILE] "<<std::endl;
    exit(0);
  }

  //Init
  int appargc = 0;
  char **appargv = NULL;
  EventDisplay *ed;
  if(argc==1){
    ed = new EventDisplay();
  }
  else {
    ed = new EventDisplay(argv[1]);
  }
  ed->Open();

  return 0;

}
