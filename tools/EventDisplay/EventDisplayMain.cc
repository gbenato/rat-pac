#include<iostream>
#include<string>
#include "EventDisplay.hh"

int main(){

  //Init
  int appargc = 0;
  char **appargv = NULL;
  EventDisplay *ed = new EventDisplay();
  ed->Open();
  
  return 0;
  
}
