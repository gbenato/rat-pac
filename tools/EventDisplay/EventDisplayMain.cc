#include<iostream>
#include<string>
#include<TApplication.h>
#include "EventDisplay.hh"

int main(){

  //Init
  int appargc = 0;
  char **appargv = NULL;
  TApplication dummy_app("App", &appargc, appargv);
  EventDisplay *ed = new EventDisplay();
  ed->Open();
  
  dummy_app.Run();
  return 0;
  
}
