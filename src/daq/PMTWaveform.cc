#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <RAT/PMTWaveform.hh>
#include <RAT/PMTPulse.hh>
#include <RAT/Log.hh>
#include <RAT/ListHelp.hh>
#include <CLHEP/Random/RandGauss.h>

namespace RAT {

PMTWaveform::PMTWaveform()
{
}

PMTWaveform::~PMTWaveform()
{
//   deepdelete_vector(fPulse);
}

inline bool Cmp_PMTPulse_TimeAscending(const PMTPulse *a,
				       const PMTPulse *b)
{
  double atime = a->GetPulseStartTime();
  double btime = b->GetPulseStartTime();
  return atime < btime;
}


// void PMTWaveform::GenerateElectronicNoise(double fNoiseAmpl)
// {
  
//   //Sort pulses in time order
//   std::sort(fPulse.begin(),fPulse.end(),Cmp_PMTPulse_TimeAscending);
//   int nsteps = (int)(fPulse.back()->GetPulseEndTime() - fPulse.front()->GetPulseStartTime())/fStepTime + 1;
//   double PulseDuty=0.0;
//   for(std::vector<PMTPulse*>::iterator ipulse = fPulse.begin();ipulse<fPulse.end();ipulse++)
//     PulseDuty += (*ipulse)->GetPulseEndTime() - (*ipulse)->GetPulseStartTime();
//   float NoiseAmpl = fNoiseAmpl/sqrt(PulseDuty/fStepTime);
//   //  int nsteps = (int)fEventTime/fStepTime;
//   fNoise.resize(nsteps);
//   for(int istep=0; istep<nsteps ;istep++){
//     fNoise[istep] = NoiseAmpl*CLHEP::RandGauss::shoot();
//   }
  
// }

  
float PMTWaveform::GetHeight(double time)
{
    float height = 0.;
    unsigned int i = 0;
    while (i<fPulse.size() && fPulse[i]->GetPulseStartTime()<=time){
        height+=fPulse[i]->GetPulseHeight(time);
        i++;
    }
    //    if(time>=fEventTime) time = fEventTime - 0.001;
    int istep = round(time/fStepTime);
    //    height+=fNoise[istep]; //add electronic noise
    //std::cout<<height<<" "<<fNoise[istep]<<" "<<istep<<std::endl;
    return height;

}

int PMTWaveform::GetNext(double time)
{
    unsigned int i = 0;
    while (i<fPulse.size() && fPulse[i]->GetPulseStartTime()<=time){
        i++;
    }
    if(i==fPulse.size())return 0;
    if(fPulse[i]->GetPulseStartTime()<=time){return 0;}
    else{return i;}
}

void PMTWaveform::SetGraph()
{

  //Sort pulses in time order
  std::sort(fPulse.begin(),fPulse.end(),Cmp_PMTPulse_TimeAscending);

  double height=0.;
  double StepTime=0.01;

  double time=0.;
  int nsteps = (int)fSamplingTime/StepTime;
  for(int istep=0; istep<=nsteps ;istep++){
    time = istep*StepTime;
    height = GetHeight(time);
    gwaveform.SetPoint(istep+1,time,height);
    //std::cout<<" istep "<<istep<<" time "<<time<<" height "<<height<<std::endl;
  }
}

} // namespace RAT
