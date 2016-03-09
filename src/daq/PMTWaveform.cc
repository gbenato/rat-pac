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

double PMTWaveform::GetHeight(double currenttime)
{
    float height = 0.;
    unsigned int i = 0;
    while (i<fPulse.size() && fPulse[i]->GetPulseStartTime()<=currenttime){
			height+=fPulse[i]->GetPulseHeight(currenttime);
			i++;
    }
    return height;

}

double PMTWaveform::GetCharge(double time1, double time2)
{
    double integral = 0.;
		for(int ipulse=0; ipulse<fPulse.size(); ipulse++){
			integral += fPulse[ipulse]->Integrate(time1,time2);
		}

    return integral;

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

  double time=0.;
  int nsteps = (int)fSamplingTime/fStepTime;
  for(int istep=0; istep<=nsteps ;istep++){
    time = istep*fStepTime;
    height = GetHeight(time);
    gwaveform.SetPoint(istep+1,time,height);
    //std::cout<<" istep "<<istep<<" time "<<time<<" height "<<height<<std::endl;
  }
}

} // namespace RAT
