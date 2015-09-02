#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <RAT/PMTPulse.hh>
#include <RAT/Log.hh>


namespace RAT {


PMTPulse::PMTPulse()
{
}

PMTPulse::~PMTPulse()
{
}

double SquarePMTPulse::GetPulseHeight(double time)
{
    double height;
    if (time-fStartTime < fPulseWidth && time >= fStartTime){
      height = fPulseCharge; //units are now DAC counts = ADC counts
    }
    else {
       height = 0;
    }
    return height;
}

void SquarePMTPulse::SetPulseStartTime(double time)
{
    fStartTime = time;
    fEndTime = fStartTime+fPulseWidth;
}

double SquarePMTPulse::Integrate(double time1, double time2)
{
    double totalQ=0.0;
    double height1=GetPulseHeight(time1);
    double height2=GetPulseHeight(time2);
    if (time1>fEndTime || time2<fStartTime){
        totalQ=0.0;
    }
    else if ((time1<fStartTime) && (time2>fEndTime)){
        totalQ=GetPulseCharge();
    }
    else{
      if (height1==height2 && (height1!=0)){
        totalQ = GetPulseHeight(time2)*(time2-time1)/(fPulseWidth);
      }
      else if (height1==0 && height2!=0){
        totalQ = height2*(time2-fStartTime)/(fPulseWidth);
      }
      else if (height1!=0 && height2==0){
        totalQ = height1*(fEndTime-time1)/(fPulseWidth);
      }
      else if (height1==0 && height2==0){
        totalQ = 0.0; //can only be completely off pulse
      }
    }
    return totalQ;
}

double RealPMTPulse::GetPulseHeight(double time)
{
    double height;
    //    double norm=fPulseCharge*exp(fPulseMean); //Orebi Gann normalization
    double delta_t = (time-fStartTime);

    // if (delta_t > 1.0) {
    //   //      height = fPulseOffset - (1.*norm/(delta_t))*exp(-0.5*pow(((log(delta_t)-fPulseMean)/fPulseWidth),2));
    //   height = fPulseOffset - (fPulseCharge/(delta_t*fPulseWidth*sqrt(2*3.14159)))*exp(-0.5*pow(((log(delta_t)-fPulseMean)/fPulseWidth),2));
    // }
    // else{
    //     height = -1.0E-33; //non-zero at start time
    // }

    if (delta_t > 0.0) {
      height = fPulseOffset - (fPulseCharge/(delta_t*fPulseWidth*sqrt(2*3.14159)))*exp(-0.5*pow(((log(delta_t)-fPulseMean)/fPulseWidth),2));
    }
    else{
      height = -1.0E-50; //non-zero at start time
    }

    //    height *= -1.0;  //we use positive-going PMT pulses...
    return height;
}

void RealPMTPulse::SetPulseStartTime(double time)
{
    fStartTime = time;
    double EndTime = fStartTime + exp(fPulseMean - fPulseWidth*fPulseWidth);
    while (GetPulseHeight(EndTime)<fPulseMin){
      EndTime+=fStepTime;
    }
    fEndTime = EndTime;
}

double RealPMTPulse::Integrate(double time1, double time2)
{
    double totalQ=0;
    //    if (time1>fEndTime || time2<fStartTime){
    if (time2<fStartTime){
        totalQ=0.0;
    }
    else if ((time1<fStartTime) && (time2>fEndTime)){
      totalQ=GetPulseCharge();
    }
    else{
      double currenttime=time1;
      while (currenttime < time2){
          totalQ+=GetPulseHeight(currenttime)*fStepTime;
          currenttime+=fStepTime;
      }
    }
    return totalQ;
}

} // namespace RAT
