#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <RAT/PMTPulse.hh>
#include <RAT/Log.hh>
#include <TMath.h>


namespace RAT {


PMTPulse::PMTPulse()
{
}

PMTPulse::~PMTPulse()
{
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
      //      height = fPulseOffset - (fPulseCharge/(delta_t*fPulseWidth*sqrt(2*3.14159)))*exp(-0.5*pow(((log(delta_t)-fPulseMean)/fPulseWidth),2));
      height = fPulseOffset - (fPulseCharge*TMath::LogNormal(delta_t, fPulseWidth, 0., fPulseMean) );
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
    double EndTime = fStartTime + fPulseMean * exp(- fPulseWidth*fPulseWidth);
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
      //Integrate lognormal
      double delta_t1 = (time1-fStartTime);
      double delta_t2 = (time2-fStartTime);

      if (delta_t1 > 0. && delta_t2 > 0.) {

        double xa = log(delta_t1/fPulseMean)/(sqrt(2.)*fPulseWidth);
        double xb = log(delta_t2/fPulseMean)/(sqrt(2.)*fPulseWidth);
        totalQ = 1./2.*(TMath::Erf(xb) - TMath::Erf(xa));
        totalQ *= fPulseCharge;

      }
      else{
        totalQ = 1.0E-50; //non-zero at start time
      }

      //      std::cout<<" RealPMTPulse::Integrate: time1 "<<time1<<" time2 "<<time2<<" charge "<<totalQ<<std::endl;

    }
    return (-1.)*totalQ; //Negative pulses
}

} // namespace RAT
