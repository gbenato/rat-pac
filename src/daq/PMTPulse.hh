////////////////////////////////////////////////////////////////////
#ifndef __RAT_PMTPulse__
#define __RAT_PMTPulse__

#include <vector>

namespace RAT {

class PMTPulse {
public:

  PMTPulse();
  virtual ~PMTPulse();

//  virtual void SetPulseStartTime(double time){fStartTime = time;};
  virtual void SetPulseWidth(double width){fPulseWidth=width;};
  virtual void SetPulseCharge(double charge)=0;
  virtual void SetPulseMean(double mean){fPulseMean = mean;};
  virtual void SetPulseOffset(double offset){fPulseOffset=offset;};
  virtual void SetPulseStartTime(double time)=0;
  virtual void SetStepTime(double step){fStepTime=step;};
  virtual void SetPulseMin(double min){fPulseMin=min;};

//  virtual double GetPulseStartTime()=0;
  virtual double GetPulseStartTime() const {return fStartTime;};
  virtual double GetPulseWidth(){return fPulseWidth;};
  virtual double GetPulseCharge()=0;
  virtual double GetPulseHeight(double time)=0;
  virtual double GetPulseEndTime()=0;
  virtual double Integrate(double time1, double time2)=0;

protected:

  double fStartTime;
  double fPulseWidth;
  double fPulseMean;
  double fPulseOffset;
  double fStepTime;
  double fPulseMin;

};

class RealPMTPulse : public PMTPulse {

public:
  RealPMTPulse(){};
  virtual ~RealPMTPulse(){};

  virtual void SetPulseCharge(double _fPulseCharge){fPulseCharge = _fPulseCharge;};
  virtual void SetPulseStartTime(double _fStartTime);

  virtual double GetPulseCharge(){return fPulseCharge;};
  virtual double GetPulseHeight(double time);
  virtual double GetPulseEndTime(){return fEndTime;};
  virtual double Integrate(double time1, double time2);

private:

//  double fStartTime;
  double fPulseCharge;
  double fEndTime;

};

} // namespace RAT

#endif
