#ifndef __RAT_PMTWaveform__
#define __RAT_PMTWaveform__

#include <vector>
#include <RAT/PMTPulse.hh>
#include <TObject.h>
#include <TGraph.h>

namespace RAT {

class PMTWaveform {
public:

  PMTWaveform();
  virtual ~PMTWaveform();
  //  virtual void GenerateElectronicNoise(double);
  virtual double GetHeight(double time);
  virtual double GetCharge(double time1, double time2);
  virtual int GetNext(double time);
  virtual void SetGraph(); //Set the graph with the pulses that are stored in the moment you call this method
  virtual void SetStepTime(double step){fStepTime=step;};
  virtual double GetStepTime(){return fStepTime;};
  virtual void SetSamplingWindow(double samplingtime){fSamplingTime=samplingtime;};
  virtual TGraph GetGraph(){return gwaveform;};

  std::vector<PMTPulse*> fPulse;

protected:

  TGraph gwaveform;
  double fStepTime;
  double fSamplingTime;
  std::vector<double> fNoise;

};

} // namespace RAT

#endif
