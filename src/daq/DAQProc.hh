#ifndef __RAT_DAQProc__
#define __RAT_DAQProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DB.hh>
#include <CLHEP/Random/RandGeneral.h>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/Digitizer.hh>
#include <RAT/DS/RunStore.hh>

namespace RAT {


class DAQProc : public Processor {
public:
  DAQProc();
  virtual ~DAQProc() { };
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual void SetS(std::string param, std::string value);

protected:

  DS::Run *run;
  RAT::DS::PMTInfo *pmtInfo;
  bool setRun;
  DBLinkPtr fLdaq;
  int fEventCounter;

  std::vector<double> fSPECharge;
  double fSamplingTime;
  double fIntTime;
  float fPulseWidth;
  float fPulseOffset;
  float fPulseTimeStep;
  float fPulseMin;
  float fNoiseAmpl;
  double fGDelay;
  double fTriggerThreshold;
  int fPulseType;
  float fPulseMean;
  double fTriggerDelay;
  double fTriggerJitter;

  Digitizer *fDigitizerV1730;
  Digitizer *fDigitizerV1742;

  std::string fTriggerType;

};


} // namespace RAT

#endif
