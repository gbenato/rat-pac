#ifndef __RAT_DAQProc__
#define __RAT_DAQProc__

#include <RAT/Processor.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DB.hh>
#include <CLHEP/Random/RandGeneral.h>

namespace RAT {


class DAQProc : public Processor {
public:
  DAQProc();
  virtual ~DAQProc() { };
  virtual Processor::Result DSEvent(DS::Root *ds);

protected:
  int fEventCounter;

  std::vector<double> fSPECharge;
  double fSamplingTimeDB; ///< sampling time in ns --- this is the size of a PMT time window
  double fIntTimeDB; ///< integration time in ns
  float fPulseWidthDB; ///< width of a PMT pulse in ns
  float fPulseOffsetDB; ///< offset of a PMT pulse in mV
  float fStepTimeDB; ///< stepping time for discrimination
  float fPulseMinDB; ///< Minimum pulse height to consider
  float fNoiseAmplDB; ///< width of noise in adc counts
  double fGDelayDB; ///< time before discriminator fires that sampling gate opens
  double fTriggerThresholdDB; ///< time before discriminator fires that sampling gate opens
  int fPulseTypeDB; ///< Pulse type: 0=square pulses, 1=real pulses
  float fPulseMeanDB; ///< mean of a PMT pulse in ns (only real pulses)

  DBLinkPtr fLdaq;

};


} // namespace RAT

#endif
