#ifndef __RAT_DAQProc__
#define __RAT_DAQProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DB.hh>
#include <CLHEP/Random/RandGeneral.h>
#include <RAT/Digitizer.hh>

namespace RAT {


class DAQProc : public Processor {
public:
  DAQProc();
  virtual ~DAQProc() { };
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual void SetS(std::string param, std::string value);

protected:
  int fEventCounter;

  std::vector<double> fSPECharge;
  double fSamplingTime; ///< sampling time in ns --- this is the size of a PMT time window
  double fIntTime; ///< integration time in ns
  float fPulseWidth; ///< width of a PMT pulse in ns
  float fPulseOffset; ///< offset of a PMT pulse in mV
  float fPulseTimeStep; ///< stepping time for analogue pulse
  float fPulseMin; ///< Minimum pulse height to consider
  float fNoiseAmpl; ///< width of noise in adc counts
  double fGDelay; ///< time before discriminator fires that sampling gate opens
  double fTriggerThreshold; ///< time before discriminator fires that sampling gate opens
  int fPulseType; ///< Pulse type: 0=square pulses, 1=real pulses
  float fPulseMean; ///< mean of a PMT pulse in ns (only real pulses)

  DBLinkPtr fLdaq;

  //Digitizer
  double fOffSet;
  double fVHigh;
  double fVLow;
  double fResistance;
  int fNBits;

  Digitizer *fDigitizer;

  std::string fTriggerType;

};


} // namespace RAT

#endif
