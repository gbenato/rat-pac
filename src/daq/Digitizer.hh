////////////////////////////////////////////////////////////////////
/// \class RAT::Digitizer
///
/// \brief   Digitizer
///
/// \author Javier Caravaca <jcaravaca@berkeley.edu>
///
/// REVISION HISTORY:\n
///     1 Feb 2015: Initial commit
///
/// \details
/// This class provides full support for a CAEN digitizer.
/// It digitizes PMTWaveforms, check thresholds crossing,
/// integrate charge, calculate front-end times ...
////////////////////////////////////////////////////////////////////
#ifndef __RAT_Digitizer__
#define __RAT_Digitizer__

#include <map>
#include <RAT/PMTWaveform.hh>

namespace RAT {
  
  class Digitizer {
  public:
    
    Digitizer();
    virtual ~Digitizer();

    virtual void SetStepTime(double _steptime){fStepTime=_steptime;};
    virtual void SetVHigh(double _vhigh){fVhigh=_vhigh;};
    virtual void SetVLow(double _vlow){fVlow=_vlow;};
    virtual void SetOffSet(double _offset){fOffset=_offset;};
    virtual void SetNBits(int _nbits){fNBits=_nbits;};
    virtual void SetResistance(double _res){fResistance=_res;};
    virtual void SetNoiseAmplitude(double _fnoise){fNoiseAmpl=_fnoise;};
    virtual void SetSamplingWindow(double _fwindow){fSamplingWindow=_fwindow;};
    virtual void SetSampleDelay(double _fdelay){fSampleDelay=_fdelay;};
    virtual void SetThreshold(double);
    virtual void Clear();
    
    virtual void AddChannel(int,PMTWaveform);
    virtual int GetNSamples(int ich){return fDigitWaveForm[ich].size();};
    virtual void GenerateElectronicNoise(int,PMTWaveform);
    //    virtual void DigitizeWaveForm(PMTWaveform);
    virtual std::vector<int> GetDigitizedWaveform(int ich){return fDigitWaveForm[ich];};
    virtual std::vector<int> SampleWaveform(std::vector<int>, int);
    virtual int GoToEndOfSample(int);
    virtual double IntegrateCharge(std::vector<int>);
    virtual double GetDigitizedThreshold(){return fDigitizedThreshold;};
    virtual double GetTimeAtPeak(int,int);
    virtual double GetTimeAtThreshold(int,int);

  protected:
    
    double fStepTime; //Time resolution in ns
    int fNBits; //N bits of the digitizer
    double fVhigh; //Higher voltage
    double fVlow; //Lower voltage
    double fOffset; //Digitizer offset
    double fResistance; //Circuit resistance
    double fNoiseAmpl; //Electronic noise amplitud
    double fDigitizedThreshold; //Trigger threshold in ADC counts
    int fSamplingWindow; //Width of the sampling windows in ns
    int fSampleDelay; //Samples before crossing threshold that we will store
    std::map< int, std::vector<double> > fNoise; //Electronic noise non-digitized for each channel
    std::map< int, std::vector<int> > fDigitWaveForm; //Digitized waveform for each channel
    
  };
    
}

#endif
