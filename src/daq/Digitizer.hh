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
#include <RAT/DB.hh>
#include <RAT/Log.hh>

namespace RAT {

  class Digitizer {
  public:

    Digitizer();
    virtual ~Digitizer();

    virtual void SetThreshold(double);
    virtual void Clear();

    virtual void AddChannel(int,PMTWaveform);
    virtual int GetNSamples(int ich){return fDigitWaveForm[ich].size();};
    virtual void GenerateElectronicNoise(int,PMTWaveform);
    virtual std::vector<double> GetAnalogueWaveform(int ich){return fAnalogueWaveForm[ich];};
    virtual std::vector<UShort_t> GetDigitizedWaveform(int ich){return fDigitWaveForm[ich];};
    virtual std::vector<UShort_t> SampleWaveform(std::vector<UShort_t>, int);
    virtual int GoToEndOfSample(int);
    virtual double GetDigitizedThreshold(){return fDigitizedThreshold;};

    //Getters
    virtual double GetStepTime(){return fStepTime;};
    virtual int GetNBits(){return fNBits;};
    virtual double GetVHigh(){return fVhigh;};
    virtual double GetVLow(){return fVlow;};
    virtual double GetVOffSet(){return fOffset;};
    virtual double GetResistance(){return fResistance;};
    virtual double GetNoiseAmpl(){return fNoiseAmpl;};
    virtual double GetDigitThres(){return fDigitizedThreshold;};
    virtual int GetSamplingWindow(){return fSamplingWindow;};
    virtual int GetSampleDelay(){return fSampleDelay;};

  protected:

    DBLinkPtr fLdaq;

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
    std::map< int, std::vector<double> > fNoise; //Channel:Electronic noise non-digitized
    std::map< int, std::vector<double> > fAnalogueWaveForm; //Channel:Real waveform for each channel
    std::map< int, std::vector<UShort_t> > fDigitWaveForm; //Channel:Digitized waveform for each channel

  };

}

#endif
