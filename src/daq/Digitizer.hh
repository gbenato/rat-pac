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

    Digitizer(){};
    virtual ~Digitizer(){};
    Digitizer(std::string);

    virtual void SetDigitizerType(std::string);
    virtual void SetThreshold(double);
    virtual void Clear();

    virtual void AddChannel(int,PMTWaveform);
    virtual int GetNSamples(int ich){return fDigitWaveForm[ich].size();};
    // virtual void GenerateElectronicNoise(int,PMTWaveform);
    virtual std::vector<double> GetAnalogueWaveform(int ich){return fAnalogueWaveForm[ich];};
    virtual std::vector<UShort_t> GetDigitizedWaveform(int ich){return fDigitWaveForm[ich];};
    virtual std::vector<UShort_t> SampleWaveform(std::vector<UShort_t>, int );
    virtual std::vector<UShort_t> SampleWaveform(std::vector<UShort_t>);
    virtual int GoToEndOfSample(int);

    //Getters
    virtual double GetThreshold(){return fThreshold;};
    virtual int GetDigitizedThreshold(){return fDigitizedThreshold;};
    virtual double GetSamplingRate(){return fSamplingRate;};
    virtual double GetTimeResolution();
    virtual int GetNBits(){return fNBits;};
    virtual double GetVHigh(){return fVhigh;};
    virtual double GetVLow(){return fVlow;};
    virtual double GetVOffSet(){return fOffset;};
    virtual double GetResistance(){return fResistance;};
    virtual double GetNoiseAmpl(){return fNoiseAmpl;};
    virtual double GetDigitThres(){return fDigitizedThreshold;};
    virtual int GetNSamples(){return fNSamples;};
    virtual double GetTimeWindow();
    virtual double GetSampleDelay(){return fSampleDelay;};

    //Methods
    virtual double GetTimeAtSample(int sample);
    virtual int GetSampleAtTime(double time);
    virtual std::vector<double> GetWaveformTime(int pmtID);

  protected:

    DBLinkPtr fLdaq;

    std::string fDigitName; //Digitizer type
    double fSamplingRate; //Time resolution in ns
    int fNBits; //N bits of the digitizer
    double fVhigh; //Higher voltage
    double fVlow; //Lower voltage
    double fOffset; //Digitizer offset
    double fResistance; //Circuit resistance
    double fNoiseAmpl; //Electronic noise amplitud
    double fThreshold; //Threshold in volts
    int fDigitizedThreshold; //Trigger threshold in ADC counts
    int fNSamples; //Total number of samples per digitized trace
    double fSampleDelay; //Time delay before threshold for the sampling window
    std::map< int, std::vector<double> > fNoise; //Channel:Electronic noise non-digitized
    std::map< int, std::vector<double> > fAnalogueWaveForm; //Channel:Real waveform for each channel
    std::map< int, std::vector<UShort_t> > fDigitWaveForm; //Channel:Digitized waveform for each channel

  };

}

#endif
