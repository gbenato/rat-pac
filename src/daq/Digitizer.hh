#ifndef __RAT_Digitizer__
#define __RAT_Digitizer__

#include <RAT/DS/PMTWaveForm.hh>

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
    virtual void Clear(){fDigitWaveForm.clear(); fNoise.clear();};
    
    virtual int GetNSamples(){return fDigitWaveForm.size();};
    virtual void GenerateElectronicNoise(DS::PMTWaveform);
    virtual void DigitizeWaveForm(DS::PMTWaveform);
    virtual std::vector<int> GetDigitizedWaveform(){return fDigitWaveForm;};
    virtual std::vector<int> SampleWaveform(std::vector<int>, int*);
    virtual double GetChargeForSample(std::vector<int>);
    
  protected:
    
    double fStepTime; //Time resolution in ns
    int fNBits; //N bits of the digitizer
    double fVhigh; //Higher voltage
    double fVlow; //Lower voltage
    double fOffset; //Digitizer offset
    double fResistance; //Circuit resistance
    double fNoiseAmpl; //Electronic noise amplitud
    int fSamplingWindow; //Width of the sampling windows in ns
    int fSampleDelay; //Samples before crossing threshold that we will store
    std::vector<double> fNoise; //Electronic noise non-digitized
    std::vector<int> fDigitWaveForm; //Digitized waveform
    
  };
    
}

#endif
