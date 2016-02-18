#include <iostream>
#include <RAT/Digitizer.hh>
#include <CLHEP/Random/RandGauss.h>

namespace RAT {

  Digitizer::Digitizer(std::string digitName){
    SetDigitizerType(digitName);
  }

  void Digitizer::SetDigitizerType(std::string digitName){

    std::cout<<" Adding digitizer "<<digitName<<std::endl;

    //Set configuration according to DAQ.ratdb
    fDigitName = digitName;
    fLdaq = DB::Get()->GetLink("DIGITIZER",fDigitName);
    //Sampling Rate
    fSamplingRate = fLdaq->GetD("sampling_rate");
    //OffsetfSampleDelay
    fOffset = fLdaq->GetD("offset");
    //High voltage
    fVhigh = fLdaq->GetD("volt_high");
    //Low voltage
    fVlow = fLdaq->GetD("volt_low");
    //Circuit resistance
    fResistance = fLdaq->GetD("resistance");
    //N bits
    fNBits = fLdaq->GetI("nbits");
    //width of noise in adc counts
    fNoiseAmpl = fLdaq->GetD("noise_amplitude");
    //PMT trigger thresholds
    SetThreshold(fLdaq->GetD("trigger_threshold"));
    //sampling time in ns
    fNSamples = fLdaq->GetD("nsamples");
    //time before discriminator fires that sampling gate opens
    fSampleDelay = fLdaq->GetD("gate_delay");

    detail << dformat("  Digitizer: Hi Freq. Channel Noise: ................ %6.2f adc counts\n", fNoiseAmpl);
    detail << dformat("  Digitizer: PMT Trigger threshold: ................. %5.1f mV\n", fLdaq->GetD("trigger_threshold"));
    detail << dformat("  Digitizer: PMT Trigger threshold (Digit): ......... %5.1f samples\n", fDigitizedThreshold);
    detail << dformat("  Digitizer: Sampling Rate: ......................... %6.2f ns\n", fSamplingRate);
    detail << dformat("  Digitizer: Total Number of Samples: ............... %6.2f ns\n", fNSamples);
    detail << dformat("  Digitizer: Gate Delay: ............................ %5.1f ns\n", fSampleDelay);
    detail << dformat("  Digitizer: Voltage offset: ........................ %6.2f mV\n", fOffset);
    detail << dformat("  Digitizer: Voltage High: .......................... %6.2f mV\n", fVhigh);
    detail << dformat("  Digitizer: Voltage Low: ........................... %6.2f mV\n", fVlow);
    detail << dformat("  Digitizer: Resistance: ............................ %6.2f mV\n", fResistance);

  }

  //Add channel to digitizer and inmdediatly digitize analogue waveform
  void Digitizer::AddChannel(int ichannel, PMTWaveform pmtwf){

    //    double starttime = fSampleDelay;
    double starttime = 0;
    double endtime = starttime + GetTimeWindow();
    double timeres = GetTimeResolution();
    int nsamples = GetNSamples();
    //Ensure enough samples in the sampling window
    nsamples *= 2;
    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double adcpervolt = nADCs/(fVhigh - fVlow);
    double charge = 0., tcharge = 0.;

    //First get analogue waveform for drawing purpose and debugging
    for(double itime = 0; itime<GetTimeWindow(); itime += pmtwf.GetStepTime()){
      //      charge = pmtwf.GetCharge(itime, itime + pmtwf.GetStepTime());
      charge = pmtwf.GetHeight(itime);
      fAnalogueWaveForm[ichannel].push_back(charge);
      //      std::cout<<" Digitizer::AddChannel: time "<<itime<<" step "<<pmtwf.GetStepTime()<<" charge "<<charge<<std::endl;
    }

    //Second digitize analogue
    double currenttime = starttime;
    double volt = 0.;
    int adcs = 0.;
    double tadcs = 0.;
    for(int isample = 0; isample<nsamples; isample++){

      //      charge = pmtwf.GetCharge(currenttime-timeres, currenttime);
      charge = pmtwf.GetHeight(currenttime)*timeres;
      volt = charge*fResistance/timeres; //convert to voltage
      volt = volt + fNoiseAmpl*CLHEP::RandGauss::shoot(); //add electronic noise
      adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
      tcharge += charge; //not used, just for sanity check

      //Manage voltage saturation
      if(adcs<0) adcs = 0;
      else if(adcs>=nADCs) adcs = nADCs - 1;

      // tadcs += (volt - fVlow + fOffset)*adcpervolt;
      tadcs += adcs;

      //Save sample
      fDigitWaveForm[ichannel].push_back(adcs);
      // std::cout<<isample<<" "<<volt<<" "<<adcs<<" "<<pmtwf.GetHeight(currenttime)<<" "<<fNoise[ichannel][isample]<<std::endl;

      //Step on time
      currenttime+=timeres;
    }

  }

  //Retrieves a chunk of the digitized waveform in a sampling
  //window defined by the user by fSampleDelay and GetTimeWindow()
  //[init_sample-fSampleDelay, thres_sample + GetTimeWindow() ]
  std::vector<UShort_t> Digitizer::SampleWaveform(std::vector<UShort_t> completewaveform, int init_sample){

    int nsamples_delay = (int)fSampleDelay/GetTimeResolution();
    int start_sample = init_sample + nsamples_delay;
    int end_sample = start_sample + GetNSamples();

    //Ensure we always have enough samples in the pedestal window
    while(start_sample<0) {
      //Insert noisy samples at the beginning of the window
      int nADCs = 1 << fNBits; //Calculate the number of adc counts
      double adcpervolt = nADCs/(fVhigh - fVlow);
      double volt = fNoiseAmpl*CLHEP::RandGauss::shoot();
      int adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
      //Manage voltage saturation
      if(adcs<0) adcs = 0;
      else if(adcs>=nADCs) adcs = nADCs - 1;

      completewaveform.insert(completewaveform.begin(),adcs);
      start_sample++;
      end_sample++;
    }

    std::vector<UShort_t> sampledwaveform;
    int current_sample = start_sample;
    while(current_sample<end_sample){
      sampledwaveform.push_back(completewaveform[current_sample]);
      current_sample++;
    }

    return sampledwaveform;

  }

  //If init_sample it's not specified -> Sample WF from the defined time delay
  std::vector<UShort_t> Digitizer::SampleWaveform(std::vector<UShort_t> completewaveform){
    return SampleWaveform(completewaveform, 0);
  }

  //Moves the sampling point towards the end of the sampling window defined by the
  //user
  int Digitizer::GoToEndOfSample(int start_sample){
    int end_sample = start_sample+(int)GetTimeWindow()/GetTimeResolution();
    return end_sample;
  }

  //Parse threshold in volts and store it in ADC counts
  void Digitizer::SetThreshold(double threhold_volts){

    fThreshold = threhold_volts;
    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double adcpervolt = nADCs/(fVhigh - fVlow);
    fDigitizedThreshold = round((threhold_volts - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
    //    std::cout<<"SetThreshold "<<fDigitizedThreshold<<std::endl;

  }

  //Get Time Resolution in ns
  double Digitizer::GetTimeResolution(){
    return 1./fSamplingRate;
  }

  double Digitizer::GetTimeWindow(){
    return fNSamples*GetTimeResolution();
  }

  double Digitizer::GetTimeAtSample(int sample){
    return (double)sample*GetTimeResolution() + fSampleDelay;
  }

  int Digitizer::GetSampleAtTime(double utime){
    return floor((utime - fSampleDelay)/GetTimeResolution());
  }

  void Digitizer::Clear(){

    for(std::map<int, std::vector<double> >::iterator it = fAnalogueWaveForm.begin(); it!=fAnalogueWaveForm.end(); it++){
      it->second.clear();
    }

    for(std::map<int, std::vector<UShort_t> >::iterator it = fDigitWaveForm.begin(); it!=fDigitWaveForm.end(); it++){
      it->second.clear();
    }

    for(std::map<int, std::vector<double> >::iterator it = fNoise.begin(); it!=fNoise.end(); it++){
      it->second.clear();
    }

    //Nicer but not compatible with C++<11
    // for(auto &it : fDigitWaveForm){
    //   it.second.clear();
    // }
    // for(auto &it : fNoise){
    //   it.second.clear();
    // }

  }
}
