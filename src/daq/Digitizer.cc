#include <iostream>
#include <RAT/Digitizer.hh>
#include <CLHEP/Random/RandGauss.h>

namespace RAT {

  Digitizer::Digitizer(){

    //Set configuration according to DAQ.ratdb
    fLdaq = DB::Get()->GetLink("DAQ");
    //time resolution
    fStepTime = fLdaq->GetD("step_time");
    //Offset
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
    this->SetThreshold(fLdaq->GetD("trigger_threshold"));
    //sampling time in ns
    fSamplingWindow = fLdaq->GetD("sampling_time");
    //time before discriminator fires that sampling gate opens
    fSampleDelay = fLdaq->GetD("gate_delay");

    detail << dformat("  Digitizer: Hi Freq. Channel Noise: ................ %6.2f adc counts\n", fNoiseAmpl);
    detail << dformat("  Digitizer: PMT Trigger threshold: ................. %5.1f mV\n", fLdaq->GetD("trigger_threshold"));
    detail << dformat("  Digitizer: PMT Trigger threshold (Digit): ......... %5.1f samples\n", fDigitizedThreshold);
    detail << dformat("  Digitizer: Stepping Time: ......................... %6.2f ns\n", fStepTime);
    detail << dformat("  Digitizer: Total Sample Time: ..................... %6.2f ns\n", fSamplingWindow);
    detail << dformat("  Digitizer: Gate Delay: ............................ %5.1f ns\n", fSampleDelay);
    detail << dformat("  Digitizer: Voltage offset: ........................ %6.2f mV\n", fOffset);
    detail << dformat("  Digitizer: Voltage High: .......................... %6.2f mV\n", fVhigh);
    detail << dformat("  Digitizer: Voltage Low: ........................... %6.2f mV\n", fVlow);
    detail << dformat("  Digitizer: Resistance: ............................ %6.2f mV\n", fResistance);

  }

  Digitizer::~Digitizer(){}


  void Digitizer::GenerateElectronicNoise(int ichannel, PMTWaveform pmtwf) {

    //*** deprecated
    //Sort pulses in time order
    //    std::sort(fPulse.begin(),fPulse.end(),Cmp_PMTPulse_TimeAscending);
    // double starttime = pmtwf.fPulse.front()->GetPulseStartTime();
    // double endtime = pmtwf.fPulse.back()->GetPulseEndTime();
    // if(starttime-endtime < fSamplingWindow) endtime = starttime + fSamplingWindow;
    //    int nsamples = (endtime - starttime)/fStepTime;
    // double PulseDuty=0.0;
    // for(int ipulse = 0; ipulse<pmtwf.fPulse.size(); ipulse++)
    //   PulseDuty +=  pmtwf.fPulse[ipulse]->GetPulseEndTime() - pmtwf.fPulse[ipulse]->GetPulseStartTime();
    // float NoiseAmpl = fNoiseAmpl/sqrt(PulseDuty/fStepTime);
    //    std::cout<<starttime<<" "<<endtime<<" "<<nsamples<<" "<<PulseDuty<<" "<<NoiseAmpl<<std::endl;
    //*********

    //Sort pulses in time order
    //    std::sort(fPulse.begin(),fPulse.end(),Cmp_PMTPulse_TimeAscending);
    // double starttime = pmtwf.fPulse.front()->GetPulseStartTime();
    // double endtime = pmtwf.fPulse.back()->GetPulseEndTime();
    // if(starttime-endtime < fSamplingWindow) endtime = starttime + fSamplingWindow;
    int nsamples = fSamplingWindow/fStepTime;
    fNoise[ichannel].resize(nsamples);
    float NoiseAmpl = fNoiseAmpl;
    for(int istep=0; istep<nsamples ;istep++){
      //      fNoise[ichannel][istep] = 0.;
      fNoise[ichannel][istep] = NoiseAmpl*CLHEP::RandGauss::shoot();
    }

  }

  //*** deprecated
  // void Digitizer::DigitizeWaveForm(DS::PMTWaveform pmtwf){

  //   GenerateElectronicNoise(pmtwf); //Fill fNoise vector

  //   double starttime = pmtwf.fPulse.front()->GetPulseStartTime();
  //   double endtime = pmtwf.fPulse.back()->GetPulseEndTime();
  //   if(starttime-endtime < fSamplingWindow) endtime = starttime + fSamplingWindow;
  //   int nsamples = (endtime - starttime)/fStepTime;
  //   int nADCs = 1 << fNBits; //Calculate the number of adc counts
  //   double adcpervolt = nADCs/(fVhigh - fVlow);
  //   double charge = 0.;

  //   //    std::cout<<starttime<<" "<<endtime<<" "<<nsamples<<" "<<nADCs<<" "<<adcpervolt<<std::endl;

  //   double currenttime = starttime;
  //   double volt = 0.;
  //   int adcs = 0.;
  //   for(int isample = 0; isample<nsamples; isample++){
  //     charge = pmtwf.GetHeight(currenttime)*fStepTime;
  //     volt = charge*fResistance/fStepTime; //convert to voltage
  //     volt = volt + fNoise[isample]; //add electronic noise
  //     adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
  //     //      charge += (pmtwf.GetHeight(currenttime)+fNoise[isample])*fStepTime; //not used, just for sanity check

  //     //Manage voltage saturation
  //     if(adcs<0) adcs = 0;
  //     else if(adcs>=nADCs) adcs = nADCs - 1;

  //     //Save sample
  //     fDigitWaveForm.push_back(adcs);
  //     //      std::cout<<isample<<" "<<volt<<" "<<adcs<<" "<<pmtwf.GetHeight(currenttime)<<" "<<fNoise[isample]<<std::endl;

  //     //Step on time
  //     currenttime+=fStepTime;
  //   }

  //   //    std::cout<<"Analogue integrated charge "<<charge<<std::endl;


  // }
  //**********



  //Add channel to digitizer and inmdediatly digitize analogue waveform
  void Digitizer::AddChannel(int ichannel, PMTWaveform pmtwf){

    GenerateElectronicNoise(ichannel,pmtwf); //Fill fNoise vector

    //    double starttime = pmtwf.fPulse.front()->GetPulseStartTime();
    //    double endtime = pmtwf.fPulse.back()->GetPulseEndTime();
    //    if(endtime-starttime < fSamplingWindow) endtime = starttime + fSamplingWindow;
    double starttime = -fSampleDelay;
    double endtime = starttime + fSamplingWindow;
    int nsamples = fSamplingWindow/fStepTime;
    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double adcpervolt = nADCs/(fVhigh - fVlow);
    double charge = 0.;

    //    std::cout<<starttime<<" "<<endtime<<" "<<nsamples<<" "<<nADCs<<" "<<adcpervolt<<std::endl;

    double currenttime = starttime;
    double volt = 0.;
    int adcs = 0.;
    for(int isample = 0; isample<nsamples; isample++){
      charge = pmtwf.GetHeight(currenttime)*fStepTime;
      volt = charge*fResistance/fStepTime; //convert to voltage
      volt = volt + fNoise[ichannel][isample]; //add electronic noise
      adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
      //      charge += (pmtwf.GetHeight(currenttime)+fNoise[isample])*fStepTime; //not used, just for sanity check

      //Manage voltage saturation
      if(adcs<0) adcs = 0;
      else if(adcs>=nADCs) adcs = nADCs - 1;

      //Save sample
      fAnalogueWaveForm[ichannel].push_back(pmtwf.GetHeight(currenttime));
      fDigitWaveForm[ichannel].push_back(adcs);
      //      std::cout<<isample<<" "<<volt<<" "<<adcs<<" "<<pmtwf.GetHeight(currenttime)<<" "<<fNoise[ichannel][isample]<<std::endl;

      //Step on time
      currenttime+=fStepTime;
    }

    //    std::cout<<"Analogue integrated charge "<<charge<<std::endl;

  }


  //Retrieves a chunk of the digitized waveform in a sampling
  //window defined by the user by fSampleDelay and fSamplingWindow
  //[init_sample-fSampleDelay, thres_sample+fSamplingWindow]
  std::vector<UShort_t> Digitizer::SampleWaveform(std::vector<UShort_t> completewaveform, int init_sample){

    int sample_delay = (int)fSampleDelay/fStepTime;
    int start_sample = init_sample - sample_delay;
    if(start_sample<0) start_sample = 0; //FIXME: this changes size of sample_delay
    int end_sample = init_sample+(int)fSamplingWindow/fStepTime;
    if(end_sample>completewaveform.size()-1) end_sample = completewaveform.size() - 1;
    std::vector<UShort_t> sampledwaveform;

    while(start_sample<=end_sample){
      sampledwaveform.push_back(completewaveform[start_sample]);
      start_sample++;
    }

    return sampledwaveform;

  }

  //Moves the sampling point towards the end of the sampling window defined by the
  //user
  int Digitizer::GoToEndOfSample(int start_sample){

    int end_sample = start_sample+(int)fSamplingWindow/fStepTime;
    return end_sample;

  }

  //Parse threshold in volts and store it in ADC counts
  void Digitizer::SetThreshold(double threhold_volts){

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double adcpervolt = nADCs/(fVhigh - fVlow);
    fDigitizedThreshold = round((threhold_volts - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
    //    std::cout<<"SetThreshold "<<fDigitizedThreshold<<std::endl;

  }

  // double Digitizer::IntegrateCharge(std::vector<int> digitizedwaveform){
  //
  //   int nADCs = 1 << fNBits; //Calculate the number of adc counts
  //   double voltsperadc = (fVhigh - fVlow)/(double)nADCs;
  //   double charge = 0.;
  //   int start_sample = 0;
  //   while(start_sample<digitizedwaveform.size()){
  //     charge += ((double)digitizedwaveform[start_sample]*voltsperadc + fVlow)*fStepTime/fResistance; //ADC to charge
  //     start_sample++;
  //   }
  //   //    std::cout<<" Digitized integrated charge "<<charge<<std::endl;
  //
  //   charge=std::abs(charge); //pulses are negative so covert to positive charge
  //   return charge;
  // }

  // //Calculates the time at which the peak of the digitized waveform occurs
  // double Digitizer::GetTimeAtPeak(int pmtID, int init_sample){
  //
  //   //Retrieve a piece of the waveform within the sampling window
  //   std::vector<int> sampledwf = this->SampleWaveform(this->GetDigitizedWaveform(pmtID), init_sample);
  //   //Sample waveform to look for the maximum
  //   double sampleatpeak = -9999.;
  //   int voltatstep = 9999.;
  //   for(int isample=0; isample<sampledwf.size(); isample++){
  //     if(sampledwf[isample]<voltatstep){
	// //	std::cout<<"GetTimeAtPeak "<<pmtID<<" "<<isample<<" "<<sampledwf[isample]<<std::endl;
	// voltatstep = sampledwf[isample];
	// sampleatpeak = isample;
  //     }
  //   }
  //
  //   return sampleatpeak*fStepTime;
  // }
  //
  // //Calculates the time at which a digitized waveform crosses threshold with respect
  // //to init_sample (that ideally is the time at which the trigger PMT crosses threshold)
  // double Digitizer::GetTimeAtThreshold(int pmtID, int trigger_sample){
  //
  //   //Retrieve a piece of the waveform within the sampling window
  //   std::vector<int> sampledwf = this->SampleWaveform(this->GetDigitizedWaveform(pmtID), trigger_sample);
  //   //Sample waveform to look for the threshold crossing
  //   int sampleatthres = 0;
  //   for(int isample=0; isample<sampledwf.size(); isample++){
  //     if(sampledwf[isample]<=this->GetDigitizedThreshold()){
	// std::cout<<"GetTimeAtThreshold "<<pmtID<<" "<<isample<<" "<<sampledwf[isample]<<std::endl;
	// sampleatthres = isample;
	// break;
  //     }
  //   }
  //
  //   return sampleatthres*fStepTime;
  // }


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
