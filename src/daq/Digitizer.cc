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


  // void Digitizer::GenerateElectronicNoise(int ichannel, PMTWaveform pmtwf) {
  //
  //   //*** deprecated
  //   //Sort pulses in time order
  //   //    std::sort(fPulse.begin(),fPulse.end(),Cmp_PMTPulse_TimeAscending);
  //   // double starttime = pmtwf.fPulse.front()->GetPulseStartTime();
  //   // double endtime = pmtwf.fPulse.back()->GetPulseEndTime();
  //   // if(starttime-endtime < fSamplingWindow) endtime = starttime + fSamplingWindow;
  //   //    int nsamples = (endtime - starttime)/fStepTime;
  //   // double PulseDuty=0.0;
  //   // for(int ipulse = 0; ipulse<pmtwf.fPulse.size(); ipulse++)
  //   //   PulseDuty +=  pmtwf.fPulse[ipulse]->GetPulseEndTime() - pmtwf.fPulse[ipulse]->GetPulseStartTime();
  //   // float NoiseAmpl = fNoiseAmpl/sqrt(PulseDuty/fStepTime);
  //   //    std::cout<<starttime<<" "<<endtime<<" "<<nsamples<<" "<<PulseDuty<<" "<<NoiseAmpl<<std::endl;
  //   //*********
  //
  //   int nsamples = fSamplingWindow/fStepTime;
  //   fNoise[ichannel].resize(nsamples);
  //   for(int istep=0; istep<nsamples ;istep++){
  //     fNoise[ichannel][istep] = fNoiseAmpl*CLHEP::RandGauss::shoot();
  //   }
  //
  // }

  //Add channel to digitizer and inmdediatly digitize analogue waveform
  void Digitizer::AddChannel(int ichannel, PMTWaveform pmtwf){

    double starttime = -fSampleDelay;
    double endtime = starttime + fSamplingWindow;
    int nsamples = fSamplingWindow/fStepTime;
    //Ensure enough samples in the sampling window
    nsamples *= 2;
    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double adcpervolt = nADCs/(fVhigh - fVlow);
    double charge = 0., tcharge = 0.;

    //    std::cout<<starttime<<" "<<endtime<<" "<<nsamples<<" "<<nADCs<<" "<<adcpervolt<<std::endl;

    double currenttime = starttime;
    double volt = 0.;
    int adcs = 0.;
    double tadcs = 0.;

    //First get analogue waveform for drawing purpose and debugging
    for(double itime = 0; itime<fSamplingWindow; itime += pmtwf.GetStepTime()){

      charge = pmtwf.GetCharge(itime, itime + pmtwf.GetStepTime());
      fAnalogueWaveForm[ichannel].push_back(charge);

    }


    for(int isample = 0; isample<nsamples; isample++){

      charge = pmtwf.GetCharge(currenttime-fStepTime, currenttime);

      volt = charge*fResistance/fStepTime; //convert to voltage
      // volt = volt + fNoise[ichannel][isample]; //add electronic noise
      volt = volt + fNoiseAmpl*CLHEP::RandGauss::shoot(); //add electronic noise
      adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC

      // std::cout<<" Round "<<(volt - fVlow + fOffset)*adcpervolt<<" "<<adcs<<std::endl;

      //adcs = (volt - fVlow + fOffset)*adcpervolt; //digitize: V->ADC
      tcharge += charge; //not used, just for sanity check
      //Manage voltage saturation
      if(adcs<0) adcs = 0;
      else if(adcs>=nADCs) adcs = nADCs - 1;

      // tadcs += (volt - fVlow + fOffset)*adcpervolt;
      tadcs += adcs;

      //Save sample
      // fAnalogueWaveForm[ichannel].push_back(charge);
      fDigitWaveForm[ichannel].push_back(adcs);
      // std::cout<<isample<<" "<<volt<<" "<<adcs<<" "<<pmtwf.GetHeight(currenttime)<<" "<<fNoise[ichannel][isample]<<std::endl;

      //Step on time
      currenttime+=fStepTime;
    }

    //Ensure we always have enough samples in the pedestal window
    //Insert samples with pure noise at the begining
    //(FIXME: might allow to have darks hit)
    int nsamples_delay = (int)fSampleDelay/fStepTime;
    for (size_t isample = 0; isample < nsamples_delay; isample++) {
      volt = fNoiseAmpl*CLHEP::RandGauss::shoot();
      adcs = round((volt - fVlow + fOffset)*adcpervolt); //digitize: V->ADC
      //Manage voltage saturation
      if(adcs<0) adcs = 0;
      else if(adcs>=nADCs) adcs = nADCs - 1;

      fDigitWaveForm[ichannel].insert(fDigitWaveForm[ichannel].begin(),adcs);
    }

    std::cout<<" Digit charge in channel "<<ichannel<<": "<<tcharge<<" "<<(tadcs/adcpervolt + fVlow*nsamples - fOffset*nsamples)*fStepTime/fResistance<<std::endl;

  }


  //Retrieves a chunk of the digitized waveform in a sampling
  //window defined by the user by fSampleDelay and fSamplingWindow
  //[init_sample-fSampleDelay, thres_sample+fSamplingWindow]
  std::vector<UShort_t> Digitizer::SampleWaveform(std::vector<UShort_t> completewaveform, int init_sample){

    //FIXME: might ensure we extract always the good number of samples
    // std::cout<<" 0 - SampleWaveform "<<completewaveform.size()<<std::endl;

    int nsamples_delay = (int)fSampleDelay/fStepTime;
    int start_sample = init_sample - nsamples_delay;
    int current_sample = start_sample;

//     //Ensure we always have enough samples in the pedestal window
//     if(start_sample<0) {
//       //Insert zeros at the begining
//       completewaveform.insert(completewaveform.begin(),-start_sample,completewaveform[0]); //FIXME
//       start_sample=0;
//     }
//
// // std::cout<<" 1 - SampleWaveform "<<completewaveform.size()<<" "<<start_sample<<std::endl;
//
    int end_sample = start_sample+(int)fSamplingWindow/fStepTime;
//     //Ensure we always have maximum number of samples
//     if(end_sample>completewaveform.size()){
//       //Insert zeros at the end
//       completewaveform.insert(completewaveform.end(),end_sample-completewaveform.size(),completewaveform[0]); //FIXME
//     }
     std::vector<UShort_t> sampledwaveform;

//    std::cout<<" 2 - SampleWaveform "<<completewaveform.size()<<" "<<end_sample<<std::endl;

    // std::cout<<" 1 -  init_sample "<<init_sample<<" start_sample "<<start_sample<<" "<<current_sample<<" "<<end_sample<<" "<<nsamples_delay<<std::endl;

    while(current_sample<end_sample){
      sampledwaveform.push_back(completewaveform[current_sample]);
//      std::cout<<" 3 - SampleWaveform "<<current_sample<<"/"<<end_sample<<" "<<completewaveform[current_sample]<<std::endl;
      current_sample++;
    }

    // std::cout<<" SampleWaveform "<<completewaveform.size()<<" "<<sampledwaveform.size()<<std::endl;

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
