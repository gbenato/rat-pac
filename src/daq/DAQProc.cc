#include <vector>
#include <RAT/DAQProc.hh>
#include <RAT/DB.hh>
#include <G4ThreeVector.hh>
#include <RAT/DetectorConstruction.hh>
#include <RAT/PMTPulse.hh>
#include <RAT/DS/PMTWaveform.hh>
#include <CLHEP/Random/RandGauss.h>

using namespace std;

namespace RAT {
  
  inline bool Cmp_PMTPulse_TimeAscending(const PMTPulse *a,
					 const PMTPulse *b)
  {
    double atime = a->GetPulseStartTime();
    double btime = b->GetPulseStartTime();
    return atime < btime;
  }
  
  DAQProc::DAQProc() : Processor("daq") {
    fLdaq = DB::Get()->GetLink("DAQ");

    //sampling time in ns --- this is the size of a PMT time window
    fSamplingTimeDB = fLdaq->GetD("sampling_time");
    //integration time in ns
    fIntTimeDB = fLdaq->GetD("int_time");
    //width of a PMT pulse in ns
    fPulseWidthDB = fLdaq->GetD("pulse_width");
    //offset of a PMT pulse in mV
    fPulseOffsetDB = fLdaq->GetD("pulse_offset");
    //Minimum pulse height to consider
    fPulseMinDB = fLdaq->GetD("pulse_min");
    //time before discriminator fires that sampling gate opens
    fGDelayDB = fLdaq->GetD("gate_delay");
    //width of noise in adc counts
    fNoiseAmplDB = fLdaq->GetD("noise_amplitude");
    //PMT thresholds
    fTriggerThresholdDB = fLdaq->GetD("trigger_threshold");
    //Pulse type: 0=square pulses, 1=real pulses
    fPulseTypeDB = fLdaq->GetI("pulse_type");
    //mean of a PMT pulse in ns
    fPulseMeanDB = fLdaq->GetD("pulse_mean");

    //Digitizer
    //stepping time for discrimination
    fStepTimeDB = fLdaq->GetD("step_time");
    //Offset
    fOffSetDB = fLdaq->GetD("offset");
    //High voltage
    fVHigh = fLdaq->GetD("volt_high");    
    //Low voltage
    fVLow = fLdaq->GetD("volt_low");
    //Circuit resistance
    fResistance = fLdaq->GetD("resistance");
    //N bits
    fNBits = fLdaq->GetI("nbits");
    
    
    detail << "DAQProc: DAQ constants loaded" << newline;
    detail << "  PMT Pulse type: " << (fPulseTypeDB==0 ? "square" : "realistic") << newline;
    detail << dformat("  PMT Pulse Width: ....................... %5.1f ns\n", fPulseWidthDB);
    detail << dformat("  PMT Pulse Offset: ...................... %5.1f ADC Counts\n", fPulseOffsetDB);
    detail << dformat("  Min PMT Pulse Height: .................. %5.1f mV\n", fPulseMinDB);
    detail << dformat("  PMT Channel Integration Time: .... %6.2f ns\n", fIntTimeDB);
    detail << dformat("  PMT Channel Total Sample Time: ......... %6.2f ns\n", fSamplingTimeDB);
    detail << dformat("  PMT Channel Stepping Time: ............. %6.2f ns\n", fStepTimeDB);
    detail << dformat("  Channel Gate Delay: .................... %5.1f ns\n", fGDelayDB);
    detail << dformat("  Hi Freq. Channel Noise: ................ %6.2f adc counts\n", fNoiseAmplDB);
    detail << dformat("  PMT Trigger threshold: ......................... %5.1f mV\n", fTriggerThresholdDB);
    detail << dformat("  PMT Pulse Mean: ........................ %5.1f\n", fPulseMeanDB);

    fEventCounter = 0;
  }

  
  void DAQProc::SetS(std::string param, std::string value)
  {
    if(param=="trigger")
      fTriggerType = value;

    if(fTriggerType!="allpmts" && fTriggerType!="triggerpmt"){
      std::cerr<<"DAQ: "<<fTriggerType<<" option unknown... EXIT "<<std::endl;
      exit(0);
    }
  }

  
  Processor::Result DAQProc::DSEvent(DS::Root *ds) {
    //This processor build waveforms for each PMT in the MC generated event, sample them and
    //store each sampled piece as a new event

    DS::MC *mc = ds->GetMC();
    if(ds->ExistEV()) {  // there is already a EV branch present 
      ds->PruneEV();     // remove it, otherwise we'll have multiple detector events
                         // in this physics event
                         // we really should warn the user what is taking place
    }


    fDigitizer.SetNBits(fNBits);
    fDigitizer.SetStepTime(fStepTimeDB);
    fDigitizer.SetOffSet(fOffSetDB);
    fDigitizer.SetVHigh(fVHigh);
    fDigitizer.SetVLow(fVLow);
    fDigitizer.SetResistance(fResistance);
    fDigitizer.SetNoiseAmplitude(fNoiseAmplDB);
    fDigitizer.SetSamplingWindow(fSamplingTimeDB);
    fDigitizer.SetSampleDelay((int)fGDelayDB);
    fDigitizer.SetThreshold(fTriggerThresholdDB);
    //Loop through the PMTs in the MC generated event
    //    std::map< int, std::vector<int> > DigitizedWaveforms; //ID-Waveform map
    for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

      DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      
      //For each PMT loop over hit photons and create a waveform for each of them
      DS::PMTWaveform pmtwf;
      pmtwf.SetStepTime(fStepTimeDB);
      double TimePhoton;
      //      double PulseDuty=0.0;

      for (size_t iph=0; iph < mcpmt->GetMCPhotonCount(); iph++) {
	
	DS::MCPhoton *mcphotoelectron = mcpmt->GetMCPhoton(iph);
	TimePhoton = mcphotoelectron->GetFrontEndTime();
	//TimePhoton = mcphotoelectron->GetHitTime(); //fixme: not sure about this time...
	
	//Produce pulses and add them to the waveform
	PMTPulse *pmtpulse;
	if (fPulseTypeDB==0){
            pmtpulse = new SquarePMTPulse; //square PMT pulses
	}
	else{
	  pmtpulse = new RealPMTPulse; //real PMT pulses shape
	}
	
	pmtpulse->SetPulseMean(fPulseMeanDB);
	pmtpulse->SetStepTime(fStepTimeDB);
	pmtpulse->SetPulseMin(fPulseMinDB);
	pmtpulse->SetPulseCharge(mcphotoelectron->GetCharge());
	pmtpulse->SetPulseWidth(fPulseWidthDB);
	pmtpulse->SetPulseOffset(fPulseOffsetDB);
	pmtpulse->SetPulseStartTime(TimePhoton); //also sets end time according to the pulse width and the pulse mean
	pmtwf.fPulse.push_back(pmtpulse);
	//	PulseDuty += pmtpulse->GetPulseEndTime() - pmtpulse->GetPulseStartTime();

      } // end mcphotoelectron loop: all pulses produces for this PMT
      
      //Sort pulses in time order
      std::sort(pmtwf.fPulse.begin(),pmtwf.fPulse.end(),Cmp_PMTPulse_TimeAscending);

      //At this point the PMT waveform is defined for the whole event for this PMT, so save it
      //only for drawing purposes
      mcpmt->SetWaveform(pmtwf);

      //Digitize waveform (electronic noise is added by the digitizer) and save it
      //in MCPMT object only for drawing purpose (in the future we might want to do
      //the charge integration off-line)
      fDigitizer.AddChannel(mcpmt->GetID(),pmtwf);
      mcpmt->SetDigitizedWaveform(fDigitizer.GetDigitizedWaveform(mcpmt->GetID()));
      //      DigitizedWaveform[mcpmt->GetID()] = fDigitizer.GetDigitizedWaveform();
      
    } //end pmt loop


    
    //FROM HERE THE TRIGGER PROCESSOR SHOULD TAKE OVER!
    //1) Check trigger condition and divide waveforms in chunks
    //2) Build events containing digitized waveform samples and integrated charges
    

    
    //So far the trigger condition correspond to any PMT crossing threshold. In the
    //future we'll change it to only the trigger PMT crosses threshold. When trigger
    //condition fulfilled create a new PMT in the event. If no PMT, no new event.

    //First trigger type: store all PMTs crossing threshold
    DS::EV *ev = ds->AddNewEV(); //Remove it if no PMT cross threshold
    ev->SetID(fEventCounter);
    if(fTriggerType=="allpmts"){
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
	int pmtID = mc->GetMCPMT(imcpmt)->GetID();
	//Sample digitized waveform and look for triggers
	int nsamples = fDigitizer.GetNSamples(pmtID);
	std::vector<int> DigitizedWaveform = fDigitizer.GetDigitizedWaveform(pmtID);
	for(int isample=0; isample<nsamples; isample++){
	  if (DigitizedWaveform[isample]<fDigitizer.GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)
	    DS::PMT* pmt = ev->AddNewPMT();
	    pmt->SetID(pmtID);
	    pmt->SetTime(isample*fStepTimeDB); //fixme: think about this time...
	    pmt->SetWaveform(fDigitizer.SampleWaveform(DigitizedWaveform,isample)); //it is defined by the sample that crosses threshold
	    pmt->SetCharge(fDigitizer.IntegrateCharge(DigitizedWaveform));
	    isample = fDigitizer.GoToEndOfSample(isample); //go forward towards the end of the sampling window
	  }//end if above trigger
	}//end sampling
	DigitizedWaveform.clear(); //prune for next round of PMTs
      }//end PMT loop
    fDigitizer.Clear();
    }
    //Second trigger type: when trigger PMT detects a hit above threshold store
    //hits in ALL the PMTs
    else if(fTriggerType=="triggerpmt"){
      //Identify the trigger PMT
      int triggerID=-1;
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
	DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
	if(mcpmt->GetType() == 0) triggerID = mcpmt->GetID();
      }
      //If trigger PMT has been hit, sample its waveform and check if it crosses
      //threshold
      if(triggerID>-1){
	int nsamples = fDigitizer.GetNSamples(triggerID);
	std::vector<int> DigitizedTriggerWaveform = fDigitizer.GetDigitizedWaveform(triggerID);
	for(int isample=0; isample<nsamples; isample++){

	  //	  std::cout<<" SAMPLE "<<isample<<" "<<DigitizedTriggerWaveform[isample]<<" "<<fDigitizer.GetDigitizedThreshold()<<std::endl;
	  
	  if (DigitizedTriggerWaveform[isample]<fDigitizer.GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)
	    //Read ALL PMTs
	    for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
	      int pmtID = mc->GetMCPMT(imcpmt)->GetID();
	      DS::PMT* pmt = ev->AddNewPMT();
	      pmt->SetID(pmtID);
	      pmt->SetWaveform(fDigitizer.SampleWaveform(fDigitizer.GetDigitizedWaveform(pmtID), isample));
	      pmt->SetCharge(fDigitizer.IntegrateCharge(fDigitizer.GetDigitizedWaveform(pmtID)));
	      pmt->SetTime(fDigitizer.GetPeakTime(pmtID,isample)); //sample the waveform and find the peak position
	      //	      pmt->SetTime(isample*fStepTimeDB);
	    } //end reading PMTs
	    isample = fDigitizer.GoToEndOfSample(isample); //go forward towards the end of the sampling window
	  } //end if: trigger above threshold
	}//end sampling
	DigitizedTriggerWaveform.clear(); //prune for next round of PMTs
      } //end if hit trigger PMT
      fDigitizer.Clear();
    } //end if second type of trigger
	
    //If got at least one PMT above threshold move forward one event so it is not
    //overwritten by the next one
    if(ev->GetPMTCount()>0){
      fEventCounter++;
    }




    //FIXME
    // else{
    //   std::cout<<"Prune event"<<fEventCounter<<" "<<ds->GetEVCount()<<std::endl;
    //   ds->PruneEV(fEventCounter);
    // }
    
    /*      
    //Sample the PMT waveform to look for photoelectron hits and create a new MCPMTSample
      //for every one of them, regarless they cross threshold. A flag is raised if the threshold
      //is crossed
      double TimeNow = 0.;
      for(int isample=0; isample<digitizer->GetNSamples(); isample++){

	TimeNow = digitizer->GetTimeForSample(isample);
	if (digitizer->GetHeightForSample()>fTriggerThresholdDB){ //hit above threshold!

	  bool IsAboveThreshold = true;
	  double HitStartTime = TimeNow;
	  
	  int lsample = isample + (int)fSamplingTimeDB/fSampleStep; //Last sample in the sampling window
	  
	  //Set PMT observables and save it as a new PMT object in the event
	  DS::PMT* mcpmtsample = mc->AddNewMCPMTSampled();
	  mcpmtsample->SetID(mcpmt->GetID());
	  mcpmtsample->SetCharge(IntegratedCharge);
	  mcpmtsample->SetTime(HitStartTime);
	  mcpmtsample->SetAboveThreshold(IsAboveThreshold);
	  mcpmtsample->SetDigitizedWaveForm(digitizer->GetDigitizedChunk(isample,lsample));

	  //Continue until de last sample in the sampling window to start over
	  isample = lsample;
	  
	} //end if above threshold

	  
      }

    

      















      while (pmtwf.fPulse.size()>0){
	
	double TimeNow = pmtwf.fPulse[0]->GetPulseStartTime();
	double LastPulseTime = pmtwf.fPulse[pmtwf.fPulse.size()-1]->GetPulseEndTime();
	int NextPulse=0;
	double wfheight = 0.0;
	//	float NoiseAmpl = fNoiseAmplDB/sqrt(PulseDuty/fStepTimeDB);
	//	float qfuzz=0.0;

	//Check if the waveform crosses the threshold
	while( wfheight < fTriggerThresholdDB && TimeNow < LastPulseTime){
	  wfheight =  pmtwf.GetHeight(TimeNow); //height of the waveform at this step
	  TimeNow += fStepTimeDB;

	  //	  qfuzz = wfheight+NoiseAmpl*CLHEP::RandGauss::shoot();	  
	  // move forward in time by step or if we're at baseline move ahead to next pulse
	  // if (wfheight==0.0){
	  //   NextPulse = pmtwf.GetNext(TimeNow);
	  //   if (NextPulse>0)
	  //     {TimeNow = pmtwf.fPulse[NextPulse]->GetPulseStartTime();}
	  //   else
	  //     {TimeNow= LastPulseTime+1.;}
	  // }
	  // else{
	  //   TimeNow += fStepTimeDB;
	  // }
	}
	
	//If this PMT crosses the threshold set the hit time as the time when it crosses the
	//threshold and set the flag to true. If doesn't, set the hit time to the starting time
	//of the pulse.
	double HitStartTime = pmtwf.fPulse[0]->GetPulseStartTime();
	bool IsAboveThreshold = false;
	if (TimeNow < LastPulseTime){ //means that we have crossed threshold in the time window
	  IsAboveThreshold = true;
	  HitStartTime = TimeNow;
	}
	
	//Integrate charge and create a new PMT in the event
	double IntegratedCharge = 0.;
	double TimeStartSample = HitStartTime - fGDelayDB; //start before several steps before thresholds
	double TimeEndSample = TimeStartSample+fSamplingTimeDB;
	double TimeEndIntegration = TimeStartSample+fIntTimeDB;
	
	unsigned int ipulse=0;
	while (ipulse < pmtwf.fPulse.size() && pmtwf.fPulse[ipulse]->GetPulseStartTime()<TimeEndSample){
	  IntegratedCharge+=pmtwf.fPulse[ipulse]->Integrate(TimeStartSample,TimeEndIntegration); //Fixme: create a function in PMTWaveForm that performs the integration
	  ipulse++;
	}

	//Go until the last pulse in the sampling window to start over again
	while (ipulse < pmtwf.fPulse.size() && pmtwf.fPulse[ipulse]->GetPulseStartTime()<TimeEndSample ) {
	  ipulse++;
	}
	
	//Set PMT observables and save it as a new PMT object in the event
	DS::PMT* mcpmtsample = mc->AddNewMCPMTSampled();
	mcpmtsample->SetID(mcpmt->GetID());
	mcpmtsample->SetCharge(IntegratedCharge);
	mcpmtsample->SetTime(HitStartTime);
	mcpmtsample->SetAboveThreshold(IsAboveThreshold);

	//Remove all the pulses whithin the sampling window
	pmtwf.fPulse.erase(pmtwf.fPulse.begin(),pmtwf.fPulse.begin()+ipulse);
	
      } //end pulse sampling
      
      //clean up just in case
      for (unsigned int i = 0; i<pmtwf.fPulse.size();i++){
        delete pmtwf.fPulse[i];
      }
      
    } //end PMT loop


    */


    return Processor::OK;

  } //DAQProc::DSEvent
  
} // namespace RAT
