#include <vector>
#include <RAT/DAQProc.hh>
#include <RAT/DB.hh>
#include <G4ThreeVector.hh>
#include <RAT/DetectorConstruction.hh>
#include <RAT/PMTPulse.hh>
#include <RAT/PMTWaveform.hh>
#include <RAT/DS/DAQHeader.hh>
#include <RAT/DS/RunStore.hh>

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
      fDigitizer = new Digitizer();

      //Pulses and waveforms features
      //sampling time in ns --- this is the size of a PMT time window
      fSamplingTimeDB = fLdaq->GetD("sampling_time");
      //width of a PMT pulse in ns
      fPulseWidthDB = fLdaq->GetD("pulse_width");
      //offset of a PMT pulse in mV
      fPulseOffsetDB = fLdaq->GetD("pulse_offset");
      //Minimum pulse height to consider
      fPulseMinDB = fLdaq->GetD("pulse_min");
      //Pulse type: 0=square pulses, 1=real pulses
      fPulseTypeDB = fLdaq->GetI("pulse_type");
      //mean of a PMT pulse in ns
      fPulseMeanDB = fLdaq->GetD("pulse_mean");
      //pulse step time
      fPulseStepTimeDB = fLdaq->GetD("pulse_step_time");

      detail << "DAQProc: DAQ constants loaded" << newline;
      detail << "  PMT Pulse type: " << (fPulseTypeDB==0 ? "square" : "realistic") << newline;
      detail << dformat("  PMT Pulse Mean: ........................ %5.1f\n", fPulseMeanDB);
      detail << dformat("  PMT Pulse Width: ....................... %5.1f ns\n", fPulseWidthDB);
      detail << dformat("  PMT Pulse Offset: ...................... %5.1f ADC Counts\n", fPulseOffsetDB);
      detail << dformat("  Min PMT Pulse Height: .................. %5.1f mV\n", fPulseMinDB);

      fEventCounter = 0;
    }


    void DAQProc::SetS(std::string param, std::string value)
    {
      if(param=="trigger")
      fTriggerType = value;

      if(fTriggerType!="allpmts" && fTriggerType!="triggerpmt" && fTriggerType!="simpledaq"){
        std::cerr<<"DAQ: "<<fTriggerType<<" option unknown... EXIT "<<std::endl;
        exit(0);
      }
    }


    Processor::Result DAQProc::DSEvent(DS::Root *ds) {
      //This processor build waveforms for each PMT in the MC generated event, sample them and
      //store each sampled piece as a new event

      //Set DAQ header in run structure
      DS::DAQHeader *daqHeader = new DS::DAQHeader();
      daqHeader->SetAttribute("DAQ_NAME","VIRTUAL_DIGITIZER");
      daqHeader->SetAttribute("NBITS",fDigitizer->GetNBits());
      daqHeader->SetAttribute("TIME_RES",fDigitizer->GetStepTime());
      daqHeader->SetAttribute("V_OFFSET",fDigitizer->GetVOffSet());
      daqHeader->SetAttribute("V_HIGH",fDigitizer->GetVHigh());
      daqHeader->SetAttribute("V_LOW",fDigitizer->GetVLow());
      daqHeader->SetAttribute("RESISTANCE",fDigitizer->GetResistance());

      DS::Run *run = DS::RunStore::GetRun(ds);
      run->SetID(1);
      run->SetType(0x00001111);
      run->SetStartTime(1440638077);
      run->SetDAQHeader(daqHeader);


      DS::MC *mc = ds->GetMC();

      //Loop through the PMTs in the MC generated event
      for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

        DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);

        //For each PMT loop over hit photons and create a waveform for each of them
        PMTWaveform pmtwf;
        pmtwf.SetStepTime(fPulseStepTimeDB);
        pmtwf.SetSamplingWindow(fSamplingTimeDB);
        //      double PulseDuty=0.0;

        for (size_t iph=0; iph < mcpmt->GetMCPhotonCount(); iph++) {

          DS::MCPhoton *mcphotoelectron = mcpmt->GetMCPhoton(iph);
          double TimePhoton = mcphotoelectron->GetFrontEndTime();

          PMTPulse *pmtpulse;
          if (fPulseTypeDB==0){
            pmtpulse = new SquarePMTPulse; //square PMT pulses
          }
          else{
            pmtpulse = new RealPMTPulse; //real PMT pulses shape
          }

          pmtpulse->SetPulseMean(fPulseMeanDB);
          pmtpulse->SetStepTime(fPulseStepTimeDB);
          pmtpulse->SetPulseMin(fPulseMinDB);
          pmtpulse->SetPulseCharge(mcphotoelectron->GetCharge());
          pmtpulse->SetPulseWidth(fPulseWidthDB);
          pmtpulse->SetPulseOffset(fPulseOffsetDB);
          pmtpulse->SetPulseStartTime(TimePhoton); //also sets end time according to the pulse width and the pulse mean
          pmtwf.fPulse.push_back(pmtpulse);
          //	PulseDuty += pmtpulse->GetPulseEndTime() - pmtpulse->GetPulseStartTime();

        } // end mcphotoelectron loop: all pulses produced for this PMT

        //Sort pulses in time order
        std::sort(pmtwf.fPulse.begin(),pmtwf.fPulse.end(),Cmp_PMTPulse_TimeAscending);

        //Digitize waveform -- electronic noise is added by the digitizer. Save analogue
        //waveform in MCPMT object only for drawing purpose and save the digitized one
        //for posterior analysis

//        std::cout<<" CHARGE "<<mcpmt->GetCharge()<<" "<<-pmtwf.GetCharge(0.,200.)<<std::endl;
        mcpmt->SetWFCharge(pmtwf.GetCharge(0.,200.));//for debugging
        // fDigitizer->AddChannel(mcpmt->GetID(),pmtwf);
        // mcpmt->SetWaveform(fDigitizer->GetAnalogueWaveform(mcpmt->GetID()));
        // mcpmt->SetDigitizedWaveform(fDigitizer->GetDigitizedWaveform(mcpmt->GetID()));

      } //end pmt loop



      //////////////////////////////////////////////////////////
      //FROM HERE THE TRIGGER PROCESSOR SHOULD TAKE OVER!
      //1) Check trigger condition and divide waveforms in chunks
      //2) Build events containing digitized waveform samples and integrated charges
      //////////////////////////////////////////////////////////

      //Switch among triggers
      //simpledaq for debugging
      if(fTriggerType=="simpledaq"){

        DS::EV *ev = ds->AddNewEV(); //Remove it if no PMT cross threshold
        ev->SetID(fEventCounter);
        fEventCounter++; //simpledaq

        //simpleDAQ
        double totalQ = 0.0;
        double calibQ = 0.0;
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
          DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
          int pmtID = mcpmt->GetID();

          if (mcpmt->GetMCPhotonCount() > 0) {
            // Need at least one photon to trigger
            DS::PMT* pmt = ev->AddNewPMT();
            pmt->SetID(pmtID);

            // Create one sample, hit time is determined by first hit,
            // "infinite" charge integration time
            // WARNING: gets multiphoton effect right, but not walk correction
            // Write directly to calibrated waveform branch

            double time = mcpmt->GetMCPhoton(0)->GetFrontEndTime();
            double charge = 0;

            for (int i=0; i < mcpmt->GetMCPhotonCount(); i++)  {
              if (time > mcpmt->GetMCPhoton(i)->GetHitTime())
              time = mcpmt->GetMCPhoton(i)->GetHitTime();
              charge += mcpmt->GetMCPhoton(i)->GetCharge();
            }

            //pmt->SetCalibratedCharge(charge);
            totalQ += charge;

            //charge *= fSPECharge[pmtID] * 1e12; /* convert to pC */
            pmt->SetTime(time);
            pmt->SetCharge(charge);
            calibQ += charge;
          }
        }

        ev->SetTotalCharge(totalQ);

      }
      //First trigger type: store all PMTs crossing threshold
      else if(fTriggerType=="allpmts"){

        bool eventAdded = false;
        DS::EV *ev;

        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
          int pmtID = mc->GetMCPMT(imcpmt)->GetID();
          //Sample digitized waveform and look for triggers
          std::vector<UShort_t> DigitizedWaveform = fDigitizer->GetDigitizedWaveform(pmtID);
          for(int isample=0; isample<DigitizedWaveform.size(); isample++){
            if (DigitizedWaveform[isample]<=fDigitizer->GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)

              if(!eventAdded){
                ev = ds->AddNewEV();
                ev->SetID(fEventCounter);
                fEventCounter++;
                eventAdded = true;
              }

              DS::PMT* pmt = ev->AddNewPMT();
              pmt->SetID(pmtID);
              pmt->SetWaveform(fDigitizer->SampleWaveform(DigitizedWaveform,isample)); //it is defined by the sample that crosses threshold
              isample = fDigitizer->GoToEndOfSample(isample); //go forward towards the end of the sampling window
            }//end if above trigger
          }//end sampling
          DigitizedWaveform.clear(); //prune for next round of PMTs
        }//end PMT loop
        fDigitizer->Clear();

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
        if(triggerID>-1){ //means the trigger PMT exists in the event and hence we have a hit
          std::vector<UShort_t> DigitizedTriggerWaveform = fDigitizer->GetDigitizedWaveform(triggerID);
          //Sample the digitized waveform to look for triggers
          for(int isample=0; isample<DigitizedTriggerWaveform.size(); isample++){

//            std::cout<<" SAMPLE "<<isample<<"/"<<DigitizedTriggerWaveform.size()<<" "<<DigitizedTriggerWaveform[isample]<<" "<<fDigitizer->GetDigitizedThreshold()<<std::endl;

            if (DigitizedTriggerWaveform[isample]<=fDigitizer->GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)

              // std::cout<<" Above threshold "<<isample<<"/"<<DigitizedTriggerWaveform.size()<<" "<<DigitizedTriggerWaveform[isample]<<" "<<fDigitizer->GetDigitizedThreshold()<<std::endl;

              //Create a new event
              DS::EV *ev = ds->AddNewEV();
              ev->SetID(fEventCounter);
              fEventCounter++;

              //Read ALL PMTs
              for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
                int pmtID = mc->GetMCPMT(imcpmt)->GetID();
                DS::PMT* pmt = ev->AddNewPMT();
                pmt->SetID(pmtID);
                pmt->SetWaveform(fDigitizer->SampleWaveform(fDigitizer->GetDigitizedWaveform(pmtID), isample));
              } //end reading PMTs

              isample = fDigitizer->GoToEndOfSample(isample); //go forward towards the end of the sampling window
            } //end if: trigger above threshold

          }//end sampling

          DigitizedTriggerWaveform.clear(); //prune for next round of PMTs
        } //end if hit trigger PMT

        fDigitizer->Clear();
      } //end if second type of trigger

      // //If got at least one PMT above threshold move forward one event so it is not
      // //overwritten by the next one
      // if(ev->GetPMTCount()>0){
      //   fEventCounter++;
      // }

























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
