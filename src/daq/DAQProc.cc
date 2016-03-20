#include <vector>
#include <RAT/DAQProc.hh>
#include <RAT/DB.hh>
#include <G4ThreeVector.hh>
#include <RAT/DetectorConstruction.hh>
#include <RAT/PMTPulse.hh>
#include <RAT/PMTWaveform.hh>
#include <RAT/DS/DAQHeader.hh>

#include <CLHEP/Random/RandGauss.h>

using namespace std;

namespace RAT {

  inline bool Cmp_PMTPulse_TimeAscending(const PMTPulse *a, const PMTPulse *b)
    {
      double atime = a->GetPulseStartTime();
      double btime = b->GetPulseStartTime();
      return atime < btime;
    }

    DAQProc::DAQProc() : Processor("daq") {

      fLdaq = DB::Get()->GetLink("DAQ");
      fDigitizerV1730 = new Digitizer("V1730");
      fDigitizerV1742 = new Digitizer("V1742");

      //Pulses and waveforms
      fSamplingTime = fLdaq->GetD("sampling_time");
      fPulseWidth = fLdaq->GetD("pulse_width");
      fPulseOffset = fLdaq->GetD("pulse_offset");
      fPulseMin = fLdaq->GetD("pulse_min");
      fPulseMean = fLdaq->GetD("pulse_mean");
      fPulseTimeStep = fLdaq->GetD("pulse_time_step");

      //Trigger system
      fTriggerDelay = fLdaq->GetD("trigger_delay");
      fTriggerJitter = fLdaq->GetD("trigger_jitter");


      detail << "DAQProc: DAQ constants loaded" << newline;
      detail << dformat("  PMT Pulse Mean: ........................ %5.1f\n", fPulseMean);
      detail << dformat("  PMT Pulse Width: ....................... %5.1f ns\n", fPulseWidth);
      detail << dformat("  PMT Pulse Offset: ...................... %5.1f ADC Counts\n", fPulseOffset);
      detail << dformat("  Min PMT Pulse Height: .................. %5.1f mV\n", fPulseMin);

      fEventCounter = 0;
      setRun = false;
    }


    void DAQProc::SetS(std::string param, std::string value)
    {
      std::cout<<" DAQProc::SetS: "<<param<<" "<<value<<std::endl;

      if(param=="trigger") fTriggerType = value;

      if(fTriggerType!="allpmts" && fTriggerType!="triggerpmt" && fTriggerType!="simple" && fTriggerType!="external"){
        std::cerr<<"DAQ: "<<fTriggerType<<" option unknown... EXIT "<<std::endl;
        exit(0);
      }
    }

    Processor::Result DAQProc::DSEvent(DS::Root *ds) {
      // This processor creates pulses per generated MC PE and store them in waveforms that are digitized
      // Then, it checks trigger conditions and store digitized samples of the waveforms in EVs

      //Set DAQ headers in run structure
      DS::DAQHeader *daqHeaderV1730 = new DS::DAQHeader();
      daqHeaderV1730->SetAttribute("DAQ_NAME","V1730");
      daqHeaderV1730->SetAttribute("NBITS",fDigitizerV1730->GetNBits());
      daqHeaderV1730->SetAttribute("TIME_RES",fDigitizerV1730->GetTimeResolution());
      daqHeaderV1730->SetAttribute("TIME_DELAY",fDigitizerV1730->GetSampleDelay());
      daqHeaderV1730->SetAttribute("V_OFFSET",fDigitizerV1730->GetVOffSet());
      daqHeaderV1730->SetAttribute("V_HIGH",fDigitizerV1730->GetVHigh());
      daqHeaderV1730->SetAttribute("V_LOW",fDigitizerV1730->GetVLow());
      daqHeaderV1730->SetAttribute("RESISTANCE",fDigitizerV1730->GetResistance());
      DS::DAQHeader *daqHeaderV1742= new DS::DAQHeader();
      daqHeaderV1742->SetAttribute("DAQ_NAME","V1742");
      daqHeaderV1742->SetAttribute("NBITS",fDigitizerV1742->GetNBits());
      daqHeaderV1742->SetAttribute("TIME_RES",fDigitizerV1742->GetTimeResolution());
      daqHeaderV1742->SetAttribute("TIME_DELAY",fDigitizerV1742->GetSampleDelay());
      daqHeaderV1742->SetAttribute("V_OFFSET",fDigitizerV1742->GetVOffSet());
      daqHeaderV1742->SetAttribute("V_HIGH",fDigitizerV1742->GetVHigh());
      daqHeaderV1742->SetAttribute("V_LOW",fDigitizerV1742->GetVLow());
      daqHeaderV1742->SetAttribute("RESISTANCE",fDigitizerV1742->GetResistance());

      //Get PMTInfo and set Run (only first event)
      if(!setRun){
        run = DS::RunStore::GetRun(ds);
        if(run == NULL) {
          std::cout<<" Run not Found "<<std::endl;
          exit(0);
        }
        run->SetID(1);
        run->SetType(0x00001111);
        run->SetStartTime(1440638077);
        run->SetDAQHeader(daqHeaderV1730,"V1730");
        run->SetDAQHeader(daqHeaderV1742,"V1742");
        pmtInfo = run->GetPMTInfo();
        //daqHeader->PrintAttributes();

        setRun = true;
      }

      DS::MC *mc = ds->GetMC();

      // Generate pulses per PE in MCPMT
      // Only for !simple triggering cases
      if(fTriggerType!="simple"){
        //Loop through the PMTs in the MC generated event

        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

          DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
          int pmtID = mcpmt->GetID();
          int pmtType = pmtInfo->GetType(pmtID);

          //Loop over PEs and create a pulse for each one
          PMTWaveform pmtwf;
          pmtwf.SetStepTime(fPulseTimeStep);
          pmtwf.SetSamplingWindow(fSamplingTime);
          //      double PulseDuty=0.0;

          for (int iph=0; iph < mcpmt->GetMCPhotonCount(); iph++) {

            DS::MCPhoton *mcphotoelectron = mcpmt->GetMCPhoton(iph);
            double TimePhoton = mcphotoelectron->GetFrontEndTime();

            PMTPulse *pmtpulse;
            pmtpulse = new RealPMTPulse; //real PMT pulses shape

            pmtpulse->SetPulseMean(fPulseMean);
            pmtpulse->SetStepTime(fPulseTimeStep);
            pmtpulse->SetPulseMin(fPulseMin);
            pmtpulse->SetPulseCharge(mcphotoelectron->GetCharge());
            pmtpulse->SetPulseWidth(fPulseWidth);
            pmtpulse->SetPulseOffset(fPulseOffset);
            pmtpulse->SetPulseStartTime(TimePhoton); //also sets end time according to the pulse width and the pulse mean
            pmtwf.fPulse.push_back(pmtpulse);
            //	PulseDuty += pmtpulse->GetPulseEndTime() - pmtpulse->GetPulseStartTime();

          } // end PE loop

          //Sort pulses in time order
          std::sort(pmtwf.fPulse.begin(),pmtwf.fPulse.end(),Cmp_PMTPulse_TimeAscending);

          //Digitize waveform -- electronic noise is added by the digitizer. Save analogue
          //waveform in MCPMT object only for drawing purpose and save the digitized one
          //for posterior analysis
          //          std::cout<<" CHARGE "<<mcpmt->GetCharge()<<" "<<-pmtwf.GetCharge(0.,200.)<<std::endl;
          mcpmt->SetWFCharge(pmtwf.GetCharge(0.,200.)); //for debugging
          if(pmtType == 1 || pmtType == 3 || pmtType == 0){
            fDigitizerV1742->AddChannel(pmtID,pmtwf);
            //mcpmt->SetWaveform(fDigitizerV1742->GetAnalogueWaveform(pmtID));
            //mcpmt->SetDigitizedWaveform(fDigitizerV1742->GetDigitizedWaveform(pmtID));
          } else if(pmtType == 2){
            fDigitizerV1730->AddChannel(pmtID,pmtwf);
            //mcpmt->SetWaveform(fDigitizerV1730->GetAnalogueWaveform(pmtID));
            //mcpmt->SetDigitizedWaveform(fDigitizerV1730->GetDigitizedWaveform(pmtID));
          }

        } //end pmt loop

        //Simulate pulse in muon tags PMT
        if(fTriggerType=="external"){

          PMTWaveform pmtwf_top;
          PMTWaveform pmtwf_bottom;
          pmtwf_top.SetStepTime(fPulseTimeStep);
          pmtwf_top.SetSamplingWindow(fSamplingTime);
          pmtwf_bottom.SetStepTime(fPulseTimeStep);
          pmtwf_bottom.SetSamplingWindow(fSamplingTime);

          PMTPulse *pmtpulse;
          pmtpulse = new RealPMTPulse; //real PMT pulses shape

          //Apply external trigger jitter
          double triggerTime = fTriggerDelay + CLHEP::RandGauss::shoot()*fTriggerJitter;
          double muon_tof = 10./30. + 0.3; // cm/(cm/ns) + cable delay

          pmtpulse->SetPulseMean(fPulseMean);
          pmtpulse->SetStepTime(fPulseTimeStep);
          pmtpulse->SetPulseMin(fPulseMin);
          pmtpulse->SetPulseCharge(4000.);
          pmtpulse->SetPulseWidth(fPulseWidth);
          pmtpulse->SetPulseOffset(fPulseOffset);
          pmtpulse->SetPulseStartTime(triggerTime); //also sets end time according to the pulse width and the pulse mean
          pmtwf_top.fPulse.push_back(pmtpulse);
          fDigitizerV1742->AddChannel(6,pmtwf_top);

          pmtpulse->SetPulseStartTime(triggerTime+muon_tof);
          pmtwf_bottom.fPulse.push_back(pmtpulse);
          fDigitizerV1742->AddChannel(7,pmtwf_bottom);

        }

      }

      //////////////////////////////////////////////////////////
      //FROM HERE THE TRIGGER PROCESSOR SHOULD TAKE OVER!
      //////////////////////////////////////////////////////////

      // Simple trigger: like simpledaq but triggering off one user-defined trigger PMT
      // identified by type==0
      if(fTriggerType=="simple"){

        DS::EV *ev; //will add a new event if trigger PMT fires

        double totalQ = 0.0;
        //Get trigger PMT id
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {

          // std::cout<<" simple "<<imcpmt<<"/"<<mc->GetMCPMTCount()<<std::endl;

          DS::MCPMT *mctriggerpmt = mc->GetMCPMT(imcpmt);
          int triggerID = mctriggerpmt->GetID();
          if(pmtInfo->GetType(triggerID) != 0) continue; //didn't found the trigger PMT yet
          //Found trigger PMT!
          if (mctriggerpmt->GetMCPhotonCount() == 0) break;
          //Trigger PMT has PEs!
          ev = ds->AddNewEV();
          ev->SetID(fEventCounter);

          for (int jmcpmt=0; jmcpmt < mc->GetMCPMTCount(); jmcpmt++) {
            DS::MCPMT *mcpmt = mc->GetMCPMT(jmcpmt);
            int pmtID = mcpmt->GetID();

            DS::PMT* pmt = ev->AddNewPMT();
            pmt->SetID(pmtID);
            pmt->SetType(pmtInfo->GetType(pmtID));

            double time = mcpmt->GetMCPhoton(0)->GetFrontEndTime();
            double charge = 0;

            for (int ipe=0; ipe < mcpmt->GetMCPhotonCount(); ipe++)  {
              if (time > mcpmt->GetMCPhoton(ipe)->GetHitTime())
              time = mcpmt->GetMCPhoton(ipe)->GetHitTime();
              charge += mcpmt->GetMCPhoton(ipe)->GetCharge();
            }

            totalQ += charge;

            pmt->SetTime(time);
            pmt->SetCharge(charge);

          } //end adding PMTs to EV

          ev->SetTotalCharge(totalQ);

          //if got here it found trigger PMT -> break now
          fEventCounter++;
          break;

        } //end look for the trigger PMT

      }
      // Fire all PMTs: if any PMT crosses threshold, its trace is recorded. It doesn't provide global
      // timing as the reference time is the time at which a given PMT crosses threshold
      else if(fTriggerType=="allpmts"){

        bool eventAdded = false;
        DS::EV *ev;

        RAT::Digitizer *digitizer;
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

          int pmtID = mc->GetMCPMT(imcpmt)->GetID();
          int pmtType = pmtInfo->GetType(pmtID);

          if(pmtType == 1 || pmtType == 3 || pmtType == 0) digitizer = fDigitizerV1742;
          else if(pmtType == 2) digitizer = fDigitizerV1730;

          //Sample digitized waveform and look for triggers
          std::vector<UShort_t> DigitizedWaveform = digitizer->GetDigitizedWaveform(pmtID);
          for(int isample=0; isample<DigitizedWaveform.size(); isample++){

            if (DigitizedWaveform[isample]<=digitizer->GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)

              if(!eventAdded){
                ev = ds->AddNewEV();
                ev->SetID(fEventCounter);
                fEventCounter++;
                eventAdded = true;
              }

              DS::PMT* pmt = ev->AddNewPMT();
              pmt->SetID(pmtID);
              pmt->SetWaveform(digitizer->SampleWaveform(DigitizedWaveform,isample)); //it is defined by the sample that crosses threshold
              pmt->SetWaveformTime(digitizer->GetWaveformTime(pmtID));

              isample = digitizer->GoToEndOfSample(isample); //go forward towards the end of the sampling window

            }//end if above trigger
          }//end sampling
          DigitizedWaveform.clear(); //prune for next round of PMTs
        }//end PMT loop
        fDigitizerV1730->Clear(); //Clear waveforms for the next round of hits
        fDigitizerV1742->Clear(); //Clear waveforms for the next round of hits

      }
      // External trigger: triggers all the channel at a time provided a condition. Like 'allpmts'
      // but with a common global timing
      // Used for cosmic triggers
      else if(fTriggerType=="external"){

        DS::EV *ev = ds->AddNewEV();
        ev->SetID(fEventCounter);

        RAT::Digitizer *digitizer;
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

          int pmtID = mc->GetMCPMT(imcpmt)->GetID();
          int pmtType = pmtInfo->GetType(pmtID);

          if(pmtType == 1 || pmtType == 3 || pmtType == 0) digitizer = fDigitizerV1742;
          else if(pmtType == 2) digitizer = fDigitizerV1730;

          DS::PMT* pmt = ev->AddNewPMT();
          pmt->SetID(pmtID);
          pmt->SetWaveform(digitizer->SampleWaveform(digitizer->GetDigitizedWaveform(pmtID),0));
          pmt->SetWaveformTime(digitizer->GetWaveformTime(pmtID));

        }//end PMT loop

        //Add muon tags FIXME: loop over digitizer channels
        DS::PMT* pmt = ev->AddNewPMT();
        pmt->SetID(6);
        pmt->SetWaveform(fDigitizerV1742->SampleWaveform(fDigitizerV1742->GetDigitizedWaveform(6),0));
        pmt->SetWaveformTime(fDigitizerV1742->GetWaveformTime(6));
        pmt = ev->AddNewPMT();
        pmt->SetID(7);
        pmt->SetWaveform(fDigitizerV1742->SampleWaveform(fDigitizerV1742->GetDigitizedWaveform(7),0));
        pmt->SetWaveformTime(fDigitizerV1742->GetWaveformTime(7));

        fDigitizerV1730->Clear(); //Clear waveforms for the next round of hits
        fDigitizerV1742->Clear(); //Clear waveforms for the next round of hits
        fEventCounter++;

      }
      // Trigger PMT: checks if special 'trigger' PMT (tagged as type==0) goes above threshold.
      // If so, records the traces of every channel
      // Used for source triggers
      else if(fTriggerType=="triggerpmt"){

        //Identify the trigger PMT
        int triggerID=-1;
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
          DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
          if(pmtInfo->GetType(mcpmt->GetID()) == 0) triggerID = mcpmt->GetID();
        }

        //If trigger PMT has been hit, sample its waveform and check if it crosses
        //threshold
        if(triggerID>-1){ //means the trigger PMT exists in the event and hence we have a hit

          std::vector<UShort_t> DigitizedTriggerWaveform = fDigitizerV1742->GetDigitizedWaveform(triggerID);
          //Sample the digitized waveform to look for triggers
          for(int isample=0; isample<DigitizedTriggerWaveform.size(); isample++){
            // std::cout<<" SAMPLE "<<isample<<"/"<<DigitizedTriggerWaveform.size()<<" "<<DigitizedTriggerWaveform[isample]<<" "<<digitizer->GetDigitizedThreshold()<<std::endl;

            if (DigitizedTriggerWaveform[isample]<=fDigitizerV1742->GetDigitizedThreshold()){ //hit above threshold! (remember the pulses are negative)

              double threstime = fDigitizerV1742->GetTimeAtSample(isample);
              // std::cout<<" Above threshold "<<isample<<"/"<<DigitizedTriggerWaveform.size()<<" "<<DigitizedTriggerWaveform[isample]<<" "<<fDigitizer->GetDigitizedThreshold()<<std::endl;
              // std::cout<<" Above threshold "<<threstime<<" "<<fDigitizerV1742->GetSampleAtTime(threstime)<<std::endl;

              //Create a new event
              DS::EV *ev = ds->AddNewEV();

              ev->SetID(fEventCounter);
              fEventCounter++;

              RAT::Digitizer *digitizer;
              //Read ALL PMTs
              for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
                int pmtID = mc->GetMCPMT(imcpmt)->GetID();
                int pmtType = pmtInfo->GetType(pmtID);

                if(pmtType == 1 || pmtType ==3 || pmtType == 0) digitizer = fDigitizerV1742;
                else if(pmtType == 2 || pmtType == 4) digitizer = fDigitizerV1730;

                DS::PMT* pmt = ev->AddNewPMT();
                pmt->SetID(pmtID);
                int thressample = digitizer->GetSampleAtTime(threstime);
                pmt->SetWaveform(digitizer->SampleWaveform(digitizer->GetDigitizedWaveform(pmtID), thressample));
                pmt->SetWaveformTime(digitizer->GetWaveformTime(pmtID));
              } //end reading PMTs

              isample = fDigitizerV1742->GoToEndOfSample(isample); //go forward towards the end of the sampling window
            } //end if: trigger above threshold

          }//end sampling

          DigitizedTriggerWaveform.clear(); //prune for next round of PMTs
        } //end if hit trigger PMT

        fDigitizerV1730->Clear(); //Clear waveforms for the next round of hits
        fDigitizerV1742->Clear(); //Clear waveforms for the next round of hits

      } //end if triggerpmt type of trigger

      return Processor::OK;

    } //DAQProc::DSEvent

  } // namespace RAT
