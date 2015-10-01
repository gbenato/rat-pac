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
      fSamplingTime = fLdaq->GetD("sampling_time");
      //width of a PMT pulse in ns
      fPulseWidth = fLdaq->GetD("pulse_width");
      //offset of a PMT pulse in mV
      fPulseOffset = fLdaq->GetD("pulse_offset");
      //Minimum pulse height to consider
      fPulseMin = fLdaq->GetD("pulse_min");
      //Pulse type: 0=square pulses, 1=real pulses
      fPulseType = fLdaq->GetI("pulse_type");
      //mean of a PMT pulse in ns
      fPulseMean = fLdaq->GetD("pulse_mean");
      //pulse step time
      fPulseTimeStep = fLdaq->GetD("pulse_time_step");

      detail << "DAQProc: DAQ constants loaded" << newline;
      detail << "  PMT Pulse type: " << (fPulseType==0 ? "square" : "realistic") << newline;
      detail << dformat("  PMT Pulse Mean: ........................ %5.1f\n", fPulseMean);
      detail << dformat("  PMT Pulse Width: ....................... %5.1f ns\n", fPulseWidth);
      detail << dformat("  PMT Pulse Offset: ...................... %5.1f ADC Counts\n", fPulseOffset);
      detail << dformat("  Min PMT Pulse Height: .................. %5.1f mV\n", fPulseMin);

      fEventCounter = 0;
    }


    void DAQProc::SetS(std::string param, std::string value)
    {
      if(param=="trigger")
      fTriggerType = value;

      if(fTriggerType!="allpmts" && fTriggerType!="triggerpmt" && fTriggerType!="simple" && fTriggerType!="external"){
        std::cerr<<"DAQ: "<<fTriggerType<<" option unknown... EXIT "<<std::endl;
        exit(0);
      }
    }

    Processor::Result DAQProc::DSEvent(DS::Root *ds) {
      // This processor creates pulses per generated MC PE and store them in waveforms that are digitized
      // Then, it checks trigger conditions and store digitized samples of the waveforms in EVs

      //Set DAQ header in run structure
      DS::DAQHeader *daqHeader = new DS::DAQHeader();
      daqHeader->SetAttribute("DAQ_NAME","VIRTUAL_DIGITIZER");
      daqHeader->SetAttribute("NBITS",fDigitizer->GetNBits());
      daqHeader->SetAttribute("TIME_RES",fDigitizer->GetTimeStep());
      daqHeader->SetAttribute("V_OFFSET",fDigitizer->GetVOffSet());
      daqHeader->SetAttribute("V_HIGH",fDigitizer->GetVHigh());
      daqHeader->SetAttribute("V_LOW",fDigitizer->GetVLow());
      daqHeader->SetAttribute("RESISTANCE",fDigitizer->GetResistance());

      DS::Run *run = DS::RunStore::GetRun(ds);
      run->SetID(1);
      run->SetType(0x00001111);
      run->SetStartTime(1440638077);
      run->SetDAQHeader(daqHeader);
//      RAT::DS::PMTInfo *pmtInfo = run->GetPMTInfo();


      DS::MC *mc = ds->GetMC();

      // Generate pulses per PE in MCPMT
      if(fTriggerType!="simple"){
        //Loop through the PMTs in the MC generated event
        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

          DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);

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
          mcpmt->SetWFCharge(pmtwf.GetCharge(0.,200.));//for debugging
          fDigitizer->AddChannel(mcpmt->GetID(),pmtwf);
          mcpmt->SetWaveform(fDigitizer->GetAnalogueWaveform(mcpmt->GetID()));
          mcpmt->SetDigitizedWaveform(fDigitizer->GetDigitizedWaveform(mcpmt->GetID()));

        } //end pmt loop
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
          if(mctriggerpmt->GetType() != 0) continue; //didn't found the trigger PMT yet
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
            pmt->SetType(mcpmt->GetType());

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
      // Fire all PMTs: if a PMT crosses threshold, its trace is recorded. It doesn't provide global
      // timing as the reference time is the time at which a given PMT crosses threshold
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
              pmt->SetType(mc->GetMCPMT(imcpmt)->GetType());
              pmt->SetWaveform(fDigitizer->SampleWaveform(DigitizedWaveform,isample)); //it is defined by the sample that crosses threshold
              isample = fDigitizer->GoToEndOfSample(isample); //go forward towards the end of the sampling window
            }//end if above trigger
          }//end sampling
          DigitizedWaveform.clear(); //prune for next round of PMTs
        }//end PMT loop
        fDigitizer->Clear();

      }
      // External trigger: triggers all the channel at a time provided a condition. Like 'allpmts'
      // but with a common global timing
      // Used for cosmic triggers
      else if(fTriggerType=="external"){

        DS::EV *ev = ds->AddNewEV();
        ev->SetID(fEventCounter);

        for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){

          int pmtID = mc->GetMCPMT(imcpmt)->GetID();
          std::vector<UShort_t> DigitizedWaveform = fDigitizer->GetDigitizedWaveform(pmtID);
          DS::PMT* pmt = ev->AddNewPMT();
          pmt->SetID(pmtID);
          pmt->SetType(mc->GetMCPMT(imcpmt)->GetType());
          pmt->SetWaveform(fDigitizer->SampleWaveform(DigitizedWaveform)); //sample from the beggining of the signal window

          DigitizedWaveform.clear(); //prune for next round of PMTs

        }//end PMT loop
        fDigitizer->Clear();
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
                pmt->SetType(mc->GetMCPMT(imcpmt)->GetType());
                pmt->SetWaveform(fDigitizer->SampleWaveform(fDigitizer->GetDigitizedWaveform(pmtID), isample));
              } //end reading PMTs

              isample = fDigitizer->GoToEndOfSample(isample); //go forward towards the end of the sampling window
            } //end if: trigger above threshold

          }//end sampling

          DigitizedTriggerWaveform.clear(); //prune for next round of PMTs
        } //end if hit trigger PMT

        fDigitizer->Clear();
      } //end if second type of trigger

      return Processor::OK;

    } //DAQProc::DSEvent

  } // namespace RAT
