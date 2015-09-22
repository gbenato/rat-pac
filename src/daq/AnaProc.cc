#include <vector>
#include <RAT/AnaProc.hh>
#include <RAT/DS/PMT.hh>

using namespace std;

namespace RAT {

  AnaProc::AnaProc() : Processor("analysis") {

    fLdaq = DB::Get()->GetLink("DAQ");
    ped_start = fLdaq->GetD("ped_start");
    ped_end = fLdaq->GetD("ped_end");
    max_spread = fLdaq->GetD("max_spread");
    int_start = fLdaq->GetD("int_start");
    int_end = fLdaq->GetD("int_end");

    gotDAQHeader = false;

  }

  Processor::Result AnaProc::DSEvent(DS::Root *ds) {
    //This processor retrieves the digitized waveforms of the triggered events EV
    //and performs some calculations like charge or time

    //Get DAQHeader (only first event)
    if(!gotDAQHeader){
      run = DS::RunStore::GetRun(ds);
      daqHeader = run->GetDAQHeader();
      //daqHeader->PrintAttributes();
      gotDAQHeader = true;
    }

    //Event loop: each MCEvent may have several real events
    for(int iev=0; iev<ds->GetEVCount(); iev++){
      DS::EV *ev = ds->GetEV(iev);

      for (int ipmt=0; ipmt < ev->GetPMTCount(); ipmt++){
        DS::PMT* pmt = ev->GetPMT(ipmt);
        std::vector<UShort_t> dWaveform = pmt->GetWaveform();
        pmt->SetTime(GetTimeAtPeak(dWaveform));
        pmt->SetCharge(IntegrateCharge(dWaveform));
      }//end PMT loop
    }

    return Processor::OK;

  } //AnaProc::DSEvent


  //Calculates the time at which the peak of the digitized waveform occurs
  double AnaProc::GetTimeAtPeak(std::vector<UShort_t> dWaveform){

    double fStepTime = daqHeader->GetDoubleAttribute("TIME_RES");

    //Sample waveform to look for the maximum
    double sampleatpeak = -9999.;
    int voltatstep = 9999.;
    for(int isample=0; isample<dWaveform.size(); isample++){
      if(dWaveform[isample]<voltatstep){
        //	std::cout<<"GetTimeAtPeak "<<pmtID<<" "<<isample<<" "<<dWaveform[isample]<<std::endl;
        voltatstep = dWaveform[isample];
        sampleatpeak = isample;
      }
    }

    return sampleatpeak*fStepTime;
  }


  double AnaProc::IntegrateCharge(std::vector<UShort_t> dWaveform){

    double fStepTime = daqHeader->GetDoubleAttribute("TIME_RES");
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor(ped_start/fStepTime);
    int s_ped_end = floor(ped_end/fStepTime);
    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = dWaveform[isample]*voltsperadc + fVLow - fVOffSet;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(ped_end-ped_start);
    // if (ped_max - ped_min > mv_max_spread) continue;

    //Integrate charge and sustract pedestal
    int s_int_start = floor(int_start/fStepTime);
    int s_int_end = floor(int_end/fStepTime);
    double charge = 0.;
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      charge += ((double)dWaveform[isample]*voltsperadc + fVLow - fVOffSet - ped_mean)*fStepTime/fResistance; //ADC to charge
    //  std::cout<<" charge "<<isample<<"/"<<dWaveform.size()<<" "<<charge<<std::endl;
    }

    charge=std::abs(charge); //pulses are negative so covert to positive charge
    // std::cout<<" Integrated charge "<<charge<<" pedestal "<<ped_mean<<" "<<s_int_start<<" "<<s_int_end<<std::endl;

    return charge;
  }




} // namespace RAT
