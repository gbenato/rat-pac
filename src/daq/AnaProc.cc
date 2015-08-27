#include <vector>
#include <RAT/AnaProc.hh>
#include <RAT/DS/PMT.hh>

using namespace std;

namespace RAT {

  AnaProc::AnaProc() : Processor("analysis") {

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
    double charge = 0.;
    int start_sample = 0;
    while(start_sample<dWaveform.size()){
      charge += ((double)dWaveform[start_sample]*voltsperadc + fVLow - fVOffSet)*fStepTime/fResistance; //ADC to charge
      start_sample++;
    }
    //    std::cout<<" Digitized integrated charge "<<charge<<std::endl;

    charge=std::abs(charge); //pulses are negative so covert to positive charge
    return charge;
  }




} // namespace RAT
