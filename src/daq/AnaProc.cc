#include <vector>
#include <RAT/AnaProc.hh>
#include <RAT/Digitizer.hh>
#include <RAT/DS/PMT.hh>

using namespace std;

namespace RAT {

  AnaProc::AnaProc() : Processor("analysis") {

    fLAnalysis = DB::Get()->GetLink("ANALYSIS");
    ped_start = fLAnalysis->GetD("ped_start");
    ped_end = fLAnalysis->GetD("ped_end");
    max_spread = fLAnalysis->GetD("max_spread");
    int_start = fLAnalysis->GetD("int_start");
    int_end = fLAnalysis->GetD("int_end");
    peak_window = fLAnalysis->GetD("peak_window");
    peak_qthres = fLAnalysis->GetD("peak_qthres");
    ped_max_fluc = fLAnalysis->GetD("ped_max_fluc");

    gotDAQHeader = false;

  }

  Processor::Result AnaProc::DSEvent(DS::Root *ds) {
    //This processor retrieves the digitized waveforms of the triggered events EV
    //and performs some calculations to get the charge and time

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
        pmt->SetTime(GetTimeAtThreshold(dWaveform));
        //        pmt->SetTime(GetTimeAtPeak(dWaveform));
        pmt->SetCharge(IntegrateCharge(dWaveform));
      }//end PMT loop
    }

    return Processor::OK;

  } //AnaProc::DSEvent

  //Calculates the time at which the waveform crosses threshold
  double AnaProc::GetTimeAtThreshold(std::vector<UShort_t> dWaveform){

    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fTimeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor(ped_start/fTimeStep);
    int s_ped_end = floor(ped_end/fTimeStep);
    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = dWaveform[isample]*voltsperadc + fVLow - fVOffSet;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > ped_max_fluc) return -9999.;

    //Correct by pedestals
    int s_int_start = floor(int_start/fTimeStep);
    int s_int_end = floor(int_end/fTimeStep);
    vector<double> values;
    values.resize(s_int_end-s_int_start);
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      values[isample-s_int_start] = dWaveform[isample]*voltsperadc + fVLow - fVOffSet - ped_mean;
    }

    Digitizer *digitizer = new Digitizer();
    double Vthres = digitizer->GetDigitizedThreshold();
    int sthres = 0;
    for(UShort_t isample=0; isample<values.size(); isample++){
//      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<" "<<values[isample]<<" "<<Vthres<<std::endl;
      if(values[isample]<Vthres) {
        sthres = isample;
        break;
      }
    }
    if(sthres == 0) {
      std::cout<<" AnaProc::GetTimeAtThreshold: waveform did not cross threshold "<<std::endl;
      return -9999;
    }
    else{
      return sthres*fTimeStep - fTimeDelay;
    }

  }

  //Calculates the time at which the peak of the digitized waveform occurs
  double AnaProc::GetTimeAtPeak(std::vector<UShort_t> dWaveform){

    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double timeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor(ped_start/fTimeStep);
    int s_ped_end = floor(ped_end/fTimeStep);
    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = dWaveform[isample]*voltsperadc + fVLow - fVOffSet;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > ped_max_fluc) return -9999.;

    //Integrate in a given window to compare vs threshold
    const size_t halfwindow = floor(peak_window/fTimeStep/2.0);
    double vsum = 0.0;
    int s_int_start = floor(int_start/fTimeStep);
    int s_int_end = floor(int_end/fTimeStep);
    vector<double> values,integral;
    values.resize(s_int_end-s_int_start);
    integral.resize(s_int_end-s_int_start);
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      integral[isample-s_int_start] = -vsum*fTimeStep/fResistance;
      vsum += (values[isample-s_int_start] = dWaveform[isample]*voltsperadc + fVLow - fVOffSet - ped_mean);
    }

    //Look for a window of width 'ns_window' that crosses threshold
    //If so, get the maximum within that window and store it as a peak
    vector<size_t> peaks;
    double acc_offset = 0.0;
    for (size_t i = 0; i < s_int_end-s_int_start; i++) {
//      std::cout<<"Integral: "<<i<<" "<<integral[i]<<" "<<acc_offset<<" "<<peak_qthres<<std::endl;
      if (integral[i] - acc_offset >= peak_qthres) { //crossed threshold, search in +- halfwindow
        const size_t find_start = i-halfwindow >= 0 ? i-halfwindow : 0;
        const size_t find_end = i+halfwindow < s_int_end-s_int_start ? i+halfwindow : s_int_end-s_int_start;
        double peak = values[find_start];
        size_t peak_idx = find_start;
        //Look for the maximum
        for (size_t j = find_start+1; j < find_end; j++) {
          if (values[j] < peak) {
            peak = values[j];
            peak_idx = j;
          }
        }
        peaks.push_back(peak_idx+s_int_start);
        acc_offset = integral[peak_idx+halfwindow];
        i = peak_idx+2*halfwindow-1;
      } else if (i-halfwindow > 0 && integral[i-halfwindow] > acc_offset) {
        acc_offset = integral[i-halfwindow];
      }
    }

    if(peaks.size()>0){
//      std::cout<<"Peak at: "<<peaks[0]<<std::endl;
      return (peaks[0]*fTimeStep - timeDelay); //FIXME: return time at first peak so far
    } else{
//      std::cout<<"Peak not found"<<std::endl;
      return -9999.;
    }

    // //Sample waveform to look for the maximum
    // double sampleatpeak = -9999.;
    // int voltatstep = 9999.;
    // for(int isample=0; isample<dWaveform.size(); isample++){
    //   if(dWaveform[isample]<voltatstep){
    //     //	std::cout<<"GetTimeAtPeak "<<pmtID<<" "<<isample<<" "<<dWaveform[isample]<<std::endl;
    //     voltatstep = dWaveform[isample];
    //     sampleatpeak = isample;
    //   }
    // }
    //
    // return sampleatpeak*fTimeStep;
  }


  double AnaProc::IntegrateCharge(std::vector<UShort_t> dWaveform){

    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor(ped_start/fTimeStep);
    int s_ped_end = floor(ped_end/fTimeStep);
    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = dWaveform[isample]*voltsperadc + fVLow - fVOffSet;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > ped_max_fluc) return -9999.;

    //Integrate charge and sustract pedestal
    int s_int_start = floor(int_start/fTimeStep);
    int s_int_end = floor(int_end/fTimeStep);
    double charge = 0.;
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      charge += ((double)dWaveform[isample]*voltsperadc + fVLow - fVOffSet - ped_mean)*fTimeStep/fResistance; //ADC to charge
    //  std::cout<<" charge "<<isample<<"/"<<dWaveform.size()<<" "<<charge<<std::endl;
    }

    charge=std::abs(charge); //pulses are negative so covert to positive charge
    // std::cout<<" Integrated charge "<<charge<<" pedestal "<<ped_mean<<" "<<s_int_start<<" "<<s_int_end<<std::endl;

    return charge;
  }




} // namespace RAT
