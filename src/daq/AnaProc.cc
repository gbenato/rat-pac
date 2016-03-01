#include <vector>
#include <RAT/AnaProc.hh>
#include <RAT/Digitizer.hh>
#include <RAT/DS/PMT.hh>

using namespace std;

namespace RAT {

  AnaProc::AnaProc() : Processor("analysis") {

    fLAnalysis = DB::Get()->GetLink("ANALYSIS","V1730");
    anaV1730.ped_start = fLAnalysis->GetD("ped_start");
    anaV1730.ped_end = fLAnalysis->GetD("ped_end");
    anaV1730.max_spread = fLAnalysis->GetD("max_spread");
    anaV1730.int_start = fLAnalysis->GetD("int_start");
    anaV1730.int_end = fLAnalysis->GetD("int_end");
    anaV1730.peak_window = fLAnalysis->GetD("peak_window");
    anaV1730.peak_qthres = fLAnalysis->GetD("peak_qthres");
    anaV1730.ped_max_fluc = fLAnalysis->GetD("ped_max_fluc");
    anaV1730.time_thres = fLAnalysis->GetD("time_thres");
    anaV1730.time_thres_frac = fLAnalysis->GetD("time_thres_frac");

    fLAnalysis = DB::Get()->GetLink("ANALYSIS","V1742");
    anaV1742.ped_start = fLAnalysis->GetD("ped_start");
    anaV1742.ped_end = fLAnalysis->GetD("ped_end");
    anaV1742.max_spread = fLAnalysis->GetD("max_spread");
    anaV1742.int_start = fLAnalysis->GetD("int_start");
    anaV1742.int_end = fLAnalysis->GetD("int_end");
    anaV1742.peak_window = fLAnalysis->GetD("peak_window");
    anaV1742.peak_qthres = fLAnalysis->GetD("peak_qthres");
    anaV1742.ped_max_fluc = fLAnalysis->GetD("ped_max_fluc");
    anaV1742.time_thres = fLAnalysis->GetD("time_thres");
    anaV1742.time_thres_frac = fLAnalysis->GetD("time_thres_frac");

    gotDAQHeader = false;

  }

  AnaProc::~AnaProc()
  {
    /* Do nothing */
  }

  Processor::Result AnaProc::DSEvent(DS::Root *ds) {
    //This processor retrieves the digitized waveforms of the triggered events EV
    //and performs some calculations to get the charge and time

    //Get DAQHeader (only first event)
    if(!gotDAQHeader){
      run = DS::RunStore::GetRun(ds);
      if(run == NULL) {
        std::cout<<" Run not Found "<<std::endl;
        exit(0);
      }
      daqHeaderV1730 = run->GetDAQHeader("V1730");
      daqHeaderV1742 = run->GetDAQHeader("V1742");
      pmtInfo = run->GetPMTInfo();
      //daqHeader->PrintAttributes();
      gotDAQHeader = true;
    }

    //Event loop: each MCEvent may have several real events
    for(int iev=0; iev<ds->GetEVCount(); iev++){
      DS::EV *ev = ds->GetEV(iev);
      for (int ipmt=0; ipmt < ev->GetPMTCount(); ipmt++){
        DS::PMT* pmt = ev->GetPMT(ipmt);
        int pmtID = pmt->GetID();
        int pmtType = pmtInfo->GetType(pmtID);
        std::vector<UShort_t> dWaveform = pmt->GetWaveform();
        std::vector<double> dWaveformTime = pmt->GetWaveformTime();
        if(pmtType==1 || pmtType==3 || pmtType==0){
          //          pmt->SetTime(GetTimeAtPeak(dWaveform), daqHeaderV1742 );
          //          pmt->SetTime( GetTimeAtThreshold(dWaveform, daqHeaderV1742, anaV1742) );
          pmt->SetTime( GetTimeAtFraction(dWaveform, dWaveformTime, daqHeaderV1742, anaV1742) );
          pmt->SetCharge(IntegrateCharge(dWaveform, dWaveformTime, daqHeaderV1742, anaV1742) );
        } else if(pmtType==2 || pmtType==4){
          //        pmt->SetTime(GetTimeAtPeak(dWaveform), daqHeaderV1730 );
          //        pmt->SetTime( GetTimeAtThreshold(dWaveform, daqHeaderV1730, anaV1730) );
          pmt->SetTime( GetTimeAtFraction(dWaveform, dWaveformTime, daqHeaderV1730, anaV1730) );
          pmt->SetCharge(IntegrateCharge(dWaveform, dWaveformTime, daqHeaderV1730, anaV1730) );
        }
      }//end PMT loop

    }

    return Processor::OK;

  } //AnaProc::DSEvent

  //Calculates the time at which the waveform crosses threshold
  double AnaProc::GetTimeAtFraction(std::vector<UShort_t> dWaveform, std::vector<double> dWaveformTime, RAT::DS::DAQHeader *daqHeader, AnaParams anaParams){

    std::string daqName = daqHeader->GetStringAttribute("DAQ_NAME");
    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fTimeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    fTimeDelay = 0.;
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = 1.;//daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor((anaParams.ped_start + fTimeDelay)/fTimeStep);
    int s_ped_end = floor((anaParams.ped_end + fTimeDelay)/fTimeStep);

    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > anaParams.ped_max_fluc) {
      std::cout<<" Pedestal above fluctuations "<<ped_max<<" "<<ped_min<<" "<<anaParams.ped_max_fluc<<std::endl;
      return -9999.;
    }

    int s_int_start = floor((anaParams.int_start + fTimeDelay)/fTimeStep);
    int s_int_end = floor((anaParams.int_end + fTimeDelay)/fTimeStep);

    vector<double> dWaveformPedCorr;
    for (size_t isample = 0; isample < dWaveform.size(); isample++) {
      dWaveformPedCorr.push_back( (dWaveform[isample]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance - ped_mean );
    }

    //Find the max and define threshold given the fraction
    double high_peak = 99999.;
    for(UShort_t isample=s_int_start; isample<s_int_end; isample++){
      //      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<"/"<<values.size()<<": "<<values[isample]<<" "<<Vthres<<std::endl;
      //      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<"/"<<values.size()<<": "<<dWaveform[isample]<<std::endl;
      if(dWaveformPedCorr[isample] < high_peak) {
        high_peak = dWaveformPedCorr[isample];
      }
    }

    double Vthres = anaParams.time_thres;
    double VthresFrac = high_peak*anaParams.time_thres_frac;
    int s_af_thres = -1;
    for(UShort_t isample=s_int_start; isample<s_int_end; isample++){
      //      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<"/"<<values.size()<<": "<<values[isample]<<" "<<Vthres<<" "<<VthresFrac<<std::endl;
      if(dWaveformPedCorr[isample]<Vthres && dWaveformPedCorr[isample]<VthresFrac) {
        s_af_thres = isample;
        break;
      }
    }
    if(s_af_thres == -1) {
      //      std::cout<<" AnaProc::GetTimeAtThreshold: waveform did not cross threshold "<<std::endl;
      return -9999;
    }
    else{
      //Interpolate to get the right time at threshold
      return dWaveformTime[s_af_thres - 1] + ( (dWaveformTime[s_af_thres] - dWaveformTime[s_af_thres-1])/(dWaveformPedCorr[s_af_thres] - dWaveformPedCorr[s_af_thres -1]) ) * (VthresFrac - dWaveformPedCorr[s_af_thres-1]);
      //      return s_after_thres*fTimeStep - fTimeDelay;
    }


  }

  //Calculates the time at which the waveform crosses threshold
  double AnaProc::GetTimeAtThreshold(std::vector<UShort_t> dWaveform, std::vector<double> dWaveformTime, RAT::DS::DAQHeader *daqHeader, AnaParams anaParams){

    std::string daqName = daqHeader->GetStringAttribute("DAQ_NAME");
    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fTimeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    fTimeDelay = 0.;
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = 1.;//daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor((anaParams.ped_start + fTimeDelay)/fTimeStep);
    int s_ped_end = floor((anaParams.ped_end + fTimeDelay)/fTimeStep);

    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > anaParams.ped_max_fluc) {
      std::cout<<" Pedestal above fluctuations "<<ped_max<<" "<<ped_min<<" "<<anaParams.ped_max_fluc<<std::endl;
      return -9999.;
    }

    int s_int_start = floor((anaParams.int_start + fTimeDelay)/fTimeStep);
    int s_int_end = floor((anaParams.int_end + fTimeDelay)/fTimeStep);

    vector<double> values;
    values.resize(s_int_end - s_int_start);
    for (size_t isample = 0; isample < values.size(); isample++) {
      values[isample] = (dWaveform[isample+s_int_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance - ped_mean;
    }

    double Vthres = anaParams.time_thres;
    int sthres = 0;
    for(UShort_t isample=0; isample<values.size(); isample++){
      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<"/"<<values.size()<<": "<<values[isample]<<" "<<Vthres<<std::endl;
      //      std::cout<<" AnaProc::GetTimeAtThreshold: "<<isample<<"/"<<values.size()<<": "<<dWaveform[isample]<<std::endl;
      if(values[isample]<Vthres) {
        sthres = isample + s_int_start;
        break;
      }
    }
    if(sthres == 0) {
      //      std::cout<<" AnaProc::GetTimeAtThreshold: waveform did not cross threshold "<<std::endl;
      return -9999;
    }
    else{
      return sthres*fTimeStep - fTimeDelay;
    }

  }

  //Calculates the time at which the peak of the digitized waveform occurs
  double AnaProc::GetTimeAtPeak(std::vector<UShort_t> dWaveform, std::vector<double> dWaveformTime, RAT::DS::DAQHeader *daqHeader, AnaParams anaParams){

    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fTimeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    fTimeDelay = 0.;
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = 1.;//daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor((anaParams.ped_start + fTimeDelay)/fTimeStep);
    int s_ped_end = floor((anaParams.ped_end + fTimeDelay)/fTimeStep);

    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > anaParams.ped_max_fluc) {
      std::cout<<" Pedestal above fluctuations "<<ped_max<<" "<<ped_min<<" "<<anaParams.ped_max_fluc<<std::endl;
      return -9999.;
    }

    //Integrate in a given window to compare vs threshold
    const size_t halfwindow = floor(anaParams.peak_window/fTimeStep/2.0);
    double vsum = 0.0;
    int s_int_start = floor((anaParams.int_start + fTimeDelay)/fTimeStep);
    int s_int_end = floor((anaParams.int_end + fTimeDelay)/fTimeStep);
    vector<double> values,integral;
    values.resize(s_int_end-s_int_start);
    integral.resize(s_int_end-s_int_start);
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      integral[isample-s_int_start] = -vsum*fTimeStep/fResistance;
      vsum += (values[isample-s_int_start] = (dWaveform[isample]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance - ped_mean);
    }

    //Look for a window of width 'ns_window' that crosses threshold
    //If so, get the maximum within that window and store it as a peak
    vector<size_t> peaks;
    double acc_offset = 0.0;
    for (size_t i = 0; i < s_int_end-s_int_start; i++) {
//      std::cout<<"Integral: "<<i<<" "<<integral[i]<<" "<<acc_offset<<" "<<peak_qthres<<std::endl;
      if (integral[i] - acc_offset >= anaParams.peak_qthres) { //crossed threshold, search in +- halfwindow
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
      return (peaks[0]*fTimeStep - fTimeDelay); //FIXME: return time at first peak so far
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


  double AnaProc::IntegrateCharge(std::vector<UShort_t> dWaveform, std::vector<double> dWaveformTime, RAT::DS::DAQHeader *daqHeader, AnaParams anaParams){

    double fTimeStep = daqHeader->GetDoubleAttribute("TIME_RES");
    double fTimeDelay = daqHeader->GetDoubleAttribute("TIME_DELAY");
    fTimeDelay = 0.;
    double fVHigh = daqHeader->GetDoubleAttribute("V_HIGH");
    double fVLow = daqHeader->GetDoubleAttribute("V_LOW");
    double fVOffSet = daqHeader->GetDoubleAttribute("V_OFFSET");
    double fResistance = 1.;//daqHeader->GetDoubleAttribute("RESISTANCE");
    int fNBits = daqHeader->GetIntAttribute("NBITS");

    int nADCs = 1 << fNBits; //Calculate the number of adc counts
    double voltsperadc = (fVHigh - fVLow)/(double)nADCs;

    //Compute pedestal charge
    int s_ped_start = floor((anaParams.ped_start + fTimeDelay)/fTimeStep);
    int s_ped_end = floor((anaParams.ped_end + fTimeDelay)/fTimeStep);

    double ped_min,ped_max,ped_mean;
    ped_min = ped_max = ped_mean = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
    for (size_t isample = s_ped_start+1; isample < s_ped_end; isample++) {
      double voltage = (dWaveform[s_ped_start]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance;
      if (ped_min > voltage) ped_min = voltage;
      if (ped_max < voltage) ped_max = voltage;
      ped_mean += voltage;
    }
    ped_mean /= (double)(s_ped_end-s_ped_start);
    if (ped_max - ped_min > anaParams.ped_max_fluc) {
      std::cout<<" Pedestal above fluctuations "<<ped_max<<" "<<ped_min<<" "<<anaParams.ped_max_fluc<<std::endl;
      return -9999.;
    }

    //Integrate charge and sustract pedestal
    int s_int_start = floor((anaParams.int_start + fTimeDelay)/fTimeStep);
    int s_int_end = floor((anaParams.int_end + fTimeDelay)/fTimeStep);
    double charge = 0.;
    for (size_t isample = s_int_start; isample < s_int_end; isample++) {
      charge += ((double)dWaveform[isample]*voltsperadc + fVLow - fVOffSet)*fTimeStep/fResistance - ped_mean; //ADC to charge
      //    std::cout<<" charge "<<isample<<"/"<<dWaveform.size()<<" "<<charge<<std::endl;
    }

    charge *= -1; //pulses are negative so covert to positive charge
    //    std::cout<<" Integrated charge "<<charge<<" pedestal "<<ped_mean<<" "<<s_int_start<<" "<<s_int_end<<std::endl;

    return charge;
  }

} // namespace RAT
