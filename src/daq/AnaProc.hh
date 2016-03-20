#ifndef __RAT_AnaProc__
#define __RAT_AnaProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DS/DAQHeader.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DB.hh>

namespace RAT {

struct AnaParams {

  bool prune_wf;
  int ped_start;
  int ped_end;
  double max_spread;
  int int_start;
  int int_end;
  int qshort_ped;
  int qshort_int;

  double peak_window;
  double peak_qthres;
  double ped_max_fluc;
  double time_thres;
  double time_thres_frac;

};

class AnaProc : public Processor {
public:
  AnaProc();
  virtual ~AnaProc();
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual double GetTimeAtPeak(std::vector<UShort_t>, std::vector<double>, RAT::DS::DAQHeader*, AnaParams);
  virtual double GetTimeAtThreshold(std::vector<UShort_t>, std::vector<double>, RAT::DS::DAQHeader*, AnaParams);
  virtual double GetTimeAtFraction(std::vector<UShort_t>, std::vector<double>, RAT::DS::DAQHeader*, AnaParams);
  virtual double IntegrateCharge(std::vector<UShort_t>, std::vector<double>, RAT::DS::DAQHeader*, AnaParams);
  virtual double IntegrateQShort(std::vector<UShort_t>, std::vector<double>, RAT::DS::DAQHeader*, AnaParams);

protected:

  DS::Run *run;
  DS::DAQHeader *daqHeaderV1730;
  DS::DAQHeader *daqHeaderV1742;
  RAT::DS::PMTInfo *pmtInfo;
  bool gotDAQHeader;

  DBLinkPtr fLAnalysis;

  AnaParams anaV1730, anaV1742;

};


} // namespace RAT

#endif
