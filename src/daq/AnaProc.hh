#ifndef __RAT_AnaProc__
#define __RAT_AnaProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DS/DAQHeader.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DB.hh>

namespace RAT {


class AnaProc : public Processor {
public:
  AnaProc();
  virtual ~AnaProc();
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual double GetTimeAtPeak(std::vector<UShort_t>, RAT::DS::DAQHeader *daqHeader);
  virtual double GetTimeAtThreshold(std::vector<UShort_t>, RAT::DS::DAQHeader *daqHeader);
  virtual double IntegrateCharge(std::vector<UShort_t>, RAT::DS::DAQHeader *daqHeader);

protected:

  DS::Run *run;
  DS::DAQHeader *daqHeaderV1730;
  DS::DAQHeader *daqHeaderV1742;
  bool gotDAQHeader;

  DBLinkPtr fLAnalysis;

  double ped_start;
  double ped_end;
  double max_spread;
  double int_start;
  double int_end;
  double peak_window;
  double peak_qthres;
  double ped_max_fluc;

};


} // namespace RAT

#endif
