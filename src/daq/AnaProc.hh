#ifndef __RAT_AnaProc__
#define __RAT_AnaProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DS/DAQHeader.hh>
#include <RAT/DS/RunStore.hh>

namespace RAT {


class AnaProc : public Processor {
public:
  AnaProc();
  virtual ~AnaProc() { };
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual double GetTimeAtPeak(std::vector<UShort_t>);
  virtual double IntegrateCharge(std::vector<UShort_t>);

protected:

  DS::Run *run;
  DS::DAQHeader *daqHeader;
  bool gotDAQHeader;

};


} // namespace RAT

#endif
