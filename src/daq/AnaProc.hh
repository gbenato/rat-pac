#ifndef __RAT_AnaProc__
#define __RAT_AnaProc__

#include <string>
#include <RAT/Processor.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DB.hh>
#include <CLHEP/Random/RandGeneral.h>
#include <RAT/Digitizer.hh>

namespace RAT {


class AnaProc : public Processor {
public:
  AnaProc();
  virtual ~AnaProc() { };
  virtual Processor::Result DSEvent(DS::Root *ds);
  virtual double GetTimeAtPeak(std::vector<unsigned short int>);
  virtual double IntegrateCharge(std::vector<unsigned short int>);

protected:
  DBLinkPtr fLdaq;


};


} // namespace RAT

#endif
