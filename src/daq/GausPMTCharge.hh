///////////////////////////////////////////////////////////////////////////////
/// \class RAT::MiniCleanPMTCharge
///
/// \brief  Implementation of PMTCharge using MiniCLEAN's parameterization.
///
/// \author Benjamin Land <benland100@berkeley.edu>
///
/// REVISION HISTORY:\n
///     2015-01-07 : B Land - Added doxygen header block \n
///
/// \details See MiniCLEAN docs for details. This is the default model if no
///          other model is specified.
///
///////////////////////////////////////////////////////////////////////////////

#ifndef __RAT_GausPMTCharge__
#define __RAT_GausPMTCharge__

#include <vector>
#include <cstddef>
#include <CLHEP/Random/RandGeneral.h>
#include "TH1.h"
#include <RAT/PMTCharge.hh>

namespace RAT {

class GausPMTCharge : public PMTCharge {
public:
  GausPMTCharge(int pmtid=0);
  virtual ~GausPMTCharge();

  /** Returns charge for one photoelectron. */
  virtual double PickCharge() const;

  /** Value of charge PDF at charge q (not normalized) */
  virtual double PDF(double q) const;

protected:
  double gaus_mean;
  double gaus_sigma;
};

} // namespace RAT

#endif
