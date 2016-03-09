#include <cmath>
#include <algorithm>
#include <vector>
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/RandFlat.h>
#include <TMath.h>
#include <TH1.h>
#include <GausPMTCharge.hh>
#include <RAT/DB.hh>

namespace RAT {

GausPMTCharge::GausPMTCharge(int pmtid) {

  DBLinkPtr fLGaus = DB::Get()->GetLink("PMTGAUSCHARGE");
  gaus_mean = fLGaus->GetDArray("gaus_mean")[pmtid];
  gaus_sigma = fLGaus->GetDArray("gaus_sigma")[pmtid];

  info << "Setting up GausPMTCharge model for PMT "<<pmtid<<" "<<gaus_mean<<" "<<gaus_sigma<<"\n";

}

GausPMTCharge::~GausPMTCharge() {}

/** Returns charge for one photoelectron. */
double GausPMTCharge::PickCharge() const {
  double charge = gaus_mean + CLHEP::RandGauss::shoot()*gaus_sigma;
  if(charge < 0.) charge = 0.;
  return charge;
}

/** Value of charge PDF at charge q */
double GausPMTCharge::PDF(double q) const {
    return 1/sqrt(2.*3.14159)/gaus_sigma * exp(-1./2.*pow( (q-gaus_mean)/gaus_sigma ,2 ));
}

} // namespace RAT
