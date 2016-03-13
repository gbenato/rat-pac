#include <cmath>
#include <algorithm>
#include <vector>
#include <Randomize.hh>
#include <RAT/GausPMTTime.hh>
#include <RAT/DB.hh>
#include <RAT/Log.hh>

using namespace std;

namespace RAT {

GausPMTTime::GausPMTTime(int pmtid) {

  DBLinkPtr fLGaus = DB::Get()->GetLink("PMTGAUSTIME");

  fCableDelay = fLGaus->GetDArray("cable_delay")[pmtid];
  fTransitTime = fLGaus->GetDArray("transit_time")[pmtid];
  fTTS = fLGaus->GetDArray("tts")[pmtid];
  fJitter = fLGaus->GetDArray("jitter")[pmtid];

  info << "Setting up GausPMTTime model for PMT "<<pmtid<<" "<<fCableDelay<<" "<<fTransitTime<<" "<<fTTS<<" "<<fJitter<<"\n";

}

GausPMTTime::~GausPMTTime() {}

double GausPMTTime::PickTime(double phtime) const {
  double fetime = phtime + fCableDelay + fTransitTime + fTTS*CLHEP::RandGauss::shoot() + fJitter*CLHEP::RandGauss::shoot();
  if(fetime<0.) fetime = 0.;
  return fetime;
}

} // namespace RAT
