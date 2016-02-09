#include <RAT/DS/MCPMT.hh>

ClassImp(RAT::DS::MCPMT)

namespace RAT {
  namespace DS {
    
    Float_t MCPMT::GetCharge() const {
      Float_t charge = 0.0;
      for (unsigned int i=0; i < photon.size(); i++)
	charge += photon[i].GetCharge();
      return charge;
    }
    
    Float_t MCPMT::GetTime() const {
      Float_t time = -9999.;
      if(photon.size()>0)
	time = photon[0].GetHitTime();
      return time;
    }
    
    Float_t MCPMT::GetFrontEndTime() const {
      Float_t time = -9999.;
      if(photon.size()>0)
	time = photon[0].GetFrontEndTime();
      return time;
    }
    
  } // namespace DS
} // namespace RAT

