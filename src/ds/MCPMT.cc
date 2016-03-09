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
      Float_t time = 9999.;
      for(int iph = 0; iph<photon.size(); iph++){
	if(photon[iph].GetHitTime() < time) time = photon[iph].GetHitTime();
      }
      return time;
    }
    
    Float_t MCPMT::GetFrontEndTime() const {
      Float_t time = 9999.;
      for(int iph = 0; iph<photon.size(); iph++){
	if(photon[iph].GetFrontEndTime() < time) time = photon[iph].GetFrontEndTime();
      }
      return time;
    }
    
  } // namespace DS
} // namespace RAT

