////////////////////////////////////////////////////////////////////////
/// \class RAT::ParentPrimaryVertexInformation
/// \brief  Concrete instance of G4 Information class to pass
///  information alongside PrimaryVertex.
///
/// \author Jeanne Wilson <j.r.wilson@qmul.ac.uk> -- contact person
///  
/// \detail
///  This class includes a G4PrimaryParticle object that stores the 
///  information required about the parent particle for a vertex.
///  This parent is not propagated but we wish to store the information 
///  into the Rat DS for later use - eg neutrino parent of ES event
////////////////////////////////////////////////////////////////////////


#ifndef __RAT_ParentPrimaryVertexInformation__
#define __RAT_ParentPrimaryVertexInformation__

#include <G4PrimaryParticle.hh>
#include <G4VUserPrimaryVertexInformation.hh>

namespace RAT {

class ParentPrimaryVertexInformation : public G4VUserPrimaryVertexInformation
{
   public:
     explicit ParentPrimaryVertexInformation() {parent = 0;};
	 
	 ~ParentPrimaryVertexInformation() {};
	 void Print() const {};
	 
	 inline bool ExistParentParticle() const
	 {
	   if(parent){
	     return true;
	   }else{
	     return false;
	   }
	 };
	 
	 inline G4PrimaryParticle *GetVertexParentParticle() const 
	 {
	   return parent;
	 };
	 
	 virtual void SetVertexParentParticle(G4PrimaryParticle *p) {parent = p;};
	 
   private:
	 G4PrimaryParticle *parent;
};
} // namespace RAT
#endif
