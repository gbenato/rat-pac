#ifndef __RAT_CosmicGen__
#define __RAT_CosmicGen__

#include <GLG4Gen.hh>
#include <G4ThreeVector.hh>

#include <globals.hh>

class G4Event;
class G4ParticleDefinition;
class GLG4TimeGen;
class GLG4PosGen;

namespace RAT {

  class CosmicGen : public GLG4Gen {
  public:
    CosmicGen();
    virtual ~CosmicGen();
    virtual void GenerateEvent(G4Event *event);
    virtual void ResetTime(double offset=0.0);
    virtual bool IsRepeatable() const { return true; };

    virtual void SetState(G4String state);
    virtual G4String GetState() const;


  protected:

    G4ThreeVector *sensiVolPos; //random point inside the sensitive volume defined by the user
    G4String sensiVolName;
    G4String stateStr;
    G4ParticleDefinition *muonm, *muonp;
    GLG4TimeGen* timeGen;


  };

} // namespace RAT

#endif // RAT_CosmicGen_h
