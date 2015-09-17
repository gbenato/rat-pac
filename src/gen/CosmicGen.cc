#include <RAT/CosmicGen.hh>

#include <RAT/GLG4PosGen.hh>
#include <RAT/GLG4TimeGen.hh>
#include <RAT/Factory.hh>
#include <RAT/GLG4StringUtil.hh>

#include <TF1.h>
#include <TMath.h>

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4RandomDirection.hh"
#include <G4Event.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4ParticleDefinition.hh>
#include <G4ThreeVector.hh>
#include <G4UnitsTable.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>

#include <CLHEP/Vector/LorentzVector.h>

#include <cstring>
#include <sstream>

#undef DEBUG

namespace RAT {

  CosmicGen::CosmicGen() :
    stateStr("")
  {
    // As in the combo generator, use a default time generator if the
    // user does not supply one.
    timeGen = new GLG4TimeGen_Poisson();
    muonm   = G4MuonMinus::MuonMinus();
    muonp   = G4MuonPlus::MuonPlus();
  }

  CosmicGen::~CosmicGen()
  {
    delete timeGen;
  }

  void CosmicGen::GenerateEvent(G4Event* event)
  {

    //Get a random point in the world surface
    G4Navigator* gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4VSolid* worldSolid = gNavigator->GetWorldVolume()->GetLogicalVolume()->GetSolid();
    std::vector<G4VPhysicalVolume*> sensiVol;
    std::vector<G4VSolid*> sensiSolid;
    std::vector<G4AffineTransform> local_to_global;
    std::vector<G4AffineTransform> global_to_local;
    for(int ivol=0; ivol<sensiVolName.size(); ivol++){
      sensiVol.push_back(gNavigator->LocateGlobalPointAndSetup(*sensiVolPos[ivol],0,false));
      sensiSolid.push_back(sensiVol[ivol]->GetLogicalVolume()->GetSolid());
      if(sensiVolName[ivol] != sensiSolid[ivol]->GetName()){
        std::cout<<" Selected point ["<<sensiVolPos[ivol]->x()<<", "<<sensiVolPos[ivol]->y()<<", "<<sensiVolPos[ivol]->z()
                 <<"] is not contained in volume "<<sensiVolName[ivol]<<" but in volume "<< sensiSolid[ivol]->GetName()<<"\n"
                 <<" Change point position or volume name!"<<std::endl;
        exit(0);
      }
      local_to_global.push_back(gNavigator->GetLocalToGlobalTransform());
      global_to_local.push_back(gNavigator->GetGlobalToLocalTransform());
    }

    //Get an angle and momentum following a given distribution
    double energy = 1000.0;
    //Get cos2 random distribution
    double phi = G4UniformRand()*CLHEP::twopi;
    TF1 *theta_dist = new TF1("f","cos(x)*cos(x)",TMath::Pi()/2.,TMath::Pi());

    G4ThreeVector startPos;
    G4double dist = 0.;
    G4ThreeVector dir(0.0,0.0,0.0);
    double theta = 0.;
    do{
      startPos = worldSolid->GetPointOnSurface();
      theta = theta_dist->GetRandom();
      dir = G4ThreeVector(sin(theta)*sin(phi),
                          sin(theta)*cos(phi),
                          cos(theta));

      dist = 0.;

      for(int ivol=0; ivol<sensiVolName.size(); ivol++){
        //  A solid has no position defined in space hence 'DistanceToIn'
        //method DO NOT take into account the global coordinates and we
        //need to work in locals and transform back to globals
        global_to_local[ivol].ApplyPointTransform(startPos); // convert to local coords
        dist += sensiSolid[ivol]->DistanceToIn(startPos,dir); // return kInfinity if intersection does not exist
        local_to_global[ivol].ApplyPointTransform(startPos); // convert back to global coords
      }

    } while(dist<0 || dist >= kInfinity);

    double mass = muonm->GetPDGMass();
    double mom = sqrt(energy*energy - mass*mass);
    G4ThreeVector vmom = dir;
    vmom.setMag(mom);

    G4PrimaryVertex* vertex = new G4PrimaryVertex(startPos, NextTime());
    //Pick a random sign for the muon
    G4PrimaryParticle* particle;
    if(G4UniformRand() >= 0.5) particle = new G4PrimaryParticle(muonm, vmom.x(), vmom.y(), vmom.z());
    else particle = new G4PrimaryParticle(muonp, vmom.x(), vmom.y(), vmom.z());
    particle->SetMass(muonm->GetPDGMass());
    vertex->SetPrimary(particle);
    event->AddPrimaryVertex(vertex);

  }

  void CosmicGen::ResetTime(double offset)
  {
    double eventTime = timeGen->GenerateEventTime();
    nextTime = eventTime + offset;
#ifdef DEBUG
    G4cout << "RAT::CosmicGen::ResetTime:"
	   << " eventTime=" << G4BestUnit(eventTime,"Time")
	   << ", offset=" << G4BestUnit(offset,"Time")
	   << ", nextTime=" << G4BestUnit(nextTime,"Time")
	   << G4endl;
#endif
  }

  void CosmicGen::SetState(G4String state)
  {

#ifdef DEBUG
    G4cout << "RAT::CosmicGen::SetState called with state='"
	   << state << "'" << G4endl;
#endif

    state = util_strip_default(state);
    std::vector<std::string> allVolumes = util_split(state, ":");
    for(int ivol=0; ivol<allVolumes.size(); ivol++){
      allVolumes[ivol] = util_strip_default(allVolumes[ivol]);
      std::vector<std::string> volume = util_split(allVolumes[ivol], " ");
      sensiVolPos.push_back(new G4ThreeVector(std::stoi(volume[0]),std::stoi(volume[1]),std::stoi(volume[2])));
      sensiVolName.push_back(volume[3]);
    }

  }

  G4String CosmicGen::GetState() const
  {
    return util_dformat("%ld\t%ld\t%ld", sensiVolPos[0]->x(), sensiVolPos[0]->y(), sensiVolPos[0]->z());
  }

} // namespace RAT
