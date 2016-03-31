#include <RAT/GeoCheSSVesselFactory.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <CLHEP/Units/SystemOfUnits.h>

using namespace std;

namespace RAT {

G4VSolid *GeoCheSSVesselFactory::ConstructSolid(DBLinkPtr table)
{
  string volume_name = table->GetIndex();
  const double r_min = 0.;
  const double r_max = table->GetD("r_max");
  const double size_z = table->GetD("size_z");
  G4Tubs *cup = new G4Tubs("cup", r_min * CLHEP::mm, r_max * CLHEP::mm, (size_z + 3.18/2.) * CLHEP::mm, 0., CLHEP::twopi);
  G4Tubs *hollow = new G4Tubs("hollow", r_min * CLHEP::mm, 20.0 * CLHEP::mm, 3.81 * CLHEP::mm, 0., CLHEP::twopi);
  G4Box *flat = new G4Box("flat", 5.0 * CLHEP::mm, size_z * CLHEP::mm, size_z * CLHEP::mm);

  G4ThreeVector *trans = new G4ThreeVector(0., 0., (size_z + 3.18/2.) * CLHEP::mm);
  G4RotationMatrix *rotation = new G4RotationMatrix();
  G4Transform3D *transf = new G4Transform3D(*rotation, *trans);
  G4SubtractionSolid *subtractionVolume = new G4SubtractionSolid("subtr0", cup, hollow, *transf);

  trans = new G4ThreeVector((r_max - 5.0) * CLHEP::mm, 0., -3.18/2.);
  transf = new G4Transform3D(*rotation,*trans);
  G4UnionSolid *unionVolume = new G4UnionSolid("union0", cup, flat, *transf);
  //G4UnionSolid *unionVolume = new G4UnionSolid("union0", subtractionVolume, flat, *transf);

  trans->rotateZ(2*atan(size_z/r_max) * CLHEP::rad);
  rotation->rotateZ(2*atan(size_z/r_max) * CLHEP::rad);
  transf = new G4Transform3D(*rotation,*trans);
  unionVolume = new G4UnionSolid("union1", unionVolume, flat, *transf);

  trans->rotateZ(-4*atan(size_z/r_max) * CLHEP::rad);
  rotation->rotateZ(-4*atan(size_z/r_max) * CLHEP::rad);
  transf = new G4Transform3D(*rotation,*trans);
  unionVolume = new G4UnionSolid(volume_name, unionVolume, flat, *transf);

  return unionVolume;

}

} // namespace RAT
