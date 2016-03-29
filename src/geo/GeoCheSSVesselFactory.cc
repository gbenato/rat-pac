#include <RAT/GeoCheSSVesselFactory.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <CLHEP/Units/SystemOfUnits.h>

using namespace std;

namespace RAT {

G4VSolid *GeoCheSSVesselFactory::ConstructSolid(DBLinkPtr table)
{
  string volume_name = table->GetIndex();
  const double r_min = 0.;
  const double r_max = table->GetD("r_max");
  const double size_z = table->GetD("size_z");
  G4Tubs *cup = new G4Tubs("cup", r_min * CLHEP::mm, r_max * CLHEP::mm, size_z * CLHEP::mm, 0., CLHEP::twopi);
  G4Box *flat = new G4Box("flat", (r_max - r_min)/2. * CLHEP::mm, size_z * CLHEP::mm, size_z * CLHEP::mm);

  G4RotationMatrix rotation;
  G4ThreeVector trans(((r_max - r_min)/2. + r_min + .5)* CLHEP::mm, 0., 0.);
  G4Transform3D *transf = new G4Transform3D(rotation,trans);
  G4UnionSolid *unionVolume = new G4UnionSolid("union0", cup, flat, *transf);

  trans.rotateZ(2*atan(size_z/r_max) * CLHEP::rad);
  rotation.rotateZ(2*atan(size_z/r_max) * CLHEP::rad);
  transf = new G4Transform3D(rotation,trans);
  unionVolume = new G4UnionSolid("union1", unionVolume, flat, *transf);

  trans.rotateZ(-4*atan(size_z/r_max) * CLHEP::rad);
  rotation.rotateZ(-4*atan(size_z/r_max) * CLHEP::rad);
  transf = new G4Transform3D(rotation,trans);
  unionVolume = new G4UnionSolid(volume_name, unionVolume, flat, *transf);

  return unionVolume;

}

} // namespace RAT
