////////////////////////////////////////////////////////////////////
/// RAT::VertexGen_Decay0
/// \brief Vertex Generator - events consist of particle types, times and momenta 
///
/// \author Aleksandra Bialek <abialek@ualberta.ca> -- contact person
///
/// REVISION HISTORY:\n
///    30 Jun 2010  A.Bialek 
///       This class is used by the vertex generator to calculate 
///       momenta and time of initial particles.
///       At that moment only the double beta decays are included. 
///       All information of isotopes included in the Decay0.ratdb.
///    10 Oct 2010  A.Bialek
///       First 24 isotopes for simulation of the background&source events
///       Database file: Decay0Backg.ratdb
///    28 Oct 2010 A.Bialek
///       Additional 10 isotopes in the list.
///    09 Nov 2010 A.Bialek
///       Next 10 isotopes in the list.  
///    03 Dec 2010 A.Bialek
///       Last 14 isotopes 
///
/// \detail  List of the isotopes included in the generator:
///   Ac228   Ar39   Ar42   As79  Bi207  Bi208  Bi210  Bi212  Bi214    C14
///    Ca48  Cd113   Co60  Cs136  Cs137  Eu147  Eu152  Eu154  Gd146  Hf182
///    I126   I133   I134   I135    K40    K42   Kr81   Kr85   Mn54   Na22
///     P32 Pa234m  Pb210  Pb211  Pb212  Pb214  Po214  Ra228   Rb87  Rh106    
///   Sb125  Sb126  Sb133   Sr90  Ta182  Te133 Te133m  Te134  Th234  Tl207
///   Tl208 Xe129m Xe131m  Xe133  Xe135    Y88    Y90   Zn65   Zr96
///
/// To use generator, define in the macro:
///   /generator/add combo decay0:[POSITION]  #POSITION generator: point,fill,..
///
/// to simulate the double beta events:
///   /generator/vtx/set 2beta [ISOTOPE] [LEVEL] [MODE] [LoELim] [HiELim]
/// OR to simulate background events
///   /generator/vtx/set backg [ISOTOPE]
///
///For detail information go to User Manual, DocDB (SNO+-doc-592)
////////////////////////////////////////////////////////////////////
#ifndef __RAT_VertexGen_Decay0__
#define __RAT_VertexGen_Decay0__

#include <RAT/GLG4VertexGen.hh>
#include "RAT/DB.hh"
#include <RAT/Decay0.hh>
#include <fstream>
#include <G4Event.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>
#include <vector>
namespace RAT {

class VertexGen_Decay0 : public GLG4VertexGen {
  public:
    VertexGen_Decay0(const char *arg_dbname="decay0");
    virtual ~VertexGen_Decay0();
  
    virtual void  GeneratePrimaryVertex(G4Event *argEvent,
                                        G4ThreeVector &dx,
                                        G4double dt); 
    virtual void   SetState(G4String newValues );
    virtual        G4String GetState();
	virtual void   GetParentAZ(G4int &A1, G4int &Z1, G4int &A2, G4int &Z2);

  private:
    G4String  fCodeToName(G4int code);
    G4ParticleDefinition *fPDef;    // particle definition
    G4String fType;                 // defines or double beta decay "2beta" or 
                                    // background and source study "backg"
    G4String   fIsotope;            // parent isotope
    G4int fLevel;                   // daughter energy level
    G4int fMode;                    // decay mode
    G4float fLoE, fHiE;             // limit for energy spectrum
    
    DBLinkPtr fLdecay;

  protected:
    Decay0 *fDecay0; 
};

} // namespace RAT

#endif

