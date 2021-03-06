#ifndef __RAT_InHDF5Producer__
#define __RAT_InHDF5Producer__

#include <string>
#include <RAT/Producer.hh>
#include <RAT/DB.hh>
#include <globals.hh>

class G4UIcmdWithAString;
class G4UIcommand;

namespace RAT {


class InHDF5Producer : public Producer {
public:
  InHDF5Producer();
  InHDF5Producer(ProcBlock *block);
  virtual ~InHDF5Producer();

  virtual bool ReadEvents(G4String filename);

  // override G4UImessenger (from Producer) methods
  virtual G4String GetCurrentValue(G4UIcommand * command);
  virtual void SetNewValue(G4UIcommand * command,G4String newValue);


protected:
  void Init();
  std::string filename;

  G4UIcmdWithAString *readCmd;
  G4UIcommand *readDefaultCmd;

  DBLinkPtr fLCalibV1742;
  std::vector<int> FastChtoID;
  std::vector<int> SlowChtoID;
  std::map<int,int> IDToDAQ;

};

} // namespace RAT

#endif
