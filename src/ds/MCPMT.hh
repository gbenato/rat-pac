/** 
 * @class DS::MCPMT
 *  Data Structure: Hit PMT in Monte Carlo
 *
 *  @author Stan Seibert <sseibert@hep.upenn.edu>
 *
 *  This class represents a PMT in which at least one photoelectron
 *  was generated by an incident photon.
 */

#ifndef __RAT_DS_MCPMT__
#define __RAT_DS_MCPMT__

#include <vector>
#include <RAT/DS/MCPhoton.hh>
#include <RAT/Log.hh>
#include <RAT/DS/PMTWaveform.hh>

namespace RAT {
  namespace DS {

class MCPMT : public TObject {
public:
  MCPMT() : TObject() {}
  virtual ~MCPMT() {}

  /** ID number */
  virtual Int_t GetID() const { return id; };
  virtual void SetID(Int_t _id) { id = _id; };

  /** Charge */
  virtual Float_t GetCharge() const;

  /** PMT type */
  virtual Int_t GetType() const { return type; };
  virtual void SetType(Int_t _type) { type = _type; };

  /** List of photoelectrons created in this PMT. */
  MCPhoton* GetMCPhoton(Int_t i) { return &photon[i]; }
  Int_t GetMCPhotonCount() const { return photon.size(); }
  MCPhoton* AddNewMCPhoton() {
    photon.resize(photon.size() + 1);
    return &photon.back();
  }
  void PruneMCPhoton() { photon.resize(0); }
 
  /** PMT waveform */
  PMTWaveform* GetWaveform() { return &waveform; };
  void SetWaveform(PMTWaveform _waveform) {_waveform.SetGraph(); waveform = _waveform; };
  void AddDigitizedWaveform(std::vector<int> _digitwaveform) {fDigitWaveForm = _digitwaveform;};
  std::vector<int> GetDigitizedWaveform() {return fDigitWaveForm;};

  ClassDef(MCPMT, 1)
    
protected:
  Int_t id;
  Int_t type;
  std::vector<MCPhoton> photon;
  PMTWaveform waveform;
  std::vector<int> fDigitWaveForm;
  
};

  } // namespace DS
} // namespace RAT

#endif

