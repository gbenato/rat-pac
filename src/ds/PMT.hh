/**
 * @class PMT
 * Data Structure: PMT in triggered event
 *
 * This represents a PMT in a detector event.
 */

#ifndef __RAT_DS_PMT__
#define __RAT_DS_PMT__

#include <Rtypes.h>

namespace RAT {
  namespace DS {

class PMT : public TObject {
public:
  PMT() : TObject() {}
  virtual ~PMT() {}

  /** ID number of PMT */
  virtual void SetID(Int_t _id) { this->id = _id; }
  virtual Int_t GetID() { return id; }

  /** Total charge in waveform (pC) */
  virtual void SetCharge(Float_t _charge) { this->charge = _charge; }
  virtual Float_t GetCharge() { return charge; }

  /** Hit time in ns */
  virtual void SetTime(Float_t _time) { this->time = _time; }
  virtual Float_t GetTime() { return time; }

  /** If it was above threshold or not */
  virtual void SetAboveThreshold(bool _AboveThreshold) { this->AboveThreshold=_AboveThreshold; }
  virtual bool IsAboveThreshold() { return AboveThreshold; }

  /** Digitzed and sampled waveform */
  virtual void SetWaveform(std::vector<unsigned short int> _waveform) {fWaveform = _waveform; }
  virtual std::vector<unsigned short int> GetWaveform() { return fWaveform; }


  ClassDef(PMT, 1);

protected:
  Int_t id;
  Float_t charge;
  Float_t time;
  std::vector<unsigned short int> fWaveform;
  bool AboveThreshold;
};

  } // namespace DS
} // namespace RAT

#endif
