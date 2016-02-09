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

  /** Type of PMT */
  virtual void SetType(Int_t _type) { this->type = _type; }
  virtual Int_t GetType() { return type; }

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
   virtual void SetWaveform(std::vector<UShort_t> _waveform) {waveform = _waveform; }
   virtual std::vector<UShort_t> GetWaveform() { return waveform; }

 ClassDef(PMT, 1);

protected:
  Int_t id;
  Int_t type;
  Float_t charge;
  Float_t time;
  bool AboveThreshold;
  std::vector<UShort_t> waveform;
};

  } // namespace DS
} // namespace RAT

#endif
