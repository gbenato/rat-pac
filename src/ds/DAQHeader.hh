/**
 * @class DS::DAQHeader
 * Data Structure: DAQ attributes
 *
 * std::coutrmation about the DAQ system. To be kept as flexible as possible.
 */

#ifndef __RAT_DS_DAQHeader__
#define __RAT_DS_DAQHeader__

#include <iostream>

#include <RAT/Log.hh>

#include <TObject.h>
#include <TVector3.h>
#include <algorithm>

namespace RAT {
  namespace DS {

class DAQHeader : public TObject {
public:
  DAQHeader() : TObject() {}
  virtual ~DAQHeader() {}

  /** Attributes */
  virtual int GetIntAttribute(std::string attr){ return int_attr[attr];}
  virtual double GetDoubleAttribute(std::string attr){ return double_attr[attr];}
  virtual std::string GetStringAttribute(std::string attr){ return string_attr[attr];}

  virtual std::map<std::string,int> GetIntAttributes() const { return int_attr; }
  virtual std::map<std::string,double> GetDoubleAttributes() const { return double_attr; }
  virtual std::map<std::string,std::string> GetStringAttributes() const { return string_attr; }
  virtual void SetAttribute(std::string attr, int val) { int_attr[attr] = val; }
  virtual void SetAttribute(std::string attr, double val) { double_attr[attr] = val; }
  virtual void SetAttribute(std::string attr, std::string val) { string_attr[attr] = val; }
  virtual void PrintAttributes() {
    std::cout<<" Int attributes "<<"\n";
    for (std::map<std::string,int>::iterator attr = int_attr.begin(); attr != int_attr.end(); attr++) {
      std::cout<<"  |->"<<attr->first<<" = "<<attr->second<<"\n";
    }
    std::cout<<" Double attributes "<<"\n";
        for (std::map<std::string,double>::iterator attr = double_attr.begin(); attr != double_attr.end(); attr++) {
      std::cout<<"  |->"<<attr->first<<" = "<<attr->second<<"\n";
    }
    std::cout<<" String attributes "<<"\n";
    for (std::map<std::string,std::string>::iterator attr = string_attr.begin(); attr != string_attr.end(); attr++) {
      std::cout<<"  |->"<<attr->first<<" = "<<attr->second<<"\n";
    }
  }

  ClassDef(DAQHeader, 1)

protected:
  std::map<std::string,int> int_attr;
  std::map<std::string,double> double_attr;
  std::map<std::string,std::string> string_attr;

};

  } // namespace DS
} // namespace RAT

#endif
