#include <H5Cpp.h>
#include <RAT/InHDF5Producer.hh>
#include <RAT/ProcBlock.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/SignalHandler.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/DS/DAQHeader.hh>

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <assert.h>
#include <regex>

typedef std::map< int, std::vector< std::vector<UShort_t> > > mw_t; //channel:[event,waveform]

namespace RAT {


  InHDF5Producer::InHDF5Producer()
  {
    mainBlock = 0;
    Init();
  }

  InHDF5Producer::InHDF5Producer(ProcBlock *block)
  {
    SetMainBlock(block);
    Init();
  }

  InHDF5Producer::~InHDF5Producer()
  {
  }

  void InHDF5Producer::Init()
  {
    // Build commands
    G4UIdirectory* DebugDir = new G4UIdirectory("/rat/inhdf5/");
    DebugDir->SetGuidance("Read Events from HDF5 file");

    readCmd = new G4UIcmdWithAString("/rat/inhdf5/read", this);
    readCmd->SetGuidance("name of input file");
    readCmd->SetParameterName("filename", false);  // required

    readDefaultCmd = new G4UIcommand("/rat/inhdf5/read_default", this);
    readDefaultCmd->SetGuidance("read from IO.default_input_filename");
  }

  G4String InHDF5Producer::GetCurrentValue(G4UIcommand * /*command*/)
  {
    Log::Die("invalid inhdf5 \"get\" command");
    return G4String("You never see this.");
  }


  void InHDF5Producer::SetNewValue(G4UIcommand * command, G4String newValue)
  {
    // readCmd
    if (command == readCmd || command == readDefaultCmd) {
      std::string filename;

      if (command == readDefaultCmd) {
        DBLinkPtr lIO = DB::Get()->GetLink("IO");
        filename = lIO->GetS("default_input_filename");
      } else {
        size_t size = newValue.size();

        // Trim extraneous quotation marks to avoid confusing people
        if (size >= 2 && newValue[(size_t)0] == '\"' && newValue[size-1] == '\"')
        filename = newValue.substr(1, size-2);
        else
        filename = newValue;
      }

      if (!mainBlock)
      Log::Die("inhdf5: No main block declared! (should never happen)");
      else if (!ReadEvents(filename))
      Log::Die("inhdf5: Could not read from " + filename);

    } else
    Log::Die("invalid inhdf5 \"set\" command");
  }

  bool InHDF5Producer::ReadEvents(G4String filename)
  {

    if (!H5::H5File::isHdf5(filename.c_str()))
    {
      // Invalid HDF5 file
      std::cout<<" Invalid HDF5 file "<<filename<<std::endl;
      return 0;
    }

    H5::H5File *h5file = new H5::H5File(filename, H5F_ACC_RDONLY);
    H5::Group *h5rootgrp = new H5::Group(h5file->openGroup("/"));
    //Loop over groups (each group should be a different channel)
    //and fill the waveforms in event order and the DAQHeader.
    mw_t waveforms;
    bool daqHeader_filled = false;
    DS::DAQHeader *daqHeader = new DS::DAQHeader();
    int ngroups = h5rootgrp->getNumObjs();
    for (int i = 0; i < ngroups; i++){
      H5std_string grpname = h5rootgrp->getObjnameByIdx(i);
      int chnumber = -9999;
      if ( std::regex_match (grpname, std::regex("(ch)(.*)")) ){
        std::size_t pos = grpname.find("ch") + 2;
        std::string schnumber = grpname.substr (pos);
        try{
          chnumber = std::stoi(schnumber);
        }
        catch (...){
          Log::Die("Channel number "+ schnumber +" is neither a positive integer nor zero");
        }
      }
      else{
        Log::Die("Group name "+ grpname +" different from \"ch*\"");
      }

      H5::Group *group = new H5::Group(h5file->openGroup("/"+grpname));
      //Get DAQ attributes from first channel
      if(!daqHeader_filled){
        uint32_t bits,voffset,ns_sample;
        H5::Attribute bits_attr = group->openAttribute("bits");
        bits_attr.read(H5::PredType::NATIVE_UINT32, &bits);
        H5::Attribute voffset_attr = group->openAttribute("offset");
        voffset_attr.read(H5::PredType::NATIVE_UINT32, &voffset);
        H5::Attribute ns_sample_attr = group->openAttribute("ns_sample");
        ns_sample_attr.read(H5::PredType::NATIVE_UINT32, &ns_sample);

        daqHeader->SetAttribute("DAQ_NAME","DIGITIZER_V16");
        daqHeader->SetAttribute("NBITS",(int)bits);
        daqHeader->SetAttribute("TIME_RES",(int)ns_sample);
        daqHeader->SetAttribute("V_OFFSET",(int)voffset);
        daqHeader->SetAttribute("V_HIGH",1000);
        daqHeader->SetAttribute("V_LOW",-1000);
        daqHeader->SetAttribute("RESISTANCE",50);
        daqHeader_filled = true;
      }
      //Now retrieve the samples
      info<<"Getting samples from group " + grpname + "\n";
      H5::DataSet *dataset;
      try{
        dataset = new H5::DataSet(group->openDataSet("samples"));
      }
      catch(...){
        Log::Die("Data set \"samples\" not found in group " + grpname);
      }

      //Get dimensions of the dataset
      H5::DataSpace dataspace = dataset->getSpace();
      int rank = dataspace.getSimpleExtentNdims();
      hsize_t *dims = new hsize_t[rank];
      dataspace.getSimpleExtentDims(dims);

      //Get traces
      const int nevents = dims[0], nsamples = dims[1];
      UShort_t *data = new UShort_t[nevents*nsamples];
      dataset->read(data,H5::PredType::NATIVE_UINT16,dataspace);
      for(int iev=0; iev<nevents; iev++){
        std::vector<UShort_t> waveform;
        for(int isample=0; isample<nsamples; isample++){
          waveform.push_back(data[iev*nsamples + isample]);
        }
        waveforms[chnumber].push_back( (std::vector<UShort_t>) waveform);
      }
    }

    //Loop over events and waveforms and fill the DS
    int nevents = waveforms.begin()->second.size();
    for(int iev=0; iev<nevents; iev++){
      DS::Root* ds = new DS::Root();
      RAT::DS::EV *ev = ds->AddNewEV();
      for(mw_t::iterator iwaveform = waveforms.begin(); iwaveform != waveforms.end(); iwaveform++){
        RAT::DS::PMT *pmt = ev->AddNewPMT();
        pmt->SetID(iwaveform->first);
        pmt->SetWaveform(iwaveform->second.at(iev));
        // info<<"Waveforms "<<pmt->GetWaveform().size()<<"\n";
        // for(int isample=0; isample<iwaveform->second.at(iev).size();isample++){
        //   info<<"   "<<isample<<" "<<pmt->GetWaveform()[isample]<<" "<<iwaveform->second.at(iev)[isample]<<"\n";
        // }
      }

      //Now set the DAQHeader
      DS::Run *run = new DS::Run();
      run->SetID(1);
      run->SetType(0x00001111);
      run->SetStartTime(1440638077);
      run->SetDAQHeader(daqHeader);

      DS::RunStore::AddNewRun(run);



      mainBlock->DSEvent(ds);
    }

    return true;
  }

} // namespace RAT
