#include <H5Cpp.h>
#include <RAT/InHDF5Producer.hh>
#include <RAT/ProcBlock.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/SignalHandler.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/DS/DAQHeader.hh>
#include <RAT/PMTFactoryBase.hh>

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

    bool DEBUG = false;
    std::map<int,int> FastChtoID;
    std::map<int,int> SlowChtoID;
    FastChtoID[0]=6;
    FastChtoID[1]=999;
    FastChtoID[2]=7;
    FastChtoID[3]=999;
    FastChtoID[4]=8;
    FastChtoID[5]=999;
    FastChtoID[6]=9;
    FastChtoID[7]=999;
    FastChtoID[8]=10;
    FastChtoID[9]=999;
    FastChtoID[10]=11;
    FastChtoID[11]=999;
    FastChtoID[12]=12;
    FastChtoID[13]=999;
    FastChtoID[14]=13;
    FastChtoID[15]=999;
    FastChtoID[16]=14;
    FastChtoID[17]=999;
    FastChtoID[18]=15;
    FastChtoID[19]=999;
    FastChtoID[20]=16;
    FastChtoID[21]=999;
    FastChtoID[22]=17;
    FastChtoID[23]=999;
    FastChtoID[24]=999;
    FastChtoID[25]=999;
    FastChtoID[26]=999;
    FastChtoID[27]=999;
    FastChtoID[28]=999;
    FastChtoID[29]=999;
    FastChtoID[30]=999;
    FastChtoID[31]=999;
    SlowChtoID[0]=999;
    SlowChtoID[1]=0;
    SlowChtoID[2]=1;
    SlowChtoID[3]=2;
    SlowChtoID[4]=3;
    SlowChtoID[5]=4;
    SlowChtoID[6]=5;
    SlowChtoID[7]=999;
    SlowChtoID[8]=999;
    SlowChtoID[9]=999;
    SlowChtoID[10]=999;
    SlowChtoID[11]=999;
    SlowChtoID[12]=999;
    SlowChtoID[13]=999;
    SlowChtoID[14]=999;
    SlowChtoID[15]=999;

    if(DEBUG) info<<"Opening FAST group..... \n";

    //Get waveforms from FAST group
    H5::H5File *h5file = new H5::H5File(filename, H5F_ACC_RDONLY);
    H5::Group *h5fastgr;
    try{
      h5fastgr = new H5::Group(h5file->openGroup("/fast"));
    }
    catch(...){
      Log::Die("Group name fast or master do not exist");
    }
    //Loop over groups (each group should be a different channel)
    //and fill the waveforms in event order and the DAQHeader.
    mw_t waveforms;
    bool daqHeaderV1742_filled = false;
    DS::DAQHeader *daqHeaderV1742 = new DS::DAQHeader();
    int ngroups = h5fastgr->getNumObjs();
    //DAQ groups loop
    for (int i = 0; i < ngroups; i++){
      H5std_string grpname = h5fastgr->getObjnameByIdx(i);
      if(DEBUG) std::cout<<" |-> Group name "<<i<<" "<<grpname<<std::endl;
      int grnumber = -9999;
      if ( std::regex_match (grpname, std::regex("(gr)(.*)")) ){
        std::size_t pos = grpname.find("gr") + 2;
        std::string sgrnumber = grpname.substr(pos);
        try{
          grnumber = std::stoi(sgrnumber);
        }
        catch (...){
          Log::Die("Group number " + sgrnumber + "is neither a positive integer nor zero");
        }
      }
      else{
        if(DEBUG) std::cout<<"Group name "+ grpname +" different from \"gr*\""<<std::endl;
      }

      H5::Group *daqgroup = new H5::Group(h5file->openGroup("/fast/"+grpname));
      int nch = daqgroup->getNumObjs();
      //Channels groups loop
      for (int i = 0; i < nch; i++){
        H5std_string chname = daqgroup->getObjnameByIdx(i);
        if(DEBUG) std::cout<<"   |-> Group name "<<i<<" "<<chname<<std::endl;
        int chnumber = -9999;
        if ( std::regex_match (chname, std::regex("(ch)(.*)")) ){
          std::size_t pos = chname.find("ch") + 2;
          std::string schnumber = chname.substr (pos);
          try{
            chnumber = std::stoi(schnumber);
            chnumber = chnumber + grnumber*8; //8 channels per group
          }
          catch (...){
            Log::Die("Channel number " + schnumber + " is neither a positive integer nor zero");
          }
        }
        else{
          if(DEBUG) std::cout<<"    Channel name "+ chname +" different from \"ch*\""<<std::endl;
          continue;
        }

        H5::Group *channel = new H5::Group(h5file->openGroup("/fast/"+grpname+"/"+chname));
        //Get DAQ attributes from first channel
        if(!daqHeaderV1742_filled){
          uint32_t bits;
          double voffset,ns_sample;
          H5::Attribute bits_attr = h5fastgr->openAttribute("bits");
          bits_attr.read(H5::PredType::NATIVE_UINT32, &bits);
          H5::Attribute voffset_attr = channel->openAttribute("offset");
          voffset_attr.read(H5::PredType::NATIVE_DOUBLE, &voffset);
          H5::Attribute ns_sample_attr = h5fastgr->openAttribute("ns_sample");
          ns_sample_attr.read(H5::PredType::NATIVE_DOUBLE, &ns_sample);

          daqHeaderV1742->SetAttribute("DAQ_NAME","V1742");
          daqHeaderV1742->SetAttribute("NBITS",(int)bits);
          daqHeaderV1742->SetAttribute("TIME_RES",(double)ns_sample);
          daqHeaderV1742->SetAttribute("V_OFFSET",(double)voffset);
          daqHeaderV1742->SetAttribute("V_HIGH",1000.);
          daqHeaderV1742->SetAttribute("V_LOW",-1000.);
          daqHeaderV1742->SetAttribute("RESISTANCE",50.);
          daqHeaderV1742_filled = true;
        }
        //Now retrieve the samples
        H5::DataSet *dataset;
        try{
          dataset = new H5::DataSet(channel->openDataSet("samples"));
        }
        catch(...){
          Log::Die("Data set \"samples\" not found in channel " + chname);
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
          waveforms[FastChtoID[chnumber]].push_back( (std::vector<UShort_t>) waveform);
        }
      }//end channels loop
    }//end groups loop


    if(DEBUG) info<<"Opening MASTER group..... \n";

    //Get waveforms from MASTER group
    H5::Group *h5mastergr;
    try{
      h5mastergr = new H5::Group(h5file->openGroup("/master"));
    }
    catch(...){
      Log::Die("Group name fast or master do not exist");
    }
    //Loop over groups (each group should be a different channel)
    //and fill the waveforms in event order and the DAQHeader.
    bool daqHeaderV1730_filled = false;
    DS::DAQHeader *daqHeaderV1730 = new DS::DAQHeader();
    int nch = h5mastergr->getNumObjs();
    //Channels groups loop
    for (int i = 0; i < nch; i++){
      H5std_string chname = h5mastergr->getObjnameByIdx(i);
      if(DEBUG) std::cout<<" |-> Group name "<<i<<" "<<chname<<std::endl;
      int chnumber = -9999;
      if ( std::regex_match (chname, std::regex("(ch)(.*)")) ){
        std::size_t pos = chname.find("ch") + 2;
        std::string schnumber = chname.substr (pos);
        try{
          chnumber = std::stoi(schnumber);
        }
        catch (...){
          Log::Die("Channel number "+ schnumber +" is neither a positive integer nor zero");
        }
      }
      else{
        if(DEBUG) std::cout<<"   Group name "+ chname +" different from \"ch*\""<<std::endl;
        continue;
      }

      H5::Group *channel = new H5::Group(h5file->openGroup("/master/"+chname));
      //Get DAQ attributes from first channel
      if(!daqHeaderV1730_filled){
        uint32_t bits;
        double voffset,ns_sample;
        H5::Attribute bits_attr = h5mastergr->openAttribute("bits");
        bits_attr.read(H5::PredType::NATIVE_UINT32, &bits);
        H5::Attribute voffset_attr = channel->openAttribute("offset");
        voffset_attr.read(H5::PredType::NATIVE_DOUBLE, &voffset);
        H5::Attribute ns_sample_attr = h5mastergr->openAttribute("ns_sample");
        ns_sample_attr.read(H5::PredType::NATIVE_DOUBLE, &ns_sample);

        daqHeaderV1730->SetAttribute("DAQ_NAME","V1730");
        daqHeaderV1730->SetAttribute("NBITS",(int)bits);
        daqHeaderV1730->SetAttribute("TIME_RES",(double)ns_sample);
        daqHeaderV1730->SetAttribute("V_OFFSET",(double)voffset);
        daqHeaderV1730->SetAttribute("V_HIGH",1000.);
        daqHeaderV1730->SetAttribute("V_LOW",-1000.);
        daqHeaderV1730->SetAttribute("RESISTANCE",50.);
        daqHeaderV1730_filled = true;
      }
      //Now retrieve the samples
      H5::DataSet *dataset;
      try{
        dataset = new H5::DataSet(channel->openDataSet("samples"));
      }
      catch(...){
        Log::Die("Data set \"samples\" not found in channel " + chname);
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
        waveforms[SlowChtoID[chnumber]].push_back( (std::vector<UShort_t>) waveform);
      }
    }//end channels loop

    int nevents = waveforms.begin()->second.size(); //FIXME: deal with different number of events...

    //Loop over events and waveforms and fill the DS
    for(int iev=0; iev<nevents; iev++){

      DS::Root* ds = new DS::Root();
      ds->SetRunID(1);
      RAT::DS::EV *ev = ds->AddNewEV();
      for(mw_t::iterator iwaveform = waveforms.begin(); iwaveform != waveforms.end(); iwaveform++){
        if(DEBUG) std::cout<<" Events for CH"<<iwaveform->first<<" "<<iwaveform->second.size()<<std::endl;
        if(iwaveform->first == 999) continue;
        if(iwaveform->second.size()-1 < iev) continue; //FIXME: deal with different number of events...
        RAT::DS::PMT *pmt = ev->AddNewPMT();
        pmt->SetID(iwaveform->first);
        if(pmt->GetID() >= 6) pmt->SetType(1);
        else if(pmt->GetID() < 6) pmt->SetType(2);
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
      run->SetDAQHeader(daqHeaderV1730,"V1730");
      run->SetDAQHeader(daqHeaderV1742,"V1742");
      run->SetPMTInfo(&PMTFactoryBase::GetPMTInfo());
      DS::RunStore::AddNewRun(run);

      DS::Run *run2 = DS::RunStore::GetRun(ds);

      if(run == NULL) {
        std::cout<<" InHD5F: Run not Found "<<std::endl;
        exit(0);
      }

      if(run2 == NULL) {
        std::cout<<" InHD5F: Run2 not Found "<<std::endl;
        exit(0);
      }

      mainBlock->DSEvent(ds);

    } //end event loop

    return true;

  } //end ReadEvents

} // namespace RAT
