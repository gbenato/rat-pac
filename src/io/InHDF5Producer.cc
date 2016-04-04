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
#include <fstream>
#include <cstdint> //64bits integer

typedef std::map< int, std::vector< std::vector<UShort_t> > > mw_t; //channel:[event,waveform]
typedef std::map< int, std::vector< std::vector<double> > > mwt_t; //channel:[event,waveformtime]

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
      std::string infilename;

      if (command == readDefaultCmd) {
        DBLinkPtr lIO = DB::Get()->GetLink("IO");
        infilename = lIO->GetS("default_input_filename");
      } else {
        size_t size = newValue.size();

        // Trim extraneous quotation marks to avoid confusing people
        if (size >= 2 && newValue[(size_t)0] == '\"' && newValue[size-1] == '\"')
        infilename = newValue.substr(1, size-2);
        else
        infilename = newValue;
      }

      if (!mainBlock)
      Log::Die("inhdf5: No main block declared! (should never happen)");
      else if (!ReadEvents(infilename))
      Log::Die("inhdf5: Could not read from " + infilename);

    } else
    Log::Die("invalid inhdf5 \"set\" command");
  }

  bool InHDF5Producer::ReadEvents(G4String inputfilename)
  {

    if (!H5::H5File::isHdf5(inputfilename.c_str()))
    {
      // Invalid HDF5 file
      std::cout<<" Invalid HDF5 file "<<inputfilename<<std::endl;
      return 0;
    }

    //Fill map digitizer channel - pmtID. FIXME: eventually this should go in a ratdb table
    bool DEBUG = false;
    std::map<int,int> FastChtoID;
    std::map<int,int> SlowChtoID;
    std::map<int,int> IDToDAQ;
    FastChtoID[0]=12;
    FastChtoID[1]=999;
    FastChtoID[2]=13;
    FastChtoID[3]=999;
    FastChtoID[4]=14;
    FastChtoID[5]=999;
    FastChtoID[6]=15;
    FastChtoID[7]=999;
    FastChtoID[8]=16;
    FastChtoID[9]=999;
    FastChtoID[10]=17;
    FastChtoID[11]=999;
    FastChtoID[12]=18;
    FastChtoID[13]=999;
    FastChtoID[14]=19;
    FastChtoID[15]=999;
    FastChtoID[16]=20;
    FastChtoID[17]=999;
    FastChtoID[18]=21;
    FastChtoID[19]=999;
    FastChtoID[20]=22;
    FastChtoID[21]=999;
    FastChtoID[22]=23;
    FastChtoID[23]=999;
    FastChtoID[24]=6;
    FastChtoID[25]=999;
    FastChtoID[26]=7;
    FastChtoID[27]=999;
    FastChtoID[28]=999;
    FastChtoID[29]=999;
    FastChtoID[30]=999;
    FastChtoID[31]=999;
    SlowChtoID[0]=24;
    SlowChtoID[1]=0;
    SlowChtoID[2]=1;
    SlowChtoID[3]=2;
    SlowChtoID[4]=3;
    SlowChtoID[5]=4;
    SlowChtoID[6]=5;
    SlowChtoID[7]=999;
    SlowChtoID[8]=8;
    SlowChtoID[9]=9;
    SlowChtoID[10]=10;
    SlowChtoID[11]=11;
    SlowChtoID[12]=999;
    SlowChtoID[13]=999;
    SlowChtoID[14]=999;
    SlowChtoID[15]=999;

    for(std::map<int,int>::iterator ich = SlowChtoID.begin(); ich != SlowChtoID.end(); ich++){
      IDToDAQ[ich->second] = 0;
    }
    for(std::map<int,int>::iterator ich = FastChtoID.begin(); ich != FastChtoID.end(); ich++){
      IDToDAQ[ich->second] = 1;
    }

    if(DEBUG) info<<"Opening FAST group..... \n";

    //Get calibration from ratdb
    // fLCalibV1742 = DB::Get()->GetLink("CALIB","V1742");
    // json::Value fTimeCalibV1742 = fLCalibV1742->GetJSON("2.5GHz");

    ifstream calibfile(Form("%s/data/CALIB.ratdb", getenv("RATROOT") ) ) ;
    json::Reader reader(calibfile);
    json::Value calib;
    reader.getValue(calib);
    //    json::Value fTimeCalibV1742 = calib["1GHz"];
    //    json::Value fTimeCalibV1742 = calib["2.5GHz"];
    json::Value fTimeCalibV1742 = calib["5GHz"];

    //Get waveforms from FAST group
    H5::H5File *h5file = new H5::H5File(inputfilename, H5F_ACC_RDONLY);
    H5::Group *h5fastgr;
    try{
      h5fastgr = new H5::Group(h5file->openGroup("/fast"));
    }
    catch(...){
      Log::Die("Group name fast or master do not exist");
    }
    //Loop over groups (each group should be a different channel)
    //and fill the waveforms in event order and the DAQHeader.
    uint64_t *data_times;
    mw_t waveforms;
    mwt_t waveformTimes;
    std::map< int, uint16_t* > start_cell; //For time calibrations (PMTID:[EV:CELL])
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
          grnumber = std::atoi(sgrnumber.c_str());
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
        if(DEBUG) std::cout<<"   |-> Channel name "<<i<<" "<<chname<<std::endl;
        int chnumber = -9999;
        if ( std::regex_match (chname, std::regex("(ch)(.*)")) ){
          std::size_t pos = chname.find("ch") + 2;
          std::string schnumber = chname.substr (pos);
          try{
            chnumber = std::atoi(schnumber.c_str());
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
          daqHeaderV1742->SetAttribute("TIME_DELAY",0.0);
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
        const int nevents = dims[0], nsamples = dims[1];

        //Get initial cell index for time calibration
        H5::DataSet sidataset = daqgroup->openDataSet("start_index");
        start_cell[FastChtoID[chnumber]] = new uint16_t[nevents];
        sidataset.read(start_cell[FastChtoID[chnumber]],H5::PredType::NATIVE_UINT16);

        if(DEBUG) info<<"Got start_index \n";

        //Get traces
        UShort_t *data = new UShort_t[nevents*nsamples];
        dataset->read(data,H5::PredType::NATIVE_UINT16,dataspace);
        for(int iev=0; iev<nevents; iev++){
          std::vector<UShort_t> waveform;
          //Get calibrated time
          std::vector<double> waveformTime = fTimeCalibV1742[grpname]["cell_delay"].toVector<double>();
          double t0 = waveformTime[start_cell[FastChtoID[chnumber]][iev]];
          std::vector<double> waveformTimeShifted(waveformTime.size(),0.);
          for(int isample=0; isample<nsamples; isample++){
            if(data[iev*nsamples + isample]>65000) data[iev*nsamples + isample] = 0; //correction
            waveform.push_back(data[iev*nsamples + isample]);
            waveformTimeShifted[isample] = waveformTime[isample] - t0;
            if(waveformTimeShifted[isample]<0) waveformTimeShifted[isample] += waveformTime[nsamples-1];
          }
          waveforms[FastChtoID[chnumber]].push_back( (std::vector<UShort_t>) waveform);
          std::rotate(waveformTimeShifted.begin(),waveformTimeShifted.begin() + start_cell[FastChtoID[chnumber]][iev],waveformTimeShifted.end());
          waveformTimes[FastChtoID[chnumber]].push_back( waveformTimeShifted );
          std::vector<UShort_t>().swap(waveform);
          std::vector<double>().swap(waveformTime);
          std::vector<double>().swap(waveformTimeShifted);
        }

        delete channel;
        delete dataset;
        delete[] dims;
        delete[] data;
      }//end channels loop

      delete daqgroup;
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
          chnumber = std::atoi(schnumber.c_str());
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
        daqHeaderV1730->SetAttribute("TIME_DELAY",0.0);
        daqHeaderV1730->SetAttribute("V_OFFSET",(double)voffset);
        daqHeaderV1730->SetAttribute("V_HIGH",1000.);
        daqHeaderV1730->SetAttribute("V_LOW",-1000.);
        daqHeaderV1730->SetAttribute("RESISTANCE",50.);
        daqHeaderV1730_filled = true;

      }

      //Now retrieve the samples
      H5::DataSet *dataset;
      H5::DataSet *dataset_times;
      try{
        dataset = new H5::DataSet(channel->openDataSet("samples"));
        dataset_times = new H5::DataSet(channel->openDataSet("times"));
      }
      catch(...){
        Log::Die("Data set \"samples\" not found in channel " + chname);
      }

      //Get traces
      H5::DataSpace dataspace = dataset->getSpace();
      int rank = dataspace.getSimpleExtentNdims();
      hsize_t *dims = new hsize_t[rank];
      dataspace.getSimpleExtentDims(dims);
      const int nevents = dims[0], nsamples = dims[1];
      UShort_t *data = new UShort_t[nevents*nsamples];
      dataset->read(data,H5::PredType::NATIVE_UINT16,dataspace);

      H5::DataSpace dataspace_times = dataset_times->getSpace();
      int rank_times = dataspace_times.getSimpleExtentNdims();
      hsize_t *dims_times = new hsize_t[rank_times];
      dataspace_times.getSimpleExtentDims(dims_times);
      const int nevents_times = dims_times[0];
      data_times = new uint64_t[nevents_times];
      dataset_times->read(data_times,H5::PredType::NATIVE_UINT64,dataspace_times);

      for(int iev=0; iev<nevents; iev++){
        std::vector<UShort_t> waveform;
        for(int isample=0; isample<nsamples; isample++){
          waveform.push_back(data[iev*nsamples + isample]);
        }
        waveforms[SlowChtoID[chnumber]].push_back( (std::vector<UShort_t>) waveform);
        //Set sample times FIXME: need calibrated times
        std::vector<double> waveformTimesV1730;
        for(int isample=0; isample<nsamples; isample++){
          waveformTimesV1730.push_back( isample*daqHeaderV1730->GetDoubleAttribute("TIME_RES") );
        }
        waveformTimes[SlowChtoID[chnumber]].push_back( waveformTimesV1730 );
        std::vector<UShort_t>().swap(waveform);
        std::vector<double>().swap(waveformTimesV1730);
      }

      delete channel;
      delete dataset;
      delete[] dims;
      delete[] data;
    }//end channels loop

    delete h5mastergr;

    //Get event mapping file and opened file index
    std::stringstream sopened_file, sevent_map;

    std::size_t pos0 = inputfilename.find_last_of("/");
    std::size_t pos1 = inputfilename.find(".",pos0)+1;
    std::size_t pos2 = inputfilename.find(".",pos1+1);
    std::string filename = inputfilename.substr(0,pos1-1);

    std::string runnumber = inputfilename.substr(pos1,pos2-pos1);
    int opened_file = std::atoi(runnumber.c_str());
    std::string filename_map = filename + ".map.csv";

    TTree *event_map = new TTree("event_map","event_map");
    event_map->ReadFile(filename_map.c_str(),"event_id:master_file:master_index:fast_file:fast_index");
    float event_id, master_file, master_index, fast_file, fast_index;
    event_map->SetBranchAddress("event_id", &event_id);
    event_map->SetBranchAddress("master_file", &master_file);
    event_map->SetBranchAddress("master_index", &master_index);
    event_map->SetBranchAddress("fast_file", &fast_file);
    event_map->SetBranchAddress("fast_index", &fast_index);
    int nevents = event_map->GetEntries();

    //Loop over events and waveforms and fill the DS
    int skippedEvents = 0;
    for(int iev=0; iev<nevents; iev++){

      event_map->GetEntry(iev);
      if(fast_file != opened_file) continue;
      if(master_file != fast_file) {
        info<<" V1730and V1742 events are in different files. Skipping event "<<event_id<<"! \n";
        skippedEvents++;
        continue;
      }
      if(DEBUG) info<<" Ordering events: "<<iev<<" "<<event_id<<" "<<master_file<<" "<<master_index<<" "<<fast_file<<" "<<fast_index<<" \n";

      DS::Root* ds = new DS::Root();
      ds->SetRunID(1);
      RAT::DS::EV *ev = ds->AddNewEV();
      ev->SetID((int)event_id);
      if(master_index==0){
        ev->SetDeltaT(-9999.);
      }
      else{
        ev->SetDeltaT( (float) (data_times[(int)master_index] - data_times[(int)master_index-1]) );
      }

      if(DEBUG) std::cout<<" DeltaT "<<master_index<<" "<<data_times[(int)master_index]<<" "<<data_times[(int)master_index-1]<<" "<<ev->GetDeltaT()<<std::endl;

      for(mw_t::iterator iwaveform = waveforms.begin(); iwaveform != waveforms.end(); iwaveform++){
        if(iwaveform->first == 999) continue;
        // if(DEBUG) std::cout<<" Events for CH"<<iwaveform->first<<" "<<iwaveform->second.size()<<std::endl;
        if(DEBUG) info<<" WF "<<iwaveform->first<<"\n";

        if(IDToDAQ[iwaveform->first] == 0) {
          try{
            iwaveform->second.at((int)master_index);
          } catch (...){
            if(DEBUG) info<<" V1730 event: "<<event_id<<" "<<master_index<<" "<<" not in file "<<master_file<<"\n";
            continue;
          }
          if(DEBUG) info<<" Passed check V1730 \n";
          RAT::DS::PMT *pmt = ev->AddNewPMT();
          pmt->SetID(iwaveform->first);
          if(DEBUG) info<<" Passed set ID \n";
          pmt->SetWaveform(iwaveform->second.at((int)master_index));
          if(DEBUG) info<<" Passed Set wf \n";
          pmt->SetWaveformTime(waveformTimes[iwaveform->first].at((int)master_index));
          if(DEBUG) info<<" Passed Set wft \n";
        }
        if(IDToDAQ[iwaveform->first] == 1) {
          try{
            iwaveform->second.at((int)fast_index);
          } catch (...){
            if(DEBUG) info<<" V1742 event: "<<event_id<<" "<<fast_index<<" "<<" not in file "<<fast_file<<"\n";
          }
          if(DEBUG) info<<" Passed check V1742 \n";
          RAT::DS::PMT *pmt = ev->AddNewPMT();
          pmt->SetID(iwaveform->first);
          if(DEBUG) info<<" Passed set ID \n";
          pmt->SetWaveform(iwaveform->second.at((int)fast_index));
          if(DEBUG) info<<" Passed Set wf \n";
          pmt->SetWaveformTime(waveformTimes[iwaveform->first].at((int)fast_index));
          if(DEBUG) info<<" Passed Set wft \n";
        }
        // info<<"Waveforms "<<pmt->GetWaveform().size()<<"\n";
        // for(int isample=0; isample<iwaveform->second.at(iev).size();isample++){
        //   info<<"   "<<isample<<" "<<pmt->GetWaveform()[isample]<<" "<<iwaveform->second.at(iev)[isample]<<"\n";
        // }
        // info<<"Waveform Times "<<pmt->GetWaveformTime().size()<<"\n";
        // for(int isample=0; isample<waveformTimes[iwaveform->first].at(iev).size();isample++){
        //   info<<"  "<<"WF times: PMT "<<iwaveform->first<<" "<<isample<<" "<<pmt->GetWaveformTime()[isample]<<" "<<waveformTimes[iwaveform->first].at(iev)[isample]<<"\n";
        // }
      }

      //Now set the DAQHeader
      DS::Run *run = new DS::Run();
      run->SetID(1);
      run->SetType(0x00001111);
      run->SetStartTime(1440638077);
      run->SetDAQHeader(daqHeaderV1730,"V1730");
      run->SetDAQHeader(daqHeaderV1742,"V1742");

      const RAT::DS::PMTInfo *pmtInfo = &PMTFactoryBase::GetPMTInfo();
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

    std::cout<<" Skipped events: "<<skippedEvents<<std::endl;


    //Delete
    delete event_map;
    delete h5file;
    delete h5fastgr;
    mw_t().swap(waveforms);
    mwt_t().swap(waveformTimes);
    std::map< int, uint16_t* >().swap(start_cell);
    // delete daqHeaderV1742;
    // delete daqHeaderV1730;

    return true;

  } //end ReadEvents

} // namespace RAT
