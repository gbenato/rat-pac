//Usage with ROOT: .L DrawProcess.C;
// extractTree("../results/some_file.root", 2)
// extractEvents("../results/some_file.root")

class Particle {
    public:
        Particle();
        virtual ~Particle();
        
        void setProperties( string name, string gen_process, string gen_volume, string end_process, string end_volume, TVector3 *steps, int nSteps);
        void addChild(Particle *child);
        
        int numRemaining();
        int numChildren();
        
        Particle *getNext();
        Particle *getChildren();
        
        void dumpTree(ostream &out);
        void dumpList(ostream &out);
        
        string name;
        string gen_process;
        string end_process;
	string gen_volume;
	string end_volume;
        TVector3 *steps; //take ownership of this
        int nSteps;

    protected:
        void dumpTree(int depth, ostream &out);
        
        //don't take ownership of any of these pointers
        Particle *next;
        Particle *children;
        Particle *child_tail;
	Particle *parent;
};

Particle::Particle() : parent(NULL), name("UNSET"), gen_process("UNSET"), end_process("UNSET"), next(NULL), children(NULL), child_tail(NULL), steps(NULL), nSteps(0) { 
}

Particle::~Particle() { 
    if (steps) delete steps;
}

void Particle::setProperties(string name, string gen_process, string gen_volume, string end_process, string end_volume, TVector3 *steps, int nSteps) {
    this->name = name;
    this->gen_process = gen_process;
    this->gen_volume = gen_volume;
    this->end_process = end_process;
    this->end_volume = end_volume;
    this->steps = steps;
    this->nSteps = nSteps;
}
        
void Particle::addChild(Particle *child) {
    if (children) {
        child_tail->next = child;
    } else {
        children = child;
    }
    child_tail = child;
}

void Particle::addParent(Particle* aparent) {
    parent = aparent;
}

int Particle::numChildren() {
    if (children) return 1+children->numRemaining();
    return 0;
}

int Particle::numRemaining() {
    int remaining = 0;
    Particle *cur = this;
    while (cur->next) {
        remaining++;
        cur = cur->next;
    }
    return remaining;
}

Particle *Particle::getNext() {
    return next;
}

Particle *Particle::getChildren() {
    return children;
}

Particle *Particle::getParent(){
   return parent;
}

void Particle::dumpList(ostream &out) {
    /***
    Dump all the track information of an event. Starting from the current particle go down in particle depth (children of children) for each
    particle until that particular track is finished, then go to the next child and repeat the process.
    ***/
    out << fixed << setprecision(5);
    out << '[';
        if(parent){
          out << '\"' << parent->name << "\",";
        }
        out << '\"' << name << "\",";
        out << '\"' << gen_process << "\",";
        out << '\"' << gen_volume << "\",";
        out << '\"' << end_process << "\",";
        out << '\"' << end_volume << "\",";
        out << '(';
            for (int i = 0; i <= nSteps-1; i++) {
                out << '{' << steps[i].X() << ',' << steps[i].Y() << ',' << steps[i].Z() << "},";
            }
//            out << '{' << steps[nSteps-1].X() << ',' << steps[nSteps-1].Y() << ',' << steps[nSteps-1].Z() << "}";
        out << "),";
        out << '{';
            for (Particle *child = children; child; child = child->next) {
                child->dumpList(out); 
            }
        out << '}';
    out << ']';
    if (next) out << ",\n";
}

void Particle::dumpTree(ostream &out) {
    dumpTree(0,out);        
}

void Particle::dumpTree(int depth, ostream &out) {
    for (int i = 0; i < depth; i++) {
        out << '|';
    }
    out << name << ' '; 
    out << gen_process << ' ';
    out << end_process << ' ';
    out << '\n';
    for (Particle *child = children; child; child = child->next) {
        child->dumpTree(depth+1,out);
    }
}




void extractTree(const char *file, int event_num) {

  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");
    
    RAT::DS::Root *rds = new RAT::DS::Root();
    tree->SetBranchAddress("ds", &rds);
    
    tree->GetEntry(event_num);
    RAT::DS::MC *mc = rds->GetMC();

    int nTracks = mc->GetMCTrackCount();
    cout << "Reading event from tree with " << nTracks << " tracks" << endl;
    while(!nTracks && event_num < (tree->GetEntries()-1)){
      cout << "Nothing to display for event "<< event_num << endl;
      event_num +=1;
      tree->GetEntry(event_num)
      nTracks = mc->GetMCTrackCount();
    }
    cout << "Reading event "<< event_num << " with " << nTracks << " tracks" << endl;
    Particle *trackmap = new Particle[nTracks+1];
    for (int j = 0; j < nTracks; j++) {
        RAT::DS::MCTrack *track = mc->GetMCTrack(j);
        int tid = track->GetID();
        int pid = track->GetParentID();

        Particle *cur = &trackmap[tid];
	Particle *par = &trackmap[pid];
        trackmap[pid].addChild(cur);
        trackmap[tid].addParent(par);

        RAT::DS::MCTrackStep *first = track->GetMCTrackStep(0);
        RAT::DS::MCTrackStep *last = track->GetLastMCTrackStep();
        
        int nSteps = track->GetMCTrackStepCount();
        cout << "Got " << nSteps << " steps in MCTrack" << endl;
        TVector3 *steps = new TVector3[nSteps];
        for (int k = 0; k < nSteps; k++) {
            steps[k] = track->GetMCTrackStep(k)->GetEndpoint();
        }
//	cout << first->GetVolume() << endl;
//        cout << first->GetVolume().c_str() << endl;
//        string s = first->GetVolume();
//	cout << j << endl;
        
	cur->setProperties( track->GetParticleName(), first->GetProcess(), first->GetVolume() ,last->GetProcess(), last->GetVolume(), steps,nSteps);
    }
    // I'm not sure what the primary particle ID is. It seems that 0 is undefined and it is probably 1, but 0 still has children - How???? 
    // oh yes probably the initial value for parentid is 0 if there is no parent --> hence the electron is a child of the empty Particle...
//    trackmap[0].dumpList(cout);
    trackmap[1].dumpList(cout);
}




void extractEvents(const char *file) {
  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");
  
  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);
  
  int nEvents = tree->GetEntries();
  
  cout << '{';
  for (int i = 0; i < nEvents; i++) {
    tree->GetEntry(i);
    RAT::DS::MC *mc = rds->GetMC();
    RAT::DS::MCParticle *prim = mc->GetMCParticle(0);
    
    cout << '{';
    cout << mc->GetMCPMTCount() << ',';
    cout << prim->GetKE() << ',';
    cout << '{' << prim->GetPosition().X() << ',' << prim->GetPosition().Y() << ',' << prim->GetPosition().Z() << "},";
    cout << '{' << prim->GetMomentum().X() << ',' << prim->GetMomentum().Y() << ',' << prim->GetMomentum().Z() << "}";
    cout << '}';
    if (i != nEvents-1) cout << ',';   
  }
  cout << "}\n";
}
