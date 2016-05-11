#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"

#include "TaggerTest/include/TopVar.h"
#include "TaggerTest/include/Utility.h"

#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"



// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./TopVar TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1]; 
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);  
  const int maxevent = std::atoi(Maxevent);
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";  
  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }

  AnaFunctions::prepareForNtupleReader();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain);
  //configure top tagger
  TopTagger tt;
  //  tt.setCfgFile("TaggerTest/test/file.cfg");
  tt.setCfgFile("file.cfg");
  bool isSignal = false;
  TString sample(subsamplename);
  if(sample.Contains("Signal")) isSignal = true;
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr->getNextEvent()){
    k++;    
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr->getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
    const vector<int> &genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");

    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    double EvtWt = tr->getVar<double>("evtWeight");

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    EventWeight = EvtWt;
    Lumiscale = Lumiscale * EventWeight;
    int Njets = Utility::Njet(jetsLVec, "jetPt30Eta24");
    int Nbjets =  Utility::Nbjet(jetsLVec, recoJetsBtag_0);
    double HT = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt30Eta24Arr);
   
    bool passMET = true, passNJET = true, passBJET = true, passHT = true;
    if(metLVec.Pt() < AnaConsts::defaultMETcut) passMET = false;
    if(Utility::Njet(jetsLVec, "jetPt50Eta24") < AnaConsts::nJetsSelPt50Eta24) passNJET = false;
    if(Utility::Njet(jetsLVec, "jetPt30Eta24") < AnaConsts::nJetsSelPt30Eta24) passNJET = false;
    if(!( (AnaConsts::low_nJetsSelBtagged == -1 || Nbjets >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || Nbjets < AnaConsts::high_nJetsSelBtagged ) ) ) passBJET = false;
    if(HT < AnaConsts::defaultHTcut) passHT = false;
    

    //construct vector of constituents 
    vector<Constituent> constituents = ttUtility::packageCandidates(jetsLVec, recoJetsBtag_0);
    //run tagger
    tt.runTagger(constituents);
    //get output of tagger
    const TopTaggerResults& ttr = tt.getResults();

    //Use result for top var
    vector<TopObject*> Ntop = ttr.getTops();
    vector<TopObject> Ntopcand = ttr.getTopCandidates();    
    ntops_ = Ntop.size();
    isSignal==true? Utility::FillInt(myBaseHistgram.hNtop_Sig,ntops_,Lumiscale):Utility::FillDouble(myBaseHistgram.hNtop_Bck,ntops_,Lumiscale);
    cout<<"No. of tops: "<<ntops_<<endl;    
    cout<<"No. of top candidates: "<<Ntopcand.size()<<endl; 
   for(const TopObject* top : ttr.getTops())
      {
	double topPt = top->p().Pt();
	double topMass = top->p().M();
	double topdRMax = top->getDRmax();
	isSignal==true? Utility::FillDouble(myBaseHistgram.htopPt_Sig,topPt,Lumiscale):Utility::FillDouble(myBaseHistgram.htopPt_Bck,topPt,Lumiscale);     
	isSignal==true? Utility::FillDouble(myBaseHistgram.htopMass_Sig,topMass,Lumiscale):Utility::FillDouble(myBaseHistgram.htopMass_Bck,topMass,Lumiscale);     
	isSignal==true? Utility::FillDouble(myBaseHistgram.hdRMax_Sig,topdRMax,Lumiscale):Utility::FillDouble(myBaseHistgram.hdRMax_Bck,topdRMax,Lumiscale);     

	cout<<"Top Pt: "<<topPt<<endl;
	cout<<"Top mass: "<<top->p().M()<<endl;
	cout<<"DeltaRmax: "<<top->getDRmax()<<endl;
	cout<<"Top const: "<<top->getNConstituents()<<endl;

	for(const Constituent *topconst : top->getConstituents()) 
	  {
	    cout<<"topconst pt:"<<topconst->p().Pt()<<endl;
	    cout<<"topconst mass:"<<topconst->p().M()<<endl;
	    cout<<"topconst deltaR: "<<ROOT::Math::VectorUtil::DeltaR(top->p(), topconst->p())<<endl;
	  }
      }

     


  } // End of loop over tree entries
  
  
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  return 0;
}

