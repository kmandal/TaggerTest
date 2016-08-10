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

#include "TaggerTest/include/TopCatagory.h"
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
      std::cerr <<" ./TopCatagory TTbarInc 1000 0 1" << std::endl;
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
  tt.setCfgFile("file.cfg");
  TString sample(subsamplename);
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;

  // Loop over the events (tree entries)
  while(tr->getNextEvent()){
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
    bool topmatch = false;
    bool topmatchCand = false;

    int Njets = Utility::Njet(jetsLVec, "Pt30Eta24");
    int Nbjets =  Utility::Nbjet(jetsLVec, recoJetsBtag_0);
    double HT = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt30Eta24Arr);
   
    bool passMET = true, passNJET = true, passBJET = true, passHT = true;
    if(metLVec.Pt() < AnaConsts::defaultMETcut) passMET = false;
    if(Utility::Njet(jetsLVec, "Pt50Eta24") < AnaConsts::nJetsSelPt50Eta24) passNJET = false;
    if(Utility::Njet(jetsLVec, "Pt30Eta24") < AnaConsts::nJetsSelPt30Eta24) passNJET = false;
    if(!( (AnaConsts::low_nJetsSelBtagged == -1 || Nbjets >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || Nbjets < AnaConsts::high_nJetsSelBtagged ) ) ) passBJET = false;
    if(HT < AnaConsts::defaultHTcut) passHT = false;
    
    //cut    
    if(!(passMET && passNJET && passBJET))continue;

    //gen info
    vector<TLorentzVector> hadtopLVec = GetTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
    //    cout<<"hadtop: "<<hadtopLVec.size()<<endl;
   
    //construct vector of constituents 
    vector<Constituent> constituents = ttUtility::packageCandidates(jetsLVec, recoJetsBtag_0);
    //run tagger
    tt.runTagger(constituents);
    //get output of tagger
    const TopTaggerResults& ttr = tt.getResults();

    //Use result for top var
    vector<TopObject*> Ntop = ttr.getTops();
    vector<TopObject> NtopCand = ttr.getTopCandidates();    

    ntops_ = Ntop.size();

    //matching
    vector<TopObject*> MatchNtop;
    vector<TopObject> MatchNtopCand;
    vector<TLorentzVector> MatchGentop; 
    vector<TLorentzVector> MatchGentopWidCand; 
    if(GetMatchedTop(Ntop, MatchNtop, hadtopLVec, MatchGentop)) topmatch = true;//final top match
    if(GetMatchedTopCand(NtopCand, MatchNtopCand, hadtopLVec, MatchGentopWidCand)) topmatchCand = true;//topcand match

    matchntops_ = MatchNtop.size();
    matchntopCands_ = MatchNtopCand.size();
 
    Utility::FillInt(myBaseHistgram.hNtop,ntops_,Lumiscale);
    Utility::FillInt(myBaseHistgram.hMatchedNtop,matchntops_,Lumiscale);
    Utility::FillInt(myBaseHistgram.hMatchedNtopCand,matchntopCands_,Lumiscale);

    vector<bool> monojetmatchCand(MatchNtopCand.size(),false);
    vector<bool> dijet0matchCand(MatchNtopCand.size(),false);
    vector<bool> dijet1matchCand(MatchNtopCand.size(),false);
    vector<bool> dijet2matchCand(MatchNtopCand.size(),false);
    vector<bool> trijet0matchCand(MatchNtopCand.size(),false);
    vector<bool> trijet1matchCand(MatchNtopCand.size(),false);
    vector<bool> trijet2matchCand(MatchNtopCand.size(),false);
    vector<bool> trijet3matchCand(MatchNtopCand.size(),false);
    
    //constituent matching
    if(matchntopCands_){
      for(unsigned tc = 0; tc<MatchNtopCand.size(); tc++){
	vector<Constituent const*> topconst1 = MatchNtopCand[tc].getConstituents();
	vector<TLorentzVector>gentopdauLVec = GettopdauLVec(MatchGentopWidCand[tc], genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
	if(topconst1.size()==1) monojetmatchCand[tc]=true;
	if(topconst1.size()==2){
	  int dimatch = GetMatchedTopConst(topconst1, gentopdauLVec);
	  if(dimatch==0)dijet0matchCand[tc] = true;
	  if(dimatch==1)dijet1matchCand[tc] = true;
	  if(dimatch==2)dijet2matchCand[tc] = true;
	}
	if(topconst1.size()==3){
	  int trimatch = GetMatchedTopConst(topconst1, gentopdauLVec);
	  if(trimatch==0)trijet0matchCand[tc] = true;
	  if(trimatch==1)trijet1matchCand[tc] = true;
	  if(trimatch==2)trijet2matchCand[tc] = true;
	  if(trimatch==3)trijet3matchCand[tc] = true;
	}
      }
    }
    
    vector<bool> monojetmatch(MatchNtop.size(),false);
    vector<bool> dijet0match(MatchNtop.size(),false);
    vector<bool> dijet1match(MatchNtop.size(),false);
    vector<bool> dijet2match(MatchNtop.size(),false);
    vector<bool> trijet0match(MatchNtop.size(),false);
    vector<bool> trijet1match(MatchNtop.size(),false);
    vector<bool> trijet2match(MatchNtop.size(),false);
    vector<bool> trijet3match(MatchNtop.size(),false);

    //constituent matching
    if(matchntops_){
      for(unsigned tc = 0; tc<MatchNtop.size(); tc++){
	vector<Constituent const*> topconst1 = MatchNtop[tc]->getConstituents();
	vector<TLorentzVector>gentopdauLVec = GettopdauLVec(MatchGentop[tc], genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
	if(topconst1.size()==1) monojetmatch[tc]=true;
	if(topconst1.size()==2){
	  int dimatch = GetMatchedTopConst(topconst1, gentopdauLVec);
	  if(dimatch==0)dijet0match[tc] = true;
	  if(dimatch==1)dijet1match[tc] = true;
	  if(dimatch==2)dijet2match[tc] = true;
	}
	if(topconst1.size()==3){
	  int trimatch = GetMatchedTopConst(topconst1, gentopdauLVec);
	  if(trimatch==0)trijet0match[tc] = true;
	  if(trimatch==1)trijet1match[tc] = true;
	  if(trimatch==2)trijet2match[tc] = true;
	  if(trimatch==3)trijet3match[tc] = true;
	}
      }
    }

    if(topmatch){
      for(unsigned a =0; a<MatchNtop.size(); a++){
	Utility::FillDouble(myBaseHistgram.htopPt_matchtop,MatchNtop[a]->p().Pt(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.htopMass_matchtop,MatchNtop[a]->p().M(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.hdRMax_matchtop,MatchNtop[a]->getDRmax(),Lumiscale);
	vector<Constituent const*>con = MatchNtop[a]->getConstituents();
	if(con.size()==3 && trijet3match[a]){
	  Utility::Fill2D(myBaseHistgram.hPt_dRMax_matchtop,MatchNtop[a]->p().Pt(), MatchNtop[a]->getDRmax(),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_dRMin_matchtop, MatchNtop[a]->getDRmax(), GetdRmin_top(MatchNtop[a]),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_area_matchtop, MatchNtop[a]->getDRmax(), GetArea_top(MatchNtop[a]),Lumiscale);
	}
      }
    }
    if(!topmatch){
      for(unsigned b =0; b<Ntop.size(); b++){
	Utility::FillDouble(myBaseHistgram.htopPt_nomatchtop,Ntop[b]->p().Pt(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.htopMass_nomatchtop,Ntop[b]->p().M(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.hdRMax_nomatchtop,Ntop[b]->getDRmax(),Lumiscale);
	Utility::Fill2D(myBaseHistgram.hPt_dRMax_nomatchtop,Ntop[b]->p().Pt(), Ntop[b]->getDRmax(),Lumiscale);
      	vector<Constituent const*>con = Ntop[b]->getConstituents();
	if(con.size()==3){
	  Utility::Fill2D(myBaseHistgram.hdRMax_dRMin_nomatchtop, Ntop[b]->getDRmax(), GetdRmin_top(Ntop[b]),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_area_nomatchtop, Ntop[b]->getDRmax(), GetArea_top(Ntop[b]),Lumiscale);
	}
      }
    }
    if(topmatchCand){
      for(unsigned c =0; c<MatchNtopCand.size(); c++){
	Utility::FillDouble(myBaseHistgram.htopPt_matchtopCand,MatchNtopCand[c].p().Pt(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.htopMass_matchtopCand,MatchNtopCand[c].p().M(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.hdRMax_matchtopCand,MatchNtopCand[c].getDRmax(),Lumiscale);
	vector<Constituent const*>con = MatchNtopCand[c].getConstituents();
	if(con.size()==3 && trijet3matchCand[c]){
	  Utility::Fill2D(myBaseHistgram.hPt_dRMax_matchtopCand,MatchNtopCand[c].p().Pt(), MatchNtopCand[c].getDRmax(),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_dRMin_matchtopCand, MatchNtopCand[c].getDRmax(), GetdRmin_topcand(MatchNtopCand[c]),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_area_matchtopCand, MatchNtopCand[c].getDRmax(), GetArea_topcand(MatchNtopCand[c]),Lumiscale);
	  vector<TLorentzVector>MatchgentopdauLVec = GettopdauLVec(MatchGentopWidCand[c], genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
	  double MatchgentopdRMax = GetdRmax(MatchGentopWidCand[c],MatchgentopdauLVec);
	  Utility::Fill2D(myBaseHistgram.hPt_dRMax_matchtopGen,MatchGentopWidCand[c].Pt(), MatchgentopdRMax,Lumiscale);
	}
      }
    }
    if(!topmatchCand){
      for(unsigned d =0; d<NtopCand.size(); d++){
	if(!(NtopCand[d].p().M()>minTopCandMass && NtopCand[d].p().M()<maxTopCandMass)) continue;//mass window
	vector<Constituent const*>con = NtopCand[d].getConstituents();
	Utility::FillDouble(myBaseHistgram.htopPt_nomatchtopCand,NtopCand[d].p().Pt(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.htopMass_nomatchtopCand,NtopCand[d].p().M(),Lumiscale);
	Utility::FillDouble(myBaseHistgram.hdRMax_nomatchtopCand,NtopCand[d].getDRmax(),Lumiscale);
	if(con.size()==3){
	  Utility::Fill2D(myBaseHistgram.hPt_dRMax_nomatchtopCand,NtopCand[d].p().Pt(), NtopCand[d].getDRmax(),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_dRMin_nomatchtopCand, NtopCand[d].getDRmax(), GetdRmin_topcand(NtopCand[d]),Lumiscale);
	  Utility::Fill2D(myBaseHistgram.hdRMax_area_nomatchtopCand, NtopCand[d].getDRmax(), GetArea_topcand(NtopCand[d]),Lumiscale);
	}
      }
    }

    //Efficency
    for(unsigned ie = 0; ie<hadtopLVec.size(); ie++){
      Utility::FillDouble(myBaseHistgram.hgentopPt_den,hadtopLVec[ie].Pt(),Lumiscale);
      bool fMatch = false;
      for(unsigned je = 0; je<MatchGentop.size(); je++){
	if(hadtopLVec[ie]!=MatchGentop[je])continue;
	//if(trijet3match[je] || dijet2match[je] || monojetmatch[je]){
	  fMatch = true; break;
	  //	}
      }
      if(fMatch)Utility::FillDouble(myBaseHistgram.hgentopPt_num,hadtopLVec[ie].Pt(),Lumiscale);
    }

    //FakeRate
    if(sample.Contains("ZJetsToNuNu")){
      Utility::FillDouble(myBaseHistgram.hfakeNjet_den, Njets, Lumiscale);
      Utility::FillDouble(myBaseHistgram.hfakeHT_den, HT, Lumiscale);
      Utility::FillDouble(myBaseHistgram.hfakeMET_den, met, Lumiscale);
      if(Ntop.size()){
	Utility::FillDouble(myBaseHistgram.hfakeNjet_num, Njets, Lumiscale);
	Utility::FillDouble(myBaseHistgram.hfakeHT_num, HT, Lumiscale);
	Utility::FillDouble(myBaseHistgram.hfakeMET_num, met, Lumiscale);
      }
    }
  } // End of loop over tree entries
  
  
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  return 0;
}

