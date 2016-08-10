#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

static const int nSB = 45;//search bin number

using namespace std;

double Lumiscale = 1.0;
double EventWeight = 1.0;

int ntops_ = 0;
int matchntops_ =0;
int matchntopCands_ =0;
int maxConst = 3;
double minTopCandMass = 100;
double maxTopCandMass = 250;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hNtop;
  TH1D *hMatchedNtop;
  TH1D *hMatchedNtopCand;
  TH1D *hdR1;
  TH1D *hdR2;
  TH1D *hdR3;
  TH1D *hdRMax_matchtop;
  TH1D *hdRMax_nomatchtop;
  TH1D *hdRMax_matchtopCand;
  TH1D *hdRMax_nomatchtopCand;
  TH1D *htopMass_matchtop;
  TH1D *htopMass_nomatchtop;
  TH1D *htopMass_matchtopCand;
  TH1D *htopMass_nomatchtopCand;
  TH1D *htopPt_matchtop;
  TH1D *htopPt_nomatchtop;
  TH1D *htopPt_matchtopCand;
  TH1D *htopPt_nomatchtopCand;
  TH2D *hPt_dRMax_matchtop;
  TH2D *hPt_dRMax_nomatchtop;
  TH2D *hPt_dRMax_matchtopCand;
  TH2D *hPt_dRMax_nomatchtopCand;
  TH2D *hPt_dRMax_matchtopGen;
  TH2D *hdRMax_dRMin_matchtop;
  TH2D *hdRMax_dRMin_nomatchtop;
  TH2D *hdRMax_dRMin_matchtopCand;
  TH2D *hdRMax_dRMin_nomatchtopCand;
  TH2D *hdRMax_area_matchtop;
  TH2D *hdRMax_area_nomatchtop;
  TH2D *hdRMax_area_matchtopCand;
  TH2D *hdRMax_area_nomatchtopCand;

  TH1D *hgentopPt_den;
  TH1D *hgentopPt_num;
  TH1D *hfakeNjet_den;
  TH1D *hfakeHT_den;
  TH1D *hfakeMET_den;
  TH1D *hfakeNjet_num;
  TH1D *hfakeHT_num;
  TH1D *hfakeMET_num;

  TH1D * hTot_fake;

  std::vector<TH1*> hconstPt;
  std::vector<TH1*> hconstM;

  TString Title(TString spec, unsigned int i);

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_TopCatagory"+index+".root";
  oFile = new TFile(filename, "recreate");

  hNtop = new TH1D("hNtop", "No. of top;N_{top};Event", 5, 0, 5);
  hNtop->Sumw2();
  hMatchedNtop = new TH1D("hMatchedNtop", "No.of top;N_{top};Event", 5, 0, 5);
  hMatchedNtop->Sumw2();
  hMatchedNtopCand = new TH1D("hMatchedNtopCand", "No.of top;N_{top};Event", 5, 0, 5);
  hMatchedNtopCand->Sumw2();
  hdR2 = new TH1D("hdR2", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR2->Sumw2(); 
  hdR3 = new TH1D("hdR3", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR3->Sumw2(); 
  hdR1 = new TH1D("hdR1", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR1->Sumw2(); 

  hdRMax_matchtop = new TH1D("hdRMax_matchtop", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_matchtop->Sumw2(); 
  hdRMax_nomatchtop = new TH1D("hdRMax_nomatchtop", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_nomatchtop->Sumw2(); 
  hdRMax_matchtopCand = new TH1D("hdRMax_matchtopCand", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_matchtopCand->Sumw2(); 
  hdRMax_nomatchtopCand = new TH1D("hdRMax_nomatchtopCand", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_nomatchtopCand->Sumw2(); 
  htopMass_matchtop = new TH1D("htopMass_matchtop","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_matchtop->Sumw2();
  htopMass_nomatchtop = new TH1D("htopMass_nomatchtop","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_nomatchtop->Sumw2();
  htopMass_matchtopCand = new TH1D("htopMass_matchtopCand","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_matchtopCand->Sumw2();
  htopMass_nomatchtopCand = new TH1D("htopMass_nomatchtopCand","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_nomatchtopCand->Sumw2();
  htopPt_matchtop = new TH1D("htopPt_matchtop","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_matchtop->Sumw2();
  htopPt_nomatchtop = new TH1D("htopPt_nomatchtop","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_nomatchtop->Sumw2();
  htopPt_matchtopCand = new TH1D("htopPt_matchtopCand","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_matchtopCand->Sumw2();
  htopPt_nomatchtopCand = new TH1D("htopPt_nomatchtopCand","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_nomatchtopCand->Sumw2();
  hPt_dRMax_matchtop = new TH2D("hPt_dRMax_matchtop","Top P_{T} vs Max #Delta R;p_{T}[GeV];Max #Delta R",25, 0, 1000, 15, 0, 1.5);
  hPt_dRMax_matchtop->Sumw2();
  hPt_dRMax_nomatchtop = new TH2D("hPt_dRMax_nomatchtop","Top P_{T} vs Max #Delta R;p_{T}[GeV];Max #Delta R",25, 0, 1000, 15, 0, 1.5);
  hPt_dRMax_nomatchtop->Sumw2();
  hPt_dRMax_matchtopCand = new TH2D("hPt_dRMax_matchtopCand","Top P_{T} vs Max #Delta R;p_{T}[GeV];Max #Delta R",25, 0, 1000, 15, 0, 1.5);
  hPt_dRMax_matchtopCand->Sumw2();
  hPt_dRMax_nomatchtopCand = new TH2D("hPt_dRMax_nomatchtopCand","Top P_{T} vs Max #Delta R;p_{T}[GeV];Max #Delta R",25, 0, 1000, 15, 0, 1.5);
  hPt_dRMax_nomatchtopCand->Sumw2();
  hPt_dRMax_matchtopGen = new TH2D("hPt_dRMax_matchtopGen","Top P_{T} vs Max #Delta R;p_{T}[GeV];Max #Delta R",25, 0, 1000, 15, 0, 1.5);
  hPt_dRMax_matchtopGen->Sumw2();
  hdRMax_dRMin_matchtop = new TH2D("hdRMax_dRMin_matchtop","Max #DeltaR vs Min #DeltaR;Max #DeltaR;Min #DeltaR",50, 0, 1.5, 50, 0, 1.5);
  hdRMax_dRMin_matchtop->Sumw2();
  hdRMax_dRMin_nomatchtop = new TH2D("hdRMax_dRMin_nomatchtop","Max #DeltaR vs Min #DeltaR;Max #DeltaR;Min #DeltaR",50, 0, 1.5, 50, 0, 1.5);
  hdRMax_dRMin_nomatchtop->Sumw2();
  hdRMax_dRMin_matchtopCand = new TH2D("hdRMax_dRMin_matchtopCand","Max #DeltaR vs Min #DeltaR;Max #DeltaR;Min #DeltaR",50, 0, 1.5, 50, 0, 1.5);
  hdRMax_dRMin_matchtopCand->Sumw2();
  hdRMax_dRMin_nomatchtopCand = new TH2D("hdRMax_dRMin_nomatchtopCand","Max #DeltaR vs Min #DeltaR;Max #DeltaR;Min #DeltaR",50, 0, 1.5, 50, 0, 1.5);
  hdRMax_dRMin_nomatchtopCand->Sumw2();
  hdRMax_area_matchtop = new TH2D("hdRMax_area_matchtop","Max #DeltaR vs area;Max #DeltaR;Area",100, 0, 1.5, 100, 0, 5);
  hdRMax_area_matchtop->Sumw2();
  hdRMax_area_nomatchtop = new TH2D("hdRMax_area_nomatchtop","Max #DeltaR vs area;Max #DeltaR;Area",100, 0, 1.5, 100, 0, 5);
  hdRMax_area_nomatchtop->Sumw2();
  hdRMax_area_matchtopCand = new TH2D("hdRMax_area_matchtopCand","Max #DeltaR vs area;Max #DeltaR;Area",100, 0, 1.5, 100, 0, 5);
  hdRMax_area_matchtopCand->Sumw2();
  hdRMax_area_nomatchtopCand = new TH2D("hdRMax_area_nomatchtopCand","Max #DeltaR vs area;Max #DeltaR;Area",100, 0, 1.5, 100, 0, 5);
  hdRMax_area_nomatchtopCand->Sumw2();

  hgentopPt_den = new TH1D("hgentopPt_den","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  hgentopPt_den->Sumw2();
  hgentopPt_num = new TH1D("hgentopPt_num","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  hgentopPt_num->Sumw2();
  hfakeNjet_den = new TH1D("hfakeNjet_den","FakeRate in N_{jet} bin;N_{jet};Event",6, 4, 10);
  hfakeNjet_den->Sumw2();
  hfakeHT_den = new TH1D("hfakeHT_den","FakeRate in H_{T} bin;H_{T}[GeV];Event",25, 0, 1000);
  hfakeHT_den->Sumw2();
  hfakeMET_den = new TH1D("hfakeMET_den","FakeRate in p_{T}^{miss} bin;p_{T}^{miss}[GeV];Event",25, 0, 1000);
  hfakeMET_den->Sumw2();
  hfakeNjet_num = new TH1D("hfakeNjet_num","FakeRate in N_{jet} bin;N_{jet};Event",6, 4, 10);
  hfakeNjet_num->Sumw2();
  hfakeHT_num = new TH1D("hfakeHT_num","FakeRate in H_{T} bin;H_{T}[GeV];Event",25, 0, 1000);
  hfakeHT_num->Sumw2();
  hfakeMET_num = new TH1D("hfakeMET_num","FakeRate in p_{T}^{miss} bin;p_{T}^{miss}[GeV];Event",25, 0, 1000);
  hfakeMET_num->Sumw2();

  hTot_fake = new TH1D("hTot_fake", "Fake rate in Zinv", 3, 0, 3);
  hTot_fake->SetBit(TH1::kCanRebin);                                                                                            
  hTot_fake->Sumw2();

  for(unsigned int i = 0; i<maxConst; i++){
    hconstPt.push_back(new TH1D(Title("hconstPt",i),";Event;p_{T}",25, 0, 1000));
    hconstPt.back()->Sumw2();
    hconstM.push_back(new TH1D(Title("hconstM",i),";Event;p_{T}",25, 0, 1000));
    hconstM.back()->Sumw2();
  }

}

TString BaseHistgram::Title(TString spec, unsigned int j){
  TString title = spec;
  title+= "_";
  title+=j;
  return title;
}

bool FillChain(TChain* &chain, const char *subsample, const string condorSpec, const int& startfile, const int& filerun){
  AnaSamples::SampleSet        allSamples = condorSpec.empty()? AnaSamples::SampleSet():AnaSamples::SampleSet(condorSpec);
  AnaSamples::SampleCollection allCollections(allSamples);
  bool find = false;  
  TString subsamplename(subsample);
  chain = new TChain(allSamples[subsample].treePath.c_str());
    if(allSamples[subsample] != allSamples.null())
      {
	allSamples[subsample].addFilesToChain(chain, startfile, filerun);
	find = true;
	Lumiscale = allSamples[subsample].getWeight();  
    }
    return find;
}
vector<TLorentzVector>GetTopLVec(vector<TLorentzVector>genDecayLVec, vector<int>genDecayPdgIdVec, vector<int>genDecayIdxVec, vector<int>genDecayMomIdxVec){
  vector<TLorentzVector> tLVec;
  for(unsigned it=0; it<genDecayLVec.size(); it++){
    int pdgId = genDecayPdgIdVec.at(it);
    if(abs(pdgId)==6){
     for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ){
	  int pdgId = genDecayPdgIdVec.at(ig);
	  if(abs(pdgId)==24){
	    int flag = 0;
	    for(unsigned iq=0; iq<genDecayLVec.size(); iq++){
	      if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) {
		int pdgid = genDecayPdgIdVec.at(iq);	
		if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++;
	      }
	    }
	    if(!flag) tLVec.push_back(genDecayLVec.at(it));
	  }
	}
      }//dau. loop
    }//top cond
  }//genloop
  return tLVec;
}
bool GetMatchedTop(vector<TopObject*> NTop, vector<TopObject*>&MachedNTop, vector<TLorentzVector>topLVec, vector<TLorentzVector>& MtopLVec){
  bool match = false; 
  if(topLVec.size()==0) return match;
  double DeltaR = 0.4;
  for(unsigned nt=0; nt<NTop.size();nt++){
    double deltaRMin = 100000.;
    unsigned tid = -1;
    for(unsigned gent = 0; gent < topLVec.size(); gent++) { // Loop over objects
      const double dr = NTop[nt]->p().DeltaR(topLVec.at(gent));
      if( dr < deltaRMin ) {deltaRMin = dr; tid = gent;}
    }
    if(deltaRMin < DeltaR){
      MachedNTop.push_back(NTop[nt]);
      MtopLVec.push_back(topLVec[tid]);
      match = true;
    }
  }
  return match;
}
bool GetMatchedTopCand(vector<TopObject> NTopCand, vector<TopObject>&MachedNTopCand, vector<TLorentzVector>topLVec, vector<TLorentzVector>& MtopLVec){
  bool match = false; 
  if(topLVec.size()==0) return match;
  double DeltaR = 0.4;
  for(unsigned nt=0; nt<NTopCand.size();nt++){
    if(!(NTopCand[nt].p().M()>minTopCandMass && NTopCand[nt].p().M()<maxTopCandMass)) continue;//mass window
    double deltaRMin = 100000.;
    unsigned tid = -1;
    for(unsigned gent = 0; gent < topLVec.size(); gent++) { // Loop over objects
      const double dr = NTopCand[nt].p().DeltaR(topLVec.at(gent));
      if( dr < deltaRMin ) {deltaRMin = dr; tid = gent;}
    }
    if(deltaRMin < DeltaR){
      MachedNTopCand.push_back(NTopCand[nt]);
      MtopLVec.push_back(topLVec[tid]);
      match = true;
    }
  }
  return match;
}
int GetMatchedTopConst(vector<Constituent const *> topconst, vector<TLorentzVector>gentopdauLVec){
  int match=0;
  if(gentopdauLVec.size()==0) return match;
  double DeltaR = 0.4;
  for(unsigned nt=0; nt<topconst.size();nt++){
    double deltaRMin = 100000.;
    for(unsigned gent = 0; gent < gentopdauLVec.size(); gent++) { // Loop over objects                                                                                           
      const double dr = topconst[nt]->p().DeltaR(gentopdauLVec.at(gent));
      if( dr < deltaRMin ) deltaRMin = dr;
    }
    if(deltaRMin < DeltaR) match++;
  }
  return match;
}
vector<TLorentzVector> GettopdauLVec(TLorentzVector top, vector<TLorentzVector>genDecayLVec, vector<int>genDecayPdgIdVec, vector<int>genDecayIdxVec, vector<int>genDecayMomIdxVec){
  vector<TLorentzVector>topdauLVec;
  for(unsigned it=0; it<genDecayLVec.size(); it++){
    if(genDecayLVec[it]==top){
      for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ){
	  int pdgId = genDecayPdgIdVec.at(ig);
	  if(abs(pdgId)==5)topdauLVec.push_back(genDecayLVec[ig]);
	  if(abs(pdgId)==24){
	    for(unsigned iq=0; iq<genDecayLVec.size(); iq++){
	      if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) {
		int pdgid = genDecayPdgIdVec.at(iq);	
		if(abs(pdgid)!= 11 && abs(pdgid)!= 13 && abs(pdgid)!= 15) topdauLVec.push_back(genDecayLVec[iq]);
	      }
	    }
	  }
	}
      }//dau. loop
    }//top cand.
  }//gen loop
  return topdauLVec;
}
double GetdRmax(TLorentzVector top, vector<TLorentzVector>dauLVec){
  double dRmax = 0;
  for(unsigned dr =0; dr<dauLVec.size(); dr++){
    double dR = dauLVec[dr].DeltaR(top);
    if(dR>dRmax) dRmax = dR;
  }
  return dRmax;
}
double GetdRmin_topcand(TopObject Top){
  vector<Constituent const*> con = Top.getConstituents();
  double dRmin = 999;
  for(unsigned dr =0; dr<con.size(); dr++){
    double dR = con[dr]->p().DeltaR(Top.p());
    if(dR<dRmin) dRmin = dR;
  }
  return dRmin;
}
double GetdRmin_top(TopObject* Top){
  vector<Constituent const*> con = Top->getConstituents();
  double dRmin = 999;
  for(unsigned dr =0; dr<con.size(); dr++){
    double dR = con[dr]->p().DeltaR(Top->p());
    if(dR<dRmin) dRmin = dR;
  }
  return dRmin;
}
double GetArea_topcand(TopObject Top){
  vector<Constituent const*> con = Top.getConstituents();
  vector<double>dR(con.size(),0);
  double area = 0;
  if(con.size()<3)return area;
  dR[0] = con[0]->p().DeltaR(con[1]->p());
  dR[1] = con[0]->p().DeltaR(con[2]->p());
  dR[2] = con[1]->p().DeltaR(con[2]->p());
  double s = (dR[0]+dR[1]+dR[2])/2;
  area = sqrt(s*(s-dR[0])*(s-dR[2])*(s-dR[2]));
  return area;
}
double GetArea_top(TopObject* Top){
  vector<Constituent const*> con = Top->getConstituents();
  vector<double>dR(con.size(),0);
  double area = 0;
  if(con.size()<3)return area;
  dR[0] = con[0]->p().DeltaR(con[1]->p());
  dR[1] = con[0]->p().DeltaR(con[2]->p());
  dR[2] = con[1]->p().DeltaR(con[2]->p());
  double s = (dR[0]+dR[1]+dR[2])/2;
  area = sqrt(s*(s-dR[0])*(s-dR[2])*(s-dR[2]));
  return area;
}
