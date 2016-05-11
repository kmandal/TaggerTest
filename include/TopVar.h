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

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hNtop_Sig;
  TH1D *hNtop_Bck;
  TH1D *hdR1_Sig;
  TH1D *hdR1_Bck;
  TH1D *hdR2_Sig;
  TH1D *hdR2_Bck;
  TH1D *hdR3_Sig;
  TH1D *hdR3_Bck;
  TH1D *hdRMax_Sig;
  TH1D *hdRMax_Bck;
  TH1D *htopMass_Sig;
  TH1D *htopMass_Bck;
  TH1D *htopPt_Sig;
  TH1D *htopPt_Bck;
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_TopVar"+index+".root";
  oFile = new TFile(filename, "recreate");

  hNtop_Sig = new TH1D("hNtop_Sig", "No. of top;N_{top};Event", 5, 0, 5);
  hNtop_Sig->Sumw2();
  hNtop_Bck = new TH1D("hNtop_Bck", "No.of top;N_{top};Event", 5, 0, 5);
  hNtop_Bck->Sumw2();
  hdR1_Sig = new TH1D("hdR1_Sig", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR1_Sig->Sumw2(); 
  hdR2_Sig = new TH1D("hdR2_Sig", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR2_Sig->Sumw2(); 
  hdR3_Sig = new TH1D("hdR3_Sig", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR3_Sig->Sumw2(); 
  hdRMax_Sig = new TH1D("hdRMax_Sig", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_Sig->Sumw2(); 
  hdR1_Bck = new TH1D("hdR1_Bck", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR1_Bck->Sumw2(); 
  hdR2_Bck = new TH1D("hdR2_Bck", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR2_Bck->Sumw2(); 
  hdR3_Bck = new TH1D("hdR3_Bck", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdR3_Bck->Sumw2(); 
  hdRMax_Bck = new TH1D("hdRMax_Bck", "dR(const,top);#Delta R;Event", 15, 0, 1.5);
  hdRMax_Bck->Sumw2(); 

  htopMass_Sig = new TH1D("htopMass_Sig","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_Sig->Sumw2();
  htopMass_Bck = new TH1D("htopMass_Bck","Top Mass;Mass[GeV];Event",25, 0, 500);
  htopMass_Bck->Sumw2();
  htopPt_Sig = new TH1D("htopPt_Sig","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_Sig->Sumw2();
  htopPt_Bck = new TH1D("htopPt_Bck","Top P_{T};p_{T}[GeV];Event",25, 0, 1000);
  htopPt_Bck->Sumw2();

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
