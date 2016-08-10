#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "SusyAnaTools/Tools/customize.h"

namespace Utility {

  void FillDouble(TH1* hist, const double &a, const double &w){
    int nbin = hist->GetNbinsX();
    double low = hist->GetBinLowEdge(nbin);
    double high = hist->GetBinLowEdge(nbin + 1);
    double copy = a;
    if(copy >= high) copy = low;
    hist->Fill(copy, w);
  }

  void FillInt(TH1* hist, const int &a, const double &w){
    int nbin = hist->GetNbinsX();
    int low = hist->GetBinLowEdge(nbin);
    int high = hist->GetBinLowEdge(nbin + 1);
    int copy = a;
    if(copy >= high) copy = low;
    hist->Fill(copy, w);
  }

  void Fill2D(TH2 *hist, const double &a, const double &b, const double &w){
    int nbinx = hist->GetNbinsX();
    int nbiny = hist->GetNbinsY();
    double lowx = hist->GetXaxis()->GetBinLowEdge(nbinx);
    double highx = hist->GetXaxis()->GetBinLowEdge(nbinx + 1);
    double lowy = hist->GetYaxis()->GetBinLowEdge(nbiny);
    double highy = hist->GetYaxis()->GetBinLowEdge(nbiny + 1);
    double copyx = a;
    if(copyx >= highx) copyx = lowx;
    double copyy = b;
    if(copyy >= highy) copyy = lowy;
    hist->Fill(copyx, copyy, w);
  }


  double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){
    const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();
    const double met = metLVec.Pt(), metphi = metLVec.Phi();
    double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );
    return mt;
  }

  void BjetIndx(const std::vector<TLorentzVector> &inputJets, const std::vector<double> &inputCSVS, const double cutCSVS, const AnaConsts::AccRec& jetCutsArr, std::vector<int> &Bjetidx, const int nbjet){
    const double minAbsEta = jetCutsArr.minAbsEta, maxAbsEta = jetCutsArr.maxAbsEta, minPt = jetCutsArr.minPt, maxPt = jetCutsArr.maxPt;
    int cntNJets =0;
    for(int ij=0; ij<inputJets.size(); ij++){
      if( !AnaFunctions::jetPassCuts(inputJets[ij], jetCutsArr) ) continue;
      if( std::isnan(inputCSVS[ij]) ) continue;
      if( inputCSVS[ij] > cutCSVS ) Bjetidx.push_back(ij);
    }
    if(nbjet!=Bjetidx.size())std::cout<<"Check bjet calculation!"<<std::endl;
  }
  double MT_B(const std::vector<int> &Bjetidx, const std::vector<TLorentzVector> &inputJets, const TLorentzVector &met){
    double MTb = 0;
    std::vector<TLorentzVector> bjetsLVec;
    std::vector<double> mt;
    for(int l=0; l<Bjetidx.size();l++){
      bjetsLVec.push_back(inputJets.at(Bjetidx[l]));
      mt.push_back(calcMT(bjetsLVec.back(), met));
    }
    if(mt.size()==1)MTb = mt.at(0);
    if(mt.size()>1) MTb = mt.at(0)<mt.at(1)? mt.at(0):mt.at(1);
    return MTb;
  }

  int Njet(const std::vector<TLorentzVector> &inputJets, const std::string spec){
    int cntNJets = 0;
    if( spec.compare("Pt50Eta24") == 0) cntNJets = AnaFunctions::countJets(inputJets, AnaConsts::pt50Eta24Arr);
    if( spec.compare("Pt30Eta24") == 0) cntNJets = AnaFunctions::countJets(inputJets, AnaConsts::pt30Eta24Arr);
    if( spec.compare("Pt30") == 0)      cntNJets = AnaFunctions::countJets(inputJets, AnaConsts::pt30Arr);
    return cntNJets;
  }
  int Nbjet(const std::vector<TLorentzVector> &inputJets, const std::vector<double> &inputCSVS){
    int cntCSVS = AnaFunctions::countCSVS(inputJets, inputCSVS, AnaConsts::cutCSVS, AnaConsts::bTagArr);
    return cntCSVS;
  }
 
}
