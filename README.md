# TaggerTest

cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src/
cmsenv
git cms-init
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone -b Ana_DataNtpV9_Final_Lumi git@github.com:susy2015/SusyAnaTools.git
git clone -b master git@github.com:susy2015/TopTagger.git

git clone -b master git@github.com:kmandal/TaggerTest.git

scram b -j9

cd TaggerTest/test
make

#TopCatagorization code: TaggerTest/src/TopCatagory.cc
