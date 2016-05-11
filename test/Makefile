ODIR       = obj

TTIDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/include
TPIDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/include

TTSDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/src
TPSDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/src

SDIR = $(CMSSW_BASE)/src/SusyAnaTools/Tools
IDIR = $(CMSSW_BASE)/src/TaggerTest/src

CXX        = g++

CXXFLAGS  += -I$(CMSSW_BASE)/src -std=c++0x
## Optimization/debug flag
CXXFLAGS += -g #-O3
## Enable for maximum warning
#CXXFLAGS += -Wall -Wextra -Weffc++

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

LIBS       = $(shell root-config --glibs)

PROGRAMS = TopVar

LIBRARIES = TopTagger

all: mkobj sampPyWrap $(PROGRAMS)


mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(TTSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TTSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<


TopVar: $(ODIR)/TopVar.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -lTopTaggerTopTagger -lTopTaggerCfgParser -o $@


#TopTagger : $(addprefix $(ODIR)/, $(notdir $(patsubst %.cc, %.o, $(patsubst %.cpp, %.o, $(wildcard $(TTSDIR)/*.cpp $(TTSDIR)/*.cc $(TPSDIR)/*.cpp $(TPSDIR)/*.cc)))))
#	ar rcs lib$@.a $^
#	gcc -shared -o lib$@.so $^

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) $(LIBRARIES) core 

-include $(ODIR)/*.d
