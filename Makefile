CXX=c++
CC=gcc
CFLAGS= -std=c++17 -O2 -Wall
INSS=-I./include
INCMSSW=-I${CMSSW_BASE}/src
INCMSSW2=-I/cvmfs/cms.cern.ch/${SCRAM_ARCH}/cms/cmssw/CMSSW_10_6_26/src/
INCMSSW3=-I/cvmfs/cms.cern.ch/${SCRAM_ARCH}/external/boost/1.67.0-pafccj/include/
INCMSSW4=-I/cvmfs/cms.cern.ch/${SCRAM_ARCH}/external/clhep/2.3.4.2/include/

CFLAGS += `root-config --cflags`
LIBS += `root-config --ldflags --glibs`
LIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore -lRooStats
#LIBS += -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lJHUGenMELAMELA
LIBS += -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lRecoJetsJetProducers -lSimDataFormatsHTXS -lPhysicsToolsPatAlgos 
LIBS += -L/cvmfs/cms.cern.ch/${SCRAM_ARCH}/cms/cmssw/CMSSW_10_6_26/lib/${SCRAM_ARCH}/ -lCondFormatsJetMETObjects -lFWCoreParameterSet

zz4lOBJ=ZZ4L_Ana.o

.PHONY: clean all main test

all: ZZ4L_Ana

ZZ4L_Ana: ZZ4L_Ana.o 
		$(CXX) -o ZZ4L_Ana.exe $(zz4lOBJ) $(LIBS)

clean:
	@rm *.o *.exe

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} ${INCMSSW2} ${INCMSSW3} ${INCMSSW4}  -I. -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} ${INCMSSW2} ${INCMSSW3} ${INCMSSW4} -I. -c $<

