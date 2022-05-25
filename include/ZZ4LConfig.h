#ifndef ZZ4LConfig_h
#define ZZ4LConfig_h


using namespace std;  

bool isData;
//bool bestCandMela=true;
bool bestCandMela=false;
//bool redoEventSelection=true;
bool redoEventSelection=false;
bool redoJets=true;
//bool redoMela=true;
bool redoMela=false;
bool redoEbE=true;
bool redoRefit=true;
bool redoDifferentialObs=true;
//bool redoGEN=true;
bool redoGEN=false;
double m4lLowCut=70.0;
double mZ2Low=12.0;
double mZ2High=120.0;
double mZ1Low=40.0;
double mZ1High=120.0;
double isoCutEl=0.35;
double isoCutMu=0.35;
double leadingPtCut=20.0;
double subleadingPtCut=10.0;
double BTagCut=0.8484;
int job=-1;
int njobs=-1;
bool doM4lSkim=true;
double GENmassZZLow=50.0;
double GENmassZZHigh=160.0;

#endif
