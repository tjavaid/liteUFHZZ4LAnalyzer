// system include files
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>
#include<TApplication.h>
// user include files 
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"
//
#include "ZZ4LConfig.h"
#include "ZZ4LAnalysisTree_jes.h"
//#include "LeptonEfficiency.h"
//#include "PileupWeight.h"
//#include "EbECorrection.h"
//#include "Helper.h"
//#include "JHUGenMELA/MELA/interface/Mela.h"
//#include "JHUGenMELA/MELA/interface/TUtil.hh"
//#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
//#include "KinZfitter/KinZfitter/interface/KinZfitter.h"


#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"


using namespace std;
/////////////////////

void SetNewTree(TTree* newtree);
void ReadTree(TTree* tree, TTree* & newtree, TString filename);

TString filename;
bool debug;
std::unique_ptr<JetCorrectionUncertainty> jecunc;
bool applyJEC_ = true;
bool isMC = true;
std::string jecUncFile_;
std::vector<string> uncSources {};  
std::vector<JetCorrectionUncertainty*> splittedUncerts_;
int year=2016;
//Mela* mela = new Mela(13.0,125.0,TVar::SILENT);

//KalmanMuonCalibrator *kalmanMuonCalibrator;

//KinZfitter *kinZfitter;

int main(int argc, char *argv[])
{    
     
  //debug = false;     
  debug = true;     

  if(argc > 6)  {
      cout<<argv[0]<<" filename "<<argv[1]<<" outfile "<<argv[2]<<" isData "<<argv[3]<<endl;
      return -1;
    }

  /////////////////////

  filename = argv[1];
  TString outfilename = argv[2];

  if(atof(argv[3])>0) { isData = true; }
  else {isData = false; }

  if(atof(argv[4])>0) { 
      job=strtol(argv[4], NULL, 10); njobs=strtol(argv[5], NULL, 10);
      cout<<"job "<<job<<" of "<<njobs<<endl;
      outfilename+="_";
      outfilename+=argv[4];
  }

 
  TFile* infile = TFile::Open(filename+".root");
  TTree* tree; 
  
  tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) tree = (TTree*) infile->Get("passedEvents");
  if(!tree) tree = (TTree*) infile->Get("selectedEvents");
  if(!tree) { cout<<"ERROR could not find the tree for "<<filename<<endl; return -1;}

  TH1F* sumWeights = (TH1F*) infile->Get("Ana/sumWeights");
  sumW = 1.0;
  if(sumWeights) sumW = sumWeights->GetBinContent(1); 

  cout<<"sumW is "<<sumW<<endl;
  TH1F *th[7];
  th[0]=(TH1F*)infile->Get("Ana/nEvents");
  th[1]=(TH1F*)infile->Get("Ana/sumWeights");
  th[2]=(TH1F*)infile->Get("Ana/sumWeightsPU");
  th[3]=(TH1F*)infile->Get("Ana/nVtx");
  th[4]=(TH1F*)infile->Get("Ana/nVtx_ReWeighted");
  th[5]=(TH1F*)infile->Get("Ana/nInteractions");
  th[6]=(TH1F*)infile->Get("Ana/nInteraction_ReWeighted");

  // read tree     

  TString name = outfilename;
  TFile* tmpFile =  new TFile(name+".root","RECREATE");
  tmpFile->mkdir("Ana");
  tmpFile->cd("Ana");
  
//  TTree* newtree = new TTree("passedEvents","passedEvents");
  TTree *newtree = tree->CloneTree(0);
  if(debug)cout<<"start setting new tree "<<endl;

  SetNewTree(newtree);

  if(debug)cout<<"start reading tree "<<endl;

  ReadTree(tree, newtree, filename);

  if(debug)cout<<"end reading tree"<<endl;

  //tmpFile->cd();


  tmpFile->cd("Ana");
  for (int ii=0; ii<7; ii++)
        {
            th[ii]->Write();
        }

  newtree->Write("passedEvents",TObject::kOverwrite);

  tmpFile->Close(); 


}


void ReadTree(TTree* tree, TTree* & newtree, TString filename){
    
    ZZ4LAnalysisTree_jes::setAddresses(tree, filename);

    float npass = 0.0;
    sumweight = 0.0;
    if (redoJets) {
       // split sources
       // JEC split
       if (year == 2016)   {    
             edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
             jecUncFile_ = jecUncFile.fullPath();
             //FILE *jecUncFile_ = fopen("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt","r"); 
             uncSources.push_back("Total");
             uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2016");
             uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2016");
             uncSources.push_back("EC2"); uncSources.push_back("EC2_2016");
             uncSources.push_back("FlavorQCD");
             uncSources.push_back("HF"); uncSources.push_back("HF_2016");
             uncSources.push_back("RelativeBal");
             uncSources.push_back("RelativeSample_2016");
             }    
       if (year == 2017)   {
             edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
             jecUncFile_ = jecUncFile.fullPath();
             //TFile *jecUncFile_ = TFile::Open("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt"); //jecUncFile.fullPath();
             //FILE *jecUncFile_ = fopen("JECUncertaintySources/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt","r");
             //jecUncFile.fullPath();
             uncSources.push_back("Total");
             uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2017");
             uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2017");
             uncSources.push_back("EC2"); uncSources.push_back("EC2_2017");
             uncSources.push_back("FlavorQCD");
             uncSources.push_back("HF"); uncSources.push_back("HF_2017");
             uncSources.push_back("RelativeBal");
             uncSources.push_back("RelativeSample_2017");
             }
       if (year == 2018)   {
             edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
             jecUncFile_ = jecUncFile.fullPath();
             //TFile *jecUncFile_ = TFile::Open("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt"); //jecUncFile.fullPath();
             //FILE *jecUncFile_ = fopen("JECUncertaintySources/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt","r");
             //jecUncFile.fullPath();
             uncSources.push_back("Total");
             uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2018");
             uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2018");
             uncSources.push_back("EC2"); uncSources.push_back("EC2_2018");
             uncSources.push_back("FlavorQCD");
             uncSources.push_back("HF"); uncSources.push_back("HF_2018");
             uncSources.push_back("RelativeBal");
             uncSources.push_back("RelativeSample_2018");
             }



            cout<<"jecUncFile_:   "<<jecUncFile_<<endl;

//            bool applyJEC_ = true;
//            bool isMC = true;
            if(applyJEC_ && isMC)
                    {
		      cout<<"uncSources.size():    "<<uncSources.size()<<endl;
                      for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
                        {
                          JetCorrectorParameters corrParams = JetCorrectorParameters(jecUncFile_, uncSources[s_unc]);
                          splittedUncerts_.push_back(new JetCorrectionUncertainty(corrParams));
			  cout<<"uncSources["<<s_unc<<"]"<<uncSources[s_unc]<<endl;
			  cout<<"splittedUncerts_["<<s_unc<<"]"<<splittedUncerts_[s_unc]<<endl;
                        }
                    }
	}



///  =======================
    if(debug) cout<<"start looping"<<endl;

    int firstevt=0; int lastevt=tree->GetEntries();
    if (job>0) {
        firstevt = tree->GetEntries()*(job-1)/njobs;
        lastevt = tree->GetEntries()*(job)/njobs-1;
    }
    
    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
    
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
       
        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;

        tree->GetEntry(evt);

        if (debug) cout<<"Run: "<<Run<<" Lumi: "<<LumiSect<<" Event: "<<Event<<endl;


        if (debug) cout<<"GENmassZZ"<<GENmassZZ<<endl;
        if (debug) cout<<"passedFullSelection"<<passedFullSelection<<endl;
        if (debug) cout<<"mass4l"<<mass4l<<endl;

        // do a Skim
        if ( doM4lSkim && (GENmassZZ<GENmassZZLow || GENmassZZ>GENmassZZHigh) && !(passedFullSelection==1 && mass4l>GENmassZZLow && mass4l<GENmassZZHigh) ) continue;       
        
        if (debug) cout<<"passed m4l skim "<<endl;


        passTrig=false;
        if (isData) {
            passTrig = passedTrig;
            /*
            // double ele
            if (strstr((*triggersPassed).c_str(),"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v")) passTrig=true;
            // double mu
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) passTrig=true;
            // single ele
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele25_eta2p1_WPTight_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele27_WPTight_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele27_eta2p1_WPLoose_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele32_eta2p1_WPTight_Gsf_v")) passTrig=true;
            // single mu
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu20_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu20_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_eta2p1_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_eta2p1_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu24_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu24_v")) passTrig=true;
            // multi lepton
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_12_10_5_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) passTrig=true;
            */
        } else {
            passTrig = true;
            //pileupWeight = float(puweight(nInt));
            sumweight += pileupWeight*genWeight;
        }

        bool foundHiggsCandidate=false;        

        if (redoEventSelection) {

            if((*lep_id).size()<4 || (*lep_pt).size()<4) goto filltree;
            unsigned int Nlep = (*lep_id).size();
            if (debug) cout<<Nlep<<" leptons in total"<<endl;

            // 2 OSSF Pairs
            bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-13) Nmm = Nmm+1;
                if((*lep_id)[i]==13) Nmp = Nmp+1;
            }
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-11) Nem = Nem+1;
                if((*lep_id)[i]==11) Nep = Nep+1;
            }
            
            if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
            if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
            if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu
            
            if(!properLep_ID) goto filltree;
            
            // First, make all Z candidates including any FSR photons
            const double Zmass = 91.1876;
            int n_Zs=0;
            vector<int> Z_lepindex1;
            vector<int> Z_lepindex2;
            vector<float> Z_pt, Z_eta, Z_phi, Z_mass;
            
            for(unsigned int i=0; i<Nlep; i++){
                for(unsigned int j=i+1; j<Nlep; j++){
                    
                    // same flavor opposite charge
                    if(((*lep_id)[i]+(*lep_id)[j])!=0) continue;
                    
                    TLorentzVector li, lj;
                    li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                    lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
                    
                    TLorentzVector lifsr, ljfsr;
                    lifsr.SetPtEtaPhiM((*lepFSR_pt)[i],(*lepFSR_eta)[i],(*lepFSR_phi)[i],(*lepFSR_mass)[i]);
                    ljfsr.SetPtEtaPhiM((*lepFSR_pt)[j],(*lepFSR_eta)[j],(*lepFSR_phi)[j],(*lepFSR_mass)[j]);
                    
                    TLorentzVector liljfsr = lifsr+ljfsr;
                    
                    if (debug) {
                        cout<<"OSSF pair: i="<<i<<" id1="<<(*lep_id)[i]<<" j="<<j<<" id2="<<(*lep_id)[j]<<" pt1: "
                            <<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
                    }
                    
                    TLorentzVector Z, Z_noFSR;
                    Z = lifsr+ljfsr;
                    Z_noFSR = li+lj;
                    
                    if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;
                    
                    if (Z.M()>0.0) {
                        n_Zs++;
                        Z_pt.push_back(Z.Pt());
                        Z_eta.push_back(Z.Eta());
                        Z_phi.push_back(Z.Phi());
                        Z_mass.push_back(Z.M());
                        Z_lepindex1.push_back(i);
                        Z_lepindex2.push_back(j);
                        if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
                    }
                    
                } // lep i
            } // lep j
            
            // Consider all ZZ candidates
            TLorentzVector Z1Vec, Z2Vec, HVec;
            double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
            double maxZ2SumPt_SR=0.0; double maxZ2SumPt_CR=0.0;
            double max_D_bkg_kin_SR=0.0; double max_D_bkg_kin_CR=0.0;
            bool foundSRCandidate=false;
            
            passedZ4lSelection=false;
            passedFullSelection=false;
            
            //vector<int> lep_Hindex, Z_Hindex;
            vector<int> Z_Hindex;
            for (int i=0; i<4; i++) {
                if (i<2) Z_Hindex.push_back(-1);
                //lep_Hindex.push_back(-1);
                lep_Hindex[i]=-1;
            }
            
            for (int i=0; i<n_Zs; i++) {
                for (int j=i+1; j<n_Zs; j++) {
                    
                    int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
                    int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            
                    
                    if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
                    
                    TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
                    lep_i1.SetPtEtaPhiM((*lepFSR_pt)[i1],(*lepFSR_eta)[i1],(*lepFSR_phi)[i1],(*lepFSR_mass)[i1]);
                    lep_i2.SetPtEtaPhiM((*lepFSR_pt)[i2],(*lepFSR_eta)[i2],(*lepFSR_phi)[i2],(*lepFSR_mass)[i2]);
                    lep_j1.SetPtEtaPhiM((*lepFSR_pt)[j1],(*lepFSR_eta)[j1],(*lepFSR_phi)[j1],(*lepFSR_mass)[j1]);
                    lep_j2.SetPtEtaPhiM((*lepFSR_pt)[j2],(*lepFSR_eta)[j2],(*lepFSR_phi)[j2],(*lepFSR_mass)[j2]);
                    
                    TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
                    lep_i1_nofsr.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
                    lep_i2_nofsr.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
                    lep_j1_nofsr.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);
                    lep_j2_nofsr.SetPtEtaPhiM((*lep_pt)[j2],(*lep_eta)[j2],(*lep_phi)[j2],(*lep_mass)[j2]);
                    
                    TLorentzVector Zi, Zj;
                    Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
                    Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
                    
                    if (debug) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}
                    
                    TLorentzVector Z1, Z2;
                    int Z1index, Z2index;
                    int Z1_lepindex[2] = {0,0};
                    int Z2_lepindex[2] = {0,0};
                    double Z1DeltaM, Z2SumPt;
                    
                    if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                        Z1index = i; Z2index = j;
                        Z1 = Zi; Z2 = Zj;                 
                        if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                        else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                        if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                        else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                        Z1DeltaM = abs(Zi.M()-Zmass); 
                        Z2SumPt = lep_j1.Pt()+lep_j2.Pt();
                    }
                    else { 
                        Z1index = j; Z2index = i;
                        Z1 = Zj; Z2 = Zi; 
                        if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                        else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                        if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                        else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                        Z1DeltaM = abs(Zj.M()-Zmass); 
                        Z2SumPt = lep_i1.Pt()+lep_i2.Pt();
                    }
                    
                    // Check isolation cut (without FSR ) for Z1 leptons
                    if ((*lep_RelIsoNoFSR)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue;
                    if ((*lep_RelIsoNoFSR)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue;
                    // Check tight ID cut for Z1 leptons
                    if (!((*lep_tightId)[Z1_lepindex[0]])) continue;
                    if (!((*lep_tightId)[Z1_lepindex[1]])) continue;
                    
                    // Check Leading and Subleading pt Cut
                    vector<double> allPt;
                    allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
                    allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
                    std::sort(allPt.begin(), allPt.end());
                    if (debug) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut
                                   <<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
                    if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
                    
                    // Check dR(li,lj)>0.02 for any i,j
                    vector<double> alldR;
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_i2_nofsr));
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_j1_nofsr));
                    alldR.push_back(lep_i1_nofsr.DeltaR(lep_j2_nofsr));
                    alldR.push_back(lep_i2_nofsr.DeltaR(lep_j1_nofsr));
                    alldR.push_back(lep_i2_nofsr.DeltaR(lep_j2_nofsr));
                    alldR.push_back(lep_j1_nofsr.DeltaR(lep_j2_nofsr));
                    if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
                    if (*min_element(alldR.begin(),alldR.end())<0.02) continue;
                    
                    // Check M(l+,l-)>4.0 GeV for any OS pair
                    // Do not include FSR photons
                    vector<double> allM;
                    TLorentzVector i1i2;
                    i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
                    TLorentzVector j1j2;
                    j1j2 = (lep_j1_nofsr)+(lep_j2_nofsr); allM.push_back(j1j2.M());            
                    
                    if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
                        TLorentzVector i1j1;
                        i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
                        TLorentzVector i2j2;
                        i2j2 = (lep_i2_nofsr)+(lep_j2_nofsr); allM.push_back(i2j2.M());
                    } else {
                        TLorentzVector i1j2;
                        i1j2 = (lep_i1_nofsr)+(lep_j2_nofsr); allM.push_back(i1j2.M());
                        TLorentzVector i2j1;
                        i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
                    }
                    if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
                    if (*min_element(allM.begin(),allM.end())<4.0) { continue;}
                    
                    // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
                    // only for 4mu or 4e ZZ candidates
                    bool passSmartCut=true;
                    if ( abs((*lep_id)[i1])==abs((*lep_id)[j1])) {
                        TLorentzVector Za, Zb;
                        if ((*lep_id)[i1]==(*lep_id)[j1]) {                  
                            Za = (lep_i1)+(lep_j2);
                            Zb = (lep_i2)+(lep_j1);                    
                        } else {
                            Za = (lep_i1)+(lep_j1);
                            Zb = (lep_i2)+(lep_j2);
                        }                
                        if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                            if (debug) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "
                                           <<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                            if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2Low ) passSmartCut=false;
                        }
                        else {
                            if (debug) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "
                                           <<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                            if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2Low ) passSmartCut=false;
                        }
                        
                    }
                    if (!passSmartCut) continue;
                    
                    if (debug) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
                    if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) || (Z2.M() < mZ2Low) || (Z2.M() > mZ2High) ) continue;
                    if (debug) cout<<" pass Z mass cuts"<<endl;
                    
                    
                    // Signal region if Z2 leptons are both tight ID Iso
                    bool signalRegion=true;
                    if ((*lep_RelIsoNoFSR)[Z2_lepindex[0]]>((abs((*lep_id)[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                    if ((*lep_RelIsoNoFSR)[Z2_lepindex[1]]>((abs((*lep_id)[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                    if (!((*lep_tightId)[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
                    if (!((*lep_tightId)[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID          
                    
                    if (debug) cout<<"signalRegion? "<<signalRegion<<endl;
                    
                    // Check if this candidate has the highest D_bkg_kin
                    vector<TLorentzVector> P4s;
                    P4s.clear();
                    vector<int> tmpIDs;
                    tmpIDs.clear();
                    
                    if (Z1_lepindex[0] == i1) {
                        P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                        if (Z2_lepindex[0] == j1) {
                            P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                        } else {
                            P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                        }
                    } else if (Z1_lepindex[0] == i2) {
                        P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                        if (Z2_lepindex[0] == j1) {
                            P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                        } else {
                            P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                        }
                    } else if (Z1_lepindex[0] == j1) {
                        P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                        if (Z2_lepindex[0] == i1) {
                            P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                        } else {
                            P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                        }
                    } else if (Z1_lepindex[0] == j2) {
                        P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                        if (Z2_lepindex[0] == i1) {
                            P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                        } else {
                            P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                        }
                    }
                    
                    tmpIDs.push_back((*lep_id)[Z1_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z1_lepindex[1]]);
                    tmpIDs.push_back((*lep_id)[Z2_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z2_lepindex[1]]);

/*
                    SimpleParticleCollection_t daughters;
                    daughters.push_back(SimpleParticle_t(tmpIDs[0],P4s[0]));
                    daughters.push_back(SimpleParticle_t(tmpIDs[1],P4s[1]));
                    daughters.push_back(SimpleParticle_t(tmpIDs[2],P4s[2]));
                    daughters.push_back(SimpleParticle_t(tmpIDs[3],P4s[3]));

                    SimpleParticleCollection_t associated;
                    mela->setInputEvent(&daughters, &associated, 0, 0);
                    mela->setCurrentCandidateFromIndex(0);

                    float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                    mela->computeP(me_0plus_JHU_tmp, true);            
                    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
                    mela->computeP(me_qqZZ_MCFM_tmp, true);
                    float D_bkg_kin_tmp = me_0plus_JHU_tmp/(me_0plus_JHU_tmp+me_qqZZ_MCFM_tmp);

                    mela->resetInputEvent(); 
  */                  
//                    if (debug) cout<<"good ZZ candidate, D_bkg_kin: "<<D_bkg_kin_tmp<<" max D_bkg_kin SR: "
  //                                 <<max_D_bkg_kin_SR<<" max D_bkg_kin CR: "<<max_D_bkg_kin_CR<<endl;
                    
                    bool same4l=false;
                    bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
                    for(int l = 0; l < 4; l++){
                        if (lep_Hindex[l]==Z1_lepindex[0]) foundZ11 = true;
                        if (lep_Hindex[l]==Z1_lepindex[1]) foundZ12 = true;
                        if (lep_Hindex[l]==Z2_lepindex[0]) foundZ21 = true;
                        if (lep_Hindex[l]==Z2_lepindex[1]) foundZ22 = true;
                    }
                    same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);
                    
/*                    if (signalRegion) { // Signal Region has priority
                        
                        if (!foundSRCandidate) same4l=false;
                        

                        if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR))) 
                             || (!bestCandMela && Z1DeltaM<=minZ1DeltaM_SR) ) { 
                        //if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR)) {                 
                            
                            max_D_bkg_kin_SR = D_bkg_kin_tmp;
                            minZ1DeltaM_SR = Z1DeltaM;
                            
                            if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_SR) continue;

                            Z_Hindex[0] = Z1index;
                            lep_Hindex[0] = Z1_lepindex[0];
                            lep_Hindex[1] = Z1_lepindex[1];
                        
                            maxZ2SumPt_SR = Z2SumPt;    
                            Z_Hindex[1] = Z2index;
                            lep_Hindex[2] = Z2_lepindex[0];
                            lep_Hindex[3] = Z2_lepindex[1];
                            
                            Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                            massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                            
                            if (debug) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
                            if (HVec.M()>m4lLowCut)  {
                                foundHiggsCandidate=true;                    
                                foundSRCandidate=true;
                            }
                        }
                    } else if (!foundSRCandidate) { // Control regions get second priority
                        

                        if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR)))
                             || (!bestCandMela && Z1DeltaM<=minZ1DeltaM_CR) ) {                 
                            //if ( (!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR) ) {

                            max_D_bkg_kin_CR = D_bkg_kin_tmp;
                            minZ1DeltaM_CR = Z1DeltaM;
                    
                            if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_CR) continue;
                            
                            Z_Hindex[0] = Z1index;
                            lep_Hindex[0] = Z1_lepindex[0];
                            lep_Hindex[1] = Z1_lepindex[1];
                            
                            maxZ2SumPt_CR = Z2SumPt;
                            Z_Hindex[1] = Z2index;
                            lep_Hindex[2] = Z2_lepindex[0];
                            lep_Hindex[3] = Z2_lepindex[1];
                        
                            Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                            massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); 
                            mass4l = HVec.M(); pT4l = HVec.Pt(); 
                            
                            if (debug) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                            if (HVec.M()>m4lLowCut) foundHiggsCandidate=true;                    
                        }
                    }
                    
                    if (debug) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                                   <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
*/                    
                } // Zj
            } // Zi
        }
    filltree:

 //       mass4l=-1.0; pT4l=-1.0; y4l=-9999.0;
//        massZ1=-1.0; massZ2=-1.0;
        
//        cosTheta1=-9999.0; cosTheta2=-9999.0;
//        cosThetaStar=-9999.0; Phi=-9999.0; Phi1=-9999.0;
        
        if ( (redoEventSelection&&foundHiggsCandidate) || passedZ4lSelection || passedZXCRSelection) {
            

            if (debug) cout <<"AGAIN"<<endl;
            if (debug) cout<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                           <<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
                    

            TLorentzVector Lep1, Lep2, Lep3, Lep4;
            TLorentzVector Jet1, Jet2, Jet1_2p5, Jet2_2p5;            
            TLorentzVector Jet1_jesup, Jet2_jesup, Jet1_2p5_jesup, Jet2_2p5_jesup;            
            TLorentzVector Jet1_jesdn, Jet2_jesdn, Jet1_2p5_jesdn, Jet2_2p5_jesdn;            
            TLorentzVector nullFourVector(0, 0, 0, 0);                 
            Lep1.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[0]],(*lepFSR_eta)[lep_Hindex[0]],(*lepFSR_phi)[lep_Hindex[0]],(*lepFSR_mass)[lep_Hindex[0]]);
            Lep2.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[1]],(*lepFSR_eta)[lep_Hindex[1]],(*lepFSR_phi)[lep_Hindex[1]],(*lepFSR_mass)[lep_Hindex[1]]);
            Lep3.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[2]],(*lepFSR_eta)[lep_Hindex[2]],(*lepFSR_phi)[lep_Hindex[2]],(*lepFSR_mass)[lep_Hindex[2]]);
            Lep4.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[3]],(*lepFSR_eta)[lep_Hindex[3]],(*lepFSR_phi)[lep_Hindex[3]],(*lepFSR_mass)[lep_Hindex[3]]);

            if (redoEventSelection) {                
                nZXCRFailedLeptons=0;
                for(unsigned int i = 0; i <= 3; i++) {
                    if (!(abs((*lep_id)[lep_Hindex[i]])==11 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutEl)) &&
                        !(abs((*lep_id)[lep_Hindex[i]])==13 && ((*lep_tightId)[lep_Hindex[i]] && (*lep_RelIsoNoFSR)[lep_Hindex[i]]<isoCutMu))){ 
                        nZXCRFailedLeptons++; 
                    }
                }
                
                if (debug) cout << nZXCRFailedLeptons<<" failing leptons in higgs candidate"<<endl;
                
                if (nZXCRFailedLeptons>0) { // at least one lepton has failed 
                    passedZ4lZXCRSelection = true;
                    if ((Lep1+Lep2).M() > mZ2Low && passedTrig) passedZXCRSelection = true;
                } else { //  signal region candidate                    
                    passedZ4lSelection = true;
                    if((Lep1+Lep2).M() > mZ2Low && passedTrig) passedFullSelection = true;
                }
                if (passedZ4lSelection) npass+=1;
                
            }
/*
            int jet1index=-1, jet2index=-1;
            float jet1pt=0.0, jet2pt=0.0;
            int jet1index2p5=-1, jet2index2p5=-1;
            float jet1pt2p5=0.0, jet2pt2p5=0.0;

            int jet1index_jesup=-1, jet2index_jesup=-1;
            float jet1pt_jesup=0.0, jet2pt_jesup=0.0;
            int jet1index2p5_jesup=-1, jet2index2p5_jesup=-1;
            float jet1pt2p5_jesup=0.0, jet2pt2p5_jesup=0.0;

            int jet1index_jesdn=-1, jet2index_jesdn=-1;
            float jet1pt_jesdn=0.0, jet2pt_jesdn=0.0;
            int jet1index2p5_jesdn=-1, jet2index2p5_jesdn=-1;
            float jet1pt2p5_jesdn=0.0, jet2pt2p5_jesdn=0.0;
*/
            if (redoJets) {

		// split sources
		jes_unc_split_Total=-999; 
		jes_unc_split_Abs=-999; 
		jes_unc_split_Abs_year=-999; 
		jes_unc_split_BBEC1=-999; 
		jes_unc_split_BBEC1_year=-999;
		jes_unc_split_EC2=-999; 
		jes_unc_split_EC2_year=-999; 
		jes_unc_split_FlavQCD=-999; 
		jes_unc_split_HF=-999; 
		jes_unc_split_HF_year=-999;
		jes_unc_split_RelBal=-999; 
		jes_unc_split_RelSample_year=-999;
	// up
		pt_jesup_split_Total=-999; 
		pt_jesup_split_Abs=-999; 
		pt_jesup_split_Abs_year=-999; 
		pt_jesup_split_BBEC1=-999; 
		pt_jesup_split_BBEC1_year=-999;
		pt_jesup_split_EC2=-999; 
		pt_jesup_split_EC2_year=-999; 
		pt_jesup_split_FlavQCD=-999; 
		pt_jesup_split_HF=-999; 
		pt_jesup_split_HF_year=-999;
		pt_jesup_split_RelBal=-999; 
		pt_jesup_split_RelSample_year=-999;
	// dn
		pt_jesdn_split_Total=-999; 
		pt_jesdn_split_Abs=-999; 
		pt_jesdn_split_Abs_year=-999; 
		pt_jesdn_split_BBEC1=-999; 
		pt_jesdn_split_BBEC1_year=-999;
		pt_jesdn_split_EC2=-999; 
		pt_jesdn_split_EC2_year=-999; 
		pt_jesdn_split_FlavQCD=-999; 
		pt_jesdn_split_HF=-999; 
		pt_jesdn_split_HF_year=-999;
		pt_jesdn_split_RelBal=-999; 
		pt_jesdn_split_RelSample_year=-999;

                njets_pt30_eta4p7=0;
                njets_pt30_eta2p5=0;
                nbjets_pt30_eta4p7=0;

                njets_pt30_eta4p7_jesup=0;
                njets_pt30_eta2p5_jesup=0;

                njets_pt30_eta4p7_jesdn=0;
                njets_pt30_eta2p5_jesdn=0;

                // Nominal
                jet_iscleanH4l->clear();
//
    	        vector<float> jes_unc_split {};
                vector<float> pt_jesup_split {};
                vector<float> pt_jesdn_split {};
                float singleContr_jes_unc = 0;


                for( unsigned int k = 0; k<(*jet_pt).size(); k++) {
		                   
                    if ((*jet_pt)[k]<30.0 || abs((*jet_eta)[k])>4.7) continue;

                    TLorentzVector thisJet;
                    thisJet.SetPtEtaPhiM((*jet_pt)[k],(*jet_eta)[k],(*jet_phi)[k],(*jet_mass)[k]);


		    bool applyJEC_ = true;  // try

   	            if(applyJEC_ && isMC)
		      {
			cout<<" uncSources.size():    "<< uncSources.size()<<endl;
			for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
			  {
			    singleContr_jes_unc = 0;
			    splittedUncerts_[s_unc]->setJetEta(thisJet.Eta());
			    splittedUncerts_[s_unc]->setJetPt(thisJet.Pt());
			    singleContr_jes_unc = splittedUncerts_[s_unc]->getUncertainty(true); 
			    jes_unc_split.push_back(singleContr_jes_unc);
			    pt_jesup_split.push_back( thisJet.Pt() * (1.0 + singleContr_jes_unc));
			    pt_jesdn_split.push_back( thisJet.Pt() * (1.0 - singleContr_jes_unc));
			    cout<<"thisJet.Pt():   "<<thisJet.Pt()<<endl;
			    cout<<"Uncertainty source is:  uncSources["<<s_unc<<"]  "<<uncSources[s_unc]<<endl;
			    cout<<"jes_unc_split["<<s_unc<<"]  "<<jes_unc_split[s_unc]<<endl;
			    cout<<"pt_jesup_split["<<s_unc<<"]  "<<pt_jesup_split[s_unc]<<endl;
			    cout<<"pt_jesdn_split["<<s_unc<<"]  "<<pt_jesdn_split[s_unc]<<endl;
			  }
		      }
		    else
		      {
			for(unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
			  {
			    jes_unc_split.push_back(-999.);
			    pt_jesup_split.push_back(-999.);
			    pt_jesdn_split.push_back(-999.);
			  }
		      }
//filling variables
	//nominal
                jes_unc_split_Total=jes_unc_split[0]; 
                jes_unc_split_Abs=jes_unc_split[1]; 
                jes_unc_split_Abs_year=jes_unc_split[2]; 
                jes_unc_split_BBEC1=jes_unc_split[3]; 
                jes_unc_split_BBEC1_year=jes_unc_split[4];
                jes_unc_split_EC2=jes_unc_split[5]; 
                jes_unc_split_EC2_year=jes_unc_split[6]; 
                jes_unc_split_FlavQCD=jes_unc_split[7]; 
                jes_unc_split_HF=jes_unc_split[8]; 
                jes_unc_split_HF_year=jes_unc_split[9];
                jes_unc_split_RelBal=jes_unc_split[10]; 
                jes_unc_split_RelSample_year=jes_unc_split[11];
	//up
                pt_jesup_split_Total=pt_jesup_split[0]; 
                pt_jesup_split_Abs=pt_jesup_split[1]; 
                pt_jesup_split_Abs_year=pt_jesup_split[2]; 
                pt_jesup_split_BBEC1=pt_jesup_split[3]; 
                pt_jesup_split_BBEC1_year=pt_jesup_split[4];
                pt_jesup_split_EC2=pt_jesup_split[5]; 
                pt_jesup_split_EC2_year=pt_jesup_split[6]; 
                pt_jesup_split_FlavQCD=pt_jesup_split[7]; 
                pt_jesup_split_HF=pt_jesup_split[8]; 
                pt_jesup_split_HF_year=pt_jesup_split[9];
                pt_jesup_split_RelBal=pt_jesup_split[10]; 
                pt_jesup_split_RelSample_year=pt_jesup_split[11];
	//dn	
                pt_jesdn_split_Total=pt_jesdn_split[0]; 
                pt_jesdn_split_Abs=pt_jesdn_split[1]; 
                pt_jesdn_split_Abs_year=pt_jesdn_split[2]; 
                pt_jesdn_split_BBEC1=pt_jesdn_split[3]; 
                pt_jesdn_split_BBEC1_year=pt_jesdn_split[4];
                pt_jesdn_split_EC2=pt_jesdn_split[5]; 
                pt_jesdn_split_EC2_year=pt_jesdn_split[6]; 
                pt_jesdn_split_FlavQCD=pt_jesdn_split[7]; 
                pt_jesdn_split_HF=pt_jesdn_split[8]; 
                pt_jesdn_split_HF_year=pt_jesdn_split[9];
                pt_jesdn_split_RelBal=pt_jesdn_split[10]; 
                pt_jesdn_split_RelSample_year=pt_jesdn_split[11];



                    
/*                    bool isclean_H4l=true;
                    
                    for (unsigned int i=0; i<(*lep_pt).size(); i++) {
                        bool passed_idiso=true;
                        if (abs((*lep_id)[i])==13 && (*lep_RelIsoNoFSR)[i]>isoCutMu) passed_idiso=false;
                        if (abs((*lep_id)[i])==11 && (*lep_RelIsoNoFSR)[i]>isoCutEl) passed_idiso=false;
                        if (!(*lep_tightId)[i]) passed_idiso=false;
                        bool candlep=false;
                        for (unsigned int l = 0; l <= 3; l++) {
                            if ((int)i==lep_Hindex[l]) candlep=true;
                        }
                        if (!(passed_idiso || candlep)) continue;
                        TLorentzVector thisLep;
                        thisLep.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                        if (thisLep.DeltaR(thisJet)<0.4) isclean_H4l=false;

                    }

                    for(unsigned int i=0; i<(*fsrPhotons_pt).size(); i++) {
            
                        // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
                        if (!(*lep_tightId)[(*fsrPhotons_lepindex)[i]]) continue;
                        double RelIsoNoFSR=(*lep_RelIsoNoFSR)[(*fsrPhotons_lepindex)[i]];
                        if (RelIsoNoFSR>((abs((*lep_id)[(*fsrPhotons_lepindex)[i]])==11) ? isoCutEl : isoCutMu)) continue;
                        
                        TLorentzVector thisPho;
                        thisPho.SetPtEtaPhiM((*fsrPhotons_pt)[i],(*fsrPhotons_eta)[i],(*fsrPhotons_phi)[i],0.0);
                        if (thisPho.DeltaR(thisJet)<0.4) isclean_H4l = false;
                    }

                    if (isclean_H4l) {
                        njets_pt30_eta4p7+=1;  
                        if ((*jet_csvv2)[k]>BTagCut) nbjets_pt30_eta4p7++;
                        jet_iscleanH4l->push_back(k);
                        if (thisJet.Pt()>jet1pt) {
                            jet2pt=jet1pt; jet2index=jet1index;
                            jet1pt=thisJet.Pt(); jet1index=k;
                        } else if (thisJet.Pt()>jet2pt) {
                            jet2pt=thisJet.Pt(); jet2index=k;
                        }
                        if (abs((*jet_eta)[k])<2.5) {
                            njets_pt30_eta2p5+=1;
                            if (thisJet.Pt()>jet1pt2p5) {
                                jet2pt2p5=jet1pt2p5; jet2index2p5=jet1index2p5;
                                jet1pt2p5=thisJet.Pt(); jet1index2p5=k;
                            } else if (thisJet.Pt()>jet2pt2p5) {
                                jet2pt2p5=thisJet.Pt(); jet2index2p5=k;
                            }

                        }
                            
                    }
                    
*/
                }  // end jets loop (nominal)

//                if (debug) cout<<njets_pt30_eta4p7<<" jets"<<endl;
/*                
                // JES up
                jet_jesup_iscleanH4l->clear();
                for( unsigned int k = 0; k<(*jet_jesup_pt).size(); k++) {
                    
                    if ((*jet_jesup_pt)[k]<30.0 || abs((*jet_jesup_eta)[k])>4.7) continue;

                    TLorentzVector thisJet_jesup;
                    thisJet_jesup.SetPtEtaPhiM((*jet_jesup_pt)[k],(*jet_jesup_eta)[k],(*jet_jesup_phi)[k],(*jet_jesup_mass)[k]);
                    
                    bool isclean_H4l=true;
                    
                    for (unsigned int i=0; i<(*lep_pt).size(); i++) {
                        bool passed_idiso=true;
                        if (abs((*lep_id)[i])==13 && (*lep_RelIsoNoFSR)[i]>isoCutMu) passed_idiso=false;
                        if (abs((*lep_id)[i])==11 && (*lep_RelIsoNoFSR)[i]>isoCutEl) passed_idiso=false;
                        if (!(*lep_tightId)[i]) passed_idiso=false;
                        bool candlep=false;
                        for (unsigned int l = 0; l <= 3; l++) {
                            if ((int)i==lep_Hindex[l]) candlep=true;
                        }
                        if (!(passed_idiso || candlep)) continue;
                        TLorentzVector thisLep;
                        thisLep.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                        if (thisLep.DeltaR(thisJet_jesup)<0.4) isclean_H4l=false;

                    }

                    for(unsigned int i=0; i<(*fsrPhotons_pt).size(); i++) {
            
                        // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
                        if (!(*lep_tightId)[(*fsrPhotons_lepindex)[i]]) continue;
                        double RelIsoNoFSR=(*lep_RelIsoNoFSR)[(*fsrPhotons_lepindex)[i]];
                        if (RelIsoNoFSR>((abs((*lep_id)[(*fsrPhotons_lepindex)[i]])==11) ? isoCutEl : isoCutMu)) continue;
                        
                        TLorentzVector thisPho;
                        thisPho.SetPtEtaPhiM((*fsrPhotons_pt)[i],(*fsrPhotons_eta)[i],(*fsrPhotons_phi)[i],0.0);
                        if (thisPho.DeltaR(thisJet_jesup)<0.4) isclean_H4l = false;
                    }

                    if (isclean_H4l) {
                        njets_pt30_eta4p7_jesup+=1;  
                        jet_jesup_iscleanH4l->push_back(k);
                        if (thisJet_jesup.Pt()>jet1pt_jesup) {
                            jet2pt_jesup=jet1pt_jesup; jet2index_jesup=jet1index_jesup;
                            jet1pt_jesup=thisJet_jesup.Pt(); jet1index_jesup=k;
                        } else if (thisJet_jesup.Pt()>jet2pt_jesup) {
                            jet2pt_jesup=thisJet_jesup.Pt(); jet2index_jesup=k;
                        }
                        if (abs((*jet_jesup_eta)[k])<2.5) {
                            njets_pt30_eta2p5_jesup+=1;
                            if (thisJet_jesup.Pt()>jet1pt2p5_jesup) {
                                jet2pt2p5_jesup=jet1pt2p5_jesup; jet2index2p5_jesup=jet1index2p5_jesup;
                                jet1pt2p5_jesup=thisJet_jesup.Pt(); jet1index2p5_jesup=k;
                            } else if (thisJet_jesup.Pt()>jet2pt2p5_jesup) {
                                jet2pt2p5_jesup=thisJet_jesup.Pt(); jet2index2p5_jesup=k;
                            }

                        }
                            
                    }
                    
                }
                if (debug) cout<<njets_pt30_eta4p7_jesup<<" jets (jesup)"<<endl;

                // JES dn
                jet_jesdn_iscleanH4l->clear();
                for( unsigned int k = 0; k<(*jet_jesdn_pt).size(); k++) {
                    
                    if ((*jet_jesdn_pt)[k]<30.0 || abs((*jet_jesdn_eta)[k])>4.7) continue;

                    TLorentzVector thisJet_jesdn;
                    thisJet_jesdn.SetPtEtaPhiM((*jet_jesdn_pt)[k],(*jet_jesdn_eta)[k],(*jet_jesdn_phi)[k],(*jet_jesdn_mass)[k]);
                    
                    bool isclean_H4l=true;
                    
                    for (unsigned int i=0; i<(*lep_pt).size(); i++) {
                        bool passed_idiso=true;
                        if (abs((*lep_id)[i])==13 && (*lep_RelIsoNoFSR)[i]>isoCutMu) passed_idiso=false;
                        if (abs((*lep_id)[i])==11 && (*lep_RelIsoNoFSR)[i]>isoCutEl) passed_idiso=false;
                        if (!(*lep_tightId)[i]) passed_idiso=false;
                        bool candlep=false;
                        for (unsigned int l = 0; l <= 3; l++) {
                            if ((int)i==lep_Hindex[l]) candlep=true;
                        }
                        if (!(passed_idiso || candlep)) continue;
                        TLorentzVector thisLep;
                        thisLep.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                        if (thisLep.DeltaR(thisJet_jesdn)<0.4) isclean_H4l=false;

                    }

                    for(unsigned int i=0; i<(*fsrPhotons_pt).size(); i++) {
            
                        // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
                        if (!(*lep_tightId)[(*fsrPhotons_lepindex)[i]]) continue;
                        double RelIsoNoFSR=(*lep_RelIsoNoFSR)[(*fsrPhotons_lepindex)[i]];
                        if (RelIsoNoFSR>((abs((*lep_id)[(*fsrPhotons_lepindex)[i]])==11) ? isoCutEl : isoCutMu)) continue;
                        
                        TLorentzVector thisPho;
                        thisPho.SetPtEtaPhiM((*fsrPhotons_pt)[i],(*fsrPhotons_eta)[i],(*fsrPhotons_phi)[i],0.0);
                        if (thisPho.DeltaR(thisJet_jesdn)<0.4) isclean_H4l = false;
                    }

                    if (isclean_H4l) {
                        njets_pt30_eta4p7_jesdn+=1;  
                        jet_jesdn_iscleanH4l->push_back(k);
                        if (thisJet_jesdn.Pt()>jet1pt_jesdn) {
                            jet2pt_jesdn=jet1pt_jesdn; jet2index_jesdn=jet1index_jesdn;
                            jet1pt_jesdn=thisJet_jesdn.Pt(); jet1index_jesdn=k;
                        } else if (thisJet_jesdn.Pt()>jet2pt_jesdn) {
                            jet2pt_jesdn=thisJet_jesdn.Pt(); jet2index_jesdn=k;
                        }
                        if (abs((*jet_jesdn_eta)[k])<2.5) {
                            njets_pt30_eta2p5_jesdn+=1;
                            if (thisJet_jesdn.Pt()>jet1pt2p5_jesdn) {
                                jet2pt2p5_jesdn=jet1pt2p5_jesdn; jet2index2p5_jesdn=jet1index2p5_jesdn;
                                jet1pt2p5_jesdn=thisJet_jesdn.Pt(); jet1index2p5_jesdn=k;
                            } else if (thisJet_jesdn.Pt()>jet2pt2p5_jesdn) {
                                jet2pt2p5_jesdn=thisJet_jesdn.Pt(); jet2index2p5_jesdn=k;
                            }

                        }
                            
                    }
                    
                }
                if (debug) cout<<njets_pt30_eta4p7_jesdn<<" jets (jes dn)"<<endl;
                
*/
            } // else {
//                jet1index=(*jet_iscleanH4l)[0]; jet2index=(*jet_iscleanH4l)[1];
//            }
/*
            idL1 = (*lep_id)[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
            idL2 = (*lep_id)[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
            idL3 = (*lep_id)[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
            idL4 = (*lep_id)[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta();

            vector<TLorentzVector> P4s; vector<int> tmpIDs;             
            P4s.push_back(Lep1); P4s.push_back(Lep2);
            P4s.push_back(Lep3); P4s.push_back(Lep4);
            tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
            tmpIDs.push_back(idL3); tmpIDs.push_back(idL4);

            TLorentzVector higgs_undec = Lep1+Lep2+Lep3+Lep4;            

            TLorentzVector LepNoFSR1, LepNoFSR2, LepNoFSR3, LepNoFSR4;
            LepNoFSR1.SetPtEtaPhiM((*lep_pt)[lep_Hindex[0]],(*lep_eta)[lep_Hindex[0]],(*lep_phi)[lep_Hindex[0]],(*lep_mass)[lep_Hindex[0]]);
            LepNoFSR2.SetPtEtaPhiM((*lep_pt)[lep_Hindex[1]],(*lep_eta)[lep_Hindex[1]],(*lep_phi)[lep_Hindex[1]],(*lep_mass)[lep_Hindex[1]]);
            LepNoFSR3.SetPtEtaPhiM((*lep_pt)[lep_Hindex[2]],(*lep_eta)[lep_Hindex[2]],(*lep_phi)[lep_Hindex[2]],(*lep_mass)[lep_Hindex[2]]);
            LepNoFSR4.SetPtEtaPhiM((*lep_pt)[lep_Hindex[3]],(*lep_eta)[lep_Hindex[3]],(*lep_phi)[lep_Hindex[3]],(*lep_mass)[lep_Hindex[3]]);

            vector<TLorentzVector> P4sNoFSR; 
            P4sNoFSR.push_back(LepNoFSR1); P4sNoFSR.push_back(LepNoFSR2);
            P4sNoFSR.push_back(LepNoFSR3); P4sNoFSR.push_back(LepNoFSR4);

            vector<double> lep_pterr_correction;
            for (int i=0; i<4; i++) {
                float pti = (*lep_pt)[lep_Hindex[i]];
                float pterri = (*lep_pterr)[lep_Hindex[i]];
                float absetai = abs((*lep_eta)[lep_Hindex[i]]);
                int idi = abs((*lep_id)[lep_Hindex[i]]);               
                bool ecaldriveni = (*lep_ecalDriven)[lep_Hindex[i]];               
                if (redoEbE) lep_pterr_correction.push_back(EbeCorrection(idi,pti,pterri,absetai,ecaldriveni,isData));
                else lep_pterr_correction.push_back(1.0);
            }
            vector<float> ptErrs;
            ptErrs.push_back(lep_pterr_correction[0]*(*lep_pterr)[lep_Hindex[0]]); ptErrs.push_back(lep_pterr_correction[1]*(*lep_pterr)[lep_Hindex[1]]);
            ptErrs.push_back(lep_pterr_correction[2]*(*lep_pterr)[lep_Hindex[2]]); ptErrs.push_back(lep_pterr_correction[3]*(*lep_pterr)[lep_Hindex[3]]);

            map<unsigned int, TLorentzVector> selectedFsrMap;
            vector<float> phoPtErrs;
            for(unsigned int i = 0; i<4; i++) {
                selectedFsrMap[i] = nullFourVector;
                phoPtErrs.push_back(0.0);
            }

            for(unsigned int i = 0; i<4; i++) {
                int index = lep_Hindex[i];
                for (unsigned int j = 0; j<(*fsrPhotons_pt).size();j++) {
                    if ( (*fsrPhotons_lepindex)[j]==index ) {
                        TLorentzVector pho;
                        pho.SetPtEtaPhiM((*fsrPhotons_pt)[j],(*fsrPhotons_eta)[j],(*fsrPhotons_phi)[j],0.0);
                        selectedFsrMap[i] = pho;
                        phoPtErrs[i] = (*fsrPhotons_pterr)[j];
                    }
                }
            }
*/

/*            if (redoRefit) {
            
                kinZfitter->Setup(P4sNoFSR, tmpIDs, ptErrs, selectedFsrMap, phoPtErrs);
                kinZfitter->KinRefitZ();
                mass4lREFIT = (float)kinZfitter->GetRefitM4l();
                mass4lErrREFIT = (float)kinZfitter->GetRefitM4lErrFullCov();
                mass4lErr = (float)kinZfitter->GetM4lErr();
                massZ1REFIT = (float)kinZfitter->GetRefitMZ1(); 
                massZ2REFIT = (float)kinZfitter->GetRefitMZ2(); 
            }
*/
//            massZ1 = (Lep1+Lep2).M(); massZ2 = (Lep3+Lep4).M(); 
//            mass4l = higgs_undec.M(); pT4l = higgs_undec.Pt();
//            y4l = higgs_undec.Rapidity();

            if (!isData) {
                isH4l = ( (((*lep_genindex)[passedFullSelection*lep_Hindex[0]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[0]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[0]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[0]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[1]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[1]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[1]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[1]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[2]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[2]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[2]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[2]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[3]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[3]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[3]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[3]])]==23) );
            }
/*
            TUtil::computeAngles(Lep1, idL1, Lep2, idL2, Lep3, idL3, Lep4, idL4, \
                                 cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);

            if (abs(idL1)==11 && abs(idL3)==11) {mass4e=mass4l; mass4mu=-1.0; mass2e2mu=-1.0;}
            else if (abs(idL1)==13 && abs(idL3)==13) {mass4e=-1.0; mass4mu=mass4l; mass2e2mu=-1.0;}
            else if (abs(idL1)!=abs(idL3)) {mass4e=-1.0; mass4mu=-1.0; mass2e2mu=mass4l;}

            if (debug) cout<<"mass4l "<<mass4l<<" mass4lErr "<<mass4lErr<<" mass4lREFIT "<<mass4lREFIT<<" mass4lErrREFIT "<<mass4lErrREFIT<<endl;

            if (redoMela) {

                SimpleParticleCollection_t daughters;
                SimpleParticleCollection_t associated;               

                if (foundHiggsCandidate) {
                    daughters.push_back(SimpleParticle_t(idL1, Lep1));
                    daughters.push_back(SimpleParticle_t(idL2, Lep2));
                    daughters.push_back(SimpleParticle_t(idL3, Lep3));
                    daughters.push_back(SimpleParticle_t(idL4, Lep4));
                } else {
                    daughters.push_back(SimpleParticle_t((*lep_id)[0], Lep1));
                    daughters.push_back(SimpleParticle_t((*lep_id)[1], Lep2));
                    daughters.push_back(SimpleParticle_t((*lep_id)[2], Lep3));
                    daughters.push_back(SimpleParticle_t((*lep_id)[3], Lep4));
                }




                // Nominal
                pTj1=-1.0; pTj2=-1.0;                
                qgj1=-1.0; qgj2=-1.0;
                csvj1=-10.0; csvj2=-10.0;
                pTj1_2p5=-1.0; pTj2=-1.0;
                mj1j2=-1.0; dEtaj1j2=-1.0;

                yj1=-1.0; yj2=-1.0;
                dPhiHj1=-1.0; dyHj1=-1.0;
                mj1j2=-1.0; dEtaj1j2=-1.0;
                dPhij1j2=-1.0; dPhiHj1j2=-1.0;
                yj1_2p5=-1.0; yj2_2p5=-1.0;
                
                dPhiHj1_2p5=-1.0; dyHj1_2p5=-1.0;
                mj1j2_2p5=-1.0; dEtaj1j2_2p5=-1.0;
                dPhij1j2_2p5=-1.0; dPhiHj1j2_2p5=-1.0;
                
                if (njets_pt30_eta4p7 > 0) {
                    Jet1.SetPtEtaPhiM((*jet_pt)[jet1index],(*jet_eta)[jet1index],(*jet_phi)[jet1index],(*jet_mass)[jet1index]);
//                    associated.push_back(SimpleParticle_t(0, Jet1));
                    pTj1=Jet1.Pt(); 
                    etaj1=Jet1.Eta();
                    yj1=Jet1.Rapidity();
                    qgj1=(*jet_QGTagger)[jet1index]; 
                    csvj1=(*jet_csvv2)[jet1index]; 
                    dPhiHj1=TMath::Abs(higgs_undec.Phi()-Jet1.Phi());
                    dyHj1=TMath::Abs(y4l-yj1);
                }
                if (njets_pt30_eta4p7 > 1) {
                    Jet2.SetPtEtaPhiM((*jet_pt)[jet2index],(*jet_eta)[jet2index],(*jet_phi)[jet2index],(*jet_mass)[jet2index]);
//                    associated.push_back(SimpleParticle_t(0, Jet2));
                    pTj2=Jet2.Pt();
                    etaj2=Jet2.Eta();
                    yj2=Jet2.Rapidity();
                    qgj2=(*jet_QGTagger)[jet2index];
                    csvj2=(*jet_csvv2)[jet2index];
                    mj1j2=(Jet1+Jet2).M();                    
                    dEtaj1j2=TMath::Abs(Jet1.Eta()-Jet2.Eta());
                    dPhij1j2=TMath::Abs(Jet1.Phi()-Jet2.Phi());
                    dPhiHj1j2=TMath::Abs(higgs_undec.Phi()-(Jet1+Jet2).Phi());
                }
                if (njets_pt30_eta2p5 > 0) {
                    Jet1_2p5.SetPtEtaPhiM((*jet_pt)[jet1index2p5],(*jet_eta)[jet1index2p5],(*jet_phi)[jet1index2p5],(*jet_mass)[jet1index2p5]);
                    pTj1_2p5=Jet1_2p5.Pt();
                    yj1_2p5=Jet1_2p5.Rapidity();
                    dPhiHj1_2p5=TMath::Abs(higgs_undec.Phi()-Jet1_2p5.Phi());
                    dyHj1_2p5=TMath::Abs(y4l-yj1_2p5);
                }
                if (njets_pt30_eta2p5 > 1) {
                    Jet2_2p5.SetPtEtaPhiM((*jet_pt)[jet2index2p5],(*jet_eta)[jet2index2p5],(*jet_phi)[jet2index2p5],(*jet_mass)[jet2index2p5]);
                    pTj2_2p5=Jet2_2p5.Pt();
                    yj2_2p5=Jet2_2p5.Rapidity();
                    mj1j2_2p5=(Jet1_2p5+Jet2_2p5).M();                    
                    dEtaj1j2_2p5=TMath::Abs(Jet1_2p5.Eta()-Jet2_2p5.Eta());
                    dPhij1j2_2p5=TMath::Abs(Jet1_2p5.Phi()-Jet2_2p5.Phi());
                    dPhiHj1j2_2p5=TMath::Abs(higgs_undec.Phi()-(Jet1_2p5+Jet2_2p5).Phi());
                }


                // JES up
                pTj1_jesup=-1.0; pTj2_jesup=-1.0;                
                pTj1_2p5_jesup=-1.0; pTj2_jesup=-1.0;
                mj1j2_jesup=-1.0; dEtaj1j2_jesup=-1.0;

                yj1_jesup=-1.0; yj2_jesup=-1.0;
                dPhiHj1_jesup=-1.0; dyHj1_jesup=-1.0;
                mj1j2_jesup=-1.0; dEtaj1j2_jesup=-1.0;
                dPhij1j2_jesup=-1.0; dPhiHj1j2_jesup=-1.0;
                yj1_2p5_jesup=-1.0; yj2_2p5_jesup=-1.0;
                
                dPhiHj1_2p5_jesup=-1.0; dyHj1_2p5_jesup=-1.0;
                mj1j2_2p5_jesup=-1.0; dEtaj1j2_2p5_jesup=-1.0;
                dPhij1j2_2p5_jesup=-1.0; dPhiHj1j2_2p5_jesup=-1.0;
                
                if (njets_pt30_eta4p7_jesup > 0) {
                    Jet1_jesup.SetPtEtaPhiM((*jet_jesup_pt)[jet1index_jesup],(*jet_jesup_eta)[jet1index_jesup],(*jet_jesup_phi)[jet1index_jesup],(*jet_jesup_mass)[jet1index_jesup]);
                    pTj1_jesup=Jet1_jesup.Pt(); 
                    etaj1_jesup=Jet1_jesup.Eta();
                    yj1_jesup=Jet1_jesup.Rapidity();
                    dPhiHj1_jesup=TMath::Abs(higgs_undec.Phi()-Jet1_jesup.Phi());
                    dyHj1_jesup=TMath::Abs(y4l-yj1_jesup);
                }
                if (njets_pt30_eta4p7_jesup > 1) {
                    Jet2_jesup.SetPtEtaPhiM((*jet_jesup_pt)[jet2index_jesup],(*jet_jesup_eta)[jet2index_jesup],(*jet_jesup_phi)[jet2index_jesup],(*jet_jesup_mass)[jet2index_jesup]);
                    pTj2_jesup=Jet2_jesup.Pt();
                    etaj2_jesup=Jet2_jesup.Eta();
                    yj2_jesup=Jet2_jesup.Rapidity();
                    dEtaj1j2_jesup=TMath::Abs(Jet1_jesup.Eta()-Jet2_jesup.Eta());
                    dPhij1j2_jesup=TMath::Abs(Jet1_jesup.Phi()-Jet2_jesup.Phi());
                    dPhiHj1j2_jesup=TMath::Abs(higgs_undec.Phi()-(Jet1_jesup+Jet2_jesup).Phi());
                }
                if (njets_pt30_eta2p5_jesup > 0) {
                    Jet1_2p5_jesup.SetPtEtaPhiM((*jet_jesup_pt)[jet1index2p5_jesup],(*jet_jesup_eta)[jet1index2p5_jesup],(*jet_jesup_phi)[jet1index2p5_jesup],(*jet_jesup_mass)[jet1index2p5_jesup]);
                    pTj1_2p5_jesup=Jet1_2p5_jesup.Pt();
                    yj1_2p5_jesup=Jet1_2p5.Rapidity();
                    dPhiHj1_2p5_jesup=TMath::Abs(higgs_undec.Phi()-Jet1_2p5_jesup.Phi());
                    dyHj1_2p5_jesup=TMath::Abs(y4l-yj1_2p5_jesup);
                }
                if (njets_pt30_eta2p5_jesup > 1) {
                    Jet2_2p5_jesup.SetPtEtaPhiM((*jet_jesup_pt)[jet2index2p5_jesup],(*jet_jesup_eta)[jet2index2p5_jesup],(*jet_jesup_phi)[jet2index2p5_jesup],(*jet_jesup_mass)[jet2index2p5_jesup]);
                    pTj2_2p5_jesup=Jet2_2p5_jesup.Pt();
                    yj2_2p5_jesup=Jet2_2p5_jesup.Rapidity();
                    mj1j2_2p5_jesup=(Jet1_2p5_jesup+Jet2_2p5_jesup).M();                    
                    dEtaj1j2_2p5_jesup=TMath::Abs(Jet1_2p5_jesup.Eta()-Jet2_2p5_jesup.Eta());
                    dPhij1j2_2p5_jesup=TMath::Abs(Jet1_2p5_jesup.Phi()-Jet2_2p5_jesup.Phi());
                    dPhiHj1j2_2p5_jesup=TMath::Abs(higgs_undec.Phi()-(Jet1_2p5_jesup+Jet2_2p5_jesup).Phi());
                }

                // JES dn
                pTj1_jesdn=-1.0; pTj2_jesdn=-1.0;                
                pTj1_2p5_jesdn=-1.0; pTj2_jesdn=-1.0;
                mj1j2_jesdn=-1.0; dEtaj1j2_jesdn=-1.0;

                yj1_jesdn=-1.0; yj2_jesdn=-1.0;
                dPhiHj1_jesdn=-1.0; dyHj1_jesdn=-1.0;
                mj1j2_jesdn=-1.0; dEtaj1j2_jesdn=-1.0;
                dPhij1j2_jesdn=-1.0; dPhiHj1j2_jesdn=-1.0;
                yj1_2p5_jesdn=-1.0; yj2_2p5_jesdn=-1.0;
                
                dPhiHj1_2p5_jesdn=-1.0; dyHj1_2p5_jesdn=-1.0;
                mj1j2_2p5_jesdn=-1.0; dEtaj1j2_2p5_jesdn=-1.0;
                dPhij1j2_2p5_jesdn=-1.0; dPhiHj1j2_2p5_jesdn=-1.0;
                
                if (njets_pt30_eta4p7_jesdn > 0) {
                    Jet1_jesdn.SetPtEtaPhiM((*jet_jesdn_pt)[jet1index_jesdn],(*jet_jesdn_eta)[jet1index_jesdn],(*jet_jesdn_phi)[jet1index_jesdn],(*jet_jesdn_mass)[jet1index_jesdn]);
                    pTj1_jesdn=Jet1_jesdn.Pt(); 
                    etaj1_jesdn=Jet1_jesdn.Eta();
                    yj1_jesdn=Jet1_jesdn.Rapidity();
                    dPhiHj1_jesdn=TMath::Abs(higgs_undec.Phi()-Jet1_jesdn.Phi());
                    dyHj1_jesdn=TMath::Abs(y4l-yj1_jesdn);
                }
                if (njets_pt30_eta4p7_jesdn > 1) {
                    Jet2_jesdn.SetPtEtaPhiM((*jet_jesdn_pt)[jet2index_jesdn],(*jet_jesdn_eta)[jet2index_jesdn],(*jet_jesdn_phi)[jet2index_jesdn],(*jet_jesdn_mass)[jet2index_jesdn]);
                    pTj2_jesdn=Jet2_jesdn.Pt();
                    etaj2_jesdn=Jet2_jesdn.Eta();
                    yj2_jesdn=TMath::Abs(Jet2_jesdn.Rapidity());
                    dEtaj1j2_jesdn=TMath::Abs(Jet1_jesdn.Eta()-Jet2_jesdn.Eta());
                    dPhij1j2_jesdn=TMath::Abs(Jet1_jesdn.Phi()-Jet2_jesdn.Phi());
                    dPhiHj1j2_jesdn=TMath::Abs(higgs_undec.Phi()-(Jet1_jesdn+Jet2_jesdn).Phi());
                }
                if (njets_pt30_eta2p5_jesdn > 0) {
                    Jet1_2p5_jesdn.SetPtEtaPhiM((*jet_jesdn_pt)[jet1index2p5_jesdn],(*jet_jesdn_eta)[jet1index2p5_jesdn],(*jet_jesdn_phi)[jet1index2p5_jesdn],(*jet_jesdn_mass)[jet1index2p5_jesdn]);
                    pTj1_2p5_jesdn=Jet1_2p5_jesdn.Pt();
                    yj1_2p5_jesdn=TMath::Abs(Jet1_2p5.Rapidity());
                    dPhiHj1_2p5_jesdn=TMath::Abs(higgs_undec.Phi()-Jet1_2p5_jesdn.Phi());
                    dyHj1_2p5_jesdn=TMath::Abs(y4l-yj1_2p5_jesdn);
                }
                if (njets_pt30_eta2p5_jesdn > 1) {
                    Jet2_2p5_jesdn.SetPtEtaPhiM((*jet_jesdn_pt)[jet2index2p5_jesdn],(*jet_jesdn_eta)[jet2index2p5_jesdn],(*jet_jesdn_phi)[jet2index2p5_jesdn],(*jet_jesdn_mass)[jet2index2p5_jesdn]);
                    pTj2_2p5_jesdn=Jet2_2p5_jesdn.Pt();
                    yj2_2p5_jesdn=TMath::Abs(Jet2_2p5_jesdn.Rapidity());
                    mj1j2_2p5_jesdn=(Jet1_2p5_jesdn+Jet2_2p5_jesdn).M();                    
                    dEtaj1j2_2p5_jesdn=TMath::Abs(Jet1_2p5_jesdn.Eta()-Jet2_2p5_jesdn.Eta());
                    dPhij1j2_2p5_jesdn=TMath::Abs(Jet1_2p5_jesdn.Phi()-Jet2_2p5_jesdn.Phi());
                    dPhiHj1j2_2p5_jesdn=TMath::Abs(higgs_undec.Phi()-(Jet1_2p5_jesdn+Jet2_2p5_jesdn).Phi());
                }

                
                pTj1_VBF=-1.0; dPhij1j2_VBF=-1.0; dPhiHj1j2_VBF=-1.0;
                if (njets_pt30_eta4p7>1 && dEtaj1j2>3.5 && mj1j2>200.0) {
                    pTj1_VBF=pTj1; dPhij1j2_VBF=dPhij1j2; dPhiHj1j2_VBF=dPhiHj1j2;
                }
                pTj1_VBF_jesup=-1.0; dPhij1j2_VBF_jesup=-1.0; dPhiHj1j2_VBF_jesup=-1.0;
                if (njets_pt30_eta4p7_jesup>1 && dEtaj1j2_jesup>3.5 && mj1j2_jesup>200.0) {
                    pTj1_VBF_jesup=pTj1_jesup; dPhij1j2_VBF_jesup=dPhij1j2_jesup; dPhiHj1j2_VBF_jesup=dPhiHj1j2_jesup;
                }
                pTj1_VBF_jesdn=-1.0; dPhij1j2_VBF_jesdn=-1.0; dPhiHj1j2_VBF_jesdn=-1.0;
                if (njets_pt30_eta4p7_jesdn>1 && dEtaj1j2_jesdn>3.5 && mj1j2_jesdn>200.0) {
                    pTj1_VBF_jesdn=pTj1_jesdn; dPhij1j2_VBF_jesdn=dPhij1j2_jesdn; dPhiHj1j2_VBF_jesdn=dPhiHj1j2_jesdn;
                }

                
                
                mela->setInputEvent(&daughters, &associated, 0, 0);
                mela->setCurrentCandidateFromIndex(0);
                
                float me_0plus_JHU=0.0;
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                mela->computeP(me_0plus_JHU, true);
                
                float me_qqZZ_MCFM=0.0;
                mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
                mela->computeP(me_qqZZ_MCFM, true);
                
                float p0plus_m4l=0.0;
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
                mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);
                

                float p0minus_VAJHU=0.0;
                mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
                mela->computeP(p0minus_VAJHU, true);

                float pg1g4_VAJHU=0.0;
                mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
                (mela->selfDHggcoupl)[0][0][0]=1.;
                (mela->selfDHzzcoupl)[0][0][0]=1.;
                (mela->selfDHzzcoupl)[0][3][0]=1.;
                mela->computeP(pg1g4_VAJHU, true);
                pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;

                float bkg_m4l=0.0;
                mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
                mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
                

                D_bkg_kin = me_0plus_JHU/(me_0plus_JHU+me_qqZZ_MCFM*getDbkgkinConstant(idL1*idL2*idL3*idL3,mass4l));
                D_bkg = me_0plus_JHU*p0plus_m4l/(me_0plus_JHU*p0plus_m4l+me_qqZZ_MCFM*bkg_m4l*getDbkgConstant(idL1*idL2*idL3*idL4,mass4l)); 
                D_g4 = me_0plus_JHU/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); // D_0-
                D_g1g4 = pg1g4_VAJHU*2.521/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); // D_CP, 2.521 since g1=1 and g4=1 is used

                pvbf_VAJHU=-9999.0; phjj_VAJHU=-9999.0; phj_VAJHU=-9999.0;
                pwh_hadronic_VAJHU=-9999.0; pzh_hadronic_VAJHU=-9999.0;
                pAux_vbf_VAJHU=-9999.0;
                if (njets_pt30_eta4p7>=2){
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                    mela->computeProdP(pvbf_VAJHU, true);
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
                    mela->computeProdP(phjj_VAJHU, true);
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                    mela->computeProdP(pwh_hadronic_VAJHU, true);
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                    mela->computeProdP(pzh_hadronic_VAJHU, true);
                    D_VBF = pvbf_VAJHU/(pvbf_VAJHU+phjj_VAJHU*getDVBF2jetsConstant(mass4l) ); // VBF(2j) vs. gg->H+2j
                    D_HadWH = pwh_hadronic_VAJHU/(pwh_hadronic_VAJHU+1e-3*phjj_VAJHU ); // W(->2j)H vs. gg->H+2j
                    D_HadZH = pzh_hadronic_VAJHU/(pzh_hadronic_VAJHU+1e-4*phjj_VAJHU ); // Z(->2j)H vs. gg->H+2j
                } else {
                    D_VBF = -1.0; D_HadWH = -1.0; D_HadZH = -1.0;
                }
                
                if (njets_pt30_eta4p7==1) {
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
                    mela->computeProdP(phj_VAJHU, true);
                    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
                    mela->computeProdP(pvbf_VAJHU, true); // Un-integrated ME
                    mela->getPAux(pAux_vbf_VAJHU); // = Integrated / un-integrated
                    D_VBF1j = pvbf_VAJHU*pAux_vbf_VAJHU/(pvbf_VAJHU*pAux_vbf_VAJHU+phj_VAJHU*getDVBF1jetConstant(mass4l)); // VBF(1j) vs. gg->H+1j
                } else {
                    D_VBF1j = -1.0;
                }

                if (njets_pt30_eta4p7>=2) {
                    pvbf_VAJHU=TMath::Log(pvbf_VAJHU); 
                    phjj_VAJHU=TMath::Log(phjj_VAJHU); 
                    pwh_hadronic_VAJHU=TMath::Log(pwh_hadronic_VAJHU); 
                    pzh_hadronic_VAJHU=TMath::Log(pzh_hadronic_VAJHU);
                }
                if (njets_pt30_eta4p7==1) {
                    phj_VAJHU=TMath::Log(phj_VAJHU);
                    pAux_vbf_VAJHU=TMath::Log(pAux_vbf_VAJHU);
                }
                
                if (njets_pt30_eta4p7>=2) {
                    float jetqgl0 =(*jet_QGTagger)[jet1index];
                    float jetqgl1 =(*jet_QGTagger)[jet2index];
                    if(jetqgl0<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                        TRandom3 rand;
                        rand.SetSeed(abs(static_cast<int>(sin((*jet_phi)[jet1index])*100000)));
                        jetqgl0 = rand.Uniform();
                    }
                    if(jetqgl1<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                        TRandom3 rand;
                        rand.SetSeed(abs(static_cast<int>(sin((*jet_phi)[jet2index])*100000)));
                        jetqgl1 = rand.Uniform();
                    }
                    float jetPgOverPq0 = 1./jetqgl0- 1.;
                    float jetPgOverPq1 = 1./jetqgl1- 1.;
                    D_VBF_QG = 1./(1.+ (1./D_VBF - 1.) * pow(jetPgOverPq0*jetPgOverPq1, 1./3.));
                    D_HadWH_QG = 1./(1.+ (1./D_HadWH - 1.) * pow(jetPgOverPq0*jetPgOverPq1, 1./3.));
                    D_HadZH_QG = 1./(1. + (1./D_HadZH - 1.) * pow(jetPgOverPq0*jetPgOverPq1,1./3.));
                } else {
                    D_VBF_QG = -1.0; D_HadWH_QG = -1.0; D_HadZH_QG = -1.0;
                }
                
                if (njets_pt30_eta4p7==1) {
                    float jetqgl0 =(*jet_QGTagger)[jet1index];
                    if(jetqgl0<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
                        TRandom3 rand;
                        rand.SetSeed(abs(static_cast<int>(sin((*jet_phi)[jet1index])*100000)));
                        jetqgl0 = rand.Uniform();
                    }
                    float jetPgOverPq0 = 1./jetqgl0- 1.;
                    D_VBF1j_QG = 1/(1+ (1./D_VBF1j - 1.) * pow(jetPgOverPq0,1./3.));
                } else {
                    D_VBF1j_QG = -1.0;
                }
                
                mela->resetInputEvent(); 
                

                nisoleptons=0;
                worstiso=9999.0;

                sumplus=0; summinus=0;
                sumflavor=0;
                int noverlapping=0;
                for(unsigned int i = 0; i < (*lep_pt).size(); i++) {

                    bool goodi=false;
                    if ((int)i==lep_Hindex[0] || (int)i==lep_Hindex[1] || (int)i==lep_Hindex[2] || (int)i==lep_Hindex[3]) {
                        goodi=true;
                    } else {
                        if ( (abs((*lep_id)[i])==11 && (*lep_tightId)[i]==1 && (*lep_RelIsoNoFSR)[i]<isoCutEl) || 
                             (abs((*lep_id)[i])==13 && (*lep_tightId)[i]==1 && (*lep_RelIsoNoFSR)[i]<isoCutMu) ) {
                            goodi=true;
                        }
                    }
                    if (!goodi) continue;
                                                
                    nisoleptons++;  

                    if ((*lep_RelIsoNoFSR)[i]<worstiso) worstiso=(*lep_RelIsoNoFSR)[i];

                    if ((*lep_id)[i]>0) sumplus++;
                    if ((*lep_id)[i]<0) summinus++;
                    sumflavor+=(*lep_id)[i];

                    TLorentzVector lepi;
                    lepi.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);

                    float minDr=9999.0;
                    for(unsigned int j = 0; j < (*lep_pt).size(); j++) {
                        
                        if (i==j) continue;

                        bool goodj=false;
                        if ((int)j==lep_Hindex[0] || (int)j==lep_Hindex[1] || (int)j==lep_Hindex[2] || (int)j==lep_Hindex[3]) {
                            goodj=true;
                        } else {
                            if ( (abs((*lep_id)[j])==11 && (*lep_tightId)[j]==1 && (*lep_RelIsoNoFSR)[j]<isoCutEl) || 
                                 (abs((*lep_id)[j])==13 && (*lep_tightId)[j]==1 && (*lep_RelIsoNoFSR)[j]<isoCutMu) ) {
                                goodj=true;
                            }
                        }
                        if (!goodj) continue;
                           
                        TLorentzVector lepj;
                        lepj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
                        float thisdR = lepi.DeltaR(lepj);
                        if (thisdR<minDr) minDr=thisdR;
                    }
                    if (minDr<0.02) noverlapping+=1;
                }
                nisoleptons-=(noverlapping/2);

                // Event Categories 
                if (nisoleptons==4 && (((njets_pt30_eta4p7==2||njets_pt30_eta4p7==3)&&nbjets_pt30_eta4p7<2)||(njets_pt30_eta4p7>=4&&nbjets_pt30_eta4p7==0)) && D_VBF>(1.043-460./(mass4l+634.))) {EventCat=2;}
                else if (nisoleptons==4 && ((njets_pt30_eta4p7==2||njets_pt30_eta4p7==3)||(njets_pt30_eta4p7>=4&&nbjets_pt30_eta4p7==0)) && (D_HadWH>0.951 || D_HadZH>0.9937)) {EventCat=4;}
                else if (njets_pt30_eta4p7<=3 && nbjets_pt30_eta4p7==0 && (nisoleptons==5 || (nisoleptons>=6&&sumplus>=3&&summinus>=3&&sumflavor==0))) {EventCat=3;}
                else if (njets_pt30_eta4p7==0 && nisoleptons>=5) {EventCat=3;}
                else if (njets_pt30_eta4p7>=4 && nbjets_pt30_eta4p7>0) {EventCat=5;}
                else if (nisoleptons>=5) {EventCat=5;}
                else if (nisoleptons==4 && njets_pt30_eta4p7<=1 && met>100.0) {EventCat=6;}
                else if (nisoleptons==4 && njets_pt30_eta4p7==1 && D_VBF1j>0.697) {EventCat=1;}
                else {EventCat=0;}


            }
            
            if (redoDifferentialObs) {
                
            }
*/
            //cout<<"idL1: "<<idL1<<" passedFullSelection "<<passedFullSelection<<" passedZ4lSelection: "<<passedZ4lSelection<<endl;                

            if(debug) cout<<"fill tree"<<endl;
            if(debug) cout<<endl;            
            if (!redoGEN) newtree->Fill();
            
        }

/*
        if (!isData && redoGEN) {

            GENmass4l=-1.0; GENpT4l=-1.0; GENy4l=-9999.0; 
            GENmassZ1=-1.0; GENmassZ2=-1.0;
            GENmass4e=-1.0; GENmass4mu=-1.0; GENmass2e2mu=-1.0;

            GENcosTheta1=-9999.0; GENcosTheta2=-9999.0;
            GENcosThetaStar=-9999.0; GENPhi=-9999.0; GENPhi1=-9999.0;

            GENnjets_pt30_eta4p7=0; GENnjets_pt30_eta2p5=0;

            GENpTj1=-1.0; GENyj1=-1.0;
            GENpTj2=-1.0; GENyj2=-1.0; 
            
            GENdPhiHj1=-1.0; GENdyHj1=-1.0; 
            GENmj1j2=-1.0; GENdEtaj1j2=-1.0; 
            GENdPhij1j2=-1.0; GENdPhiHj1j2=-1.0; 
            
            GENpTj1_2p5=-1.0; GENyj1_2p5=-1.0;
            GENpTj2_2p5=-1.0; GENyj2_2p5=-1.0;
            
            GENdPhiHj1_2p5=-1.0; GENdyHj1_2p5=-1.0; 
            GENmj1j2_2p5=-1.0; GENdEtaj1j2_2p5=-1.0; 
            GENdPhij1j2_2p5=-1.0; GENdPhiHj1j2_2p5=-1.0;

            GENfinalState=0;

            if ( (*GENZ_MomId).size()>2 ) {
                if ( (((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[1]==25 && (*GENZ_DaughtersId)[0]==11 && (*GENZ_DaughtersId)[1]==11) || ((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[2]==25 && (*GENZ_DaughtersId)[0]==11 && (*GENZ_DaughtersId)[2]==11) || ((*GENZ_MomId)[1]==25 && (*GENZ_MomId)[2]==25 && (*GENZ_DaughtersId)[1]==11 && (*GENZ_DaughtersId)[2]==11)) ) GENfinalState=1;
                else if ( (((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[1]==25 && (*GENZ_DaughtersId)[0]==13 && (*GENZ_DaughtersId)[1]==13) || ((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[2]==25 && (*GENZ_DaughtersId)[0]==13 && (*GENZ_DaughtersId)[2]==13) || ((*GENZ_MomId)[1]==25 && (*GENZ_MomId)[2]==25 && (*GENZ_DaughtersId)[1]==13 && (*GENZ_DaughtersId)[2]==13)) ) GENfinalState=2;                
                else if ( (((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[1]==25 && ((*GENZ_DaughtersId)[0]==11 || (*GENZ_DaughtersId)[0]==13) && ((*GENZ_DaughtersId)[1]==11 || (*GENZ_DaughtersId)[1]==13)) || ((*GENZ_MomId)[0]==25 && (*GENZ_MomId)[2]==25 && ((*GENZ_DaughtersId)[0]==11 || (*GENZ_DaughtersId)[2]==13) && ((*GENZ_DaughtersId)[2]==11 || (*GENZ_DaughtersId)[2]==13)) || ((*GENZ_MomId)[1]==25 && (*GENZ_MomId)[2]==25 && ((*GENZ_DaughtersId)[1]==11 || (*GENZ_DaughtersId)[1]==13) && ((*GENZ_DaughtersId)[2]==11 || (*GENZ_DaughtersId)[2]==13))) ) GENfinalState=3;
            } else {
                if ((abs((*GENZ_DaughtersId)[0])==11 && abs((*GENZ_DaughtersId)[1])==11) ) GENfinalState=1;
                else if ( ((*GENZ_DaughtersId)[0]==13 && (*GENZ_DaughtersId)[1]==13) ) GENfinalState=2;
                else if ( (((*GENZ_DaughtersId)[0]==11 && (*GENZ_DaughtersId)[1]==13) || ((*GENZ_DaughtersId)[0]==13 && (*GENZ_DaughtersId)[1]==11)) ) GENfinalState=3;
            }

            if (GENfinalState==0) {
                int nel=0; int nmu=0; int ntau=0;
                for (int l=0; l<(int)(*GENlep_id).size(); l++) {
                    if (abs((*GENlep_id)[l])==11) nel++;
                    if (abs((*GENlep_id)[l])==13) nmu++;
                    if (abs((*GENlep_id)[l])==15) ntau++;
                }
                if (nel>=4 && ntau==0) GENfinalState=1;
                if (nmu>=4 && ntau==0) GENfinalState=2;
                if (nel>=2 && nmu>=2 && ntau==0) GENfinalState=3;
            }

            int nzl=2; int nzq=0; int nznu=0;
            if ((*GENZ_DaughtersId).size()>2) {
                nzl=0;
                for (int z=0; z<(int)(*GENZ_DaughtersId).size(); z++) {
                    if ((*GENZ_DaughtersId)[z]==11) nzl+=1;
                    if ((*GENZ_DaughtersId)[z]==13) nzl+=1;
                    if ((*GENZ_DaughtersId)[z]==15) nzl+=1;
                    
                    if ((*GENZ_DaughtersId)[z]==12) nznu+=1;
                    if ((*GENZ_DaughtersId)[z]==14) nznu+=1;
                    if ((*GENZ_DaughtersId)[z]==16) nznu+=1;
                    
                    if ((*GENZ_DaughtersId)[z]<10) nzq+=1;
                }
            }

            
            GENEventCat=-2;
            if (stage0cat==10 || stage0cat==20 || stage0cat==22 || stage0cat==30 || stage0cat==40 || stage0cat==50 || stage0cat==60 ) GENEventCat=-1;  
            else if (stage0cat==11 || (stage0cat==21 && njets_pt30_eta4p7==0) ) GENEventCat=0; // ggH
            else if (stage0cat==21 && njets_pt30_eta4p7>=2) GENEventCat=2; // VBF 2jet reco
            else if (stage0cat==21 && njets_pt30_eta4p7==1) GENEventCat=1; // VBF 1jet reco
            else if (stage0cat==23 || (nzq==1) ) GENEventCat=4; // VHhad
            else if (stage0cat==31 || (stage0cat==41 && nzl==3) ) GENEventCat=3; // VHlep
            else if (stage0cat==41 && nznu==1) GENEventCat=6; // ZHmet
            else if (stage0cat==61) GENEventCat=5; // ttH

            if (passedFiducialSelection) {

                TLorentzVector Lep1, Lep2, Lep3, Lep4,  GENJet1, GENJet2, GENJet1_2p5, GENJet2_2p5;
                TLorentzVector nullFourVector(0, 0, 0, 0);
                Lep1.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[0]],(*GENlep_eta)[GENlep_Hindex[0]],(*GENlep_phi)[GENlep_Hindex[0]],(*GENlep_mass)[GENlep_Hindex[0]]);
                Lep2.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[1]],(*GENlep_eta)[GENlep_Hindex[1]],(*GENlep_phi)[GENlep_Hindex[1]],(*GENlep_mass)[GENlep_Hindex[1]]);
                Lep3.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[2]],(*GENlep_eta)[GENlep_Hindex[2]],(*GENlep_phi)[GENlep_Hindex[2]],(*GENlep_mass)[GENlep_Hindex[2]]);
                Lep4.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[3]],(*GENlep_eta)[GENlep_Hindex[3]],(*GENlep_phi)[GENlep_Hindex[3]],(*GENlep_mass)[GENlep_Hindex[3]]);

                TLorentzVector GENH = Lep1+Lep2+Lep3+Lep4;
                GENmassZ1 = (Lep1+Lep2).M(); GENmassZ2 = (Lep3+Lep4).M();
                GENmass4l = GENH.M(); GENpT4l=GENH.Pt(); GENy4l = GENH.Rapidity();                 

                int GENidL1 = (*GENlep_id)[GENlep_Hindex[0]];
                int GENidL2 = (*GENlep_id)[GENlep_Hindex[1]];
                int GENidL3 = (*GENlep_id)[GENlep_Hindex[2]];
                int GENidL4 = (*GENlep_id)[GENlep_Hindex[3]];
                
                if (abs(GENidL1)==11 and abs(GENidL3)==11) GENmass4e=GENmass4l;
                if (abs(GENidL1)==13 and abs(GENidL3)==13) GENmass4mu=GENmass4l;
                if (abs(GENidL1)!=abs(GENidL3)) GENmass2e2mu=GENmass4l;
                
                GENisH4l = ( ((*GENlep_MomMomId)[GENlep_Hindex[0]]==25 && (*GENlep_MomMomId)[GENlep_Hindex[1]]==25 && (*GENlep_MomMomId)[GENlep_Hindex[2]]==25 && (*GENlep_MomMomId)[GENlep_Hindex[3]]==25) );


//                TUtil::computeAngles(Lep1, GENidL1, Lep2, GENidL2, Lep3, GENidL3, Lep4, GENidL4, \
//                                     GENcosThetaStar,GENcosTheta1,GENcosTheta2,GENPhi,GENPhi1);


                int GENjet1index=0; int GENjet2index=0; int GENjet1index_2p5=0; int GENjet2index_2p5=0;
            

                for (unsigned int k=0; k<(*GENjet_pt).size(); k++) {
                    
                    TLorentzVector thisGENJet;
                    thisGENJet.SetPtEtaPhiM((*GENjet_pt)[k],(*GENjet_eta)[k],(*GENjet_phi)[k],(*GENjet_mass)[k]);
                    
                    if ((*GENjet_pt)[k]<30.0 || abs((*GENjet_eta)[k])>4.7) continue;
                     
                    GENnjets_pt30_eta4p7+=1;
                    
                    if (thisGENJet.Pt()>GENpTj1) {
                        GENpTj2=GENpTj1; GENjet2index=GENjet1index;
                        GENpTj1=thisGENJet.Pt(); GENjet1index=k;
                    } else if (thisGENJet.Pt()>GENpTj2) {
                        GENpTj2=thisGENJet.Pt(); GENjet2index=k;
                    }
                    if (abs((*GENjet_eta)[k])<2.5) {
                        GENnjets_pt30_eta2p5+=1;
                        if (thisGENJet.Pt()>GENpTj1_2p5) {
                            GENpTj2_2p5=GENpTj1_2p5; GENjet2index_2p5=GENjet1index_2p5;
                            GENpTj1_2p5=thisGENJet.Pt(); GENjet1index_2p5=k;
                        } else if (thisGENJet.Pt()>GENpTj2_2p5) {
                            GENpTj2_2p5=thisGENJet.Pt(); GENjet2index_2p5=k;
                        }
                        
                    }

                }

                if (GENnjets_pt30_eta4p7 > 0) {
                    GENJet1.SetPtEtaPhiM((*GENjet_pt)[GENjet1index],(*GENjet_eta)[GENjet1index],(*GENjet_phi)[GENjet1index],(*GENjet_mass)[GENjet1index]);
                    GENpTj1=GENJet1.Pt();
                    GENyj1=TMath::Abs(GENJet1.Rapidity());
                    GENdPhiHj1=TMath::Abs(GENH.Phi()-GENJet1.Phi());
                    GENdyHj1=TMath::Abs(GENy4l-GENyj1);
                }
                if (GENnjets_pt30_eta4p7 > 1) {
                    GENJet2.SetPtEtaPhiM((*GENjet_pt)[GENjet2index],(*GENjet_eta)[GENjet2index],(*GENjet_phi)[GENjet2index],(*GENjet_mass)[GENjet2index]);
                    GENpTj2=GENJet2.Pt();
                    GENyj2=TMath::Abs(GENJet2.Rapidity());
                    GENmj1j2=(GENJet1+GENJet2).M();
                    GENdEtaj1j2=TMath::Abs(GENJet1.Eta()-GENJet2.Eta());
                    GENdPhij1j2=TMath::Abs(GENJet1.Phi()-GENJet2.Phi());
                    GENdPhiHj1j2=TMath::Abs(GENH.Phi()-(GENJet1+GENJet2).Phi());
                }
                if (GENnjets_pt30_eta2p5 > 0) {
                    GENJet1_2p5.SetPtEtaPhiM((*GENjet_pt)[GENjet1index_2p5],(*GENjet_eta)[GENjet1index_2p5],(*GENjet_phi)[GENjet1index_2p5],(*GENjet_mass)[GENjet1index_2p5]);
                    GENpTj1_2p5=GENJet1_2p5.Pt();
                    GENyj1_2p5=TMath::Abs(GENJet1.Rapidity());
                    GENdPhiHj1_2p5=TMath::Abs(GENH.Phi()-GENJet1_2p5.Phi());
                    GENdyHj1_2p5=TMath::Abs(GENy4l-GENyj1_2p5);
                }
                if (GENnjets_pt30_eta2p5 > 1) {
                    GENJet2_2p5.SetPtEtaPhiM((*GENjet_pt)[GENjet2index_2p5],(*GENjet_eta)[GENjet2index_2p5],(*GENjet_phi)[GENjet2index_2p5],(*GENjet_mass)[GENjet2index_2p5]);
                    GENpTj2_2p5=GENJet2_2p5.Pt();
                    GENyj2_2p5=TMath::Abs(GENJet2.Rapidity());
                    GENmj1j2_2p5=(GENJet1_2p5+GENJet2_2p5).M();
                    GENdEtaj1j2_2p5=TMath::Abs(GENJet1_2p5.Eta()-GENJet2_2p5.Eta());
                    GENdPhij1j2_2p5=TMath::Abs(GENJet1_2p5.Phi()-GENJet2_2p5.Phi());
                    GENdPhiHj1j2_2p5=TMath::Abs(GENH.Phi()-(GENJet1_2p5+GENJet2_2p5).Phi());
                }
                
                GENpTj1_VBF=-1.0; GENdPhij1j2_VBF=-1.0; GENdPhiHj1j2_VBF=-1.0;
                if (GENnjets_pt30_eta4p7>1 && GENdEtaj1j2>3.5 && GENmj1j2>200.0) {
                    GENpTj1_VBF=GENpTj1; GENdPhij1j2_VBF=GENdPhij1j2; GENdPhiHj1j2_VBF=GENdPhiHj1j2;
                }
                
                
            }


        }  
        if (redoGEN) newtree->Fill();
*/
        
       
    }
    cout<<"sumweight: "<<sumweight<<endl;
    cout<<"npass: "<<npass<<endl;
    
}

void SetNewTree(TTree* newtree){

/*    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    newtree->Branch("nVtx",&nVtx,"nVtx/I");
    newtree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    newtree->Branch("passedFullSelection",&passedFullSelection,"passedFullSelection/O");
    newtree->Branch("passedFiducialSelection",&passedFiducialSelection,"passedFiducialSelection/O");
    newtree->Branch("passedZ4lSelection",&passedZ4lSelection,"passedZ4lSelection/O");
    newtree->Branch("passedZXCRSelection",&passedZXCRSelection,"passedZXCRSelection/O");
    newtree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons,"nZXCRFailedLeptons/I");
    newtree->Branch("finalState",&finalState,"finalState/I");    
    newtree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F"); */
/*    newtree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    newtree->Branch("genWeight",&genWeight,"genWeight/F");
    newtree->Branch("nnloWeights",&nnloWeights);
    newtree->Branch("qcdWeights",&qcdWeights);
    newtree->Branch("pdfENVup",&pdfENVup,"pdfENVup/F");
    newtree->Branch("pdfENVdown",&pdfENVdown,"pdfENVdown/F");
    newtree->Branch("pdfRMSup",&pdfRMSup,"pdfRMSup/F");
    newtree->Branch("crossSection",&crossSection,"crossSection/F");
    newtree->Branch("sumW",&sumW,"sumW/F");
    newtree->Branch("k_qqZZ_qcd_M",&k_qqZZ_qcd_M,"k_qqZZ_qcd_M/F");
    newtree->Branch("k_qqZZ_ewk",&k_qqZZ_ewk,"k_qqZZ_ewk/F");
    newtree->Branch("k_ggZZ",&k_ggZZ,"k_ggZZ/F");

    newtree->Branch("pTL1",&pTL1,"pTL1/F");
    newtree->Branch("pTL2",&pTL2,"pTL2/F");
    newtree->Branch("pTL3",&pTL3,"pTL3/F");
    newtree->Branch("pTL4",&pTL4,"pTL4/F");
    newtree->Branch("idL1",&idL1,"idL1/I");
    newtree->Branch("idL2",&idL2,"idL2/I");
    newtree->Branch("idL3",&idL3,"idL3/I");
    newtree->Branch("idL4",&idL4,"idL4/I");
    newtree->Branch("etaL1",&etaL1,"etaL1/F");
    newtree->Branch("etaL2",&etaL2,"etaL2/F");
    newtree->Branch("etaL3",&etaL3,"etaL3/F");
    newtree->Branch("etaL4",&etaL4,"etaL4/F");

    newtree->Branch("isH4l",&isH4l,"isH4l/O");
    newtree->Branch("mass4l",&mass4l,"mass4l/F");
    newtree->Branch("mass4lErr",&mass4lErr,"mass4lErr/F");
    newtree->Branch("mass4lREFIT",&mass4lREFIT,"mass4lREFIT/F");
    newtree->Branch("mass4lErrREFIT",&mass4lErrREFIT,"mass4lErrREFIT/F");
    newtree->Branch("mass4mu",&mass4mu,"mass4mu/F");
    newtree->Branch("mass4e",&mass4e,"mass4e/F");
    newtree->Branch("mass2e2mu",&mass2e2mu,"mass2e2mu/F");
    newtree->Branch("pT4l",&pT4l,"pT4l/F");
    newtree->Branch("y4l",&y4l,"y4l/F");
    newtree->Branch("cosTheta1",&cosTheta1,"cosTheta1/F");
    newtree->Branch("cosTheta2",&cosTheta2,"cosTheta2/F");
    newtree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    newtree->Branch("Phi",&Phi,"Phi/F");
    newtree->Branch("Phi1",&Phi1,"Phi1/F");
    newtree->Branch("massZ1",&massZ1,"massZ1/F");
    newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
    newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
    newtree->Branch("njets_pt30_eta4p7_jesup",&njets_pt30_eta4p7_jesup,"njets_pt30_eta4p7_jesup/I");
    newtree->Branch("njets_pt30_eta4p7_jesdn",&njets_pt30_eta4p7_jesdn,"njets_pt30_eta4p7_jesdn/I");
    newtree->Branch("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7,"nbjets_pt30_eta4p7/I");
    newtree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5,"njets_pt30_eta2p5/I");
    newtree->Branch("njets_pt30_eta2p5_jesup",&njets_pt30_eta2p5_jesup,"njets_pt30_eta2p5_jesup/I");
    newtree->Branch("njets_pt30_eta2p5_jesdn",&njets_pt30_eta2p5_jesdn,"njets_pt30_eta2p5_jesdn/I");
    newtree->Branch("nisoleptons",&nisoleptons,"nisoleptons/I");
    newtree->Branch("worstiso",&worstiso,"worstiso/F");
    newtree->Branch("sumplus",&sumplus,"sumplus/I");
    newtree->Branch("summinus",&summinus,"summinus/I");
    newtree->Branch("sumflavor",&sumflavor,"sumflavor/I");
    newtree->Branch("met",&met,"met/F"); 
    newtree->Branch("met_phi",&met_phi,"met_phi/F"); 
*/

    newtree->Branch("isH4l",&isH4l,"isH4l/O");
    // Nominal

    //newtree->Branch("jes_unc_split_Total", jes_unc_split[0]);

    newtree->Branch("jes_unc_split_Total",&jes_unc_split_Total,"jes_unc_split_Total/F");
    newtree->Branch("jes_unc_split_Abs",&jes_unc_split_Abs,"jes_unc_split_Abs/F");
    newtree->Branch("jes_unc_split_Abs_year",&jes_unc_split_Abs_year,"jes_unc_split_Abs_year/F");
    newtree->Branch("jes_unc_split_BBEC1",&jes_unc_split_BBEC1,"jes_unc_split_BBEC1/F");
    newtree->Branch("jes_unc_split_BBEC1_year",&jes_unc_split_BBEC1_year,"jes_unc_split_BBEC1_year/F");
    newtree->Branch("jes_unc_split_EC2", &jes_unc_split_EC2,"jes_unc_split_EC2/F");
    newtree->Branch("jes_unc_split_EC2_year", &jes_unc_split_EC2_year,"jes_unc_split_EC2_year/F");
    newtree->Branch("jes_unc_split_FlavQCD", &jes_unc_split_FlavQCD,"jes_unc_split_FlavQCD/F");
    newtree->Branch("jes_unc_split_HF",&jes_unc_split_HF,"jes_unc_split_HF/F"); 
    newtree->Branch("jes_unc_split_HF_year", &jes_unc_split_HF_year,"jes_unc_split_HF_year/F");
    newtree->Branch("jes_unc_split_RelBal", &jes_unc_split_RelBal, "jes_unc_split_RelBal/F");
    newtree->Branch("jes_unc_split_RelSample_year", &jes_unc_split_RelSample_year, "jes_unc_split_RelSample_year/F");
    // JES splitted up/dn
    newtree->Branch("pt_jesup_split_Total",&pt_jesup_split_Total,"pt_jesup_split_Total/F");
    newtree->Branch("pt_jesup_split_Abs",&pt_jesup_split_Abs,"pt_jesup_split_Abs/F");
    newtree->Branch("pt_jesup_split_Abs_year",&pt_jesup_split_Abs_year,"pt_jesup_split_Abs_year/F");
    newtree->Branch("pt_jesup_split_BBEC1",&pt_jesup_split_BBEC1,"pt_jesup_split_BBEC1/F");
    newtree->Branch("pt_jesup_split_BBEC1_year",&pt_jesup_split_BBEC1_year,"pt_jesup_split_BBEC1_year/F");
    newtree->Branch("pt_jesup_split_EC2", &pt_jesup_split_EC2,"pt_jesup_split_EC2/F");
    newtree->Branch("pt_jesup_split_EC2_year", &pt_jesup_split_EC2_year,"pt_jesup_split_EC2_year/F");
    newtree->Branch("pt_jesup_split_FlavQCD", &pt_jesup_split_FlavQCD,"pt_jesup_split_FlavQCD/F");
    newtree->Branch("pt_jesup_split_HF",&pt_jesup_split_HF,"pt_jesup_split_HF/F"); 
    newtree->Branch("pt_jesup_split_HF_year", &pt_jesup_split_HF_year,"pt_jesup_split_HF_year/F");
    newtree->Branch("pt_jesup_split_RelBal", &pt_jesup_split_RelBal, "pt_jesup_split_RelBal/F");
    newtree->Branch("pt_jesup_split_RelSample_year", &pt_jesup_split_RelSample_year, "pt_jesup_split_RelSample_year/F");

    newtree->Branch("pt_jesdn_split_Total",&pt_jesdn_split_Total,"pt_jesdn_split_Total/F");
    newtree->Branch("pt_jesdn_split_Abs",&pt_jesdn_split_Abs,"pt_jesdn_split_Abs/F");
    newtree->Branch("pt_jesdn_split_Abs_year",&pt_jesdn_split_Abs_year,"pt_jesdn_split_Abs_year/F");
    newtree->Branch("pt_jesdn_split_BBEC1",&pt_jesdn_split_BBEC1,"pt_jesdn_split_BBEC1/F");
    newtree->Branch("pt_jesdn_split_BBEC1_year",&pt_jesdn_split_BBEC1_year,"pt_jesdn_split_BBEC1_year/F");
    newtree->Branch("pt_jesdn_split_EC2", &pt_jesdn_split_EC2,"pt_jesdn_split_EC2/F");
    newtree->Branch("pt_jesdn_split_EC2_year", &pt_jesdn_split_EC2_year,"pt_jesdn_split_EC2_year/F");
    newtree->Branch("pt_jesdn_split_FlavQCD", &pt_jesdn_split_FlavQCD,"pt_jesdn_split_FlavQCD/F");
    newtree->Branch("pt_jesdn_split_HF",&pt_jesdn_split_HF,"pt_jesdn_split_HF/F"); 
    newtree->Branch("pt_jesdn_split_HF_year", &pt_jesdn_split_HF_year,"pt_jesdn_split_HF_year/F");
    newtree->Branch("pt_jesdn_split_RelBal", &pt_jesdn_split_RelBal, "pt_jesdn_split_RelBal/F");
    newtree->Branch("pt_jesdn_split_RelSample_year", &pt_jesdn_split_RelSample_year, "pt_jesdn_split_RelSample_year/F");

    // jet observables
/*
    newtree->Branch("pTj1",&pTj1,"pTj1/F");
    newtree->Branch("pTj1_VBF",&pTj1_VBF,"pTj1_VBF/F");
    newtree->Branch("etaj1",&etaj1,"etaj1/F");
    newtree->Branch("qgj1",&qgj1,"qgj1/F");
    newtree->Branch("pTj2",&pTj2,"pTj2/F");
    newtree->Branch("etaj2",&etaj2,"etaj2/F");
    newtree->Branch("qgj2",&qgj2,"qgj1/F");
    newtree->Branch("csvj1",&csvj1,"csvj1/F");
    newtree->Branch("csvj2",&csvj2,"csvj2/F");

    newtree->Branch("yj1",&yj1,"yj1/F");
    newtree->Branch("yj2",&yj2,"yj2/F");

    newtree->Branch("dPhiHj1",&dPhiHj1,"dPhiHj1/F"); 
    newtree->Branch("dyHj1",&dyHj1,"dyHj1/F");
    newtree->Branch("mj1j2",&mj1j2,"mj1j2/F"); 
    newtree->Branch("dEtaj1j2",&dEtaj1j2,"dEtaj1j2/F");
    newtree->Branch("dPhij1j2",&dPhij1j2,"dPhij1j2/F"); 
    newtree->Branch("dPhij1j2_VBF",&dPhij1j2_VBF,"dPhij1j2_VBF/F"); 
    newtree->Branch("dPhiHj1j2",&dPhiHj1j2,"dPhiHj1j2/F");
    newtree->Branch("dPhiHj1j2_VBF",&dPhiHj1j2_VBF,"dPhiHj1j2_VBF/F"); 

    newtree->Branch("pTj1_2p5",&pTj1_2p5,"pTj1_2p5/F"); 
    newtree->Branch("yj1_2p5",&yj1_2p5,"yj1_2p5/F");
    newtree->Branch("pTj2_2p5",&pTj2_2p5,"pTj2_2p5/F"); 
    newtree->Branch("yj2_2p5",&yj2_2p5,"yj2_2p5/F");

    newtree->Branch("dPhiHj1_2p5",&dPhiHj1_2p5,"dPhiHj1_2p5/F"); 
    newtree->Branch("dyHj1_2p5",&dyHj1_2p5,"dyHj1_2p5/F");
    newtree->Branch("mj1j2_2p5",&mj1j2_2p5,"mj1j2_2p5/F"); 
    newtree->Branch("dEtaj1j2_2p5",&dEtaj1j2_2p5,"dEtaj1j2_2p5/F");
    newtree->Branch("dPhij1j2_2p5",&dPhij1j2_2p5,"dPhij1j2_2p5/F"); 
    newtree->Branch("dPhiHj1j2_2p5",&dPhiHj1j2_2p5,"dPhiHj1j2_2p5/F");
    // JES up
    newtree->Branch("pTj1_jesup",&pTj1_jesup,"pTj1_jesup/F");
    newtree->Branch("pTj1_VBF_jesup",&pTj1_VBF_jesup,"pTj1_VBF_jesup/F");
    newtree->Branch("etaj1_jesup",&etaj1_jesup,"etaj1_jesup/F");
    newtree->Branch("pTj2_jesup",&pTj2_jesup,"pTj2_jesup/F");
    newtree->Branch("etaj2_jesup",&etaj2_jesup,"etaj2_jesup/F");

    newtree->Branch("yj1_jesup",&yj1_jesup,"yj1_jesup/F");
    newtree->Branch("yj2_jesup",&yj2_jesup,"yj2_jesup/F");

    newtree->Branch("dPhiHj1_jesup",&dPhiHj1_jesup,"dPhiHj1_jesup/F"); 
    newtree->Branch("dyHj1_jesup",&dyHj1_jesup,"dyHj1_jesup/F");
    newtree->Branch("mj1j2_jesup",&mj1j2_jesup,"mj1j2_jesup/F"); 
    newtree->Branch("dEtaj1j2_jesup",&dEtaj1j2_jesup,"dEtaj1j2_jesup/F");
    newtree->Branch("dPhij1j2_jesup",&dPhij1j2_jesup,"dPhij1j2_jesup/F"); 
    newtree->Branch("dPhiHj1j2_jesup",&dPhiHj1j2_jesup,"dPhiHj1j2_jesup/F");
    newtree->Branch("dPhij1j2_VBF_jesup",&dPhij1j2_VBF_jesup,"dPhij1j2_VBF_jesup/F"); 
    newtree->Branch("dPhiHj1j2_VBF_jesup",&dPhiHj1j2_VBF_jesup,"dPhiHj1j2_VBF_jesup/F");

    newtree->Branch("pTj1_2p5_jesup",&pTj1_2p5_jesup,"pTj1_2p5_jesup/F"); 
    newtree->Branch("yj1_2p5_jesup",&yj1_2p5_jesup,"yj1_2p5_jesup/F");
    newtree->Branch("pTj2_2p5_jesup",&pTj2_2p5_jesup,"pTj2_2p5_jesup/F"); 
    newtree->Branch("yj2_2p5_jesup",&yj2_2p5_jesup,"yj2_2p5_jesup/F");

    newtree->Branch("dPhiHj1_2p5_jesup",&dPhiHj1_2p5_jesup,"dPhiHj1_2p5_jesup/F"); 
    newtree->Branch("dyHj1_2p5_jesup",&dyHj1_2p5_jesup,"dyHj1_2p5_jesup/F");
    newtree->Branch("mj1j2_2p5_jesup",&mj1j2_2p5_jesup,"mj1j2_2p5_jesup/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup",&dEtaj1j2_2p5_jesup,"dEtaj1j2_2p5_jesup/F");
    newtree->Branch("dPhij1j2_2p5_jesup",&dPhij1j2_2p5_jesup,"dPhij1j2_2p5_jesup/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup",&dPhiHj1j2_2p5_jesup,"dPhiHj1j2_2p5_jesup/F");
    // JES dn
    newtree->Branch("pTj1_jesdn",&pTj1_jesdn,"pTj1_jesdn/F");
    newtree->Branch("pTj1_VBF_jesdn",&pTj1_VBF_jesdn,"pTj1_VBF_jesdn/F");
    newtree->Branch("etaj1_jesdn",&etaj1_jesdn,"etaj1_jesdn/F");
    newtree->Branch("pTj2_jesdn",&pTj2_jesdn,"pTj2_jesdn/F");
    newtree->Branch("etaj2_jesdn",&etaj2_jesdn,"etaj2_jesdn/F");

    newtree->Branch("yj1_jesdn",&yj1_jesdn,"yj1_jesdn/F");
    newtree->Branch("yj2_jesdn",&yj2_jesdn,"yj2_jesdn/F");

    newtree->Branch("dPhiHj1_jesdn",&dPhiHj1_jesdn,"dPhiHj1_jesdn/F"); 
    newtree->Branch("dyHj1_jesdn",&dyHj1_jesdn,"dyHj1_jesdn/F");
    newtree->Branch("mj1j2_jesdn",&mj1j2_jesdn,"mj1j2_jesdn/F"); 
    newtree->Branch("dEtaj1j2_jesdn",&dEtaj1j2_jesdn,"dEtaj1j2_jesdn/F");
    newtree->Branch("dPhij1j2_jesdn",&dPhij1j2_jesdn,"dPhij1j2_jesdn/F"); 
    newtree->Branch("dPhiHj1j2_jesdn",&dPhiHj1j2_jesdn,"dPhiHj1j2_jesdn/F");
    newtree->Branch("dPhij1j2_VBF_jesdn",&dPhij1j2_VBF_jesdn,"dPhij1j2_VBF_jesdn/F"); 
    newtree->Branch("dPhiHj1j2_VBF_jesdn",&dPhiHj1j2_VBF_jesdn,"dPhiHj1j2_VBF_jesdn/F");

    newtree->Branch("pTj1_2p5_jesdn",&pTj1_2p5_jesdn,"pTj1_2p5_jesdn/F"); 
    newtree->Branch("yj1_2p5_jesdn",&yj1_2p5_jesdn,"yj1_2p5_jesdn/F");
    newtree->Branch("pTj2_2p5_jesdn",&pTj2_2p5_jesdn,"pTj2_2p5_jesdn/F"); 
    newtree->Branch("yj2_2p5_jesdn",&yj2_2p5_jesdn,"yj2_2p5_jesdn/F");

    newtree->Branch("dPhiHj1_2p5_jesdn",&dPhiHj1_2p5_jesdn,"dPhiHj1_2p5_jesdn/F"); 
    newtree->Branch("dyHj1_2p5_jesdn",&dyHj1_2p5_jesdn,"dyHj1_2p5_jesdn/F");
    newtree->Branch("mj1j2_2p5_jesdn",&mj1j2_2p5_jesdn,"mj1j2_2p5_jesdn/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn",&dEtaj1j2_2p5_jesdn,"dEtaj1j2_2p5_jesdn/F");
    newtree->Branch("dPhij1j2_2p5_jesdn",&dPhij1j2_2p5_jesdn,"dPhij1j2_2p5_jesdn/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn",&dPhiHj1j2_2p5_jesdn,"dPhiHj1j2_2p5_jesdn/F");

    newtree->Branch("D_bkg_kin", &D_bkg_kin, "D_bkg_kin/F");
    newtree->Branch("D_bkg", &D_bkg, "D_bkg/F");
    newtree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/F");
    newtree->Branch("D_g4", &D_g4, "D_g4/F");
    newtree->Branch("D_VBF", &D_VBF, "D_VBF/F");
    newtree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/F");
    newtree->Branch("D_HadWH",&D_HadWH,"D_HadWH/F");
    newtree->Branch("D_HadZH",&D_HadZH,"D_HadZH/F");
    newtree->Branch("D_VBF_QG",&D_VBF_QG,"D_VBF_QG/F");
    newtree->Branch("D_VBF1j_QG",&D_VBF1j_QG,"D_VBF1j_QG/F");
    newtree->Branch("D_HadWH_QG",&D_HadWH_QG,"D_HadWH_QG/F");
    newtree->Branch("D_HadZH_QG",&D_HadZH_QG,"D_HadZH_QG/F");
    newtree->Branch("EventCat",&EventCat,"EventCat/I");
    newtree->Branch("GENEventCat",&GENEventCat,"GENEventCat/I");
    newtree->Branch("stage1cat",&stage1cat,"stage1cat/I");
    newtree->Branch("phj_VAJHU",&phj_VAJHU,"phj_VAJHU/F");
    newtree->Branch("phjj_VAJHU",&phjj_VAJHU,"phjj_VAJHU/F");
    newtree->Branch("pvbf_VAJHU",&pvbf_VAJHU,"pvbf_VAJHU/F");
    newtree->Branch("pAux_vbf_VAJHU",&pAux_vbf_VAJHU,"pAux_vbf_VAJHU/F");
    newtree->Branch("pwh_hadronic_VAJHU",&pwh_hadronic_VAJHU,"pwh_hadronic_VAJHU/F");
    newtree->Branch("pzh_hadronic_VAJHU",&pzh_hadronic_VAJHU,"pzh_hadronic_VAJHU/F");

    newtree->Branch("GENisH4l",&GENisH4l,"GENisH4l/O");
    newtree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");
    newtree->Branch("GENmass4l",&GENmass4l,"GENmass4l/F");
    newtree->Branch("GENmass4e",&GENmass4e,"GENmass4e/F");
    newtree->Branch("GENmass4mu",&GENmass4mu,"GENmass4mu/F");
    newtree->Branch("GENmass2e2mu",&GENmass2e2mu,"GENmass2e2mu/F");
    newtree->Branch("GENpT4l",&GENpT4l,"GENpT4l/F");
    newtree->Branch("GENy4l",&GENy4l,"GENy4l/F");
    newtree->Branch("GENmassZ1",&GENmassZ1,"GENmassZ1/F");
    newtree->Branch("GENmassZ2",&GENmassZ2,"GENmassZ2/F"); 

    newtree->Branch("GENcosTheta1",&GENcosTheta1,"GENcosTheta1/F");
    newtree->Branch("GENcosTheta2",&GENcosTheta2,"GENcosTheta2/F");
    newtree->Branch("GENcosThetaStar",&GENcosThetaStar,"GENcosThetaStar/F");
    newtree->Branch("GENPhi",&GENPhi,"GENPhi/F");
    newtree->Branch("GENPhi1",&GENPhi1,"GENPhi1/F");

    newtree->Branch("GENnjets_pt30_eta4p7",&GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/I"); 
    newtree->Branch("GENnjets_pt30_eta2p5",&GENnjets_pt30_eta2p5,"GENnjets_pt30_eta2p5/I");

    newtree->Branch("GENpTj1",&GENpTj1,"GENpTj1/F"); 
    newtree->Branch("GENpTj1_VBF",&GENpTj1_VBF,"GENpTj1_VBF/F"); 
    newtree->Branch("GENyj1",&GENyj1,"GENyj1/F");
    newtree->Branch("GENpTj2",&GENpTj2,"GENpTj2/F"); 
    newtree->Branch("GENyj2",&GENyj2,"GENyj2/F");

    newtree->Branch("GENdPhiHj1",&GENdPhiHj1,"GENdPhiHj1/F"); 
    newtree->Branch("GENdyHj1",&GENdyHj1,"GENdyHj1/F");
    newtree->Branch("GENmj1j2",&GENmj1j2,"GENmj1j2/F"); 
    newtree->Branch("GENdEtaj1j2",&GENdEtaj1j2,"GENdEtaj1j2/F");
    newtree->Branch("GENdPhij1j2",&GENdPhij1j2,"GENdPhij1j2/F"); 
    newtree->Branch("GENdPhij1j2_VBF",&GENdPhij1j2_VBF,"GENdPhij1j2_VBF/F"); 
    newtree->Branch("GENdPhiHj1j2",&GENdPhiHj1j2,"GENdPhiHj1j2/F");
    newtree->Branch("GENdPhiHj1j2_VBF",&GENdPhiHj1j2_VBF,"GENdPhiHj1j2_VBF/F");

    newtree->Branch("GENpTj1_2p5",&GENpTj1_2p5,"GENpTj1_2p5/F"); 
    newtree->Branch("GENyj1_2p5",&GENyj1_2p5,"GENyj1_2p5/F");
    newtree->Branch("GENpTj2_2p5",&GENpTj2_2p5,"GENpTj2_2p5/F"); 
    newtree->Branch("GENyj2_2p5",&GENyj2_2p5,"GENyj2_2p5/F");

    newtree->Branch("GENdPhiHj1_2p5",&GENdPhiHj1_2p5,"GENdPhiHj1_2p5/F"); 
    newtree->Branch("GENdyHj1_2p5",&GENdyHj1_2p5,"GENdyHj1_2p5/F");
    newtree->Branch("GENmj1j2_2p5",&GENmj1j2_2p5,"GENmj1j2_2p5/F"); 
    newtree->Branch("GENdEtaj1j2_2p5",&GENdEtaj1j2_2p5,"GENdEtaj1j2_2p5/F");
    newtree->Branch("GENdPhij1j2_2p5",&GENdPhij1j2_2p5,"GENdPhij1j2_2p5/F"); 
    newtree->Branch("GENdPhiHj1j2_2p5",&GENdPhiHj1j2_2p5,"GENdPhiHj1j2_2p5/F");
*/
}

