#ifndef ZZ4LAnalysisTree_jes_h
#define ZZ4LAnalysisTree_jes_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TPaveLabel.h"
#include "TLegend.h"

#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TRandom3.h"

#include <algorithm>

#include "TColor.h"

using namespace std;  

////
// output tree
  
bool passedFullSelection, passedZ4lSelection;
bool passedZXCRSelection, passedZ4lZXCRSelection;
bool passedFiducialSelection;
int nZXCRFailedLeptons;
int finalState;
bool isH4l;

bool passTrig;
float pTL1, etaL1;
float pTL2, etaL2;
float pTL3, etaL3;
float pTL4, etaL4;

int idL1, idL2, idL3, idL4;

float mass4l, mass4lErr;
float mass4lREFIT, mass4lErrREFIT;
float massZ1REFIT, massZ2REFIT;
float mass4mu, mass4e, mass2e2mu;

float pT4l;

float massZ1, massZ2;

float y4l;
float cosTheta1; float cosTheta2;
float cosThetaStar; float Phi; float Phi1;

// Nominal
int njets_pt30_eta4p7;
int njets_pt30_eta2p5;
float pTj1, pTj1_VBF, etaj1;
float pTj2, etaj2;
float qgj1, qgj2;
float pTj1_2p5, pTj2_2p5;
float csvj1, csvj2;

float yj1; float yj2;

float dPhiHj1; float dyHj1;
float mj1j2; float dEtaj1j2;
float dPhij1j2; float dPhiHj1j2;
float dPhij1j2_VBF; float dPhiHj1j2_VBF;

float yj1_2p5; float yj2_2p5;

float dPhiHj1_2p5; float dyHj1_2p5;
float mj1j2_2p5; float dEtaj1j2_2p5;
float dPhij1j2_2p5; float dPhiHj1j2_2p5;
// JES up
int njets_pt30_eta4p7_jesup;
int njets_pt30_eta2p5_jesup;
float pTj1_jesup, pTj1_VBF_jesup, etaj1_jesup;
float pTj2_jesup, etaj2_jesup;
float pTj1_2p5_jesup, pTj2_2p5_jesup;

float yj1_jesup; float yj2_jesup;

float dPhiHj1_jesup; float dyHj1_jesup;
float mj1j2_jesup; float dEtaj1j2_jesup;
float dPhij1j2_jesup; float dPhiHj1j2_jesup;
float dPhij1j2_VBF_jesup; float dPhiHj1j2_VBF_jesup;

float yj1_2p5_jesup; float yj2_2p5_jesup;

float dPhiHj1_2p5_jesup; float dyHj1_2p5_jesup;
float mj1j2_2p5_jesup; float dEtaj1j2_2p5_jesup;
float dPhij1j2_2p5_jesup; float dPhiHj1j2_2p5_jesup;
// JES dn
int njets_pt30_eta4p7_jesdn;
int njets_pt30_eta2p5_jesdn;
float pTj1_jesdn, pTj1_VBF_jesdn, etaj1_jesdn;
float pTj2_jesdn, etaj2_jesdn;
float pTj1_2p5_jesdn, pTj2_2p5_jesdn;

float yj1_jesdn; float yj2_jesdn;

float dPhiHj1_jesdn; float dyHj1_jesdn;
float mj1j2_jesdn; float dEtaj1j2_jesdn;
float dPhij1j2_jesdn; float dPhiHj1j2_jesdn;
float dPhij1j2_VBF_jesdn; float dPhiHj1j2_VBF_jesdn;

float yj1_2p5_jesdn; float yj2_2p5_jesdn;

float dPhiHj1_2p5_jesdn; float dyHj1_2p5_jesdn;
float mj1j2_2p5_jesdn; float dEtaj1j2_2p5_jesdn;
float dPhij1j2_2p5_jesdn; float dPhiHj1j2_2p5_jesdn;

float D_bkg_kin;
float D_bkg;
float Dgg10_VAMCFM;
float D_g4;
float D_g1g4;
float D_VBF;
float D_VBF1j;
float D_HadWH;
float D_HadZH;
float D_VBF_QG;
float D_VBF1j_QG;
float D_HadWH_QG;
float D_HadZH_QG;
float phj_VAJHU;
float phjj_VAJHU;
float pvbf_VAJHU;
float pAux_vbf_VAJHU;
float pwh_hadronic_VAJHU;
float pzh_hadronic_VAJHU;

int EventCat;
int GENEventCat;
int stage1cat;
int GENfinalState;
int stage0cat;
int nisoleptons, nbjets_pt30_eta4p7;
float worstiso;
float met, met_phi;
int sumplus; int summinus;
int sumflavor;
 
float GENmassZZ;
float GENmass4l; float GENpT4l; float GENy4l;
float GENmassZ1; float GENmassZ2;
float GENmass4e, GENmass4mu, GENmass2e2mu;

float GENcosTheta1; float GENcosTheta2;
float GENcosThetaStar; float GENPhi; float GENPhi1;

int GENnjets_pt30_eta4p7;
int GENnjets_pt30_eta2p5;

float GENpTj1; float GENpTj1_VBF; float GENyj1;
float GENpTj2; float GENyj2;

float GENdPhiHj1; float GENdyHj1;
float GENmj1j2; float GENdEtaj1j2;
float GENdPhij1j2; float GENdPhiHj1j2;
float GENdPhij1j2_VBF; float GENdPhiHj1j2_VBF;

float GENpTj1_2p5; float GENyj1_2p5;
float GENpTj2_2p5; float GENyj2_2p5;

float GENdPhiHj1_2p5; float GENdyHj1_2p5;
float GENmj1j2_2p5; float GENdEtaj1j2_2p5;
float GENdPhij1j2_2p5; float GENdPhiHj1j2_2p5;
//   JES variables
/*vector<float> jes_unc_split {};
vector<float> pt_jesup_split {};
vector<float> pt_jesdn_split {};

float singleContr_jes_unc = 0;
*/
float jes_unc_split_Total; 
float jes_unc_split_Abs; 
float jes_unc_split_Abs_year; 
float jes_unc_split_BBEC1; 
float jes_unc_split_BBEC1_year;
float jes_unc_split_EC2; 
float jes_unc_split_EC2_year; 
float jes_unc_split_FlavQCD; 
float jes_unc_split_HF; 
float jes_unc_split_HF_year;
float jes_unc_split_RelBal; 
float jes_unc_split_RelSample_year; 
// up
float pt_jesup_split_Total; 
float pt_jesup_split_Abs; 
float pt_jesup_split_Abs_year; 
float pt_jesup_split_BBEC1; 
float pt_jesup_split_BBEC1_year;
float pt_jesup_split_EC2; 
float pt_jesup_split_EC2_year; 
float pt_jesup_split_FlavQCD; 
float pt_jesup_split_HF; 
float pt_jesup_split_HF_year;
float pt_jesup_split_RelBal; 
float pt_jesup_split_RelSample_year; 
//dn
float pt_jesdn_split_Total; 
float pt_jesdn_split_Abs; 
float pt_jesdn_split_Abs_year; 
float pt_jesdn_split_BBEC1; 
float pt_jesdn_split_BBEC1_year;
float pt_jesdn_split_EC2; 
float pt_jesdn_split_EC2_year; 
float pt_jesdn_split_FlavQCD; 
float pt_jesdn_split_HF; 
float pt_jesdn_split_HF_year;
float pt_jesdn_split_RelBal; 
float pt_jesdn_split_RelSample_year; 

// input tree variables
std::string *triggersPassed;
ULong64_t Run, LumiSect, Event;
bool passedTrig;

float dataMCWeight, genWeight, pileupWeight, crossSection, sumweight;
float k_qqZZ_qcd_M,k_qqZZ_ewk,k_ggZZ;
float sumW;
int nVtx, nInt; //nPV

std::vector<float>* lep_mass;
std::vector<float> *lep_pt; std::vector<float> *lep_eta; std::vector<float> *lep_phi;
std::vector<float>* lepFSR_mass;
std::vector<float> *lepFSR_pt; std::vector<float> *lepFSR_eta; std::vector<float> *lepFSR_phi;

std::vector<int> *lep_tightId;
std::vector<int> *lep_ecalDriven;
std::vector<int> *lep_id;
int lep_Hindex[4];
std::vector<int> *lep_genindex;
std::vector<float> *lep_RelIso;
std::vector<float> *lep_RelIsoNoFSR;
std::vector<float> *lep_pterr;
std::vector<float> *lep_dataMC;

std::vector<float> *jet_mass;
std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
std::vector<int> *jet_iscleanH4l; std::vector<float> *jet_QGTagger; std::vector<float> *jet_csvv2;

std::vector<float> *jet_jesup_mass; std::vector<float> *jet_jesup_pt; std::vector<float> *jet_jesup_eta; std::vector<float> *jet_jesup_phi;
std::vector<int> *jet_jesup_iscleanH4l; 

std::vector<float> *jet_jesdn_mass; std::vector<float> *jet_jesdn_pt; std::vector<float> *jet_jesdn_eta; std::vector<float> *jet_jesdn_phi;
std::vector<int> *jet_jesdn_iscleanH4l; 

std::vector<int> *fsrPhotons_lepindex;
std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
std::vector<float> *fsrPhotons_pterr;

std::vector<float>* nnloWeights;
std::vector<float>* qcdWeights;
float pdfENVup, pdfENVdown, pdfRMSup;

std::vector<float>* GENlep_mass;
std::vector<float> *GENlep_pt; std::vector<float> *GENlep_eta; std::vector<float> *GENlep_phi;
std::vector<int> *GENlep_id;
std::vector<int> *GENlep_MomId;
std::vector<int> *GENlep_MomMomId;
int GENlep_Hindex[4];
bool GENisH4l;

std::vector<float> *GENjet_mass;
std::vector<float> *GENjet_pt; std::vector<float> *GENjet_eta; std::vector<float> *GENjet_phi;

std::vector<int> *GENZ_DaughtersId;
std::vector<int> *GENZ_MomId;

//namespace ZZ4LAnalysisTree {
namespace ZZ4LAnalysisTree_jes {

    void setAddresses(TTree* tree, TString filename){
        
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("Run",1);
        tree->SetBranchStatus("LumiSect",1);
        tree->SetBranchStatus("Event",1);
        tree->SetBranchStatus("nVtx",1);
        tree->SetBranchStatus("nInt",1);
        tree->SetBranchStatus("genWeight",1);
        tree->SetBranchStatus("crossSection",1);
        tree->SetBranchStatus("pileupWeight",1);
        tree->SetBranchStatus("passedTrig",1);
        tree->SetBranchStatus("triggersPassed",1);
        tree->SetBranchStatus("lep_id",1);
        tree->SetBranchStatus("lep_genindex",1);
        tree->SetBranchStatus("lep_tightId",1);
        tree->SetBranchStatus("lep_ecalDriven",1);
        tree->SetBranchStatus("lep_pt",1);
        tree->SetBranchStatus("lep_pterr",1);
        tree->SetBranchStatus("lep_eta",1);
        tree->SetBranchStatus("lep_phi",1);
        tree->SetBranchStatus("lep_mass",1);
        tree->SetBranchStatus("lep_RelIso",1);
        tree->SetBranchStatus("lep_RelIsoNoFSR",1);
        tree->SetBranchStatus("lepFSR_pt",1);
        tree->SetBranchStatus("lepFSR_eta",1);
        tree->SetBranchStatus("lepFSR_phi",1);
        tree->SetBranchStatus("lepFSR_mass",1);        
        tree->SetBranchStatus("jet_pt",1);
        tree->SetBranchStatus("jet_eta",1);
        tree->SetBranchStatus("jet_phi",1);
        tree->SetBranchStatus("jet_mass",1);
        tree->SetBranchStatus("jet_iscleanH4l",1);
        tree->SetBranchStatus("jet_QGTagger",1);
        tree->SetBranchStatus("jet_csvv2",1);
        tree->SetBranchStatus("jet_jesup_pt",1);
        tree->SetBranchStatus("jet_jesup_eta",1);
        tree->SetBranchStatus("jet_jesup_phi",1);
        tree->SetBranchStatus("jet_jesup_mass",1);
        tree->SetBranchStatus("jet_jesup_iscleanH4l",1);
        tree->SetBranchStatus("jet_jesdn_pt",1);
        tree->SetBranchStatus("jet_jesdn_eta",1);
        tree->SetBranchStatus("jet_jesdn_phi",1);
        tree->SetBranchStatus("jet_jesdn_mass",1);
        tree->SetBranchStatus("jet_jesdn_iscleanH4l",1);
        tree->SetBranchStatus("Phi",1);
        tree->SetBranchStatus("GENPhi",1);

        tree->SetBranchStatus("pTj1",1);
        tree->SetBranchStatus("pTj1_jesup",1);
        tree->SetBranchStatus("pTj1_jesdn",1);
        tree->SetBranchStatus("pTj2",1);
        tree->SetBranchStatus("pTj2_jesup",1);
        tree->SetBranchStatus("pTj2_jesdn",1);
        tree->SetBranchStatus("mj1j2",1);
        tree->SetBranchStatus("mj1j2_jesup",1);
        tree->SetBranchStatus("mj1j2_jesdn",1);
        tree->SetBranchStatus("dEtaj1j2",1);
        tree->SetBranchStatus("dEtaj1j2_jesup",1);
        tree->SetBranchStatus("dEtaj1j2_jesdn",1);
        tree->SetBranchStatus("dPhij1j2",1);
        tree->SetBranchStatus("dPhij1j2_jesup",1);
        tree->SetBranchStatus("dPhij1j2_jesdn",1);
        tree->SetBranchStatus("pT4lj",1);
        tree->SetBranchStatus("pT4lj_jesup",1);
        tree->SetBranchStatus("pT4lj_jesdn",1);
        tree->SetBranchStatus("mass4lj",1);
        tree->SetBranchStatus("mass4lj_jesup",1);
        tree->SetBranchStatus("mass4lj_jesdn",1);
        tree->SetBranchStatus("pT4ljj",1);
        tree->SetBranchStatus("pT4ljj_jesup",1);
        tree->SetBranchStatus("pT4ljj_jesdn",1);
        tree->SetBranchStatus("mass4ljj",1);
        tree->SetBranchStatus("mass4ljj_jesup",1);
        tree->SetBranchStatus("mass4ljj_jesdn",1);
//
        tree->SetBranchStatus("pTj2_2p5",1);
        tree->SetBranchStatus("pTj2_2p5_jesup",1);
        tree->SetBranchStatus("pTj2_2p5_jesdn",1);
        tree->SetBranchStatus("mj1j2_2p5",1);
        tree->SetBranchStatus("mj1j2_2p5_jesup",1);
        tree->SetBranchStatus("mj1j2_2p5_jesdn",1);
        tree->SetBranchStatus("dEtaj1j2_2p5",1);
        tree->SetBranchStatus("dEtaj1j2_2p5_jesup",1);
        tree->SetBranchStatus("dEtaj1j2_2p5_jesdn",1);
        tree->SetBranchStatus("dPhij1j2_2p5",1);
        tree->SetBranchStatus("dPhij1j2_2p5_jesup",1);
        tree->SetBranchStatus("dPhij1j2_2p5_jesdn",1);
        tree->SetBranchStatus("pTj1_2p5",1);
        tree->SetBranchStatus("etaj1_2p5",1);
        tree->SetBranchStatus("phij1_2p5",1);
        tree->SetBranchStatus("mj1_2p5",1);
        //tree->SetBranchStatus("pTj2_2p5",1);
        tree->SetBranchStatus("etaj2_2p5",1);
        tree->SetBranchStatus("phij2_2p5",1);
        tree->SetBranchStatus("mj2_2p5",1);
        tree->SetBranchStatus("pTj1_2p5_jesup",1);
        tree->SetBranchStatus("etaj1_2p5_jesup",1);
        tree->SetBranchStatus("phij1_2p5_jesup",1);
        tree->SetBranchStatus("mj1_2p5_jesup",1);
        //tree->SetBranchStatus("pTj2_2p5_jesup",1);
        tree->SetBranchStatus("etaj2_2p5_jesup",1);
        tree->SetBranchStatus("phij2_2p5_jesup",1);
        tree->SetBranchStatus("mj2_2p5_jesup",1);
        tree->SetBranchStatus("pTj1_2p5_jesdn",1);
        tree->SetBranchStatus("etaj1_2p5_jesdn",1);
        tree->SetBranchStatus("phij1_2p5_jesdn",1);
        tree->SetBranchStatus("mj1_2p5_jesdn",1);
        //tree->SetBranchStatus("pTj2_2p5_jesdn",1);
        tree->SetBranchStatus("etaj2_2p5_jesdn",1);
        tree->SetBranchStatus("phij2_2p5_jesdn",1);
        tree->SetBranchStatus("mj2_2p5_jesdn",1);
// GEN
        tree->SetBranchStatus("GENpTj1",1);
        tree->SetBranchStatus("GENpTj2",1);
        tree->SetBranchStatus("GENmj1j2",1);
        tree->SetBranchStatus("GENdEtaj1j2",1);
        tree->SetBranchStatus("GENdPhij1j2",1);
        tree->SetBranchStatus("GENpT4lj",1);
        tree->SetBranchStatus("GENmass4lj",1);
        tree->SetBranchStatus("GENpT4ljj",1);
        tree->SetBranchStatus("GENmass4ljj",1);
//
        tree->SetBranchStatus("GENpTj2_2p5",1);
        tree->SetBranchStatus("GENmj1j2_2p5",1);
        tree->SetBranchStatus("GENdEtaj1j2_2p5",1);
        tree->SetBranchStatus("GENdPhij1j2_2p5",1);
        tree->SetBranchStatus("GENpTj1_2p5",1);
        tree->SetBranchStatus("GENetaj1_2p5",1);
        tree->SetBranchStatus("GENphij1_2p5",1);
        tree->SetBranchStatus("GENmj1_2p5",1);
        //tree->SetBranchStatus("pTj2_2p5",1);
        tree->SetBranchStatus("GENetaj2_2p5",1);
        tree->SetBranchStatus("GENphij2_2p5",1);
        tree->SetBranchStatus("GENmj2_2p5",1);


        tree->SetBranchStatus("fsrPhotons_pt",1);
        tree->SetBranchStatus("fsrPhotons_pterr",1);
        tree->SetBranchStatus("fsrPhotons_eta",1);
        tree->SetBranchStatus("fsrPhotons_phi",1);
        tree->SetBranchStatus("fsrPhotons_lepindex",1);
        tree->SetBranchStatus("met",1);
        tree->SetBranchStatus("met_phi",1);
        tree->SetBranchStatus("GENZ_DaughtersId",1);
        tree->SetBranchStatus("GENZ_MomId",1);
        tree->SetBranchStatus("stage0cat",1);
        tree->SetBranchStatus("stage1cat",1);
        tree->SetBranchStatus("passedFiducialSelection",1);
        tree->SetBranchStatus("nnloWeights",1);
        tree->SetBranchStatus("qcdWeights",1);
        tree->SetBranchStatus("pdfENVup",1);
        tree->SetBranchStatus("pdfENVdown",1);        
        tree->SetBranchStatus("pdfRMSup",1);        

        tree->SetBranchStatus("GENlep_Hindex", 1);
        tree->SetBranchStatus("GENlep_id",1);
        tree->SetBranchStatus("GENlep_MomId",1);
        tree->SetBranchStatus("GENlep_MomMomId",1);
        tree->SetBranchStatus("GENlep_pt",1);
        tree->SetBranchStatus("GENlep_eta",1);
        tree->SetBranchStatus("GENlep_phi",1);
        tree->SetBranchStatus("GENlep_mass",1);

        tree->SetBranchStatus("GENjet_pt",1);
        tree->SetBranchStatus("GENjet_eta",1);
        tree->SetBranchStatus("GENjet_phi",1);
        tree->SetBranchStatus("GENjet_mass",1);

        tree->SetBranchStatus("GENmassZZ",1);

        tree->SetBranchAddress("Run",&Run);
        tree->SetBranchAddress("LumiSect",&LumiSect);
        tree->SetBranchAddress("Event",&Event);
        tree->SetBranchAddress("nVtx",&nVtx);
        tree->SetBranchAddress("nInt",&nInt);
        tree->SetBranchAddress("genWeight",&genWeight);
        tree->SetBranchAddress("crossSection",&crossSection);
        tree->SetBranchAddress("pileupWeight",&pileupWeight);
        tree->SetBranchAddress("passedTrig",&passedTrig);
        tree->SetBranchAddress("triggersPassed",&triggersPassed);       
        tree->SetBranchAddress("lep_tightId", &lep_tightId);
        tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);
        tree->SetBranchAddress("lep_id", &lep_id);
        tree->SetBranchAddress("lep_genindex",&lep_genindex);
        tree->SetBranchAddress("lep_pt",&lep_pt);
        tree->SetBranchAddress("lep_pterr",&lep_pterr);
        tree->SetBranchAddress("lep_eta",&lep_eta);
        tree->SetBranchAddress("lep_phi",&lep_phi);
        tree->SetBranchAddress("lep_mass",&lep_mass);
        tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
        tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
        tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
        tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
        tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
        tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);
        tree->SetBranchAddress("jet_pt",&jet_pt);
        tree->SetBranchAddress("jet_eta",&jet_eta);
        tree->SetBranchAddress("jet_phi",&jet_phi);
        tree->SetBranchAddress("jet_mass",&jet_mass);
        tree->SetBranchAddress("jet_iscleanH4l",&jet_iscleanH4l);
        tree->SetBranchAddress("jet_QGTagger",&jet_QGTagger);
        tree->SetBranchAddress("jet_csvv2",&jet_csvv2);
        tree->SetBranchAddress("jet_jesup_pt",&jet_jesup_pt);
        tree->SetBranchAddress("jet_jesup_eta",&jet_jesup_eta);
        tree->SetBranchAddress("jet_jesup_phi",&jet_jesup_phi);
        tree->SetBranchAddress("jet_jesup_mass",&jet_jesup_mass);
        tree->SetBranchAddress("jet_jesup_iscleanH4l",&jet_jesup_iscleanH4l);
        tree->SetBranchAddress("jet_jesdn_pt",&jet_jesdn_pt);
        tree->SetBranchAddress("jet_jesdn_eta",&jet_jesdn_eta);
        tree->SetBranchAddress("jet_jesdn_phi",&jet_jesdn_phi);
        tree->SetBranchAddress("jet_jesdn_mass",&jet_jesdn_mass);
        tree->SetBranchAddress("jet_jesdn_iscleanH4l",&jet_jesdn_iscleanH4l);
        tree->SetBranchAddress("fsrPhotons_pt",&fsrPhotons_pt);
        tree->SetBranchAddress("fsrPhotons_pterr",&fsrPhotons_pterr);
        tree->SetBranchAddress("fsrPhotons_eta",&fsrPhotons_eta);
        tree->SetBranchAddress("fsrPhotons_phi",&fsrPhotons_phi);
        tree->SetBranchAddress("fsrPhotons_lepindex",&fsrPhotons_lepindex);
        tree->SetBranchAddress("met",&met);
        tree->SetBranchAddress("met_phi",&met_phi);
        tree->SetBranchAddress("nnloWeights",&nnloWeights);
        tree->SetBranchAddress("qcdWeights",&qcdWeights);
        tree->SetBranchAddress("pdfENVup",&pdfENVup);
        tree->SetBranchAddress("pdfENVdown",&pdfENVdown);        
        tree->SetBranchAddress("pdfRMSup",&pdfRMSup);        
        tree->SetBranchAddress("GENZ_DaughtersId",&GENZ_DaughtersId);
        tree->SetBranchAddress("GENZ_MomId",&GENZ_MomId);
        tree->SetBranchAddress("stage0cat",&stage0cat);
        tree->SetBranchAddress("stage1cat",&stage1cat);
        tree->SetBranchAddress("passedFiducialSelection",&passedFiducialSelection);
        tree->SetBranchAddress("GENlep_Hindex", &GENlep_Hindex);
        tree->SetBranchAddress("GENlep_id", &GENlep_id);
        tree->SetBranchAddress("GENlep_MomId", &GENlep_MomId);
        tree->SetBranchAddress("GENlep_MomMomId", &GENlep_MomMomId);
        tree->SetBranchAddress("GENlep_pt",&GENlep_pt);
        tree->SetBranchAddress("GENlep_eta",&GENlep_eta);
        tree->SetBranchAddress("GENlep_phi",&GENlep_phi);
        tree->SetBranchAddress("GENlep_mass",&GENlep_mass);

        tree->SetBranchAddress("GENjet_pt",&GENjet_pt);
        tree->SetBranchAddress("GENjet_eta",&GENjet_eta);
        tree->SetBranchAddress("GENjet_phi",&GENjet_phi);
        tree->SetBranchAddress("GENjet_mass",&GENjet_mass);

        tree->SetBranchAddress("GENmassZZ",&GENmassZZ);

        // Event Selection

        tree->SetBranchStatus("lep_Hindex",1);
        tree->SetBranchStatus("mass4l",1);
        tree->SetBranchStatus("passedZ4lSelection",1);
        tree->SetBranchStatus("passedFullSelection",1);
        tree->SetBranchStatus("passedZXCRSelection",1);
        tree->SetBranchStatus("nZXCRFailedLeptons",1);
        tree->SetBranchStatus("finalState",1);
        tree->SetBranchStatus("dataMCWeight",1);
        tree->SetBranchStatus("k_qqZZ_qcd_M",1);
        tree->SetBranchStatus("k_qqZZ_ewk",1);
        tree->SetBranchStatus("k_ggZZ",1);
        tree->SetBranchStatus("D_bkg_kin",1);
        tree->SetBranchStatus("njets_pt30_eta4p7",1);
        tree->SetBranchStatus("njets_pt30_eta2p5",1);
        tree->SetBranchStatus("nbjets_pt30_eta4p7",1);
        tree->SetBranchStatus("EventCat",1); 
       
        tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
        tree->SetBranchAddress("mass4l",&mass4l);
        tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
        tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
        tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
        tree->SetBranchAddress("nZXCRFailedLeptons",&nZXCRFailedLeptons);
        tree->SetBranchAddress("finalState",&finalState);
        tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
        tree->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
        tree->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);
        tree->SetBranchAddress("k_ggZZ",&k_ggZZ);
        tree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
        tree->SetBranchAddress("njets_pt30_eta4p7",&njets_pt30_eta4p7);
        tree->SetBranchAddress("njets_pt30_eta2p5",&njets_pt30_eta2p5);
        tree->SetBranchAddress("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7);
        tree->SetBranchAddress("EventCat",&EventCat);

        cout <<"filename "<<filename<<endl;
        cout<<"end loading the tree"<<endl;
        
    }

}

#endif
    
