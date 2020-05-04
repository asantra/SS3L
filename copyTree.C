#include <map>
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TKey.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "CalcGenericMT2/src/MT2_ROOT.h"
#include "TMath.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>;
#endif

#define RESET   "\033[0m"
#define BLACK   "\033[30m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define DEBUG 0
#define DEBUG2 0

enum SKIMMING
{
	SR_ONLY_SKIMMING=1,
	VR_ONLY_SKIMMING=2,
	SRVR_ONLY_SKIMMING=3,
	THREELEP_ONLY_SKIMMING=4,
	KEEP_ALL
};

TString basepath("/lustre/ific.uv.es/grid/atlas/t3/asantra/SS3L_Run2/ss3l_hf/prepare/xsections//");

/// map of input branch names
std::map<TString,TString> guessBranchMapping(TTree* tree)
{
	std::map<TString,TString> bnames;

	if(tree->GetBranch("EtMiss"))
	{
		cout<<"Could detected tree origin"<<endl;
		bnames["evn"]="EventNumber";
		bnames["runn"]="RunNumber";
		bnames["MCId"]="MCId";
        bnames["ht"]="Ht";
		bnames["mcweight"]="MCWeight";
		bnames["puweight"]="PileUpWeight";
		bnames["trigSF"]="TriggerSF";
		bnames["totweight"]="TotalWeight";
		bnames["met"]="EtMiss";
		bnames["phi_met"]="EtMissPhi";
		bnames["NlepBL"]="NlepBL";
		bnames["jetBtag"]="jetBtag";
		bnames["jetTLV"]="jetTLV";
		bnames["lepTLV"]="lepTLV";
		bnames["lepCharges"]="lepCharges";	
		bnames["lepTruth"]="lepTruth";
		bnames["SSChannel"]="SSChannel"; /// ==1 mumu, ==2 ee, ==3 emu
		bnames["sysNames"]="SysNames";
		bnames["sysSF"]="SysWeights";	
		/// bnames["pdfw"]="PDFWeights";
		bnames["chflipSF"]="CFTWeights";
		/// bnames["PDGId1"]="PDGId1";
		/// bnames["PDGId2"]="PDGId2";
		bnames["GLDec1"]="GLDec1";
		bnames["GLDec2"]="GLDec2";
		
		cout<<"Setup branches"<<endl;
	} 
	else {  
		cout<<"ERROR: cannot detect input tree origin. Code will crash. Have fun."<<endl;
	}

	return bnames;
}



/// set branch to local address, creating it if not existing
void setupBranch(TTree* tree, TString name, Long64_t* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/L");
	else tree->SetBranchAddress(name,address);
}

void setupBranch(TTree* tree, TString name, char* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/O");
	else tree->SetBranchAddress(name,address);
}

void setupBranch(TTree* tree, TString name, bool* address){
  TBranch* _b=tree->GetBranch(name);
  if(_b==0) _b= tree->Branch(name,address,name+"/O");
  else tree->SetBranchAddress(name,address);
}

void setupBranch(TTree* tree, TString name, int* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/I");
	else tree->SetBranchAddress(name,address);
}

void setupBranch(TTree* tree, TString name, float* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/F");
	else tree->SetBranchAddress(name,address);
}


bool isAtlfast(Long64_t RunNumber){
  /// signal samples:
  if((RunNumber>=376776 && RunNumber<=376834) || (RunNumber>=376844 && RunNumber<=376859) || RunNumber==377110) return true; /// 2STEPWZ
  if(RunNumber>=436533 && RunNumber<=436594) return true; /// BTT
  if((RunNumber>=376733 && RunNumber<=376775) || (RunNumber>=376835 && RunNumber<=376843)) return true; /// GG_Rpvlampp331
  if(RunNumber>=377000 && RunNumber<=377071) return true; /// GG_RpvLQD
  if((RunNumber>=375818 && RunNumber<=375932) || RunNumber==376122 || (RunNumber>=377605 && RunNumber<=377640)) return true; /// GTT  
  if(RunNumber>=377094 && RunNumber<=377102) return true; /// TT2STEP
  if(RunNumber>=449781 && RunNumber<=449822) return true; /// bRPV  
  if(RunNumber>=377408 && RunNumber<=377593) return true; /// Wino offshell WZ 
  if((RunNumber>=392200 && RunNumber<=392280) || (RunNumber>=394776 && RunNumber<=394779) || (RunNumber>=396206 && RunNumber<=397314) || RunNumber==398177) return true; /// Wino onshell WZ 
  
  /// bkg samples
  if(RunNumber==412043) return true; /// SM4topsNLO
  
  return false;
}

bool isSignal(int RunNumber){
  if((RunNumber>=376776 && RunNumber<=376834) || (RunNumber>=376844 && RunNumber<=376859) || RunNumber==377110) return true; /// 2STEPWZ
  if(RunNumber>=436533 && RunNumber<=436594) return true; /// BTT
  if((RunNumber>=376733 && RunNumber<=376775) || (RunNumber>=376835 && RunNumber<=376843)) return true; /// GG_Rpvlampp331
  if(RunNumber>=377000 && RunNumber<=377071) return true; /// GG_RpvLQD
  if((RunNumber>=375818 && RunNumber<=375932) || RunNumber==376122 || (RunNumber>=377605 && RunNumber<=377640)) return true; /// GTT  
  if(RunNumber>=377094 && RunNumber<=377102) return true; /// TT2STEP
  if(RunNumber>=449781 && RunNumber<=449822) return true; /// bRPV
  if(RunNumber>=377408 && RunNumber<=377593) return true; /// Wino offshell WZ 
  if((RunNumber>=392200 && RunNumber<=392280) || (RunNumber>=394776 && RunNumber<=394779) || (RunNumber>=396206 && RunNumber<=397314) || RunNumber==398177) return true; /// Wino onshell WZ 
  return false;
}


int hasZee(std::vector<TLorentzVector>  leptons, std::vector<float> charges){

  float M[2] = {81000., 101000};
  std::vector<TLorentzVector> electrons(0);
  std::vector<float> el_charges(0);
  for(unsigned int i(0); i<leptons.size(); i++){
    if(leptons.at(i).M()<50.) {
      electrons.push_back(leptons.at(i));
      el_charges.push_back(charges.at(i));
    }
  }

  int Zee(0);
  switch( (int)electrons.size() ){
  case 2:
    if( (electrons.at(0)+electrons.at(1)).M() > M[0] && (electrons.at(0)+electrons.at(1)).M() < M[1] && el_charges.at(0)==el_charges.at(1)) Zee = 1;
    break;
  case 3:
    if( (electrons.at(0)+electrons.at(1)).M() > M[0] && (electrons.at(0)+electrons.at(1)).M() < M[1] && el_charges.at(0)==el_charges.at(1)) Zee = 1;
    if( (electrons.at(0)+electrons.at(2)).M() > M[0] && (electrons.at(0)+electrons.at(2)).M() < M[1] && el_charges.at(0)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(1)+electrons.at(2)).M() > M[0] && (electrons.at(1)+electrons.at(2)).M() < M[1] && el_charges.at(1)==el_charges.at(2)) Zee = 1;
    break;
  case 4:
    if( (electrons.at(0)+electrons.at(1)).M() > M[0] && (electrons.at(0)+electrons.at(1)).M() < M[1] && el_charges.at(0)==el_charges.at(1)) Zee = 1;
    if( (electrons.at(0)+electrons.at(2)).M() > M[0] && (electrons.at(0)+electrons.at(2)).M() < M[1] && el_charges.at(0)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(1)+electrons.at(2)).M() > M[0] && (electrons.at(1)+electrons.at(2)).M() < M[1] && el_charges.at(1)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(0)+electrons.at(3)).M() > M[0] && (electrons.at(0)+electrons.at(3)).M() < M[1] && el_charges.at(0)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(1)+electrons.at(3)).M() > M[0] && (electrons.at(1)+electrons.at(3)).M() < M[1] && el_charges.at(1)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(2)+electrons.at(3)).M() > M[0] && (electrons.at(2)+electrons.at(3)).M() < M[1] && el_charges.at(2)==el_charges.at(3)) Zee = 1;
    break;
  case 5:
    if( (electrons.at(0)+electrons.at(1)).M() > M[0] && (electrons.at(0)+electrons.at(1)).M() < M[1] && el_charges.at(0)==el_charges.at(1)) Zee = 1;
    if( (electrons.at(0)+electrons.at(2)).M() > M[0] && (electrons.at(0)+electrons.at(2)).M() < M[1] && el_charges.at(0)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(1)+electrons.at(2)).M() > M[0] && (electrons.at(1)+electrons.at(2)).M() < M[1] && el_charges.at(1)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(0)+electrons.at(3)).M() > M[0] && (electrons.at(0)+electrons.at(3)).M() < M[1] && el_charges.at(0)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(1)+electrons.at(3)).M() > M[0] && (electrons.at(1)+electrons.at(3)).M() < M[1] && el_charges.at(1)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(2)+electrons.at(3)).M() > M[0] && (electrons.at(2)+electrons.at(3)).M() < M[1] && el_charges.at(2)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(0)+electrons.at(4)).M() > M[0] && (electrons.at(0)+electrons.at(4)).M() < M[1] && el_charges.at(0)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(1)+electrons.at(4)).M() > M[0] && (electrons.at(1)+electrons.at(4)).M() < M[1] && el_charges.at(1)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(2)+electrons.at(4)).M() > M[0] && (electrons.at(2)+electrons.at(4)).M() < M[1] && el_charges.at(2)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(3)+electrons.at(4)).M() > M[0] && (electrons.at(3)+electrons.at(4)).M() < M[1] && el_charges.at(3)==el_charges.at(4)) Zee = 1;
    break;
  case 6:
    if( (electrons.at(0)+electrons.at(1)).M() > M[0] && (electrons.at(0)+electrons.at(1)).M() < M[1] && el_charges.at(0)==el_charges.at(1)) Zee = 1;
    if( (electrons.at(0)+electrons.at(2)).M() > M[0] && (electrons.at(0)+electrons.at(2)).M() < M[1] && el_charges.at(0)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(1)+electrons.at(2)).M() > M[0] && (electrons.at(1)+electrons.at(2)).M() < M[1] && el_charges.at(1)==el_charges.at(2)) Zee = 1;
    if( (electrons.at(0)+electrons.at(3)).M() > M[0] && (electrons.at(0)+electrons.at(3)).M() < M[1] && el_charges.at(0)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(1)+electrons.at(3)).M() > M[0] && (electrons.at(1)+electrons.at(3)).M() < M[1] && el_charges.at(1)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(2)+electrons.at(3)).M() > M[0] && (electrons.at(2)+electrons.at(3)).M() < M[1] && el_charges.at(2)==el_charges.at(3)) Zee = 1;
    if( (electrons.at(0)+electrons.at(4)).M() > M[0] && (electrons.at(0)+electrons.at(4)).M() < M[1] && el_charges.at(0)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(1)+electrons.at(4)).M() > M[0] && (electrons.at(1)+electrons.at(4)).M() < M[1] && el_charges.at(1)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(2)+electrons.at(4)).M() > M[0] && (electrons.at(2)+electrons.at(4)).M() < M[1] && el_charges.at(2)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(3)+electrons.at(4)).M() > M[0] && (electrons.at(3)+electrons.at(4)).M() < M[1] && el_charges.at(3)==el_charges.at(4)) Zee = 1;
    if( (electrons.at(0)+electrons.at(5)).M() > M[0] && (electrons.at(0)+electrons.at(5)).M() < M[1] && el_charges.at(0)==el_charges.at(5)) Zee = 1;
    if( (electrons.at(1)+electrons.at(5)).M() > M[0] && (electrons.at(1)+electrons.at(5)).M() < M[1] && el_charges.at(1)==el_charges.at(5)) Zee = 1;
    if( (electrons.at(2)+electrons.at(5)).M() > M[0] && (electrons.at(2)+electrons.at(5)).M() < M[1] && el_charges.at(2)==el_charges.at(5)) Zee = 1;
    if( (electrons.at(3)+electrons.at(5)).M() > M[0] && (electrons.at(3)+electrons.at(5)).M() < M[1] && el_charges.at(3)==el_charges.at(5)) Zee = 1;
    if( (electrons.at(4)+electrons.at(5)).M() > M[0] && (electrons.at(4)+electrons.at(5)).M() < M[1] && el_charges.at(4)==el_charges.at(5)) Zee = 1;
  default: break;
    
  }
  return Zee;
}

bool isResonance(std::vector<TLorentzVector> vectors, std::vector< float > charges, int channel){
	bool res(false);
	bool OS[6] = {0, 0, 0, 0, 0, 0};
	float M[6] = {0., 0., 0., 0., 0., 0.};
	float MZlow(81200.), MZup(101200.);

	if(vectors.size() < 2) return res;
	if(vectors.size() == 2 && channel==3) return res;

	if(vectors.size() == 2){
		M[0] = (vectors.at(0) + vectors.at(1)).M();
		if(charges.at(0) != charges.at(1)) OS[0] = true;
	}
	else if(vectors.size() == 3){
		M[0] = (vectors.at(0) + vectors.at(1)).M();
		M[1] = (vectors.at(1) + vectors.at(2)).M();
		M[2] = (vectors.at(0) + vectors.at(2)).M();
		if(charges.at(0) != charges.at(1)) OS[0] = true;
		if(charges.at(1) != charges.at(2)) OS[1] = true;
		if(charges.at(0) != charges.at(2)) OS[2] = true;
	}
	else{
		M[0] = (vectors.at(0) + vectors.at(1)).M();
		M[1] = (vectors.at(1) + vectors.at(2)).M();
		M[2] = (vectors.at(0) + vectors.at(2)).M();
		M[3] = (vectors.at(3) + vectors.at(1)).M();
		M[4] = (vectors.at(3) + vectors.at(2)).M();
		M[5] = (vectors.at(3) + vectors.at(0)).M();
		if(charges.at(0) != charges.at(1)) OS[0] = true;
		if(charges.at(1) != charges.at(2)) OS[1] = true;
		if(charges.at(0) != charges.at(2)) OS[2] = true;
		if(charges.at(3) != charges.at(1)) OS[3] = true;
		if(charges.at(3) != charges.at(2)) OS[4] = true;
		if(charges.at(3) != charges.at(0)) OS[5] = true;
	}

	if( (OS[0] && (M[0] > MZlow && M[0] < MZup)) 
	  || (OS[1] && (M[1] > MZlow && M[1] < MZup)) 
	  || (OS[2] && (M[2] > MZlow && M[2] < MZup)) 
	  || (OS[3] && (M[3] > MZlow && M[3] < MZup))
	  || (OS[4] && (M[4] > MZlow && M[4] < MZup))
	  || (OS[5] && (M[5] > MZlow && M[5] < MZup)))
	  {
		  res = true;
		  /// cout << M[0] << " " << M[1] << " " << M[2] << " " << M[3] << " " << M[4] << " " << M[5] << " " << vectors.size() << endl;
	  }	 
	  
	return res;
}


float getSumPtLep(std::vector<TLorentzVector> leptons){
	float sumPt(0);
	if(leptons.size()==0) return sumPt;
	for(unsigned int i(0); i<leptons.size(); i++) sumPt += (leptons.at(i)).Pt();
	return sumPt;

}
float getSumMinvLep(std::vector<TLorentzVector> leptons){
	float sumM(0);
	if(leptons.size()==0) return sumM;
	for(unsigned int i(0); i<leptons.size(); i++) sumM += (leptons.at(i)).M();
	return sumM;
}

float getSumPtJet(std::vector<TLorentzVector> jets, std::vector<float> btags, TString type){
  float sumPt(0);
  if(jets.size()==0 || btags.size()==0) return sumPt;
  for(unsigned int i(0); i<jets.size(); i++){
    float pt(0);
    if(type=="B"){ pt = btags.at(i) ? (jets.at(i)).Pt() : 0.; }
    if(type=="L"){ pt = !btags.at(i) ? (jets.at(i)).Pt() : 0.; }
    if(type=="ALL"){ pt = (jets.at(i)).Pt(); }
    if(pt<25000 && type=="ALL")  continue;
    sumPt += pt;
  }
  return sumPt;
}

float getSumPtJet_noPtCut(std::vector<TLorentzVector> jets){
  float sumPt(0);
  if(jets.size()==0) return sumPt;
  for(unsigned int i(0); i<jets.size(); i++){
    float pt = (jets.at(i)).Pt(); 
    sumPt += pt;
  }
  return sumPt;
}


float dR(float var1, float phi1, float var2, float phi2){
	float deta = TMath::Abs(var1 - var2);
	float dphi = TMath::Abs(phi1 - phi2) < TMath::Pi() ? TMath::Abs(phi1 - phi2) : 2*TMath::Pi() - TMath::Abs(phi1 - phi2);
	return TMath::Sqrt(deta*deta + dphi*dphi);
}

float dReta(TLorentzVector v1, TLorentzVector v2){
	return dR(v1.Eta(), v1.Phi(), v2.Eta(), v2.Phi());
}

float dRy(TLorentzVector v1, TLorentzVector v2){
	return dR(v1.Rapidity(), v1.Phi(), v2.Rapidity(), v2.Phi());
}

float getDeltaRLepLep(std::vector<TLorentzVector> leptons){
	if(leptons.size()==0) return -99.;
	if(leptons.size() < 2) return -99.;
	return dReta(leptons.at(0), leptons.at(1)); 
}

float getDeltaRJJ(std::vector<TLorentzVector> jets){
    if(jets.size()==0) return -99.;
    if(jets.size() < 2) return -99.;
    return dReta(jets.at(0), jets.at(1)); 
}

float getDeltaRLepJet(std::vector<TLorentzVector> leptons, std::vector<TLorentzVector> jets, unsigned int Nl){
  if (leptons.size()<2) return -99.;
  float minDR = 99.;
  for (unsigned int j=0;j<jets.size();j++){
    if (dReta(leptons.at(Nl), jets.at(j))<minDR && jets.at(j).Pt()>25000){
      minDR=dReta(leptons.at(Nl), jets.at(j));
    }
  }
  return minDR;
}


bool match(TLorentzVector v1, TLorentzVector v2){
	if(v1.Pt()<1E-4 || v2.Pt()<1E-4) return false;
	if(dReta(v1,v2)<1E-4 && TMath::Abs(v1.Pt()-v2.Pt())<1E-4) return true;
	else return false;
}


int getSSNegative(std::vector<float> *charges){
	int count=0;
	for(unsigned int i=0;i<charges->size();i++){
		if (charges->at(i)<0) count ++;
	}
	return count;
}

int has3LSSPrompt(std::vector<float> charges, std::vector<float> lepTruth){
  int count_pos=0;
  int count_neg=0;
  /// std::cout<<"**********"<<std::endl;
  for(unsigned int i=0;i<charges.size();i++){
    /// std::cout<<"charge: "<<charges.at(i)<<std::endl;
    if(charges.at(i)>0 && lepTruth.at(i)==1) count_pos++;
    if(charges.at(i)<0 && lepTruth.at(i)==1) count_neg++;
  }
  /// std::cout<<"count_pos: "<<count_pos<<"\t count_neg: "<<count_neg<<std::endl;
  if (count_pos>2 || count_neg>2) return 1;
  else return 0;
}

int is3LSSprocess(int MCId)
{
  bool isBkg = (MCId == 407321 ||
		MCId == 364245 ||
		MCId == 364247 ||
		MCId == 342284 ||
		MCId == 342285 ||
		MCId == 345940 ||
		MCId == 345941);
		
  bool isSig = (MCId>=377094 && MCId <=377102);

  return (int)(isBkg || isSig);
}

int isSSLepPtCut(std::vector<TLorentzVector> leptons,std::vector<float> charges){
  int isAbove30 = 0;
  int isSSpair = 0;
  for(unsigned int i(0); i<leptons.size(); i++){
    for(unsigned int j(0); j<leptons.size(); j++){
      if(i==j || isSSpair==1)
        continue;
      if(charges.at(i)==charges.at(j))
	isSSpair = 1;
      if(isSSpair==1 && leptons.at(i).Pt()>30000 && leptons.at(j).Pt()>30000)
        isAbove30=1;
    }
  }
  return isAbove30;
}


float GetSFOSmass(std::vector<TLorentzVector> leptons,std::vector<float> charges)
{
  float c(0);
  int nLep = leptons.size();
  switch(nLep){
  case 2:
    c = 0;
    break;
  case 3:
    if( charges.at(0) != charges.at(1) && (int)leptons.at(0).M()==(int)leptons.at(1).M() ){
      c=(leptons.at(0)+leptons.at(1)).M();
    }
    if( charges.at(0) != charges.at(2) && (int)leptons.at(0).M()==(int)leptons.at(2).M() ){
      if(TMath::Abs((leptons.at(0)+leptons.at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(0)+leptons.at(2)).M();
    }
    if( charges.at(2) != charges.at(1) && (int)leptons.at(2).M()==(int)leptons.at(1).M() ){
      if(TMath::Abs((leptons.at(2)+leptons.at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(2)+leptons.at(1)).M();
    }
    break;
  case 4:
    if( charges.at(0) != charges.at(1) && (int)leptons.at(0).M()==(int)leptons.at(1).M() ){
      if(TMath::Abs((leptons.at(0)+leptons.at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(0)+leptons.at(1)).M();
    }
    if( charges.at(0) != charges.at(2) && (int)leptons.at(0).M()==(int)leptons.at(2).M() ){
      if(TMath::Abs((leptons.at(0)+leptons.at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(0)+leptons.at(2)).M();
    }
    if( charges.at(0) != charges.at(3) && (int)leptons.at(0).M()==(int)leptons.at(3).M() ){
      if(TMath::Abs((leptons.at(0)+leptons.at(3)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(0)+leptons.at(3)).M();
    }
    if( charges.at(1) != charges.at(2) && (int)leptons.at(1).M()==(int)leptons.at(2).M() ){
      if(TMath::Abs((leptons.at(1)+leptons.at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(1)+leptons.at(2)).M();
    }
    if( charges.at(1) != charges.at(3) && (int)leptons.at(1).M()==(int)leptons.at(3).M() ){
      if(TMath::Abs((leptons.at(3)+leptons.at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(3)+leptons.at(1)).M();
    }
    if( charges.at(2) != charges.at(3) && (int)leptons.at(2).M()==(int)leptons.at(3).M() ){
      if(TMath::Abs((leptons.at(3)+leptons.at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons.at(2)+leptons.at(3)).M();
    }
    break;
  default:
    break;
  }
  return c;
}

float getMt(std::vector<TLorentzVector> leptons, float met, float metPhi)
{
	if( leptons.empty() ) return 0;

	TVector3 lepVec(0.,0.,0.);
	lepVec.SetPtEtaPhi(leptons.at(0).Pt(), 0., leptons.at(0).Phi());
	TLorentzVector Vt(0.,0.,0.,0.);
	Vt.SetPxPyPzE( met*TMath::Cos(metPhi)+lepVec.Px(),met*TMath::Sin(metPhi)+lepVec.Py(), 0., met+lepVec.Pt() );

	return (float)Vt.M();
}

TString printVec(std::vector<float> v){
	if( v.size() < 3 ) return "";
	return Form(" [%.5f|%.5f|%.5f]",v[0],v[1],v[2]);
}



//// taken from Marco

// mljj_highestPt
auto cmljj_highestPt( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    if((LepPt.size()==0) || (JetPt.size()<2))
        return (-99.0);
    else{
        TLorentzVector lep_vec;
        lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
        TLorentzVector jet_vec; 
        TLorentzVector jet_vec0; 
        jet_vec0.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
        TLorentzVector jet_vec1; 
        jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
        jet_vec = jet_vec0 + jet_vec1;
        return ( ( jet_vec + lep_vec ).M() );
    }
};




// jets closest to leptons
auto closestJet( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, size_t &minimum, size_t &secondMin ) {
    if(JetPt.size()<2 || LepPt.size()==0)
        return (-1.0);
    else{
        vector<float> drVec;
        drVec.clear();
        for(size_t i=0; i < JetPt.size(); ++i){
            float dR_each = dR(LepEta.at(0), LepPhi.at(0), JetEta.at(i), JetPhi.at(i));
            drVec.push_back(dR_each);
        }
        std::sort(drVec.begin(),drVec.end());
        for(size_t i=0; i < JetPt.size(); ++i){
            float dR_each = dR(LepEta.at(0), LepPhi.at(0), JetEta.at(i), JetPhi.at(i));
            if(dR_each == drVec.at(0)) minimum=i;
            if(dR_each == drVec.at(1)) secondMin=i;
        }
        return (1.0);
    }
};




// mljj_closest
auto cmljj_closest( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    if((LepPt.size()==0) || (JetPt.size()<2))
        return (-99.0);
    else{
        size_t minimum(99), secondMin(99);
        closestJet( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM, minimum, secondMin );
        if(minimum==99){
            std::cout << "Something wrong in finding closest jets for mljj. Exiting" << std::endl;
            exit(1);
        }
        else{
            TLorentzVector lep_vec;
            lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
            TLorentzVector jet_vec; 
            TLorentzVector jet_vec0; 
            jet_vec0.SetPtEtaPhiM( JetPt[minimum], JetEta[minimum], JetPhi[minimum], JetM[minimum] );
            TLorentzVector jet_vec1; 
            jet_vec1.SetPtEtaPhiM( JetPt[secondMin], JetEta[secondMin], JetPhi[secondMin], JetM[secondMin] );
            jet_vec = jet_vec0 + jet_vec1;
            return ( ( jet_vec + lep_vec ).M() );
        }
    }
};




// mljj_closestToW
auto cmljj_closestToW( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    if((LepPt.size()==0) || (JetPt.size()<2))
        return (-999.0);
    else{
        size_t minimum(99), secondMin(99);
        closestJet( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM, minimum, secondMin );
        if(minimum==99){
            std::cout << "Something wrong in finding closest jets for mljj. Exiting" << std::endl;
            exit(1);
        }
        else{
            TLorentzVector lep_vec;
            lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
            TLorentzVector jet_vec; 
            TLorentzVector jet_vec0; 
            jet_vec0.SetPtEtaPhiM( JetPt[minimum], JetEta[minimum], JetPhi[minimum], JetM[minimum] );
            TLorentzVector jet_vec1; 
            jet_vec1.SetPtEtaPhiM( JetPt[secondMin], JetEta[secondMin], JetPhi[secondMin], JetM[secondMin] );
            jet_vec = jet_vec0 + jet_vec1;
            float massW = 80400.0;
            if(fabs(jet_vec.M() - massW)<= 10000.0){
                return ( ( jet_vec + lep_vec ).M() );
            }
            else{
                return (-999.0);
            }
        }
    }
};


// mjj_InW
auto cmjj_InW( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, size_t &min, size_t &secMin ) {
    if((JetPt.size()<2))
        return (-999.0);
    else{
        vector<float> dRvec;
        dRvec.clear();
        if(DEBUG2) std::cout << "In cmjj_InW function" << std::endl;
        for(size_t i = 0; i < JetPt.size()-1; ++i){
            if(DEBUG2) std::cout << "In side first jet for " << i << std::endl;
            TLorentzVector jet_vec0; 
            jet_vec0.SetPtEtaPhiM( JetPt[i], JetEta[i], JetPhi[i], JetM[i] );
            for(size_t j=i+1; j < JetPt.size(); ++j){
                if(DEBUG2) std::cout << "In side second jet for " << j << std::endl;
                TLorentzVector jet_vec, jet_vec1; 
                jet_vec1.SetPtEtaPhiM( JetPt[j], JetEta[j], JetPhi[j], JetM[j] );
                jet_vec = jet_vec0 + jet_vec1;
                if(DEBUG2) std::cout << "jet_vec mass in general for " << i << " and " << j << " : " << jet_vec.M() << std::endl;
                float massW = 80400.0;
                if(fabs(jet_vec.M() - massW)<= 10000.0){
                    if(DEBUG2) std::cout << "jet_vec mass in W window for " << i << " and " << j << " : " << jet_vec.M() << std::endl;
                    float dR_each = dR(JetEta.at(i), JetPhi.at(i), JetEta.at(j), JetPhi.at(j));
                    if(DEBUG2) std::cout << "dR of jets in W window for " << i << " and " << j << " : " << dR_each << std::endl;
                    dRvec.push_back(dR_each);
                }
            }
        }
        if(dRvec.size()>0){
            if(DEBUG2) std::cout << "When jet mass is inside W window " << std::endl;
            std::sort(dRvec.begin(),dRvec.end());
            for(size_t i = 0; i < JetPt.size()-1; ++i){
                for(size_t j=i+1; j < JetPt.size(); ++j){
                    float dR_each = dR(JetEta.at(i), JetPhi.at(i), JetEta.at(j), JetPhi.at(j));
                    if(DEBUG2) std::cout << "dR of jets for " << i << " and " << j << " : " << dR_each << std::endl;
                    if(dR_each == dRvec[0]){
                        min    = i;
                        secMin = j;
                        if(DEBUG2) std::cout << "found the minimum of dR for " << i << " and " << j << ": " << dR_each <<std::endl;
                        if(DEBUG2) std::cout << "the min :" << i << " secondMin : " << j << std::endl;
                    }
                }
            }
            return(1.0);
        }
        else{
            return(-999.0);
        }
    }
};

// lepton closest to jets
auto cmljj_shiftToW( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    if(JetPt.size()<2 || LepPt.size()==0)
        return (-999.0);
    else{
        size_t firstJet(99), secondJet(99);
        float success = cmjj_InW( JetPt, JetPhi, JetEta, JetM, firstJet, secondJet );
        if(success==1.0){
            if(DEBUG2) std::cout << " Inside cmljj_shiftToW function " << std::endl;
            if(DEBUG2) std::cout << " The first jet in cmljj_shiftToW " << firstJet << " and the secondJet: " << secondJet << std::endl;
            TLorentzVector jet1, jet2, jet;
            jet1.SetPtEtaPhiM( JetPt[firstJet], JetEta[firstJet], JetPhi[firstJet], JetM[firstJet] );
            jet2.SetPtEtaPhiM( JetPt[secondJet], JetEta[secondJet], JetPhi[secondJet], JetM[secondJet] );
            jet = jet1 + jet2;
            vector<float> dRvec;
            dRvec.clear();
            for(size_t i = 0; i <LepPt.size(); ++i){
                TLorentzVector lep;
                lep.SetPtEtaPhiM( LepPt[i], LepEta[i], LepPhi[i], LepM[i] );
                if(DEBUG2) std::cout << " The deltaR of combined jet and lep for lep" << i << " : " << jet.DeltaR(lep) << std::endl;
                dRvec.push_back(jet.DeltaR(lep));
            }
            if(dRvec.size()>0){
                std::sort(dRvec.begin(),dRvec.end());
                for(size_t i = 0; i < LepPt.size(); ++i){
                    TLorentzVector lep2, total;
                    lep2.SetPtEtaPhiM( LepPt[i], LepEta[i], LepPhi[i], LepM[i] );
                    total = jet+lep2;
                    if(DEBUG2) std::cout << " Second Time - The deltaR of combined jet and lep for lep" << i << " : " << jet.DeltaR(lep2) << std::endl;
                    float jetLepDr = jet.DeltaR(lep2);
                    if(jetLepDr == dRvec[0]){
                        if(DEBUG2) std::cout << " Second Time - The deltaR of jet and lep for minimum lep " << i << " : " << jet.DeltaR(lep2) << std::endl;
                        if(DEBUG2) std::cout << " The shifted mass " << (total.M() - jet.M() + 80400.0) << std::endl;
                        return (total.M() - jet.M() + 80400.0);
                    }
                }
            }
        }
        return (-999.0);
    }
};




/// end code Marco


std::vector<float> correctChFlipSF(std::vector<float> *cfv, std::vector<TLorentzVector> *leptons, bool use){
	bool msg(false);
	std::vector<float> cfvOld(0), cfvNew(0);
	for(unsigned int i(0); i<cfv->size(); i++){ 
		cfvNew.push_back( 1 );
		cfvOld.push_back( (float)cfv->at(i) );
	}
	if( !use ) return cfvOld;
	if( cfv->at(0) > 0.001 ){ 
		if(msg){std::cout <<"correctChFlipSF \t"<<"ChargeFlipSF ="<<cfv->at(0)<<" return "<<printVec(cfvOld)<< std::endl;} 
		return cfvOld;
	}
	for(unsigned int i(0); i<leptons->size(); i++){
		if( leptons->at(i).M()<50. && leptons->at(i).Pt()<15000){ 
			if(msg){std::cout <<"correctChFlipSF \t"<<"ChargeFlipSF ="<< cfv->at(0)<<Form("\tLepton %i (M=%.3f) has Pt=%.3f",i,leptons->at(i).M(),leptons->at(i).Pt())<<" return "<<printVec(cfvNew)<< std::endl;}
			return cfvNew;
		}
	}
	if(msg){std::cout << "correctChFlipSF \t" << "ChargeFlipSF =" <<cfv->at(0)<<" but no lep pT<15 GeV found... stay with old SF" << std::endl;}
	return cfvOld;
}


void copyTree(TTree*inTree,TTree*outTree,float xsec,float runc,Long64_t MC_ID,float nMCevents,float lumi,std::vector<Long64_t>& writtenEvents,bool DoSYS,SKIMMING mySkimming)
{
	//std::cout<<"copying tree. already written "<<writtenEvents.size()<<endl;

	// input and output trees, tag/folder name, xsection, #events in original MC sample, target lumi
	//ofstream myfile;
	//std::string s = std::to_string(MC_ID);
	//std::string l = std::to_string(lumi);
	//myfile.open("Rpc3LSS1b_"+s+"_"+l+".txt");

	Long64_t evn64=0,runn64=0;
	Int_t MCId=0;
	
	Int_t PDGId1=0,PDGId2=0,GLDec1=0,GLDec2=0;
		
    float totweight=0, TotWeightNew=0, MC_campaign_weight=1, mcweight=0, puweight=0;
	
	int nBJets20=0,nJets25=0,nJets35=0,nJets40=0,nJets50;
	int nSigLep=0,NlepBL=0;
	int SSChannel=0,isSS30=0,isZee=0,is3LSS=0,SSNegative=0,is3LSSproc=0;
	bool Zevent=0;
	
	float met=0,meff=0,mt=0, MT2=0;
	float ht=0,minvAllSigL=0,Mjj=0,Mll=0,Pt_l=0, Pt_subl=0,mSFOS=0,phi_met=0;
    
    float mljj_highestPt = 0, mljj_closest=0, mljj_closestToW=0, mlW_shifted=0;
    
    float jetPt_1=0, jetPt_2=0;

	float sumPtLep(0), SumJetPt(0), SumBJetPt(0);
	float dRll(-99.), dRl1j(-99.), dRl2j(-99.), dRJJ(-99.);
	
	bool isSR=false, isVR=false;

	float elID_rescale_up=1, elID_rescale_down=1;
	std::vector<float> *SysSF=0;
	float TrigSF=0;
	std::vector<float> *ChFlipSF=0;
	std::vector<float> *PdfSF=0;
	std::vector<float> *pdfw=0;
	std::vector<string> *SysNames=0;
	std::vector<TLorentzVector> *lepTLV=0, *jetTLV=0;
	std::vector<TLorentzVector> tmp_jetTLV,tmp_lepTLV;
	std::vector<float> tmp_jetBtag,tmp_lepCharges,tmp_lepTruth;
	std::vector<float> *lepCharges=0,*lepTruth=0;
	std::vector<float> *jetBtag=0;
    std::vector<float> lepPt, lepEta, lepPhi, lepM;
    std::vector<float> jetPt, jetEta, jetPhi, jetM;

	std::map<TString,TString> bnames=guessBranchMapping(inTree);

	float wmu_nom,
		  wmu_stat_up, wmu_stat_down,
		  wmu_sys_up, wmu_sys_down,
		  wmu_bad_sys_up, wmu_bad_sys_down,
		  wmu_bad_stat_up, wmu_bad_stat_down,
		  wmu_stat_lowpt_up, wmu_stat_lowpt_down,
		  wmu_sys_lowpt_up, wmu_sys_lowpt_down,
		  wmu_trig_stat_up, wmu_trig_stat_down,
		  wmu_trig_sys_up, wmu_trig_sys_down, 
		  wmu_iso_stat_up, wmu_iso_stat_down,
		  wmu_iso_sys_up, wmu_iso_sys_down, 
		  wmu_ttva_stat_up, wmu_ttva_stat_down,
		  wmu_ttva_sys_up, wmu_ttva_sys_down,

		  wel_nom,
		  wel_cid_up, wel_cid_down,
		  wel_id_up, wel_id_down,
		  wel_reco_up, wel_reco_down,
		  wel_trig_up, wel_trig_down,
		  wel_trigEff_up, wel_trigEff_down,
		  wel_iso_up, wel_iso_down,

		  wjet_nom,
		  wjet_b_up, wjet_b_down,
		  wjet_c_up, wjet_c_down,
		  wjet_light_up, wjet_light_down,
		  wjet_extra1_up, wjet_extra1_down,
		  wjet_extra2_up, wjet_extra2_down,

		  wjet_jvt_up, wjet_jvt_down,

		  wpu_nom_sig,
		  wpu_up_sig, wpu_down_sig,

		  wpu_nom_bkg,
		  wpu_up_bkg, wpu_down_bkg,

		  wtrig_nom,
		  wtrig_up, wtrig_down,

		  wchflip_nom,
		  wchflip_up, wchflip_down;

		  ///wpdf_up, wpdf_down;


	int nentries=inTree->GetEntries();
	/// std::cout<<"inTree has "<<nentries<<std::endl;
	/// cout<<"got xsec nentries nMCevents "<<xsec<< " "<<nentries<<" "<<nMCevents<<endl;
	float lumiScaling=xsec*lumi*1/nMCevents;
	/// data
	if(xsec==0) lumiScaling=1;
	/// cout<<"final scaling "<<lumiScaling<<endl;
	float lumiScalingUP=lumiScaling;
	float lumiScalingDOWN=lumiScaling;


	/// avoid reading what we do not need
	inTree->SetBranchStatus("*",0);
	/// turn on only the branches that will be read
	for(auto branch: bnames)
		inTree->SetBranchStatus(branch.second,1);
	///cout<<"SetBranchStatus set to 1. for the needed branches"<<endl;


	/// configure input branches
	TClass* expectedClass = 0;
	EDataType expectedType = kOther_t;
	inTree->GetBranch(bnames["evn"])->GetExpectedType(expectedClass,expectedType);
	inTree->SetBranchAddress(bnames["evn"],&evn64);
	inTree->SetBranchAddress(bnames["runn"],&runn64);
    inTree->SetBranchAddress(bnames["ht"], &ht);
	inTree->SetBranchAddress(bnames["MCId"],&MCId);
	inTree->SetBranchAddress(bnames["mcweight"],&mcweight);
	inTree->SetBranchAddress(bnames["puweight"],&puweight);
	inTree->SetBranchAddress(bnames["trigSF"],&TrigSF);
	inTree->SetBranchAddress(bnames["totweight"],&totweight);	
	inTree->SetBranchAddress(bnames["met"],&met);
	inTree->SetBranchAddress(bnames["phi_met"],&phi_met);
	inTree->SetBranchAddress(bnames["NlepBL"],&NlepBL);
	inTree->SetBranchAddress(bnames["jetBtag"],&jetBtag);
	inTree->SetBranchAddress(bnames["jetTLV"],&jetTLV);	
	inTree->SetBranchAddress(bnames["lepTLV"],&lepTLV);
	inTree->SetBranchAddress(bnames["lepCharges"],&lepCharges);		
	inTree->SetBranchAddress(bnames["lepTruth"],&lepTruth);
	inTree->SetBranchAddress(bnames["SSChannel"],&SSChannel); /// ==1 mumu, ==2 ee, ==3 emu
	inTree->SetBranchAddress(bnames["sysNames"],&SysNames);
	inTree->SetBranchAddress(bnames["sysSF"],&SysSF);
	inTree->SetBranchAddress(bnames["chflipSF"],&ChFlipSF);	
	/// inTree->SetBranchAddress(bnames["PDGId1"],&PDGId1);
	/// inTree->SetBranchAddress(bnames["PDGId2"],&PDGId2);
	inTree->SetBranchAddress(bnames["GLDec1"],&GLDec1);
	inTree->SetBranchAddress(bnames["GLDec2"],&GLDec2);
	cout <<"inTree->SetBranchAddress done"<<endl;


	/// configure output branches
	setupBranch(outTree,"evn",&evn64);
	setupBranch(outTree,"runn",&runn64);
	setupBranch(outTree,"MCId",&MCId);
	setupBranch(outTree,"mcweight",&mcweight);
	setupBranch(outTree,"MC_campaign_weight",&MC_campaign_weight);
	setupBranch(outTree,"puweight",&puweight);
	setupBranch(outTree,"totweight",&TotWeightNew);
	setupBranch(outTree,"lumiScaling",&lumiScaling);
	setupBranch(outTree,"NlepBL",&NlepBL);
	setupBranch(outTree,"nSigLep",&nSigLep);	
	setupBranch(outTree,"nBJets20",&nBJets20);
	setupBranch(outTree,"nJets25",&nJets25);
	setupBranch(outTree,"nJets35",&nJets35);
	setupBranch(outTree,"nJets40",&nJets40);
	setupBranch(outTree,"nJets50",&nJets50);
	setupBranch(outTree,"sumPtLep",&sumPtLep);
	setupBranch(outTree,"SumJetPt",&SumJetPt);
	setupBranch(outTree,"SumBJetPt",&SumBJetPt);
	setupBranch(outTree,"Pt_l",&Pt_l);
	setupBranch(outTree,"Pt_subl",&Pt_subl);
    setupBranch(outTree,"mljj_highestPt", &mljj_highestPt);
    setupBranch(outTree,"mljj_closest", &mljj_closest);
    setupBranch(outTree,"mljj_closestToW", &mljj_closestToW);
    setupBranch(outTree,"mlW_shifted", &mlW_shifted);
    setupBranch(outTree,"jetPt_1",&jetPt_1);
    setupBranch(outTree,"jetPt_2",&jetPt_2);
	/// setupBranch(outTree,"SSNegative",&SSNegative);
	setupBranch(outTree,"Zevent",&Zevent);
	setupBranch(outTree,"isSS30",&isSS30);
	setupBranch(outTree,"mt",&mt);
	setupBranch(outTree,"MT2",&MT2);
	setupBranch(outTree,"ht",&ht);
	setupBranch(outTree,"Mjj",&Mjj);
    setupBranch(outTree,"Mll",&Mll);
	setupBranch(outTree,"minvAllSigL",&minvAllSigL); /// all signal leptons
	setupBranch(outTree,"mSFOS",&mSFOS);
	setupBranch(outTree,"meff",&meff);
	setupBranch(outTree,"met",&met);
	/// setupBranch(outTree,"dRll", &dRll);
	setupBranch(outTree,"dRl1j",&dRl1j);
    setupBranch(outTree,"dRJJ",&dRJJ);
	/// setupBranch(outTree,"dRl2j",&dRl2j);
	setupBranch(outTree,"isZee",&isZee);
	setupBranch(outTree,"is3LSS",&is3LSS);
	setupBranch(outTree,"is3LSSproc",&is3LSSproc);
	setupBranch(outTree,"SSChannel",&SSChannel);
	setupBranch(outTree,"isSR",&isSR);


	setupBranch(outTree,"wmu_nom",&wmu_nom);	
	setupBranch(outTree,"wel_nom",&wel_nom);
	setupBranch(outTree,"wjet_nom",&wjet_nom);		
	if(DoSYS){
		setupBranch(outTree,"wmu_stat_up",&wmu_stat_up);
		setupBranch(outTree,"wmu_stat_down",&wmu_stat_down ); 
		setupBranch(outTree,"wmu_sys_up",&wmu_sys_up);		  
		setupBranch(outTree,"wmu_sys_down",&wmu_sys_down);
		setupBranch(outTree,"wmu_bad_sys_up", &wmu_bad_sys_up);
		setupBranch(outTree,"wmu_bad_sys_down", &wmu_bad_sys_down);
		setupBranch(outTree,"wmu_bad_stat_up", &wmu_bad_stat_up);
		setupBranch(outTree,"wmu_bad_stat_down", &wmu_bad_stat_down);
		setupBranch(outTree,"wmu_stat_lowpt_up",&wmu_stat_lowpt_up);
		setupBranch(outTree,"wmu_stat_lowpt_down",&wmu_stat_lowpt_down );
		setupBranch(outTree,"wmu_sys_lowpt_up",&wmu_sys_lowpt_up);
		setupBranch(outTree,"wmu_sys_lowpt_down",&wmu_sys_lowpt_down);
		setupBranch(outTree,"wmu_trig_stat_up",&wmu_trig_stat_up);	  
		setupBranch(outTree,"wmu_trig_stat_down",&wmu_trig_stat_down);	  
		setupBranch(outTree,"wmu_trig_sys_up",&wmu_trig_sys_up);	  
		setupBranch(outTree,"wmu_trig_sys_down",&wmu_trig_sys_down);	  
		setupBranch(outTree,"wmu_iso_stat_up",&wmu_iso_stat_up);	  
		setupBranch(outTree,"wmu_iso_stat_down",&wmu_iso_stat_down);	  
		setupBranch(outTree,"wmu_iso_sys_up",&wmu_iso_sys_up);	
		setupBranch(outTree,"wmu_iso_sys_down",&wmu_iso_sys_down);	  
		setupBranch(outTree,"wmu_ttva_stat_up",&wmu_ttva_stat_up);
		setupBranch(outTree,"wmu_ttva_stat_down",&wmu_ttva_stat_down);
		setupBranch(outTree,"wmu_ttva_sys_up",&wmu_ttva_sys_up);
		setupBranch(outTree,"wmu_ttva_sys_down",&wmu_ttva_sys_down);

		setupBranch(outTree,"wel_cid_up",&wel_cid_up);
		setupBranch(outTree,"wel_cid_down",&wel_cid_down);		  
		setupBranch(outTree,"wel_id_up",&wel_id_up);		  
		setupBranch(outTree,"wel_id_down",&wel_id_down);
		setupBranch(outTree,"wel_iso_up",&wel_iso_up);
		setupBranch(outTree,"wel_iso_down",&wel_iso_down);	  
		setupBranch(outTree,"wel_reco_up",&wel_reco_up);	  
		setupBranch(outTree,"wel_reco_down",&wel_reco_down); 
		setupBranch(outTree,"wel_trig_up",&wel_trig_up);	  
		setupBranch(outTree,"wel_trig_down",&wel_trig_down);  
		setupBranch(outTree,"wel_trigEff_up",&wel_trigEff_up);	  
		setupBranch(outTree,"wel_trigEff_down",&wel_trigEff_down);  
	
		setupBranch(outTree,"wjet_b_up",&wjet_b_up);		  
		setupBranch(outTree,"wjet_b_down",&wjet_b_down);
		setupBranch(outTree,"wjet_c_up",&wjet_c_up);		  
		setupBranch(outTree,"wjet_c_down",&wjet_c_down);
		setupBranch(outTree,"wjet_light_up",&wjet_light_up ); 
		setupBranch(outTree,"wjet_light_down",&wjet_light_down);
		setupBranch(outTree,"wjet_extra1_up",&wjet_extra1_up ); 
		setupBranch(outTree,"wjet_extra1_down",&wjet_extra1_down);
		setupBranch(outTree,"wjet_extra2_up",&wjet_extra2_up ); 
		setupBranch(outTree,"wjet_extra2_down",&wjet_extra2_down);
		setupBranch(outTree,"wjet_jvt_up",&wjet_jvt_up);
		setupBranch(outTree,"wjet_jvt_down",&wjet_jvt_down);
	}

	// to fix: not all are needed
	setupBranch(outTree,"wpu_nom_sig",&wpu_nom_sig);	
	setupBranch(outTree,"wpu_up_sig",&wpu_up_sig);		  
	setupBranch(outTree,"wpu_down_sig",&wpu_down_sig);
	setupBranch(outTree,"wpu_nom_bkg",&wpu_nom_bkg);
	setupBranch(outTree,"wpu_up_bkg",&wpu_up_bkg);
	setupBranch(outTree,"wpu_down_bkg",&wpu_down_bkg);

	setupBranch(outTree,"wtrig_nom",&wtrig_nom);
	if(DoSYS){ // to fix
		setupBranch(outTree,"wtrig_up",&wtrig_up);
		setupBranch(outTree,"wtrig_down",&wtrig_down);
	}
	
	setupBranch(outTree,"wchflip_nom",&wchflip_nom);
	if(DoSYS){
		setupBranch(outTree,"wchflip_up",&wchflip_up);
		setupBranch(outTree,"wchflip_down",&wchflip_down);

		///setupBranch(outTree,"wpdf_up",&wpdf_up);		  
		///setupBranch(outTree,"wpdf_down",&wpdf_down);
	}

	/// signals
	if(runc!=0){
		lumiScalingUP=lumiScaling*(1+runc);
		lumiScalingDOWN=lumiScaling*(1-runc);
		std::cout<< "got runc !=0 "<<runc<<" "<<lumiScaling<<" "<<lumiScalingUP<<" "<<lumiScalingDOWN<<std::endl;
	}
	setupBranch(outTree,"lumiScaling_UP",&lumiScalingUP);
	setupBranch(outTree,"lumiScaling_DOWN",&lumiScalingDOWN);


	//inTree->GetEntry(0);

	//std::cout<<"Testing weights:"<<" " <<SysNames->size()<<std::endl;
	//for(auto n: *SysNames)
	//  cout<<n<<endl;

	int nocut=0,bveto=0;
	/// Rescaling for missing samples
	/**if(MC_ID==363507 || MC_ID==363509){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16e SAMPLES"<<std::endl;
	  MC_campaign_weight = 1.7444;
	}*/
	/**if(MC_ID==345706){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16a SAMPLES"<<std::endl;
          MC_campaign_weight = 1.3473;
	  }*/
	/**if(MC_ID==410219){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16d SAMPLES"<<std::endl;
          MC_campaign_weight = 1.4608;
	  }*/
	
	/// if(MCId==412063) continue; /// Filtering out new tllq sample
	  

	for (Long64_t i=0;i<nentries; i++) 
	{
		inTree->GetEntry(i);
		
		nocut = nocut+1;
		isSR=isVR=false;		
		if(MCId>=372444 && MCId<=372511 ){ if(GLDec1==5 || GLDec2==5){continue;} }
		bveto = bveto+1;
		
		///writtenEvents.push_back(evn64);
		///if(evn64!=1780910)
		///continue;
		///std::cout<<"EVENT!"<<std::endl;
		///FIXING DUPLICATING LEP/JET ENTRIES
		
		
		/// Lep BL ==> NlepBL

		/// SIG LEP
		int isGoodLep = 1;
		tmp_lepTLV.clear();
		tmp_lepCharges.clear();
		tmp_lepTruth.clear();
        lepPt.clear();
        lepPhi.clear();
        lepEta.clear();
        lepM.clear();
        
		for(unsigned iLep=0;iLep<lepTLV->size();iLep++){
		  isGoodLep = 1;
		  for(unsigned kLep=iLep+1;kLep<lepTLV->size();kLep++){
		    if(lepTLV->at(iLep).Pt()==lepTLV->at(kLep).Pt() && lepTLV->at(iLep).Eta()==lepTLV->at(kLep).Eta()){
               isGoodLep = 0;
			}
          }
		  if(isGoodLep==1){
            tmp_lepTLV.push_back(lepTLV->at(iLep));
		    tmp_lepCharges.push_back(lepCharges->at(iLep));
		    tmp_lepTruth.push_back(lepTruth->at(iLep));
            
            lepPt.push_back(lepTLV->at(iLep).Pt());
            lepPhi.push_back(lepTLV->at(iLep).Phi());
            lepEta.push_back(lepTLV->at(iLep).Eta());
            lepM.push_back(lepTLV->at(iLep).M());
            
          }
        }
		nSigLep = tmp_lepTLV.size();

		/// SIG JETS
		tmp_jetTLV.clear();
		tmp_jetBtag.clear();
        
        jetPt.clear();
        jetPhi.clear();
        jetEta.clear();
        jetM.clear();
        
		int isGoodJet = 1;
		for(unsigned iJet=0;iJet<jetTLV->size();iJet++){
		  isGoodJet = 1;
		  for(unsigned kJet=iJet+1;kJet<jetTLV->size();kJet++){
		    if(jetTLV->at(iJet).Pt()==jetTLV->at(kJet).Pt() && jetTLV->at(iJet).Eta()==jetTLV->at(kJet).Eta() && jetBtag->at(iJet)==jetBtag->at(kJet)){
		      isGoodJet = 0;}
		  }
		  if(isGoodJet==1){
		    tmp_jetTLV.push_back(jetTLV->at(iJet));
		    tmp_jetBtag.push_back(jetBtag->at(iJet));
            
            jetPt.push_back(jetTLV->at(iJet).Pt());
            jetPhi.push_back(jetTLV->at(iJet).Phi());
            jetEta.push_back(jetTLV->at(iJet).Eta());
            jetM.push_back(jetTLV->at(iJet).M());
		  }
		}
		nJets25 = 0;
		nJets35 = 0;
		nJets40 = 0;
		nJets50 = 0;
		nBJets20= 0;
		for(unsigned iJet=0;iJet<tmp_jetTLV.size();iJet++){
		  if(tmp_jetTLV.at(iJet).Pt()>25000) nJets25++;
		  if(tmp_jetTLV.at(iJet).Pt()>35000) nJets35++;
		  if(tmp_jetTLV.at(iJet).Pt()>40000) nJets40++;
		  if(tmp_jetTLV.at(iJet).Pt()>50000) nJets50++;
		  if(tmp_jetTLV.at(iJet).Pt()>20000 && tmp_jetBtag.at(iJet)>0) nBJets20++;
		}

		/// Lepton pT for Validation Regions 
		Pt_l = tmp_lepTLV.at(0).Pt();
		Pt_subl = tmp_lepTLV.at(1).Pt();

		/// sumPt
		sumPtLep  = getSumPtLep(tmp_lepTLV);
		minvAllSigL = getSumMinvLep(tmp_lepTLV);
		SumJetPt  = getSumPtJet(tmp_jetTLV,tmp_jetBtag,"ALL");
		SumBJetPt = getSumPtJet(tmp_jetTLV,tmp_jetBtag,"B");
		meff = getSumPtJet_noPtCut(tmp_jetTLV) + sumPtLep + met;
		Mjj = tmp_jetTLV.size()>1 ? (tmp_jetTLV.at(0) + tmp_jetTLV.at(1)).M() : -99.0;
        Mll = tmp_lepTLV.size()>1 ? (tmp_lepTLV.at(0) + tmp_lepTLV.at(1)).M() : -99.0;
		mt = getMt(tmp_lepTLV,met,phi_met);
		ht = getSumPtJet_noPtCut(tmp_jetTLV);
        
        /// jet pt
        if(tmp_jetTLV.size()>1){
            jetPt_1 = tmp_jetTLV.at(0).Pt();
            jetPt_2 = tmp_jetTLV.at(1).Pt();
        }
        else if(tmp_jetTLV.size()>0){
            jetPt_1 = tmp_jetTLV.at(0).Pt();
            jetPt_2 = -999;
        }
        else{
            jetPt_1 = -999;
            jetPt_2 = -999;
        }
		/// dR
		dRll  = getDeltaRLepLep(tmp_lepTLV);
		dRl1j = getDeltaRLepJet(tmp_lepTLV, tmp_jetTLV, 0);
        
        dRJJ  = getDeltaRJJ(tmp_jetTLV);
		/// dRl2j = getDeltaRLepJet(tmp_lepTLV, jetTLV, 1);
        /// variables for Judita
        if(jetPt.size()>1){
            mljj_highestPt    = cmljj_highestPt(lepPt, lepPhi, lepEta, lepM, jetPt, jetPhi, jetEta, jetM);
            mljj_closest      = cmljj_closest(lepPt, lepPhi, lepEta, lepM, jetPt, jetPhi, jetEta, jetM); 
            mljj_closestToW   = cmljj_closestToW(lepPt, lepPhi, lepEta, lepM, jetPt, jetPhi, jetEta, jetM); 
            mlW_shifted       = cmljj_shiftToW(lepPt, lepPhi, lepEta, lepM, jetPt, jetPhi, jetEta, jetM);
        }
        else{
            mljj_highestPt  = -999.;
            mljj_closest    = -999.;
            mljj_closestToW = -999.;
            mlW_shifted     = -999.;
        }
        
		/// Rescaling for missing samples
		TotWeightNew = MC_campaign_weight*totweight;

		/// MT2
		TLorentzVector metTLV = TLorentzVector( met*TMath::Cos(phi_met) , met*TMath::Sin(phi_met) , 0 , met);
        if(tmp_lepTLV.size()>1){
            ComputeMT2 mycalc = ComputeMT2(tmp_lepTLV.at(0),tmp_lepTLV.at(1),metTLV,tmp_lepTLV.at(0).M(),tmp_lepTLV.at(1).M());
            MT2 = mycalc.Compute();
        }
        else
            MT2 = -999;

		/// take charge related variables
		TString SSPos;
		/// SSNegative = getSSNegative(tmp_lepCharges);
		if(tmp_lepCharges.size()<3)  is3LSS=0;
		else is3LSS = has3LSSPrompt(tmp_lepCharges,tmp_lepTruth); 
		is3LSSproc = is3LSSprocess(MCId); 

		/// new variables for VR
		if(tmp_lepTLV.size()>=2){
		  mSFOS=GetSFOSmass(tmp_lepTLV, tmp_lepCharges);
		  isSS30 = isSSLepPtCut(tmp_lepTLV, tmp_lepCharges);
		}
		else{
		  mSFOS=0;
		  isSS30=0;
		}
		isZee=hasZee(tmp_lepTLV,tmp_lepCharges);
		Zevent = isResonance(tmp_lepTLV, tmp_lepCharges, SSChannel);

		/// MC
		/// https://twiki.cern.ch/twiki/bin/view/AtlasProtected/SUSYSystematicUncertaintiesRun2
		if(xsec!=0 && SysNames->size()!=0){
			if(DEBUG) cout<<"--------------------------------"<<endl;
			if(DEBUG) cout << "Setting "<<SysNames->at(0)<<" to wmu_nom"<<endl;
			wmu_nom=SysSF->at(0);
			if(DoSYS){
				if(DEBUG) cout << "Setting "<<SysNames->at(29)<<" to wmu_bad_stat_up"<<endl;
				wmu_bad_stat_up=SysSF->at(29);
				if(DEBUG) cout << "Setting "<<SysNames->at(28)<<" to wmu_bad_stat_down"<<endl;
				wmu_bad_stat_down=SysSF->at(28);
				if(DEBUG) cout << "Setting "<<SysNames->at(31)<<" to wmu_bad_sys_up"<<endl;
				wmu_bad_sys_up=SysSF->at(31);
				if(DEBUG) cout << "Setting "<<SysNames->at(30)<<" to wmu_bad_sys_down"<<endl;
				wmu_bad_sys_down=SysSF->at(30);
				if(DEBUG) cout << "Setting "<<SysNames->at(33)<<" to wmu_stat_up"<<endl;
				wmu_stat_up=SysSF->at(33);
				if(DEBUG) cout <<"Setting "<<SysNames->at(32)<<" to wmu_stat_down"<<endl;
				wmu_stat_down=SysSF->at(32);
				if(DEBUG) cout <<"Setting "<<SysNames->at(35)<<" to wmu_stat_lowpt_up"<<endl;
				wmu_stat_lowpt_up=SysSF->at(35);
				if(DEBUG) cout <<"Setting "<<SysNames->at(34)<<" to wmu_stat_lowpt_down"<<endl;
				wmu_stat_lowpt_down=SysSF->at(34);
				if(DEBUG) cout <<"Setting "<<SysNames->at(37)<<" to wmu_sys_up"<<endl;
				wmu_sys_up=SysSF->at(37);
				if(DEBUG) cout <<"Setting "<<SysNames->at(36)<<" to wmu_sys_down"<<endl;
				wmu_sys_down=SysSF->at(36);
				if(DEBUG) cout <<"Setting "<<SysNames->at(39)<<" to wmu_sys_lowpt_up"<<endl;
				wmu_sys_lowpt_up=SysSF->at(39);
				if(DEBUG) cout <<"Setting "<<SysNames->at(38)<<" to wmu_sys_lowpt_down"<<endl;
				wmu_sys_lowpt_down=SysSF->at(38);
				if(DEBUG) cout <<"Setting "<<SysNames->at(41)<<" to wmu_trig_stat_up"<<endl;
				wmu_trig_stat_up=SysSF->at(41);
				if(DEBUG) cout <<"Setting "<<SysNames->at(40)<<" to wmu_trig_stat_down"<<endl;
				wmu_trig_stat_down=SysSF->at(40);
				if(DEBUG) cout <<"Setting "<<SysNames->at(43)<<" to wmu_trig_sys_up"<<endl;
				wmu_trig_sys_up=SysSF->at(43);
				if(DEBUG) cout <<"Setting "<<SysNames->at(42)<<" to wmu_trig_sys_down"<<endl;
				wmu_trig_sys_down=SysSF->at(42);
				if(DEBUG) cout <<"Setting "<<SysNames->at(45)<<" to wmu_iso_stat_up"<<endl;
				wmu_iso_stat_up=SysSF->at(45);
				if(DEBUG) cout <<"Setting "<<SysNames->at(44)<<" to wmu_iso_stat_down"<<endl;
				wmu_iso_stat_down=SysSF->at(44);
				if(DEBUG) cout <<"Setting "<<SysNames->at(47)<<" to wmu_iso_sys_up"<<endl;
				wmu_iso_sys_up=SysSF->at(47);
				if(DEBUG) cout <<"Setting "<<SysNames->at(46)<<" to wmu_iso_sys_down"<<endl;
				wmu_iso_sys_down=SysSF->at(46);
				if(DEBUG) cout <<"Setting "<<SysNames->at(49)<<" to wmu_ttva_stat_up"<<endl;
				wmu_ttva_stat_up=SysSF->at(49);
				if(DEBUG) cout <<"Setting "<<SysNames->at(48)<<" to wmu_ttva_stat_down"<<endl;
				wmu_ttva_stat_down=SysSF->at(48);
				if(DEBUG) cout <<"Setting "<<SysNames->at(51)<<" to wmu_ttva_sys_up"<<endl;
				wmu_ttva_sys_up=SysSF->at(51);
				if(DEBUG) cout <<"Setting "<<SysNames->at(50)<<" to wmu_ttva_sys_down"<<endl;
				wmu_ttva_sys_down=SysSF->at(50);
			}
			if(DEBUG) cout <<"Setting "<<SysNames->at(1)<<" to wel_nom"<<endl;
			wel_nom=SysSF->at(1);
			if(DoSYS){
				if(DEBUG) cout <<"Setting "<<SysNames->at(5)<<" to wel_cid_up"<<endl;
				wel_cid_up=SysSF->at(5);
				if(DEBUG) cout <<"Setting "<<SysNames->at(4)<<" to wel_cid_down"<<endl;
				wel_cid_down=SysSF->at(4);
				if(DEBUG) cout <<"Setting "<<SysNames->at(7)<<" to wel_id_up"<<endl;
				if(DEBUG) cout <<"Setting "<<SysNames->at(6)<<" to wel_id_down"<<endl;
				elID_rescale_up=1.;
				elID_rescale_down=1.;
				for(unsigned int ilep=0;ilep<tmp_lepTLV.size();ilep++){ // this is for both electrons and muons
				  if (getDeltaRLepJet(tmp_lepTLV,tmp_jetTLV,ilep)<0.4){
					elID_rescale_up=1.2;
					elID_rescale_down=0.8;
				  }
				}
				if(SysSF->at(7)>SysSF->at(6)){
				  wel_id_up=SysSF->at(7)*elID_rescale_up;
				  wel_id_down=SysSF->at(6)*elID_rescale_down;
				}
				else
				{
				  wel_id_up=SysSF->at(7)*elID_rescale_down;
				  wel_id_down=SysSF->at(6)*elID_rescale_up;
				}
				if(DEBUG) cout <<"Setting "<<SysNames->at(9)<<" to wel_iso_up"<<endl;
				wel_iso_up = SysSF->at(9);
				if(DEBUG) cout <<"Setting "<<SysNames->at(8)<<" to wel_iso_down"<<endl;
				wel_iso_down = SysSF->at(8);
				if(DEBUG) cout <<"Setting "<<SysNames->at(11)<<" to wel_reco_up"<<endl;
				wel_reco_up=SysSF->at(11);
				if(DEBUG) cout <<"Setting "<<SysNames->at(10)<<" to wel_reco_down"<<endl;
				wel_reco_down=SysSF->at(10);
				if(DEBUG) cout <<"Setting "<<SysNames->at(13)<<" to wel_trigEff_up"<<endl;
				wel_trigEff_up=SysSF->at(13);
				if(DEBUG) cout <<"Setting "<<SysNames->at(12)<<" to wel_trigEff_down"<<endl;
				wel_trigEff_down=SysSF->at(12);                 
				if(DEBUG) cout <<"Setting "<<SysNames->at(15)<<" to wel_trig_up"<<endl;
				wel_trig_up=SysSF->at(15);
				if(DEBUG) cout <<"Setting "<<SysNames->at(14)<<" to wel_trig_down"<<endl;
				wel_trig_down=SysSF->at(14);
			}
			if(DEBUG) cout <<"Setting "<<SysNames->at(2)<<" to wjet_nom"<<endl;
			wjet_nom=SysSF->at(2);
			if(DoSYS){
				if(DEBUG) cout <<"Setting "<<SysNames->at(17)<<" to wjet_b_up"<<endl;
				wjet_b_up=SysSF->at(17);
				if(DEBUG) cout <<"Setting "<<SysNames->at(16)<<" to wjet_b_down"<<endl;
				wjet_b_down=SysSF->at(16);
				if(DEBUG) cout <<"Setting "<<SysNames->at(19)<<" to wjet_c_up"<<endl;
				wjet_c_up=SysSF->at(19);
				if(DEBUG) cout <<"Setting "<<SysNames->at(18)<<" to wjet_c_down"<<endl;
				wjet_c_down=SysSF->at(18);
				if(DEBUG) cout <<"Setting "<<SysNames->at(21)<<" to wjet_light_up"<<endl;
				wjet_light_up=SysSF->at(21);
				if(DEBUG) cout <<"Setting "<<SysNames->at(20)<<" to wjet_light_down"<<endl;
				wjet_light_down=SysSF->at(20);
				if(DEBUG) cout <<"Setting "<<SysNames->at(23)<<" to wjet_extra1_up"<<endl;
				wjet_extra1_up=SysSF->at(23);
				if(DEBUG) cout <<"Setting "<<SysNames->at(22)<<" to wjet_extra1_down"<<endl;
				wjet_extra1_down=SysSF->at(22);
				if(DEBUG) cout <<"Setting "<<SysNames->at(25)<<" to wjet_extra2_up"<<endl;
				wjet_extra2_up=SysSF->at(25);
				if(DEBUG) cout <<"Setting "<<SysNames->at(24)<<" to wjet_extra2_down"<<endl;
				wjet_extra2_down=SysSF->at(24);
				if(DEBUG) cout <<"Setting "<<SysNames->at(27)<<" to wjet_jvt_up"<<endl;
				wjet_jvt_up=SysSF->at(27);
				if(DEBUG) cout <<"Setting "<<SysNames->at(26)<<" to wjet_jvt_down"<<endl;
				wjet_jvt_down=SysSF->at(26);
			}

			if(!isSignal(MCId)){
				wpu_nom_sig=1;
				wpu_up_sig=1;
				wpu_down_sig=1;
			}
			else{
				wpu_nom_sig=SysSF->at(3);
				wpu_up_sig=SysSF->at(47);
				wpu_down_sig=SysSF->at(46);
			}
			if(isSignal(MCId)){
				wpu_nom_bkg=1;
				wpu_up_bkg=1;
				wpu_down_bkg=1;
			}
			else{
				wpu_nom_bkg=SysSF->at(3);
				wpu_up_bkg=SysSF->at(47);
				wpu_down_bkg=SysSF->at(46);
			}

			wtrig_nom=TrigSF;//NEED TO PUT BACK IN THE UP AND DOWN VARIATIONS!
			if(DoSYS){
				wtrig_up=TrigSF;
				wtrig_down=TrigSF;
			}

			//std::vector<float> ChFlip_SF = correctChFlipSF(ChFlipSF, lepTLV, true); TBC!!!
			/*wchflip_nom=ChFlip_SF.at(0);
			wchflip_up=ChFlip_SF.at(1);
			wchflip_down=ChFlip_SF.at(2);*/
			wchflip_nom=ChFlipSF->at(0); /// see why this one is 0!			
			if(DoSYS){
				wchflip_up=ChFlipSF->at(1);
				wchflip_down=ChFlipSF->at(2);

				///wpdf_up=1;
				///wpdf_down=1;
			}
		}
		
		/// Rpc2L1b
		if(nSigLep>=2 && nBJets20>=1 && nJets40>=6 && met/meff>0.25) isSR=true;
		/// Rpc2L2b
		if(nSigLep>=2 && nBJets20>=2 && nJets25>=6 && met>300000 && meff>1400000 && met/meff>0.14) isSR=true;
		/// Rpc2L0b
		if(nSigLep>=2 && nBJets20==0 && nJets40>=6 && met>200000 && meff>1000000 && met/meff>0.2) isSR=true;
		/// Rpc3LSS1b
		if(nSigLep>=3 && nBJets20>=1 && is3LSS>0 && !isZee && met/meff>0.14)  isSR=true; 
		/// Rpv2L
		if(nSigLep>=2 && nJets40>=6 && meff>2600000) isSR=true;
		/// VRWZ4j
		if(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=4 && mSFOS>81000 && mSFOS<101000 && met>50000 && met<250000 && meff>600000 && meff<1500000 ){
		  isVR=true;
		}
		/// VRWZ5j
		if(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=5 && meff>400000 && met>50000 && mSFOS>81000 && mSFOS<101000 && meff<1500000 && met<250000){
		  isVR=true;
		  }
		/// VRttV
		if(nSigLep>=2 && NlepBL>=2 && nBJets20>=1 && nJets40>=3 && meff>600000 && meff<1500000 && met<250000 && isSS30 && dRl1j>1.1 && SumBJetPt/SumJetPt>0.4 && met/meff>0.1){
		  isVR=true;
		}


		bool cout_info = false;  
		if(cout_info){
			std::cout<<"evn: "<<evn64<<std::endl;
			std::cout<<"nSigLep: "<<nSigLep<<"\t NlepBL: "<<NlepBL<<"\t nBJets20: "<<nBJets20<<"\t nJets25: "<<nJets25<<"\t nJets40: "<<nJets40
					<<"\t meff: "<<meff<<"\t met: "<<met<<"\t isSS30: "<<isSS30<<"\t dRl1j: "<<dRl1j<<"\t SumBJetPt/SumJetPt: "<<SumBJetPt/SumJetPt
					<<"\t met/meff: "<<met/meff<<"\t mSFOS: "<<mSFOS<<"\t is3LSS: "<<is3LSS<<"\t isZee: "<<isZee<<std::endl;
		}
		/// fill the tree only for the SRs and VRs
		if(mySkimming==SR_ONLY_SKIMMING && isSR==false) {if(i<3)cout << "mySkimming==SR_ONLY_SKIMMING\n"; continue;}
		if(mySkimming==VR_ONLY_SKIMMING && isVR==false) {if(i<3)cout << "mySkimming==VR_ONLY_SKIMMING\n"; continue;}
		if(mySkimming==SRVR_ONLY_SKIMMING && !(isVR || isSR)) {cout << "mySkimming==SRVR_ONLY_SKIMMING\n"; continue;}
		if(mySkimming==THREELEP_ONLY_SKIMMING && !(nSigLep==3 && NlepBL==3)) {if(i<3)cout << "mySkimming==THREELEP_ONLY_SKIMMING\n"; continue;} /// for the CRs studies -- for data & MC
		
		/// fill the tree now:
		outTree->Fill();
		if(cout_info) std::cout<<"TAKEN!"<<std::endl;
		writtenEvents.push_back(evn64);
	}
	cout<<"No cut: "<<nocut<<"\t Bveto: "<<bveto<<endl;

	outTree->Write("",TObject::kOverwrite);
	cout<<"Done.."<<endl;
}


/// open a file, get the tree and copy it to the output. attach proper lumi rescaling
void processFile(TString filein, TString fileout, TString treein, TString treeout, float xsec,float runc, Long64_t MC_ID, float lumi,  std::vector<Long64_t>& writtenEvents,bool DoSYS,SKIMMING mySkimming)
{
    std::cout<<"*** PROCESS FILE ***"<<std::endl;
	TFile* inFile=TFile::Open(filein, "READ");
	if(inFile->IsZombie() || !inFile)  std::cout<<"FILE IS IN ZOMBIE MODE!!!"<<std::endl;
	else  std::cout<<"inFile is open"<<std::endl;
	/// need to remove "_nom" for input nominal sample
	TTree* inTree=(TTree*)inFile->Get(treein.ReplaceAll("_nom",""));
	if(!inTree){
		std::cout<<"Skipping corrupted/empty file (_nom do not exist) "<<std::endl;
		inFile->Close();
		return;
	}
	else std::cout<<"inTree is open"<<std::endl;


	TFile* outFile=TFile::Open(fileout,"UPDATE");
	if(outFile->IsZombie() || !outFile) std::cout<<"FILE IS IN ZOMBIE MODE!!!"<<std::endl;
	else std::cout<<"outFile is open"<<std::endl;
	TTree* outTree=(TTree*)gROOT->FindObject(treeout);
	if(!outTree){
		cout<<"output tree not found, creating"<<endl;
		outTree=new TTree(treeout,treeout);
	}
	else cout<<RED<<"adding entries to existing tree "<<treeout<<RESET<<endl;

	TTree* control=(TTree*)inFile->Get("ControlTree");


	/// get sum of weights. simple method using TTree::Draw and TArrayD crashes for large number of events
	Int_t NFiles=0, fileentry=0,raw=0,nofilt=0;
	float xAODWeightSum=0;
	float totw=0;
	float MCWeight=0;
	int MCId=0;

	/// get filter eff for MC (data has xsec=0)
	if(xsec){
		control->SetBranchAddress("FileEntry",&fileentry);
		control->SetBranchAddress("xAODWeightSum",&xAODWeightSum);
		control->SetBranchAddress("NFiles",&NFiles);
		control->SetBranchAddress("MCId",&MCId);    
		control->SetBranchAddress("MCWeight",&MCWeight);

		for(int i=0;i< control->GetEntries();++i){
			control->GetEntry(i);
			if(fileentry==1) totw += xAODWeightSum;  
		}
		/// if((MCId>=372444 && MCId<=372524) || (MCId>=373462 && MCId<=373478)) totw = totw*0.64; /// b-filtering for Slepton grid
	}
	else totw=1;

	cout<<"Xsec: "<<xsec<<endl;
	copyTree(inTree,outTree,xsec,runc,MC_ID,totw,lumi,writtenEvents,DoSYS,mySkimming);

	outFile->Write("",TObject::kOverwrite);

	outFile->Close();
	inFile->Close();

}


/// functions to retrieve cross section from common SUSYTools/PMG files
/// this is just a grep. returns the line(s) to be parsed by the next function
TString lookForSample(TString runnumber)
{
	/// uses perl regexp. matches trailing whitespace (space OR tab). requires begin of line.
	/// take the PMG file with the x-sect, k-factors, filter eff, etc
	/// gSystem->Exec("grep -r -P \"^"+runnumber+"\\s\" /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PMGTools/PMGxsecDB_mc16.txt  > xSec.out");
    // to uncomment of the PMG file was not saved in this dir before!
	// gSystem->Exec("cp /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PMGTools/PMGxsecDB_mc16.txt "+basepath+"/.");
	
    ///////// This takes the xsec for signals with preferable etags
    /*
    Long64_t runNr = atol(runnumber);
	TString  etag = "";
    
	if(isSignal(runNr)){
		if((runNr>=376733 && runNr<=376775) || (runNr>=376835 && runNr<=376843)) etag = "e7996"; /// GG_Rpvlampp331 etag = "e7359" if samples produced by arka, matching bug on
		if(runNr>=377000 && runNr<=377071) etag = "e7153"; /// GG_RpvLQD
		if((runNr>=377605 && runNr<=377640) || runNr==376122) etag = "e7078"; /// GTT1 
		if(runNr>=375818 && runNr<=375932) etag = "e6351"; /// GTT2
		if(runNr>=377094 && runNr<=377102) etag = "e7026"; /// TT2STEP
		if(runNr>=449781 && runNr<=449822) etag = "e7585"; /// bRPV
		if(runNr>=377408 && runNr<=377549) etag = "e7126"; /// Wino offshell WZ, part 1
		if(runNr>=377578 && runNr<=377593) etag = "e7505"; /// Wino offshell WZ 
		if((runNr>=392200 && runNr<=392280) || (runNr>=394776 && runNr<=394779)) etag = "e6777"; /// Wino onshell WZ, part 1 
		if(runNr>=396206 && runNr<=396208) etag = "e7041"; /// Wino onshell WZ, part 2
		if(runNr>=397285 && runNr<=397314) etag = "e7167"; /// Wino onshell WZ, part 3
		if(runNr==398177) etag = "e7437"; /// Wino onshell WZ, part 4
  
  
        gSystem->Exec("grep -r -P \"^"+runnumber+"\\s\" /lustre/ific.uv.es/grid/atlas/t3/asantra/SS3L_Run2/ss3l_hf/prepare/xsections//PMGxsecDB_mc16.txt > xSec1.out");
		/// when checking with root: gSystem->Exec("grep -r -P \"^376749\\s\" "+ basepath+"PMGxsecDB_mc16.txt > xSec1.out");
		gSystem->Exec("tail xSec1.out");
		gSystem->Exec("grep -r -P \"\\s"+etag+"\\s*$\" xSec1.out > xSec.out");
		/// when checking with root: gSystem->Exec("grep -r -P \"\\se7359\\s*$\" xSec1.out");
	}
	
	else
	{
        gSystem->Exec("grep -m 1 -r -P \"^"+runnumber+"\\s\" /lustre/ific.uv.es/grid/atlas/t3/asantra/SS3L_Run2/ss3l_hf/prepare/xsections//PMGxsecDB_mc16.txt > xSec.out");
	}
	*/
    
    //// takes only the first occurrence of the DSID
    gSystem->Exec("grep -m 1 -r -P \"^"+runnumber+"\\s\" /lustre/ific.uv.es/grid/atlas/t3/asantra/SS3L_Run2/ss3l_hf/prepare/xsections/PMGxsecDB_mc16.txt > xSec.out");
	/// gSystem->Exec("ls");
		
	std::string line; 
	gSystem->Exec("tail xSec.out");
	std::ifstream result("xSec.out");
	std::getline(result, line);
	if(result.eof()){
		cout<<RED<<"ERROR cannot find cross section for sample "<<runnumber<<" will skip it "<<RESET<<endl;
		return "";
	}

	TString linestr(line);
	//gSystem->Exec("rm xSec1.out");
	gSystem->Exec("rm xSec.out");

	/// sanity check: there must be only one possibility
	std::getline(result, line);
	if(!result.eof()){
		cout<< "ERROR retrieving cross section: too many possibilities found. Aborting"<<endl;
		exit(1);
	}

	return linestr.ReplaceAll("\t"," ");
}

/// parse the result from grep, check that there is only one possibility
void getNormFactor(TString runnumber, float& norm, float& runc, Long64_t &MC_ID)
{
	TString grepped=lookForSample(runnumber);

	if(grepped=="") return;

	TObjArray* tokens= grepped.Tokenize(" ");
	Float_t xsec= ((TObjString*)(tokens->At(2)))->String().Atof();
	Float_t kfac= ((TObjString*)(tokens->At(4)))->String().Atof();
	Float_t feff= ((TObjString*)(tokens->At(3)))->String().Atof();

	MC_ID =  ((TObjString*)(tokens->At(0)))->String().Atof();
	runc = ((TObjString*)(tokens->At(5)))->String().Atof();

	norm = xsec*kfac*feff;

}


/// reads all files in a folder; adds the contents to outtreename in outfname
void readFolder(TString dirname, TString outfname, TString intreename, TString outtreename, std::vector<string>& sys, 
	float targetlumi_16a, float targetlumi_16d, float targetlumi_16e, bool DoSYS, SKIMMING mySkimming)
{
	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	cout << "files = " << files << endl;
	if (files) {
		bool isData=outtreename.Contains("data");
		bool isSignal=outtreename.BeginsWith("signal");
		cout << "isData = " << isData << ", isSignal = " << isSignal << endl;

		/// outer loop is on systematics
		std::cout<<"Systematics: "<<sys.size()<<std::endl;
		for(auto sysname: sys)
		{	
			std::vector<Long64_t> writtenEvents;
			TString _intreename=intreename+"_"+sysname;
			TSystemFile *file;
			TString fname;
			TIter next(files);   
			
			/// reinit the file iterator
			while ((file=(TSystemFile*)next())) 
			{
				fname = file->GetName();
				if(fname=="." || fname=="..") continue;

				TObjArray* temp_tokens=fname.Tokenize(".");
				TString temp_runnumber="";
				temp_runnumber=((TObjString*)(temp_tokens->At(2)))->String();
				
				TObjArray* tokens=temp_runnumber.Tokenize("_");

				if(tokens->GetEntries()>1 || isData){

					cout<<"processing file "<<fname<< " for systematic "<<sysname<<endl;
					TString runnumber="";
					float norm=-1,runc=-1;
					Long64_t MC_ID=-1;
					if(! isData){
						runnumber=((TObjString*)(tokens->At(1)))->String();
						getNormFactor(runnumber,norm,runc,MC_ID);
						/// cout<<" got cross normalization "<<norm<<endl;
					} else {
						/// signal "data" to downstream functions using norm=0
						norm=0; runc=0; MC_ID=0;
					}

					Long64_t mcid = atol(runnumber);
					bool Atlfast = isAtlfast(mcid);

					/// local (per file, per syst) output tree name
					TString _outtreename="signal"+runnumber;
					if(!isSignal){
						/// error on bkg MC treated differently
						runc=0;
						_outtreename=outtreename;
					} 
					else {
						/// std::cout<<"outtree name is "<<outtreename<<std::endl;
						/// signal tree names are different, i.e. they are like signalXXXX_nom
						TObjArray* tokens=outtreename.Tokenize("_");
						if(tokens->GetEntries()==2){
							_outtreename="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
						}
					}

					///cout<<"run "<<runnumber<<" got norm factors "<<norm<<" "<<runc<<endl;
					///Setting correct lumi
					float tLumi = 0;
					///**** MC16A ****//
					if(fname.Contains("_r9364_")){
					  std::cout<<"taking MC16a Sample: "<<targetlumi_16a<<std::endl;
					  tLumi = targetlumi_16a;
					}
					///**** MC16D ****//
					if(fname.Contains("_r10201_")){
					  std::cout<<"taking MC16d Sample: "<<targetlumi_16d<<std::endl;
                      tLumi= targetlumi_16d;
					}
					///**** MC16E ****//
					if(fname.Contains("_r10724_")){
					  std::cout<<"taking MC16e Sample: "<<targetlumi_16e<<std::endl;
                      tLumi= targetlumi_16e;
                    }
					std::cout<<"tLumi: "<<tLumi<<std::endl;

					/// sanity check
					if( norm!=-1 && runc!=-1 ){
					  _outtreename=_outtreename+"_"+TString(sysname);
					  processFile(dirname+"/"+fname, outfname,_intreename,_outtreename,norm,runc,MC_ID,tLumi,writtenEvents,DoSYS,mySkimming);
					} /// end sanity check

					bool doAFII(true);
					/// AFII for Fullsim samples
					if(!isData && DoSYS && !Atlfast && sysname=="nom" && doAFII){
						/// create dummy AFII syst for fullsim samples, copying the nominal one
						cout<<"The sample is fullsim, ading AFII trees..."<<endl;
						TString sysname2("EG_SCALE_AF2__1down");
						TString sysname3("EG_SCALE_AF2__1up");
						TString _outtreename2="signal"+runnumber;
						TString _outtreename3="signal"+runnumber;

						if(!isSignal){
							/// error on bg MC treated differently
							runc=0;
							_outtreename2=outtreename;
							_outtreename3=outtreename;
						} 
						else {

							if(tokens->GetEntries()==2){
								_outtreename2="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
								_outtreename3="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
							}
						}

						_outtreename2=_outtreename2+"_"+TString(sysname2);
						processFile(dirname+"/"+fname, outfname,_intreename,_outtreename2,norm,runc,MC_ID,tLumi,writtenEvents,DoSYS,mySkimming);
						_outtreename3=_outtreename3+"_"+TString(sysname3);
						processFile(dirname+"/"+fname, outfname,_intreename,_outtreename3,norm,runc,MC_ID,tLumi,writtenEvents,DoSYS,mySkimming);
					}
				}

				/// for MC, clear event record once per file/sample
				if(!isData)
					writtenEvents.clear();
			} /// loop on files in folder
			/// cout << "done: loop on files in folder"<<endl;
		} /// loop on systematics
		/// cout << "done: loop on systematics"<<endl;
	} /// loop over list of files
	cout << "done: loop over list of files\n\n"<<endl;
}



/// main script
#include <limits>   

void lumiRescale(int inLumi, int outLumi)
{
	TString inbg=TString("background.")+Form("%d", inLumi)+TString(".root");
	TString outbg=TString("background.")+Form("%d", outLumi)+TString(".root");
	TString insig=TString("signal.")+Form("%d", inLumi)+TString(".root");
	TString outsig=TString("signal.")+Form("%d", outLumi)+TString(".root");

	TFile* finbg=TFile::Open(inbg,"READ");
	TIter next(finbg->GetListOfKeys());
	TKey *key;

	while ((key = (TKey*)next())) 
	{
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TTree")) continue;

		TTree* inTree=(TTree*)key->ReadObj();
		cout<<"Read tree "<<inTree->GetName()<<endl;

		TFile* foutbg=TFile::Open(outbg,"UPDATE");
		foutbg->cd();

		TTree* outTree=new TTree(inTree->GetName()+TString("copy"),inTree->GetTitle());
		outTree->Write(inTree->GetName(),TObject::kOverwrite);
		cout<<"Wrote tree "<<inTree->GetName()<<endl;

		foutbg->Write("",TObject::kOverwrite);
		foutbg->Close();

	}


	TFile* finsig=TFile::Open(insig,"READ");
	TIter nextsig(finbg->GetListOfKeys());
	TKey *keysig;

	while ((keysig = (TKey*)nextsig())) {
		TClass *cl = gROOT->GetClass(keysig->GetClassName());
		if (!cl->InheritsFrom("TTree")) continue;

		TTree* inTree=(TTree*)keysig->ReadObj();
		cout<<"Read tree "<<inTree->GetName()<<endl;

		TFile* foutsig=TFile::Open(outsig,"UPDATE");
		foutsig->cd();
		
		TTree* outTree=new TTree(inTree->GetName()+TString("copy"),inTree->GetTitle());
		///copyTree(inTree,outTree,float(outLumi)/float(inLumi));    
		outTree->Write(inTree->GetName(),TObject::kOverwrite);
		cout<<"Wrote tree "<<inTree->GetName()<<endl;

		foutsig->Write("",TObject::kOverwrite);
		foutsig->Close();

	}

}


void doAll(float targetlumi_16a, float targetlumi_16d, float targetlumi_16e,bool DoSYS=true,SKIMMING mySkimming=KEEP_ALL)
{  
  /// only 
  const bool isVRsSRs = false;


  int targetlumi = targetlumi_16a+targetlumi_16d+targetlumi_16e;
  cout<<"targetlumi="<<targetlumi<<", DoSYS = " << DoSYS <<", mySkimming = " << mySkimming <<endl;
  gROOT->ProcessLine("#include <vector>");

  gErrorIgnoreLevel=kError;

  const char* names[] = {
    "nom", // the nominal sample
    //// turning off the systematic trees
//     "EG_RESOLUTION_ALL__1down", "EG_RESOLUTION_ALL__1up",
//     "EG_SCALE_AF2__1down", "EG_SCALE_AF2__1up",
//     "EG_SCALE_ALL__1down", "EG_SCALE_ALL__1up",
//     "JET_EtaIntercalibration_Modelling__1up", "JET_EtaIntercalibration_Modelling__1down",
//     "JET_EtaIntercalibration_NonClosure_highE__1up", "JET_EtaIntercalibration_NonClosure_highE__1down",
//     "JET_EtaIntercalibration_NonClosure_negEta__1up", "JET_EtaIntercalibration_NonClosure_negEta__1down",
//     "JET_EtaIntercalibration_NonClosure_posEta__1up", "JET_EtaIntercalibration_NonClosure_posEta__1down",
// 	"JET_Flavor_Composition__1up", "JET_Flavor_Composition__1down",
// 	"JET_Flavor_Response__1up", "JET_Flavor_Response__1down",
//     "JET_GroupedNP_1__1up", "JET_GroupedNP_1__1down",
//     "JET_GroupedNP_2__1up", "JET_GroupedNP_2__1down",
//     "JET_GroupedNP_3__1up", "JET_GroupedNP_3__1down",
//     "JET_JER_DataVsMC_MC16__1up", "JET_JER_DataVsMC_MC16__1down",
// 	"JET_JER_EffectiveNP_1__1up", "JET_JER_EffectiveNP_1__1down",
//     "JET_JER_EffectiveNP_2__1up", "JET_JER_EffectiveNP_2__1down",
//     "JET_JER_EffectiveNP_3__1up", "JET_JER_EffectiveNP_3__1down",
//     "JET_JER_EffectiveNP_4__1up", "JET_JER_EffectiveNP_4__1down",
//     "JET_JER_EffectiveNP_5__1up", "JET_JER_EffectiveNP_5__1down",
//     "JET_JER_EffectiveNP_6__1up", "JET_JER_EffectiveNP_6__1down",
//     "JET_JER_EffectiveNP_7restTerm__1up", "JET_JER_EffectiveNP_7restTerm__1down",
//     "JET_SingleParticle_HighPt__1up","JET_SingleParticle_HighPt__1down",
// 	"MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp",
//     "MET_SoftTrk_ScaleDown", "MET_SoftTrk_ScaleUp",
//     "MUON_ID__1down", "MUON_ID__1up",
//     "MUON_MS__1down", "MUON_MS__1up",
//     "MUON_SAGITTA_RESBIAS__1down", "MUON_SAGITTA_RESBIAS__1up",
//     "MUON_SAGITTA_RHO__1down", "MUON_SAGITTA_RHO__1up",
//     "MUON_SCALE__1down", "MUON_SCALE__1up"
  };

	std::vector<std::string> sysnamesAll(names, names + sizeof(names)/sizeof(names[0]));
	/// data has only nominal
	std::vector<std::string> sysnamesNom; sysnamesNom.push_back("nom");
	
	std::vector<std::string>sysnames;
	if(DoSYS==true) sysnames = sysnamesAll;
	else sysnames = sysnamesNom;
		
	cout<<"Path ="<<basepath<<endl;
	

    //// Path for DAOD_PHYS samples
    //// new ntuple area 21p2p10
    TString basepathBkg("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/Merged/");
    TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/Merged/");
    
    cout<<"PathBkg: "<<basepathBkg<<endl;
    cout<<"PathSig: "<<basepathSig<<endl;
    
    if(DEBUG)cout << "3 processed" << endl;
    
    // Data  
    //readFolder(basepath+"DATA",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnamesData,1);
    //readFolder(basepathSig+"Data",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnamesData,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);

    //// TString dirname, TString outfname, TString intreename, TString outtreename, std::vector<string>& sys, float targetlumi, float targetlumi_16a, float targetlumi_16d, float targetlumi_16e
    
    ///// for other WPs  //////////
    ////// for WP2,3,4   ///////////
    /// Backgrounds
    
//     readFolder(basepathBkg+"ttH",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttH",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"RareMore",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","RareMore",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"VH",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","VH",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"3and4tSM",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ThreeAndFourTopSM",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"ttBar",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttBar",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"TTV",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","TTV",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"Vjets",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","Vjets",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     readFolder(basepathBkg+"MultiBoson",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","MultiBoson",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
//     
//     // signals
//     readFolder(basepathSig+"SUSY",TString("signal.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    
    
    
    //// for WP1   ///////////
    /// background
    readFolder(basepathBkg+"ttH",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttH",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"RareMore",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","RareMore",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"VH",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","VH",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"3and4tSM",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ThreeAndFourTopSM",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"ttBar",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttBar",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"TTV",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","TTV",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"Vjets",TString("backgroundOther.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","Vjets",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"MultiBoson/MultiBoson",TString("backgroundMulti1.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","MultiBoson",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"MultiBoson/DiBoson",TString("backgroundMultiWithDi2.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","MultiBoson",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathBkg+"MultiBoson/DiBoson",TString("backgroundDiBoson.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","DiBoson",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    
    
    ///// signals
    readFolder(basepathSig+"SUSY/FirstHalf",TString("signal1.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    readFolder(basepathSig+"SUSY/SecondHalf",TString("signal2.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    
    
    /************
    /////// DAOD_PHYS samples
    //// new ntuple area 21p2p102
    TString basepathBkg("/lustre/ific.uv.es/grid/atlas/t3/asantra/DAODPhys_WP1/Merged/");
    TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/DAODPhys_WP1/Merged/");
    
    cout<<"PathBkg: "<<basepathBkg<<endl;
    cout<<"PathSig: "<<basepathSig<<endl;
    
    if(DEBUG)cout << "3 processed" << endl;
    
    // Data  
    //readFolder(basepath+"DATA",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnamesData,1);
    readFolder(basepathSig+"Data",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    
    // TString dirname, TString outfname, TString intreename, TString outtreename, std::vector<string>& sys, float targetlumi, float targetlumi_16a, float targetlumi_16d, float targetlumi_16e
    ///// Backgrounds
    readFolder(basepathBkg+"Bkg",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","Bkg",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    ///// Signals
    readFolder(basepathBkg+"SUSY",TString("signal.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi_16a,targetlumi_16d,targetlumi_16e, DoSYS,mySkimming);
    ************/
    
	
}
