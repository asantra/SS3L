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
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define DEBUG 0

// map of input branch names
std::map<TString,TString> guessBranchMapping(TTree* tree){

	std::map<TString,TString> bnames;

	if(tree->GetBranch("EtMiss")){

		cout<<"Detected tree origin: Freiburg"<<endl;
		bnames["meff"]="Meff";
		bnames["minv"]="Minv";
		bnames["ht"]="Ht";
		bnames["mt"]="Mt";
		bnames["met"]="EtMiss";
		bnames["nBJets20"]="Nbjet";
		bnames["nJets25"]="Njet25";
		bnames["nJets35"]="Njet35";
		bnames["nJets40"]="Njet40";
		bnames["nJets50"]="Njet50";
		bnames["SSChannel"]="SSChannel";
		bnames["puweight"]="PileUpWeight";
		bnames["totweight"]="TotalWeight";
		bnames["mcweight"]="MCWeight";
		bnames["sysSF"]="SysWeights";
		bnames["sysNames"]="SysNames";
		//bnames["trigSF"]="TriggerSF";
		bnames["chflipSF"]="CFTWeights";
		//bnames["lepCFmedium"]="lepCFmedium";
		//bnames["sysNames"]="SysNames";
		bnames["lepTLV"]="lepTLV";
		bnames["lepTLV_BL"]="lepTLV_BL";
		bnames["lepTruth"]="lepTruth";
		bnames["jetTLV"]="jetTLV";
		bnames["jetBtag"]="jetBtag";
		bnames["lepCharges"]="lepCharges";
		bnames["NlepBL"]="NlepBL";
		//bnames["pdfw"]="PDFWeights";
		bnames["runn"]="RunNumber";
		bnames["evn"]="EventNumber";
		bnames["Mjj"]="Mjj";
		bnames["Zevent"]="Zevent";
		bnames["MCId"]="MCId";
		bnames["PDGId1"]="PDGId1";
		bnames["PDGId2"]="PDGId2";
		bnames["GLDec1"]="GLDec1";
		bnames["GLDec2"]="GLDec2";
		bnames["phi_met"]="EtMissPhi";
		//bnames["SumJetPt"] = "SumJetPt";


		cout<<"Setup branches"<<endl;
	} 
	else {  
		cout<<"ERROR: cannot detect input tree origin. Code will crash. Have fun."<<endl;
	}
	if(DEBUG)cout << "guessBranchMapping function done" << endl;
	return bnames;
}

// set branch to local address, creating it if not existing

void setupBranch(TTree* tree, TString name, Long64_t* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/L");
	else tree->SetBranchAddress(name,address);
    if(DEBUG)cout << "setupBranch function done" << endl;
}

void setupBranch(TTree* tree, TString name, char* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/O");
	else tree->SetBranchAddress(name,address);
    if(DEBUG)cout << "setupBranch function done" << endl;
}

void setupBranch(TTree* tree, TString name, bool* address){
  TBranch* _b=tree->GetBranch(name);
  if(_b==0) _b= tree->Branch(name,address,name+"/O");
  else tree->SetBranchAddress(name,address);
  if(DEBUG)cout << "setupBranch function done" << endl;
}

void setupBranch(TTree* tree, TString name, int* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/I");
	else tree->SetBranchAddress(name,address);
    if(DEBUG)cout << "setupBranch function done" << endl;
}

void setupBranch(TTree* tree, TString name, float* address){
	TBranch* _b=tree->GetBranch(name);
	if(_b==0) _b= tree->Branch(name,address,name+"/F");
	else tree->SetBranchAddress(name,address);
    if(DEBUG)cout << "setupBranch function done" << endl;
}

bool isAtlfast(Long64_t RunNumber){
    if(DEBUG)cout << "isAtlfast " << RunNumber << endl;
  /*if(RunNumber==304014) return true;                          //3top
	if(RunNumber==341177) return true;                          //ttH dilep
	if(RunNumber==331270 || RunNumber==331271) return true;     //ttH semilep and hadr
	if(RunNumber>=371200 && RunNumber<=371317) return true;     //2STEP 
	if(RunNumber>=372444 && RunNumber<=372524) return true;     //GSL 
	if(RunNumber>=373462 && RunNumber<=373478) return true;     //GSL     
	if(RunNumber>=370188 && RunNumber<=370238) return true;     //Gtt 4-5 body 
	if(RunNumber>=372300 && RunNumber<=372368) return true;     //bb1step
	if(RunNumber>=403380 && RunNumber<=403481) return true;     //RPV
	if(RunNumber>=998000 && RunNumber<=999999) return true;     //RPV copies
	if(RunNumber>=370600 && RunNumber<=370601) return true;     //NUHM2
	if(RunNumber>=370603 && RunNumber<=370623) return true;     //NUHM2 
	if(RunNumber>=388230 && RunNumber<=388238) return true;     //tt2step
	if(RunNumber>=403100 && RunNumber<=403108) return true;     //RPVLQD
	if(RunNumber>=403693 && RunNumber<=403742) return true;     //RPVLQD 
	if(RunNumber>=404420 && RunNumber<=404425) return true;     //RPVLQD*/
  if((RunNumber>=376776 && RunNumber <=376834) || (RunNumber>=376844 && RunNumber <=376859) || RunNumber==377110) return true; //2STEPWZ
  if(RunNumber>=375818 && RunNumber <=375932) return true; //GTT
  if(RunNumber>=436533 && RunNumber <=436594) return true; //BTT
  if(RunNumber>=377094 && RunNumber <=377102) return true; //TT2STEP
  if((RunNumber>=376733 && RunNumber <=376775) || (RunNumber>=376835 && RunNumber <=376843)) return true; //GGRPV
  if(RunNumber==413000 || RunNumber==413001) return true; //ttW
  if(RunNumber==345940 || RunNumber==345941) return true; //ttH
  
  
  return false;
}

bool isSignal(int RunNumber){

  /*	if(RunNumber>=370100 && RunNumber<=370249) return true;  //Gtt
	if(RunNumber>=373421 && RunNumber<=373448) return true;  //Gtt
	if(RunNumber>=371200 && RunNumber<=371317) return true;  //2step
	if(RunNumber>=372444 && RunNumber<=372524) return true;  //GSL
	if(RunNumber>=373462 && RunNumber<=373478) return true;  //GSL
	if(RunNumber>=372300 && RunNumber<=372368) return true;  //bb1step
	if(RunNumber>=403380 && RunNumber<=403481) return true;  //RPV
	if(RunNumber>=998000 && RunNumber<=999999) return true;  //RPV copies
	if(RunNumber>=370600 && RunNumber<=370623) return true;  //NUHM2
	if(RunNumber>=388230 && RunNumber<=388238) return true;  //tt2step
	if(RunNumber>=403100 && RunNumber<=403176) return true;  //RPVLQD  
	if(RunNumber>=403693 && RunNumber<=403742) return true;  //RPVLQD 
	if(RunNumber>=404420 && RunNumber<=404425) return true;  //RPVLQD
  */

  if(DEBUG)cout << "isSignal " << RunNumber << endl;
  if((RunNumber>=376776 && RunNumber <=376834) || (RunNumber>=376844 && RunNumber <=376859) || RunNumber==377110) return true; //2STEPWZ
  if(RunNumber>=375818 && RunNumber <=375932) return true; //GTT
  if(RunNumber>=436533 && RunNumber <=436594) return true; //BTT
  if(RunNumber>=377094 && RunNumber <=377102) return true; //TT2STEP
  if((RunNumber>=376733 && RunNumber <=376775) || (RunNumber>=376835 && RunNumber <=376843)) return true; //GGRPV
  return false;
}

bool isttZ(int MCID){
    if(DEBUG)cout << "isttZ " << MCID << endl;
	if( MCID==410218 || MCID==410219 || MCID==410220 ) return true;
	return false;
}

bool isttW(int MCID){
    if(DEBUG)cout << "isttW " << MCID << endl;
	if( MCID==410155 ) return true;
	return false;
}

std::vector<float> ttV_SF(int MCID, int Nb, bool use){

	std::vector<float> SF(0);
	float sf(1.), sf_up(1.), sf_down(1.);
	if( use && Nb>=3 && isttW(MCID) ){ 
		sf      = 1.34;
		sf_up   = 1.65;
		sf_down = 1.03;
	}
	if( use && Nb>=3 && isttZ(MCID) ){ 
		sf      = 1.37;
		sf_up   = 1.65;
		sf_down = 1.09;
	}
	SF.push_back(sf);
	SF.push_back(sf_up);
	SF.push_back(sf_down);
 
	if(DEBUG){std::cout << "ttV_SF \t" << Form("MCId(%i) Nb=%i, SF(%.3f|%.3f|%.3f)", MCID, Nb, SF[0], SF[1], SF[2]) << std::endl;}
	return SF;
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
  if(DEBUG)cout << "hasZee " << Zee << endl;
  return Zee;
}

/*
int hasZee(std::vector<TLorentzVector>  *leptons, std::vector<float> *charges){

	float M[2] = {81000., 101000};
	std::vector<TLorentzVector> electrons(0);
	std::vector<float> el_charges(0);
	for(unsigned int i(0); i<leptons->size(); i++){ 
		if(leptons->at(i).M()<50.) {
			electrons.push_back(leptons->at(i)); 
			el_charges.push_back(charges->at(i));
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
		default: break;
	}
	return Zee;
}
*/
float getSumPtLep(std::vector<TLorentzVector> leptons){
	float sumPt(0);
	if(leptons.size()==0) return sumPt;
	for(unsigned int i(0); i<leptons.size(); i++) sumPt += (leptons.at(i)).Pt();
    if(DEBUG)cout << "getSumPtLep " << sumPt << endl;
	return sumPt;
}

float getSumPtJet(std::vector<TLorentzVector> jets, std::vector<float> btags, TString type){
  float sumPt(0);
  if(jets.size()==0 || btags.size()==0) return sumPt;
  for(unsigned int i(0); i<jets.size(); i++){
    float pt(0);
    if(type=="B"){ pt = btags.at(i) ? (jets.at(i)).Pt() : 0.; }
    if(type=="L"){ pt = !btags.at(i) ? (jets.at(i)).Pt() : 0.; }
    if(type=="ALL"){ pt = (jets.at(i)).Pt(); }
    if(pt<25000 && type=="ALL")
      continue;
    sumPt += pt;
  }
  if(DEBUG)cout << "getSumPtJet " << sumPt << endl;
  return sumPt;
}
/*
float getSumPtJet(std::vector<TLorentzVector> *jets, std::vector<float> *btags, TString type){
	float sumPt(0);
	if(!jets || !btags) return sumPt;
	for(unsigned int i(0); i<jets->size(); i++){
		float pt(0);
		if(type=="B"){ pt = btags->at(i) ? (jets->at(i)).Pt() : 0.; }
		if(type=="L"){ pt = !btags->at(i) ? (jets->at(i)).Pt() : 0.; }
		if(type=="ALL"){ pt = (jets->at(i)).Pt(); }
		if(pt<25000 && type=="ALL")
		  continue;
		sumPt += pt;
	}
	return sumPt;
}
*/
float dR(float var1, float phi1, float var2, float phi2){
	float deta = TMath::Abs(var1 - var2);
	float dphi = TMath::Abs(phi1 - phi2) < TMath::Pi() ? TMath::Abs(phi1 - phi2) : 2*TMath::Pi() - TMath::Abs(phi1 - phi2);
    if(DEBUG)cout << "in dR " << endl;
	return TMath::Sqrt(deta*deta + dphi*dphi);
}

float dReta(TLorentzVector v1, TLorentzVector v2){
	return dR(v1.Eta(), v1.Phi(), v2.Eta(), v2.Phi());
}

float getml1j(TLorentzVector v1, TLorentzVector v2){
    return (v1+v2).M();
}

float dRy(TLorentzVector v1, TLorentzVector v2){
	return dR(v1.Rapidity(), v1.Phi(), v2.Rapidity(), v2.Phi());
}

float getDeltaRLepLep(std::vector<TLorentzVector> leptons){
	if(leptons.size()==0) return -99.;
	if(leptons.size() < 2) return -99.;
	return dReta(leptons.at(0), leptons.at(1)); 
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
// float getMassLepJet(std::vector<TLorentzVector> leptons, std::vector<TLorentzVector> jets, unsigned int Nl){
//   if (leptons.size()<2) return -99.;
//   float ml1j = 99999999.;
//   float minDR = 99.;
//   for (unsigned int j=0;j<jets.size();j++){
//     if (dReta(leptons.at(Nl), jets.at(j))<minDR && jets.at(j).Pt()>25000){
//       minDR=dReta(leptons.at(Nl), jets.at(j));
//       ml1j=getml1j(leptons.at(Nl), jets.at(j));
//     }
//   }
//   return ml1j;
// }

/*
float getDeltaRLepJet(std::vector<TLorentzVector> *leptons, std::vector<TLorentzVector> *jets, unsigned int Nl){
	if (leptons->size()<2) return -99.;
	float minDR = 99.;
	for (unsigned int j=0;j<jets->size();j++){
		if (dReta(leptons->at(Nl), jets->at(j))<minDR && jets->at(j).Pt()>25000){
			minDR=dReta(leptons->at(Nl), jets->at(j));
		}
	}
	return minDR;
}
*/
bool match(TLorentzVector v1, TLorentzVector v2){
	if(v1.Pt()<1E-4 || v2.Pt()<1E-4) return false;
	if(dReta(v1,v2)<1E-4 && TMath::Abs(v1.Pt()-v2.Pt())<1E-4) return true;
	else return false;
}

int getSSNegative(std::vector<float> charges){
	int count=0;
	for(unsigned int i=0;i<charges.size();i++){
		if (charges.at(i)<0) count ++;
	}
	return count;
}

// int has3LSSPrompt(std::vector<float> charges, std::vector<float> lepTruth){
//   int count_pos=0;
//   int count_neg=0;
//   //std::cout<<"**********"<<std::endl;
//   for(unsigned int i=0;i<charges.size();i++){
//     //std::cout<<"charge: "<<charges.at(i)<<std::endl;
//     if(charges.at(i)>0 && lepTruth.at(i)==1) count_pos++;
//     if(charges.at(i)<0 && lepTruth.at(i)==1) count_neg++;
//   }
//   //std::cout<<"count_pos: "<<count_pos<<"\t count_neg: "<<count_neg<<std::endl;
//   if (count_pos>2 || count_neg>2) return 1;
//   else return 0;
// }

int has3LSSPrompt(std::vector<float> charges, std::vector<float> lepTruth){
  int count_pos=0;
  int count_neg=0;
  for(unsigned int i=0;i<charges.size();i++){
    if(charges.at(i)>0) count_pos++;
    if(charges.at(i)<0) count_neg++;
  }
  if(DEBUG)cout << "has3LSSPrompt "  << endl;
  if (count_pos>2 || count_neg>2) return 1;
  else return 0;
}

int is3LSSprocess(int MCId){
  /*bool isBkg = (MCId==407321 || 
			MCId==361623 || 
			MCId==343365 || 
			MCId==343366 || 
			MCId==342284 || 
			MCId==342285);

	bool isSig = (MCId>=388230 && MCId<=388238);
  */
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
/*
int isSSLepPtCut(std::vector<TLorentzVector> *leptons,std::vector<float> *charges){
  int isAbove30 = 0;
  int isSSpair = 0;
  for(unsigned int i(0); i<leptons->size(); i++){
    for(unsigned int j(0); j<leptons->size(); j++){
      if(i==j || isSSpair==1)
	continue;
      if(charges->at(i)==charges->at(j))
	isSSpair = 1;
      if(isSSpair==1 && leptons->at(i).Pt()>30000 && leptons->at(j).Pt()>30000)
	isAbove30=1;
    }
  }
  return isAbove30;
}
*/
float GetSFOSmass(std::vector<TLorentzVector> leptons,std::vector<float> charges){

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
/*
float GetSFOSmass(std::vector<TLorentzVector> *leptons,std::vector<float> *charges){

	float c(0);
	int nLep = leptons->size();
	switch(nLep){
		case 2:
			c = 0;
			break;
		case 3:
			if( charges->at(0) != charges->at(1) && (int)leptons->at(0).M()==(int)leptons->at(1).M() ){
				c=(leptons->at(0)+leptons->at(1)).M();
			}
			if( charges->at(0) != charges->at(2) && (int)leptons->at(0).M()==(int)leptons->at(2).M() ){
				if(TMath::Abs((leptons->at(0)+leptons->at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(0)+leptons->at(2)).M();
			}
			if( charges->at(2) != charges->at(1) && (int)leptons->at(2).M()==(int)leptons->at(1).M() ){
				if(TMath::Abs((leptons->at(2)+leptons->at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(2)+leptons->at(1)).M();
			}
			break;
		case 4:
			if( charges->at(0) != charges->at(1) && (int)leptons->at(0).M()==(int)leptons->at(1).M() ){
				if(TMath::Abs((leptons->at(0)+leptons->at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(0)+leptons->at(1)).M();
			}
			if( charges->at(0) != charges->at(2) && (int)leptons->at(0).M()==(int)leptons->at(2).M() ){
				if(TMath::Abs((leptons->at(0)+leptons->at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(0)+leptons->at(2)).M();
			}
			if( charges->at(0) != charges->at(3) && (int)leptons->at(0).M()==(int)leptons->at(3).M() ){
				if(TMath::Abs((leptons->at(0)+leptons->at(3)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(0)+leptons->at(3)).M();
			}
			if( charges->at(1) != charges->at(2) && (int)leptons->at(1).M()==(int)leptons->at(2).M() ){
				if(TMath::Abs((leptons->at(1)+leptons->at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(1)+leptons->at(2)).M();
			}
			if( charges->at(1) != charges->at(3) && (int)leptons->at(1).M()==(int)leptons->at(3).M() ){
				if(TMath::Abs((leptons->at(3)+leptons->at(1)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(3)+leptons->at(1)).M();
			}
			if( charges->at(2) != charges->at(3) && (int)leptons->at(2).M()==(int)leptons->at(3).M() ){
				if(TMath::Abs((leptons->at(3)+leptons->at(2)).M()-91190) < TMath::Abs(c-91190)) c=(leptons->at(2)+leptons->at(3)).M();
			}
			break;
		default:
			break;
	}
	return c;
}
*/


//// taken from Marco
// lepton closest to jet(s)
auto closestLep( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    if ( LepPt.size()>0 ) {
        if ( JetPt.size() == 1 ) {	
            TLorentzVector jet_vec; 
            jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
            TLorentzVector lep_vec;
            lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
            float dr = jet_vec.DeltaR( lep_vec ); // distance from leading lepton
            TLorentzVector lep_vec1;
            lep_vec1.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] ); 
            float dr1 = jet_vec.DeltaR( lep_vec1 ); // distance from sub-leading lepton
            if ( dr <= dr1 )
                return 0;
            else
                return 1;
        }
        if ( JetPt.size() >= 2 ) {	
            TLorentzVector jet_vec; 
            jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
            TLorentzVector jet_vec1; 
            jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
            TLorentzVector dijet_vec = jet_vec + jet_vec1;
            TLorentzVector lep_vec;
            lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
            float dr = dijet_vec.DeltaR( lep_vec ); // distance from leading lepton
            TLorentzVector lep_vec1;
            lep_vec1.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] ); 
            float dr1 = dijet_vec.DeltaR( lep_vec1 ); // distance from sub-leading lepton
            if ( dr <= dr1 )
                return 0;
            else
                return 1;
        }
    }
    if(DEBUG)cout << "closest Lep" << endl;
    return (-1);
};

// mlj
auto cmlj( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
    int icl = closestLep( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
    if ( icl>-1 ) { 
        TLorentzVector lep_vec;
        lep_vec.SetPtEtaPhiM( LepPt[icl], LepEta[icl], LepPhi[icl], LepM[icl] ); 
        TLorentzVector jet_vec; 
        if ( JetPt.size() == 1 )	
            jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
        else if ( JetPt.size() >= 2 ) {	
            TLorentzVector jet_vec0; 
            jet_vec0.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
            TLorentzVector jet_vec1; 
            jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
            jet_vec = jet_vec0 + jet_vec1;
        }
        return ( ( jet_vec + lep_vec ).M() );
    }
    if(DEBUG)cout << "cmlj" << endl;
    return (double)(-1);
};




// mTl2
float cmTl2( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
    if ( LepPt.size()>0 ) {
        TLorentzVector lep_vec;
        lep_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
        TVector2 v2; 
        v2.SetMagPhi( eT_miss, EtMissPhi );
        TLorentzVector met_vec;
        met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
        float mm = (lep_vec.Mt() + met_vec.Mt())*(lep_vec.Mt() + met_vec.Mt()) - (lep_vec+met_vec).Perp2();
        return ( mm>=0. ? TMath::Sqrt(mm) : TMath::Sqrt(-mm) );
    }
    else 
        return (double)(-1);
};

// mTlmin
float cmTlmin( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
    if ( LepPt.size()>0 ) {
        TLorentzVector lep1_vec;
        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
        TLorentzVector lep2_vec;
        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
        TVector2 v2; 
        v2.SetMagPhi( eT_miss, EtMissPhi );
        TLorentzVector met_vec;
        met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
        float mm1 = (lep1_vec.Mt() + met_vec.Mt())*(lep1_vec.Mt() + met_vec.Mt()) - (lep1_vec+met_vec).Perp2();
        float mt1 = mm1>=0. ? TMath::Sqrt(mm1) : TMath::Sqrt(-mm1);
        float mm2 = (lep2_vec.Mt() + met_vec.Mt())*(lep2_vec.Mt() + met_vec.Mt()) - (lep2_vec+met_vec).Perp2();
        float mt2 = mm2>=0. ? TMath::Sqrt(mm2) : TMath::Sqrt(-mm2);
        if ( mt1 <= mt2 ) return mt1;
        else		  return mt2;
    }
    else 
        return (float)(-1);
};

// mTlmax
float cmTlmax( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
    if ( LepPt.size()>0 ) {
        TLorentzVector lep1_vec;
        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
        TLorentzVector lep2_vec;
        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
        TVector2 v2; 
        v2.SetMagPhi( eT_miss, EtMissPhi );
        TLorentzVector met_vec;
        met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
        float mm1 = (lep1_vec.Mt() + met_vec.Mt())*(lep1_vec.Mt() + met_vec.Mt()) - (lep1_vec+met_vec).Perp2();
        float mt1 = mm1>=0. ? TMath::Sqrt(mm1) : TMath::Sqrt(-mm1);
        float mm2 = (lep2_vec.Mt() + met_vec.Mt())*(lep2_vec.Mt() + met_vec.Mt()) - (lep2_vec+met_vec).Perp2();
        float mt2 = mm2>=0. ? TMath::Sqrt(mm2) : TMath::Sqrt(-mm2);
        if ( mt1 >= mt2 ) return mt1;
        else	return mt2;
    }
    else 
        return (float)(-1);
};



/// ends Marco's code



float getLepTruth(std::vector<float> *lepTruth){
	int sum=0;
	for(unsigned int ilep=0; ilep<lepTruth->size();ilep++){
		sum += lepTruth->at(ilep);
	}
	if(sum>0)
		return 1.;
	else
		return 0.;
}

TString printVec(std::vector<float> v){
	if( v.size() < 3 ) return "";
	return Form(" [%.5f|%.5f|%.5f]",v[0],v[1],v[2]);
}

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

void copyTree(TTree* inTree, TTree* outTree, float xsec, float runc, Long64_t MC_ID, float nMCevents, float lumi,std::vector<Long64_t>& writtenEvents){

	//std::cout<<"copying tree. already written "<<writtenEvents.size()<<endl;

	// input and output trees, tag/folder name, xsection, #events in original MC sample, target lumi
        //ofstream myfile;
        //std::string s = std::to_string(MC_ID);
        //std::string l = std::to_string(lumi);
        //myfile.open("Rpc3LSS1b_"+s+"_"+l+".txt");

	float totweight=0, TotWeightNew=0, MC_campaign_weight=1;
	int nBJets20=0,nJets20=0,nJets25=0,nJets35=0,nJets40=0,nJets50,channel=0,nSigLep=0,NlepBL=0,Zevent=0, isSS30=0;

	Long64_t evn64=0,runn64=0;
	Int_t evn=0,runn=0,MCId=0;
	Int_t PDGId1=0,PDGId2=0,GLDec1=0,GLDec2=0,SSChannel=0;
	float met=0,meff=0,mt=0, MT2=0;
	int isZee=0,is3LSS=0,SSNegative=0,is3LSSproc=0;

	float  ht=0,minv=0,mcweight=0,puweight=0,Mjj=0, Pt_l=0, Pt_subl=0,mSFOS=0,m_ee=0,phi_met=0;
	float SSCharge=0;

	float sumPtLep(0), SumJetPt(0), SumBJetPt(0);
	float dRll(-99.), dRl1j(-99.), dRl2j(-99.), ml1j(-999.);
        float mlj(-999.), mTl2(-999.), mTlmin(-999.), mTlmax(-999.);

	float elID_rescale_up=1, elID_rescale_down=1;

	std::vector<float> *SysSF=0;
	float TrigSF=0;
	std::vector<float> *ChFlipSF=0, *lepCFmedium=0, *tmp_lepCFmedium=0;
	std::vector<float> *PdfSF=0;
	std::vector<float> *pdfw=0;
	std::vector<string> *SysNames=0;
	std::vector<TLorentzVector> *lepTLV=0, *lepTLV_BL=0, *jetTLV=0;
	std::vector<TLorentzVector> tmp_jetTLV,tmp_lepTLV,tmp_lepTLV_BL;
	std::vector<float> tmp_jetBtag,tmp_lepCharges,tmp_lepTruth;
	std::vector<float> *lepCharges=0,*lepTruth=0;
	std::vector<float> *jetBtag=0;
    std::vector<float> lepPt, lepEta, lepPhi, lepM;
    std::vector<float> jetPt, jetEta, jetPhi, jetM;

	std::map<TString,TString> bnames=guessBranchMapping(inTree);
	bool isSR=false, isVR=false;

	int nentries=inTree->GetEntries();

	std::cout<<"inTree has "<<nentries<<std::endl;

	cout<<"got xsec nentries nMCevents "<<xsec<< " "<<nentries<<" "<<nMCevents<<endl;

	float lumiScaling=xsec*lumi*1/nMCevents;

	// data
	if(xsec==0) lumiScaling=1;

	//cout<<"final scaling "<<lumiScaling<<endl;

	float lumiScalingUP=lumiScaling;
	float lumiScalingDOWN=lumiScaling;

	// avoid reading what we do not need
	inTree->SetBranchStatus("*",0);

	// turn on only the branches that will be read
	for(auto branch: bnames)
		inTree->SetBranchStatus(branch.second,1);


	TClass* expectedClass = 0;
	EDataType expectedType = kOther_t;
	inTree->GetBranch(bnames["evn"])->GetExpectedType(expectedClass,expectedType);


	// configure input branches
	inTree->SetBranchAddress(bnames["evn"],&evn64);
	inTree->SetBranchAddress(bnames["runn"],&runn64);
	inTree->SetBranchAddress(bnames["MCId"],&MCId);
	inTree->SetBranchAddress(bnames["meff"],&meff);
	inTree->SetBranchAddress(bnames["met"],&met);
	inTree->SetBranchAddress(bnames["mt"],&mt);
	//inTree->SetBranchAddress(bnames["ht"],&ht);
	inTree->SetBranchAddress(bnames["minv"],&minv);
	inTree->SetBranchAddress(bnames["nJets25"],&nJets25);
	inTree->SetBranchAddress(bnames["nJets35"],&nJets35);
	inTree->SetBranchAddress(bnames["nJets40"],&nJets40);
	inTree->SetBranchAddress(bnames["nJets50"],&nJets50);
	inTree->SetBranchAddress(bnames["nBJets20"],&nBJets20);
	inTree->SetBranchAddress(bnames["sysSF"],&SysSF);
	inTree->SetBranchAddress(bnames["sysNames"],&SysNames);
	//inTree->SetBranchAddress(bnames["trigSF"],&TrigSF);
	inTree->SetBranchAddress(bnames["chflipSF"],&ChFlipSF);
	//inTree->SetBranchAddress(bnames["lepCFmedium"],&tmp_lepCFmedium);
	//inTree->SetBranchAddress(bnames["sysNames"],&SysNames);
	inTree->SetBranchAddress(bnames["lepTLV"],&lepTLV);
	inTree->SetBranchAddress(bnames["lepTLV_BL"],&lepTLV_BL);
	inTree->SetBranchAddress(bnames["lepTruth"],&lepTruth);
	inTree->SetBranchAddress(bnames["jetTLV"],&jetTLV);
	inTree->SetBranchAddress(bnames["jetBtag"],&jetBtag);
	inTree->SetBranchAddress(bnames["lepCharges"],&lepCharges);
	inTree->SetBranchAddress(bnames["mcweight"],&mcweight);
	inTree->SetBranchAddress(bnames["totweight"],&totweight);
	inTree->SetBranchAddress(bnames["puweight"],&puweight);
	//inTree->SetBranchAddress(bnames["pdfw"],&pdfw);
	//inTree->SetBranchAddress(bnames["Zevent"],&Zevent);
	inTree->SetBranchAddress(bnames["Mjj"],&Mjj);
	inTree->SetBranchAddress(bnames["NlepBL"], &NlepBL);
	inTree->SetBranchAddress(bnames["MCId"], &MCId);
	inTree->SetBranchAddress(bnames["PDGId1"],&PDGId1);
	inTree->SetBranchAddress(bnames["PDGId2"],&PDGId2);
	inTree->SetBranchAddress(bnames["GLDec1"],&GLDec1);
	inTree->SetBranchAddress(bnames["GLDec2"],&GLDec2);
	inTree->SetBranchAddress(bnames["SSChannel"],&SSChannel);
	inTree->SetBranchAddress(bnames["phi_met"],&phi_met);
		
	float wmu_nom,
		  wmu_stat_up,
		  wmu_stat_down,
		  wmu_sys_up,
		  wmu_sys_down,
		  wmu_bad_sys_up,
		  wmu_bad_sys_down,
		  wmu_bad_stat_up,
		  wmu_bad_stat_down,
		  wmu_stat_lowpt_up,
		  wmu_stat_lowpt_down,
		  wmu_sys_lowpt_up,
		  wmu_sys_lowpt_down,
		  wmu_trig_stat_up,
		  wmu_trig_stat_down,
		  wmu_trig_sys_up,
		  wmu_trig_sys_down, 
		  wmu_iso_stat_up,
		  wmu_iso_stat_down,
		  wmu_iso_sys_up,
		  wmu_iso_sys_down, 
		  wmu_ttva_stat_up,
		  wmu_ttva_stat_down,
		  wmu_ttva_sys_up,
		  wmu_ttva_sys_down,

		  wel_nom,
		  wel_cid_up,
		  wel_cid_down,
		  wel_id_up,
		  wel_id_down,
		  wel_reco_up,
		  wel_reco_down,
		  wel_trig_up,
		  wel_trig_down,
		  wel_trigEff_up,
		  wel_trigEff_down,
		  wel_iso_up,
		  wel_iso_down,

		  wjet_nom,
		  wjet_b_up,
		  wjet_b_down,
		  wjet_c_up,
		  wjet_c_down,
		  wjet_light_up,
		  wjet_light_down,
		  wjet_extra1_up,
		  wjet_extra1_down,
		  wjet_extra2_up,
		  wjet_extra2_down,

		  wjet_jvt_up,
		  wjet_jvt_down,

		  wpu_nom_sig,
		  wpu_up_sig,
		  wpu_down_sig,

		  wpu_nom_bkg,
		  wpu_up_bkg,
		  wpu_down_bkg,

		  wtrig_nom,
		  wtrig_up,
		  wtrig_down,

		  wchflip_nom,
		  wchflip_up,
		  wchflip_down,

		  wpdf_up,
		  wpdf_down,

		  wttV_nom,
		  wttV_up,
		  wttV_down;

	// configure output branches
	setupBranch(outTree,"MCId",&MCId);
	setupBranch(outTree,"runn",&runn64);
	setupBranch(outTree,"evn",&evn64);
	setupBranch(outTree,"nBJets20",&nBJets20);
	setupBranch(outTree,"nJets20",&nJets20);
	setupBranch(outTree,"nJets25",&nJets25);
	setupBranch(outTree,"nJets35",&nJets35);
	setupBranch(outTree,"nJets40",&nJets40);
	setupBranch(outTree,"nJets50",&nJets50);
	setupBranch(outTree,"sumPtLep",&sumPtLep);
	setupBranch(outTree,"SumJetPt",&SumJetPt);
	setupBranch(outTree,"SumBJetPt",&SumBJetPt);
	setupBranch(outTree,"Pt_l",&Pt_l);
	setupBranch(outTree,"Pt_subl",&Pt_subl);
	//setupBranch(outTree,"Pt_3",&Pt_3);
	//setupBranch(outTree,"Mass_3",&Mass_3);
	setupBranch(outTree,"NlepBL",&NlepBL);
	setupBranch(outTree,"Mjj",&Mjj);
	//setupBranch(outTree,"MinvLep",&MinvLep);
	setupBranch(outTree,"SSCharge",&SSCharge);
	setupBranch(outTree,"SSNegative",&SSNegative);
	//setupBranch(outTree,"Zevent",&Zevent);
	setupBranch(outTree,"isSS30",&isSS30);
	//setupBranch(outTree,"SRveto",&SRveto);
	setupBranch(outTree,"nSigLep",&nSigLep);
	//setupBranch(outTree,"nSigLepEta",&nSigLepEta);
	setupBranch(outTree,"mt",&mt);
	setupBranch(outTree,"ht",&ht);
	setupBranch(outTree,"minv",&minv);
	setupBranch(outTree,"meff",&meff);
	setupBranch(outTree,"met",&met);
    setupBranch(outTree,"phi_met",&phi_met);
    setupBranch(outTree,"mlj",&mlj);
    setupBranch(outTree,"mTl2",&mTl2);
    setupBranch(outTree,"mTlmin",&mTlmin);
    setupBranch(outTree,"mTlmax",&mTlmax);
	setupBranch(outTree,"mcweight",&mcweight);
	setupBranch(outTree,"puweight",&puweight);
	setupBranch(outTree,"totweight",&TotWeightNew);
	setupBranch(outTree,"lumiScaling",&lumiScaling);
	setupBranch(outTree,"MC_campaign_weight",&MC_campaign_weight);
	//setupBranch(outTree,"Pt_j2",&Pt_j2);
	//setupBranch(outTree,"m3l",&m3l);
	//setupBranch(outTree,"Veto_ee",&Veto_ee);
	setupBranch(outTree,"mSFOS",&mSFOS);
	setupBranch(outTree,"m_ee",&m_ee);
	setupBranch(outTree,"isZee",&isZee);
	setupBranch(outTree,"SSChannel",&SSChannel);
	setupBranch(outTree,"is3LSS",&is3LSS);
	setupBranch(outTree,"is3LSSproc",&is3LSSproc);
	setupBranch(outTree,"MT2",&MT2);
	setupBranch(outTree,"dRll", &dRll);
	setupBranch(outTree,"dRl1j",&dRl1j);
        //setupBranch(outTree,"ml1j",&ml1j);
	setupBranch(outTree,"dRl2j",&dRl2j);
	setupBranch(outTree,"isSR",&isSR);

// 	setupBranch(outTree,"wmu_nom",&wmu_nom);	  
// 	setupBranch(outTree,"wmu_stat_up",&wmu_stat_up);
// 	setupBranch(outTree,"wmu_stat_down",&wmu_stat_down ); 
// 	setupBranch(outTree,"wmu_sys_up",&wmu_sys_up);		  
// 	setupBranch(outTree,"wmu_sys_down",&wmu_sys_down);
// 	setupBranch(outTree,"wmu_bad_sys_up", &wmu_bad_sys_up);
// 	setupBranch(outTree,"wmu_bad_sys_down", &wmu_bad_sys_down);
// 	setupBranch(outTree,"wmu_bad_stat_up", &wmu_bad_stat_up);
// 	setupBranch(outTree,"wmu_bad_stat_down", &wmu_bad_stat_down);
// 	setupBranch(outTree,"wmu_stat_lowpt_up",&wmu_stat_lowpt_up);
// 	setupBranch(outTree,"wmu_stat_lowpt_down",&wmu_stat_lowpt_down );
// 	setupBranch(outTree,"wmu_sys_lowpt_up",&wmu_sys_lowpt_up);
// 	setupBranch(outTree,"wmu_sys_lowpt_down",&wmu_sys_lowpt_down);
// 	setupBranch(outTree,"wmu_trig_stat_up",&wmu_trig_stat_up);	  
// 	setupBranch(outTree,"wmu_trig_stat_down",&wmu_trig_stat_down);	  
// 	setupBranch(outTree,"wmu_trig_sys_up",&wmu_trig_sys_up);	  
// 	setupBranch(outTree,"wmu_trig_sys_down",&wmu_trig_sys_down);	  
// 	setupBranch(outTree,"wmu_iso_stat_up",&wmu_iso_stat_up);	  
// 	setupBranch(outTree,"wmu_iso_stat_down",&wmu_iso_stat_down);	  
// 	setupBranch(outTree,"wmu_iso_sys_up",&wmu_iso_sys_up);	  
// 	setupBranch(outTree,"wmu_iso_sys_down",&wmu_iso_sys_down);	  
// 	setupBranch(outTree,"wmu_ttva_stat_up",&wmu_ttva_stat_up);
// 	setupBranch(outTree,"wmu_ttva_stat_down",&wmu_ttva_stat_down);
// 	setupBranch(outTree,"wmu_ttva_sys_up",&wmu_ttva_sys_up);
// 	setupBranch(outTree,"wmu_ttva_sys_down",&wmu_ttva_sys_down);
// 
// 	setupBranch(outTree,"wel_nom",&wel_nom);
// 	setupBranch(outTree,"wel_cid_up",&wel_cid_up);
// 	setupBranch(outTree,"wel_cid_down",&wel_cid_down);		  
// 	setupBranch(outTree,"wel_id_up",&wel_id_up);		  
// 	setupBranch(outTree,"wel_id_down",&wel_id_down);
// 	setupBranch(outTree,"wel_iso_up",&wel_iso_up);
// 	setupBranch(outTree,"wel_iso_down",&wel_iso_down);	  
// 	setupBranch(outTree,"wel_reco_up",&wel_reco_up);	  
// 	setupBranch(outTree,"wel_reco_down",&wel_reco_down); 
// 	setupBranch(outTree,"wel_trig_up",&wel_trig_up);	  
// 	setupBranch(outTree,"wel_trig_down",&wel_trig_down);  
// 	setupBranch(outTree,"wel_trigEff_up",&wel_trigEff_up);	  
// 	setupBranch(outTree,"wel_trigEff_down",&wel_trigEff_down);  
// 
// 
// 	setupBranch(outTree,"wjet_nom",&wjet_nom);		  
// 	setupBranch(outTree,"wjet_b_up",&wjet_b_up);		  
// 	setupBranch(outTree,"wjet_b_down",&wjet_b_down);
// 	setupBranch(outTree,"wjet_c_up",&wjet_c_up);		  
// 	setupBranch(outTree,"wjet_c_down",&wjet_c_down);
// 	setupBranch(outTree,"wjet_light_up",&wjet_light_up ); 
// 	setupBranch(outTree,"wjet_light_down",&wjet_light_down);
// 	setupBranch(outTree,"wjet_extra1_up",&wjet_extra1_up ); 
// 	setupBranch(outTree,"wjet_extra1_down",&wjet_extra1_down);
// 	setupBranch(outTree,"wjet_extra2_up",&wjet_extra2_up ); 
// 	setupBranch(outTree,"wjet_extra2_down",&wjet_extra2_down);
// 
// 	setupBranch(outTree,"wjet_jvt_up",&wjet_jvt_up);
// 	setupBranch(outTree,"wjet_jvt_down",&wjet_jvt_down);
// 
// 	setupBranch(outTree,"wpu_nom_sig",&wpu_nom_sig);		  
// 	setupBranch(outTree,"wpu_up_sig",&wpu_up_sig);		  
// 	setupBranch(outTree,"wpu_down_sig",&wpu_down_sig);
// 	setupBranch(outTree,"wpu_nom_bkg",&wpu_nom_bkg);
// 	setupBranch(outTree,"wpu_up_bkg",&wpu_up_bkg);
// 	setupBranch(outTree,"wpu_down_bkg",&wpu_down_bkg);
// 
// 	setupBranch(outTree,"wtrig_nom",&wtrig_nom);
// 	setupBranch(outTree,"wtrig_up",&wtrig_up);
// 	setupBranch(outTree,"wtrig_down",&wtrig_down);
// 
// 	setupBranch(outTree,"wchflip_nom",&wchflip_nom);
// 	setupBranch(outTree,"wchflip_up",&wchflip_up);
// 	setupBranch(outTree,"wchflip_down",&wchflip_down);
// 
// 	setupBranch(outTree,"wpdf_up",&wpdf_up);		  
// 	setupBranch(outTree,"wpdf_down",&wpdf_down);
// 
// 	setupBranch(outTree,"wttV_nom",  &wttV_nom);
// 	setupBranch(outTree,"wttV_up",   &wttV_up);
// 	setupBranch(outTree,"wttV_down", &wttV_down);


	// signals
	if(runc!=0){
		lumiScalingUP=lumiScaling*(1+runc);
		lumiScalingDOWN=lumiScaling*(1-runc);
		std::cout<< "got runc !=0 "<<runc<<" "<<lumiScaling<<" "<<lumiScalingUP<<" "<<lumiScalingDOWN<<std::endl;
	}

	setupBranch(outTree,"lumiScaling_UP",&lumiScalingUP);
	setupBranch(outTree,"lumiScaling_DOWN",&lumiScalingDOWN);

	//now all branches are set. will get once and for all the systematic mapping

	//inTree->GetEntry(0);

	//std::cout<<"Testing weights:"<<" " <<SysNames->size()<<std::endl;
	//for(auto n: *SysNames)
	//  cout<<n<<endl;

	int nocut=0,bveto=0;
	//bool isSR=false, isVR=false;
	//Rescaling for missing samples
	if(MC_ID==363507 || MC_ID==363509){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16e SAMPLES"<<std::endl;
	  MC_campaign_weight = 1.7444;
	}

	/*if(MC_ID==345706){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16a SAMPLES"<<std::endl;
          MC_campaign_weight = 1.3473;
	  }*/

	/*if(MC_ID==410219){
	  std::cout<<"ADDING EXTRA WEIGHT FOR MISSING MC16d SAMPLES"<<std::endl;
          MC_campaign_weight = 1.4608;
	  }*/
	
// 	if(MCId==412063) //Filtering out old tllq sample
// 	   continue;


	for (Long64_t i=0;i<nentries; i++) {
		inTree->GetEntry(i);
		nocut = nocut+1;
		isSR=false;
		isVR=false;
		if(MCId>=372444 && MCId<=372511 ){
                   if(GLDec1==5 || GLDec2==5){continue;}
                }
		bveto = bveto+1;
		writtenEvents.push_back(evn64);
		//if(evn64!=9705782 && evn64!=7686484 && evn64!=7745977 && evn64!=4674221 && evn64!=8013792 && evn64!=7467861 && evn64!=11167679)
		//if(evn64!=1780910)
		//continue;
		//std::cout<<"EVENT!"<<std::endl;
		//FIXING DUPLICATING LEP/JET ENTRIES
		//Lep BL
		tmp_lepTLV_BL.clear();
        int isGoodLep = 1;
        for(unsigned iLep=0;iLep<lepTLV_BL->size();iLep++){
           isGoodLep = 1;
           for(unsigned kLep=iLep+1;kLep<lepTLV_BL->size();kLep++){
                if((lepTLV_BL->at(iLep).Pt()==lepTLV_BL->at(kLep).Pt() && lepTLV_BL->at(iLep).Eta()==lepTLV_BL->at(kLep).Eta()))
                    isGoodLep = 0;
            }
            if(isGoodLep==1){
		        //std::cout<<"Lepton Pt: "<<lepTLV_BL->at(iLep).Pt()<<"\t Eta: "<<lepTLV_BL->at(iLep).Eta()<<std::endl;
                tmp_lepTLV_BL.push_back(lepTLV_BL->at(iLep));
            }
        }
		
		//SIG LEP
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
		    if(lepTLV->at(iLep).Pt()==lepTLV->at(kLep).Pt() && lepTLV->at(iLep).Eta()==lepTLV->at(kLep).Eta()){isGoodLep = 0;}
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

		//JETS
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
		    if(jetTLV->at(iJet).Pt()==jetTLV->at(kJet).Pt() && jetTLV->at(iJet).Eta()==jetTLV->at(kJet).Eta() && jetBtag->at(iJet)==jetBtag->at(kJet)){isGoodJet = 0;}
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
                
                nJets20 = 0;
		nJets25 = 0;
		nJets40 = 0;
		nBJets20 = 0;

		for(unsigned iJet=0;iJet<tmp_jetTLV.size();iJet++){
		  if(tmp_jetTLV.at(iJet).Pt()>20000)
                      nJets20++;
                  if(tmp_jetTLV.at(iJet).Pt()>25000)
		      nJets25++;
		  if(tmp_jetTLV.at(iJet).Pt()>40000)
                      nJets40++;
		  if(tmp_jetTLV.at(iJet).Pt()>20000 && tmp_jetBtag.at(iJet)>0)
		      nBJets20++;
		}

		//END
		nSigLep = tmp_lepTLV.size();
		NlepBL = tmp_lepTLV_BL.size();

		//Validation Regions 
		Pt_l = tmp_lepTLV.at(0).Pt();
		Pt_subl = tmp_lepTLV.at(1).Pt();

		//sumPt
		sumPtLep  = getSumPtLep(tmp_lepTLV);
		SumJetPt  = getSumPtJet(tmp_jetTLV,tmp_jetBtag,"ALL");
		SumBJetPt = getSumPtJet(tmp_jetTLV,tmp_jetBtag,"B");
	
		//dR
		dRll  = getDeltaRLepLep(tmp_lepTLV);
        
		dRl1j = getDeltaRLepJet(tmp_lepTLV, tmp_jetTLV, 0);
                
        dRl2j = getDeltaRLepJet(tmp_lepTLV, tmp_jetTLV, 1);
                
        /// new variables added from Marco
        mlj   = cmlj(lepPt, lepPhi, lepEta, lepM, jetPt, jetPhi, jetEta, jetM);
              
        mTl2  = cmTl2(lepPt, lepEta, lepPhi, lepM, met, phi_met);
                
        mTlmin = cmTlmin(lepPt, lepEta, lepPhi, lepM, met, phi_met);
                
        mTlmax = cmTlmax(lepPt, lepEta, lepPhi, lepM, met, phi_met);
                
                
        /// Marco variables end here

		//ttV rescaling //TBC
		/*std::vector<float> ttVSF = ttV_SF(MCId, nBJets20, true);
		wttV_nom  = ttVSF.at(0);
		wttV_up   = ttVSF.at(1);
		wttV_down = ttVSF.at(2);

		//Recompute TotalWeight
		TotWeightNew = totweight*wttV_nom;
		*/
		//Rescaling for missing samples
		TotWeightNew = MC_campaign_weight*totweight;

		//MT2
		TLorentzVector metTLV = TLorentzVector( met*TMath::Cos(phi_met) , met*TMath::Sin(phi_met) , 0 , met);
        if(tmp_lepTLV.size()>1){
             ComputeMT2 mycalc = ComputeMT2(tmp_lepTLV.at(0),tmp_lepTLV.at(1),metTLV,tmp_lepTLV.at(0).M(),tmp_lepTLV.at(1).M());
             MT2 = mycalc.Compute();
        }
        else
             MT2 = -999;

		TString SSPos;
		SSNegative = getSSNegative(tmp_lepCharges);
		if(tmp_lepCharges.size()<3)
		  is3LSS=0;
		else
		  is3LSS = has3LSSPrompt(tmp_lepCharges,tmp_lepTruth); 
		is3LSSproc = is3LSSprocess(MCId); 

		//new variables for VR
		if(tmp_lepTLV.size()>=2){
		  mSFOS=GetSFOSmass(tmp_lepTLV, tmp_lepCharges);
		  isSS30 = isSSLepPtCut(tmp_lepTLV, tmp_lepCharges);
		}
		else{
		  mSFOS=0;
		  isSS30=0;
		}
		
		m_ee=0;
		isZee=hasZee(tmp_lepTLV,tmp_lepCharges);

        
		// MC
// 		if(xsec!=0){
// 
// 			if(DEBUG) cout<<"--------------------------------"<<endl;
//                         if(DEBUG) cout << "Size of Sysnames " << SysNames->size() << endl;
// 
// 			if(DEBUG) cout << "Setting "<<SysNames->at(0)<<" to wmu_nom"<<endl;
// 			wmu_nom=SysSF->at(0);
// 			if(DEBUG) cout << "Setting "<<SysNames->at(29)<<" to wmu_bad_stat_up"<<endl;
// 			wmu_bad_stat_up=SysSF->at(29);
// 			if(DEBUG) cout << "Setting "<<SysNames->at(28)<<" to wmu_bad_stat_down"<<endl;
// 			wmu_bad_stat_down=SysSF->at(28);
// 			if(DEBUG) cout << "Setting "<<SysNames->at(31)<<" to wmu_bad_sys_up"<<endl;
// 			wmu_bad_sys_up=SysSF->at(31);
// 			if(DEBUG) cout << "Setting "<<SysNames->at(30)<<" to wmu_bad_sys_down"<<endl;
// 			wmu_bad_sys_down=SysSF->at(30);
// 			if(DEBUG) cout << "Setting "<<SysNames->at(33)<<" to wmu_stat_up"<<endl;
// 			wmu_stat_up=SysSF->at(33);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(32)<<" to wmu_stat_down"<<endl;
// 			wmu_stat_down=SysSF->at(32);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(35)<<" to wmu_stat_lowpt_up"<<endl;
// 			wmu_stat_lowpt_up=SysSF->at(35);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(34)<<" to wmu_stat_lowpt_down"<<endl;
// 			wmu_stat_lowpt_down=SysSF->at(34);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(37)<<" to wmu_sys_up"<<endl;
// 			wmu_sys_up=SysSF->at(37);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(36)<<" to wmu_sys_down"<<endl;
// 			wmu_sys_down=SysSF->at(36);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(39)<<" to wmu_sys_lowpt_up"<<endl;
// 			wmu_sys_lowpt_up=SysSF->at(39);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(38)<<" to wmu_sys_lowpt_down"<<endl;
// 			wmu_sys_lowpt_down=SysSF->at(38);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(41)<<" to wmu_trig_stat_up"<<endl;
// 			wmu_trig_stat_up=SysSF->at(41);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(40)<<" to wmu_trig_stat_down"<<endl;
// 			wmu_trig_stat_down=SysSF->at(40);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(43)<<" to wmu_trig_sys_up"<<endl;
// 			wmu_trig_sys_up=SysSF->at(43);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(42)<<" to wmu_trig_sys_down"<<endl;
// 			wmu_trig_sys_down=SysSF->at(42);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(45)<<" to wmu_iso_stat_up"<<endl;
// 			wmu_iso_stat_up=SysSF->at(45);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(44)<<" to wmu_iso_stat_down"<<endl;
// 			wmu_iso_stat_down=SysSF->at(44);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(47)<<" to wmu_iso_sys_up"<<endl;
// 			wmu_iso_sys_up=SysSF->at(47);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(46)<<" to wmu_iso_sys_down"<<endl;
// 			wmu_iso_sys_down=SysSF->at(46);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(49)<<" to wmu_ttva_stat_up"<<endl;
// 			wmu_ttva_stat_up=SysSF->at(49);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(48)<<" to wmu_ttva_stat_down"<<endl;
// 			wmu_ttva_stat_down=SysSF->at(48);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(51)<<" to wmu_ttva_sys_up"<<endl;
// 			wmu_ttva_sys_up=SysSF->at(51);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(50)<<" to wmu_ttva_sys_down"<<endl;
// 			wmu_ttva_sys_down=SysSF->at(50);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(1)<<" to wel_nom"<<endl;
// 			wel_nom=SysSF->at(1);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(5)<<" to wel_cid_up"<<endl;
// 			wel_cid_up=SysSF->at(5);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(4)<<" to wel_cid_down"<<endl;
// 			wel_cid_down=SysSF->at(4);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(7)<<" to wel_id_up"<<endl;
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(6)<<" to wel_id_down"<<endl;
        
// 			elID_rescale_up=1.;
// 			elID_rescale_down=1.;
// 			for(unsigned int ilep=0;ilep<tmp_lepTLV.size();ilep++){
// 			  if (getDeltaRLepJet(tmp_lepTLV,tmp_jetTLV,ilep)<0.4){
// 			    elID_rescale_up=1.2;
// 			    elID_rescale_down=0.8;
// 			  }
// 			}
// 			if(SysSF->at(7)>SysSF->at(6)){
// 			  wel_id_up=SysSF->at(7)*elID_rescale_up;
// 			  wel_id_down=SysSF->at(6)*elID_rescale_down;
// 			}
// 			else
// 			{
// 			  wel_id_up=SysSF->at(7)*elID_rescale_down;
// 			  wel_id_down=SysSF->at(6)*elID_rescale_up;
// 			}
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(9)<<" to wel_iso_up"<<endl;
// 			wel_iso_up = SysSF->at(9);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(8)<<" to wel_iso_down"<<endl;
// 			wel_iso_down = SysSF->at(8);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(11)<<" to wel_reco_up"<<endl;
// 			wel_reco_up=SysSF->at(11);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(10)<<" to wel_reco_down"<<endl;
// 			wel_reco_down=SysSF->at(10);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(13)<<" to wel_trigEff_up"<<endl;
// 			wel_trigEff_up=SysSF->at(13);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(12)<<" to wel_trigEff_down"<<endl;
// 			wel_trigEff_down=SysSF->at(12);                 
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(15)<<" to wel_trig_up"<<endl;
// 			wel_trig_up=SysSF->at(15);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(14)<<" to wel_trig_down"<<endl;
// 			wel_trig_down=SysSF->at(14);
// 
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(2)<<" to wjet_nom"<<endl;
// 			wjet_nom=SysSF->at(2);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(17)<<" to wjet_b_up"<<endl;
// 			wjet_b_up=SysSF->at(17);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(16)<<" to wjet_b_down"<<endl;
// 			wjet_b_down=SysSF->at(16);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(19)<<" to wjet_c_up"<<endl;
// 			wjet_c_up=SysSF->at(19);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(18)<<" to wjet_c_down"<<endl;
// 			wjet_c_down=SysSF->at(18);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(21)<<" to wjet_light_up"<<endl;
// 			wjet_light_up=SysSF->at(21);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(20)<<" to wjet_light_down"<<endl;
// 			wjet_light_down=SysSF->at(20);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(23)<<" to wjet_extra1_up"<<endl;
// 			wjet_extra1_up=SysSF->at(23);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(22)<<" to wjet_extra1_down"<<endl;
// 			wjet_extra1_down=SysSF->at(22);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(25)<<" to wjet_extra2_up"<<endl;
// 			wjet_extra2_up=SysSF->at(25);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(24)<<" to wjet_extra2_down"<<endl;
// 			wjet_extra2_down=SysSF->at(24);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(27)<<" to wjet_jvt_up"<<endl;
// 			wjet_jvt_up=SysSF->at(27);
// 			if(DEBUG) cout <<"Setting "<<SysNames->at(26)<<" to wjet_jvt_down"<<endl;
// 			wjet_jvt_down=SysSF->at(26);
// 
// 			if(!isSignal(MCId)){
// 				wpu_nom_sig=1;
// 				wpu_up_sig=1;
// 				wpu_down_sig=1;
// 				//wpu_up_sig=SysSF->at(47);
// 				//wpu_down_sig=SysSF->at(46);
// 			}
// 			else{
// 				wpu_nom_sig=SysSF->at(3);
// 				wpu_up_sig=SysSF->at(47);
// 				wpu_down_sig=SysSF->at(46);
// 			}
// 			if(isSignal(MCId)){
// 				wpu_nom_bkg=1;
// 				wpu_up_bkg=1;
// 				wpu_down_bkg=1;
// 				//wpu_up_bkg=SysSF->at(47);
// 				//wpu_down_bkg=SysSF->at(46);
// 			}
// 			else{
// 				wpu_nom_bkg=SysSF->at(3);
// 				wpu_up_bkg=SysSF->at(47);
// 				wpu_down_bkg=SysSF->at(46);
// 			}
// 
// 			wtrig_nom=TrigSF;//NEED TO PUT BACK IN THE UP AND DOWN VARIATIONS!
// 			wtrig_up=TrigSF;
// 			wtrig_down=TrigSF;
// 			
// 
// 			//std::vector<float> ChFlip_SF = correctChFlipSF(ChFlipSF, lepTLV, true); TBC!!!
// 			/*wchflip_nom=ChFlip_SF.at(0);
// 			wchflip_up=ChFlip_SF.at(1);
// 			wchflip_down=ChFlip_SF.at(2);*/
// 			wchflip_nom=ChFlipSF->at(0);
// 			wchflip_up=ChFlipSF->at(1);
// 			wchflip_down=ChFlipSF->at(2);
// 
// 			wpdf_up=1;
// 			wpdf_down=1;
// 		}
		//Rpc2L1b
		if(nSigLep>=2 && nBJets20>=1 && nJets40>=6 && met/meff>0.25)
		  isSR=true;
		//Rpc2L2b
		if(nSigLep>=2 && nBJets20>=2 && nJets25>=6 && met>300000 && meff>1400000 && met/meff>0.14)
		  isSR=true;
		//Rpc2L0b
		if(nSigLep>=2 && nBJets20==0 && nJets40>=6 && met>200000 && meff>1000000 && met/meff>0.2)
		  isSR=true;
		//Rpc3LSS1b
		if(nSigLep>=3 && nBJets20>=1 && is3LSS>0 && !isZee && met/meff>0.14){
		  isSR=true;
		}
		//Rpv2L
		if(nSigLep>=2 && nJets40>=6 && meff>2600000)
		  isSR=true;
		//VRWZ4j
		if(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=4 && mSFOS>81000 && mSFOS<101000 && met>50000 && met<250000 && meff>600000 && meff<1500000 ){
		  isVR=true;
		}
		//VRWZ5j
		if(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=5 && meff>400000 && met>50000 && mSFOS>81000 && mSFOS<101000 && meff<1500000 && met<250000){
		  isVR=true;
                }
		//VRttV
		if(nSigLep>=2 && NlepBL>=2 && nBJets20>=1 && nJets40>=3 && meff>600000 && meff<1500000 && met<250000 && isSS30 && dRl1j>1.1 && SumBJetPt/SumJetPt>0.4 && met/meff>0.1){
		  isVR=true;
		}
		//std::cout<<"evn: "<<evn64<<std::endl;
		//std::cout<<"nSigLep: "<<nSigLep<<"\t NlepBL: "<<NlepBL<<"\t nBJets20: "<<nBJets20<<"\t nJets40: "<<nJets40<<"\t meff: "<<meff<<"\t met: "<<met<<"\t isSS30: "<<isSS30<<"\t dRl1j: "<<dRl1j<<"\t SumBJetPt/SumJetPt: "<<SumBJetPt/SumJetPt<<"\t met/meff: "<<met/meff<<"\t nJets25: "<<nJets25<<"\t mSFOS: "<<mSFOS<<"\t is3LSS: "<<is3LSS<<"\t isZee: "<<isZee<<std::endl;
		//if(isSR || isVR){
        outTree->Fill();
		//std::cout<<"TAKEN!"<<std::endl;
		writtenEvents.push_back(evn64);
		//}
    }
	cout<<"No cut: "<<nocut<<"\t Bveto: "<<bveto<<endl;

	outTree->Write("",TObject::kOverwrite);
	cout<<"Done..."<<endl;
}


// open a file, get the tree and copy it to the output. attach proper lumi rescaling

void processFile(TString filein, TString fileout, TString treein, TString treeout, float xsec,float runc, Long64_t MC_ID, float lumi,  std::vector<Long64_t>& writtenEvents){

    std::cout<<"*** PROCESS FILE ***"<<std::endl;
	TFile* inFile=TFile::Open(filein, "READ");
	if(inFile->IsZombie() || !inFile)
	  std::cout<<"FILE IS IN ZOMBIE MODE!!!"<<std::endl;
	else
	  std::cout<<"inFile is open"<<std::endl;
	// need to remove "_nom" for input nominal sample
	TTree* inTree=(TTree*)inFile->Get(treein.ReplaceAll("_nom",""));
	
	if(!inTree){

		std::cout<<"Skipping corrupted file "<<std::endl;
		inFile->Close();
		return;

	}
	else
	  std::cout<<"inTree is open"<<std::endl;

	TFile* outFile=TFile::Open(fileout,"UPDATE");
	if(outFile->IsZombie() || !outFile)
	  std::cout<<"FILE IS IN ZOMBIE MODE!!!"<<std::endl;
	else
	  std::cout<<"outFile is open"<<std::endl;

	TTree* outTree=(TTree*)gROOT->FindObject(treeout);

	if(!outTree){
		cout<<"output tree not found, creating"<<endl;
		outTree=new TTree(treeout,treeout);
	}

	else cout<<RED<<"adding entries to existing tree "<<treeout<<RESET<<endl;

	TTree* control=(TTree*)inFile->Get("ControlTree");


	// get sum of weights. simple method using TTree::Draw and TArrayD crashes for large number of events
	Int_t NFiles=0, fileentry=0,raw=0,nofilt=0;
	float xAODWeightSum=0;
	float totw=0;
	float MCWeight=0;
	int MCId=0;

	// get filter eff for MC (data has xsec=0)

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
		if((MCId>=372444 && MCId<=372524) || (MCId>=373462 && MCId<=373478)) totw = totw*0.64; //b-filtering for Slepton grid
	}
	else totw=1;

	cout<<"Xsec: "<<xsec<<endl;
	copyTree(inTree,outTree,xsec,runc,MC_ID,totw,lumi,writtenEvents);

	outFile->Write("",TObject::kOverwrite);

	outFile->Close();
	inFile->Close();

}

// functions to retrieve cross section from common SUSYTools files

// this is just a grep. returns the line(s) to be parsed by the next function
TString lookForSample(TString runnumber){
    //uses perl regexp. matches trailing whitespace (space OR tab). requires begin of line.
    // old xsec file provided by Peter
    //gSystem->Exec("grep -r -P \"^"+runnumber+"\\s\" /publicfs/atlas/atlasnew/SUSY/users/liuy/SS3LWorkSpace/SS3L_Packages/SS3L_HF/FullRun2/prepare/xsections/mc15_13TeV/*  > xSec.out");
    
    // new xsec files from PMG tool, takes only the first appearance in case multiple xsec values listed in the file
    gSystem->Exec("grep -m 1 -r -P \"^"+runnumber+"\\s\" /afs/ific.uv.es/user/a/asantra/BTaggingWork/SS3Lep/SS3L_HF/FullRun2/prepare/xsections/PMGxsecDB_mc16.txt  > xSec.out");
    std::string line;
    std::ifstream result("xSec.out");
    std::getline(result, line);
    if(result.eof()){
        cout<<RED<<"ERROR cannot find cross section for sample "<<runnumber<<" will skip it "<<RESET<<endl;
        return "";
    }
    
    TString linestr(line);
    gSystem->Exec("rm xSec.out");
    // sanity check: there must be only one possibility
    std::getline(result, line);

    if(!result.eof()){
            cout<< "ERROR retrieving cross section: too many possibilities found. Aborting"<<endl;
            exit(1);
    }
    return linestr.ReplaceAll("\t"," ");
}

// parse the result from grep, check that there is only one possibility
void getNormFactor(TString runnumber, float& norm, float& runc, Long64_t &MC_ID){

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


// reads all files in a folder. adds the contents to outtreename in outfname
void readFolder(TString dirname, TString outfname, TString intreename, TString outtreename, std::vector<string>& sys, float targetlumi, float targetlumi_16a, float targetlumi_16d, float targetlumi_16e){

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) {


		bool isData=outtreename.Contains("data");
		bool isSignal=outtreename.BeginsWith("signal");

		// outer loop is on systematics
		std::cout<<"Systematics: "<<sys.size()<<std::endl;
		for(auto sysname: sys){	

			std::vector<Long64_t> writtenEvents;

			TString _intreename=intreename+"_"+sysname;

			TSystemFile *file;
			TString fname;
			TIter next(files);    
			// reinit the file iterator

			while ((file=(TSystemFile*)next())) {

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
						//cout<<" got cross normalization "<<norm<<endl;
					} else {
						// signal "data" to downstream functions using norm=0
						norm=0; runc=0; MC_ID=0;
					}

					bool Atlfast;
					Long64_t mcid = atol(runnumber);
					Atlfast = isAtlfast(mcid);

					// local (per file, per syst) output tree name
					TString _outtreename="signal"+runnumber;

					if(!isSignal){
						// error on bg MC treated differently
						runc=0;
						_outtreename=outtreename;
					} else {
						//std::cout<<"outtree name is "<<outtreename<<std::endl;
						// signal tree names are different, i.e. they are like signalXXXX_nom
						TObjArray* tokens=outtreename.Tokenize("_");
						if(tokens->GetEntries()==2){
							_outtreename="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
						}
					}

					//cout<<"run "<<runnumber<<" got norm factors "<<norm<<" "<<runc<<endl;
					//Setting correct lumi
					float tLumi = 0;
					//**** MC16A ****//
					if(fname.Contains("_r9364_")){
					  std::cout<<"taking MC16a Sample: "<<targetlumi_16a<<std::endl;
					  tLumi = targetlumi_16a;
					}
					//**** MC16D ****//
					if(fname.Contains("_r10201_")){
					  std::cout<<"taking MC16d Sample: "<<targetlumi_16d<<std::endl;
                                          tLumi= targetlumi_16d;
					}
					
					//**** MC16E ****//
					if(fname.Contains("_r10724_")){
					  std::cout<<"taking MC16e Sample: "<<targetlumi_16e<<std::endl;
                                          tLumi= targetlumi_16e;
                                        }

					std::cout<<"tLumi: "<<tLumi<<std::endl;

					// sanity check
					if( norm!=-1 && runc!=-1 ){
					  _outtreename=_outtreename+"_"+TString(sysname);
					  processFile(dirname+"/"+fname, outfname,_intreename,_outtreename,norm,runc,MC_ID,tLumi,writtenEvents);
					} // end sanity check

					bool doAFII(true);

					//AFII for Fullsim samples
					if(!Atlfast && sysname=="nom" && doAFII){
						// create dummy AFII syst for fullsim samples, copying the nominal one
						cout<<"The sample is fullsim, ading AFII trees..."<<endl;
						TString sysname2("JET_RelativeNonClosure_AFII__1up");
						TString sysname3("JET_RelativeNonClosure_AFII__1down");
						TString _outtreename2="signal"+runnumber;
						TString _outtreename3="signal"+runnumber;

						if(!isSignal){
							// error on bg MC treated differently
							runc=0;
							_outtreename2=outtreename;
							_outtreename3=outtreename;
						} else {

							if(tokens->GetEntries()==2){
								_outtreename2="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
								_outtreename3="signal"+runnumber+"_"+((TObjString*)(tokens->At(1)))->String();
							}
						}

						_outtreename2=_outtreename2+"_"+TString(sysname2);
						processFile(dirname+"/"+fname, outfname,_intreename,_outtreename2,norm,runc,MC_ID,tLumi,writtenEvents);
						_outtreename3=_outtreename3+"_"+TString(sysname3);
						processFile(dirname+"/"+fname, outfname,_intreename,_outtreename3,norm,runc,MC_ID,tLumi,writtenEvents);
					}

				}

				// for MC, clear event record once per file/sample
				if(!isData)
					writtenEvents.clear();

			} // loop on files in folder
		} // loop on systematics
	}

}

// main script
#include <limits>    


void lumiRescale(int inLumi, int outLumi){

	TString inbg=TString("background.")+Form("%d", inLumi)+TString(".root");
	TString outbg=TString("background.")+Form("%d", outLumi)+TString(".root");

	TString insig=TString("signal.")+Form("%d", inLumi)+TString(".root");
	TString outsig=TString("signal.")+Form("%d", outLumi)+TString(".root");

	TFile* finbg=TFile::Open(inbg,"READ");

	TIter next(finbg->GetListOfKeys());
	TKey *key;

	while ((key = (TKey*)next())) {
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

		//copyTree(inTree,outTree,float(outLumi)/float(inLumi));    

		outTree->Write(inTree->GetName(),TObject::kOverwrite);

		cout<<"Wrote tree "<<inTree->GetName()<<endl;

		foutsig->Write("",TObject::kOverwrite);
		foutsig->Close();

	}

}


void doAll(int targetlumi_16a, int targetlumi_16d, int targetlumi_16e){
    if(DEBUG)cout << "1 processed" << endl;
    int targetlumi = targetlumi_16a+targetlumi_16d+targetlumi_16e; 
    gROOT->ProcessLine("#include <vector>");

    gErrorIgnoreLevel=kError;
  
    const char* names[] = {
        "nom", // the nominal sample
//         "EG_RESOLUTION_ALL__1down",
//         "EG_RESOLUTION_ALL__1up",
//         "EG_SCALE_ALL__1down",
//         "EG_SCALE_ALL__1up",
//         //"JET_EtaIntercalibration_NonClosure__1up",
//         //"JET_EtaIntercalibration_NonClosure__1down",
//         "JET_EtaIntercalibration_NonClosure_highE__1up",
//         "JET_EtaIntercalibration_NonClosure_highE__1down",
//         "JET_EtaIntercalibration_NonClosure_negEta__1up",
//         "JET_EtaIntercalibration_NonClosure_negEta__1down",
//         "JET_EtaIntercalibration_NonClosure_posEta__1up",
//         "JET_EtaIntercalibration_NonClosure_posEta__1down",
//         "JET_GroupedNP_1__1up",
//         "JET_GroupedNP_1__1down",
//         "JET_GroupedNP_2__1up",
//         "JET_GroupedNP_2__1down",
//         "JET_GroupedNP_3__1up",
//         "JET_GroupedNP_3__1down",
//         //"JET_JER_SINGLE_NP__1up",
//         "JET_JER_DataVsMC__1up",
//         "JET_JER_DataVsMC__1down",                                                                                                      "JET_JER_EffectiveNP_1__1up",
//         "JET_JER_EffectiveNP_1__1down",
//         "JET_JER_EffectiveNP_2__1up",
//         "JET_JER_EffectiveNP_2__1down",
//         "JET_JER_EffectiveNP_3__1up",
//         "JET_JER_EffectiveNP_3__1down",
//         "JET_JER_EffectiveNP_4__1up",                                                                                                   "JET_JER_EffectiveNP_4__1down",
//         "JET_JER_EffectiveNP_5__1up",                                                                                                   "JET_JER_EffectiveNP_5__1down",
//         "JET_JER_EffectiveNP_6__1up",                                                                                                   "JET_JER_EffectiveNP_6__1down",
//         "JET_JER_EffectiveNP_7restTerm__1up",
//         "JET_JER_EffectiveNP_7restTerm__1down",
//         "JET_RelativeNonClosure_AFII__1up",
//         "JET_RelativeNonClosure_AFII__1down",
//         "MET_SoftTrk_ResoPara",
//         "MET_SoftTrk_ResoPerp",
//         "MET_SoftTrk_ScaleDown",
//         "MET_SoftTrk_ScaleUp",
//         "MUON_ID__1down",
//         "MUON_ID__1up",
//         "MUON_MS__1down",
//         "MUON_MS__1up",
//         "MUON_SAGITTA_RESBIAS__1down",
//         "MUON_SAGITTA_RESBIAS__1up",
//         "MUON_SAGITTA_RHO__1down",
//         "MUON_SAGITTA_RHO__1up",
//         "MUON_SCALE__1down",
//         "MUON_SCALE__1up"
    };

    std::vector<std::string> sysnames(names, names + sizeof(names)/sizeof(names[0]));
    // data has only nominal
    std::vector<std::string> sysnamesData; sysnamesData.push_back("nom");
    if(DEBUG)cout << "2 processed" << endl;

    ///////// old ntuples area, 21p2p70
	//TString basepathBkg("/publicfs/atlas/atlasnew/SUSY/users/liuy/SS3LWorkSpace/BackUps/HF_Moriond_2017/SS3l-FullLumi-HistFitter/DiLepton-NTUP/BKG/");
    //TString basepathBkg("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/MergedCommonSamples/");
    //TString basepathBkg("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/bRPV_MultiBosonSamples_400HFT_Rel21p2p70_MyBaselineLepton/Merged/");
    //TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/MergedCommonSamples/");
    
    /// this is just for Nicola's cutflow with WP3
    //TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP3/MergedCommonSamples/SUSYNicola/Merged/");
    
    /// this is just for Nicola's cutflow with WP1
    ///TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/MergedCommonSamples/SUSYNicolaAB21p2p70/");
    ///TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/bRPV_MultiBosonSamples_400HFT_Rel21p2p70_MyBaselineLepton/Merged/");/
    
    
    //// new ntuple area 21p2p102
    TString basepathBkg("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/Merged/");
    TString basepathSig("/lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/Merged/");
	
    cout<<"PathBkg: "<<basepathBkg<<endl;
    cout<<"PathSig: "<<basepathSig<<endl;

    if(DEBUG)cout << "3 processed" << endl;
    
    // Data  
    //readFolder(basepath+"DATA",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnamesData,1);
    //readFolder(basepathSig+"Data",TString("data.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","data",sysnamesData,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);

    
    // TString dirname, TString outfname, TString intreename, TString outtreename, std::vector<string>& sys, float targetlumi, float targetlumi_16a, float targetlumi_16d, float targetlumi_16e
    // Backgrounds
//     readFolder(basepathBkg+"ttH",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttH",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"RareMore",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","RareMore",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"VH",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","VH",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"3and4tSM",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ThreeAndFourTopSM",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"ttBar",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttBar",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"MultiBoson",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","MultiBoson",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"TTV",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","TTV",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
//     readFolder(basepathBkg+"Vjets",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","Vjets",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
    
    
    ////////// old folder  /////////
    // readFolder(basepathBkg+"ttBarMET",TString("background.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","ttBarMET",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
    ///lustre/ific.uv.es/grid/atlas/t3/asantra/bTaggingWP1/MergedCommonSamples/SUSYNicolaAB21p2p70
    // signals, for NICOLA
    // readFolder(basepathSig+"SUSYNicolaAB21p2p70",TString("signal.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
    
    /// all signals, will work on it after PMGxsec fix
    readFolder(basepathSig+"SUSY",TString("signal.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
    
    // only bRPV signals
    //readFolder(basepathSig+"BRPVSUSY",TString("signalBRPV.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
   
    
    /// problematic tree
    //readFolder(basepathSig+"Problem",TString("signal.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","signal",sysnames,targetlumi,targetlumi_16a,targetlumi_16d,targetlumi_16e);
    // just for test
    //readFolder(basepath+"TEST",TString("test.")+Form("%d", targetlumi)+TString(".root"),"DiLeptonTree_SRall","Test",sysnames,targetlumi);

}
