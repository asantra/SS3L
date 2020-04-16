
#include "TTree.h"
#include "TFile.h"
#include <iostream>
using namespace std;

TTree* harvesttree(const char* textfile=0) {
  const char* filename    = "$outfile";
  const char* description = "$description";
  TTree* tree = new TTree("tree","data from ascii file");
  Long64_t nlines(0);
  if (textfile!=0) {
    nlines = tree->ReadFile(textfile,description);
  } else if (filename!=0) {
    nlines = tree->ReadFile(filename,description);
  } else {
    cout << "WARNING: file name is empty. No tree is read." << endl;
  }
  tree->SetMarkerStyle(8);
  tree->SetMarkerSize(0.5);
  return tree;
}

void writetree() {
  TTree* tree = (TTree *)gDirectory->Get("tree");
  if (tree==0) {
    tree = harvesttree();
    if (tree==0) return;
  }
  TFile* file = TFile::Open("$outfile.root","RECREATE");
  file->cd();
  tree->Write();
  file->Close();
}

void summary_harvest_tree_description() {
  TTree* tree = harvesttree();
  gDirectory->Add(tree);
}
