c#include <iomanip>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THStack.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TF1.h>
#include <TH1F.h>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

const UInt_t nSamples 		= 25;
const UInt_t gNCATEGORIES	= 3;
const UInt_t gNCHANNELS		= 4;
const TString gCatLabel	[gNCATEGORIES] 	= {"2lSS","3lSS","Total"};
const TString gChanLabel[gNCHANNELS] 	= {"MuMu","ElEl","ElMu","All"};
const TString sample	[nSamples] 		= {
  "TTWToLNu_ext2", "TTZToLLNuNu_ext", "TTZToLLNuNu_ext2", "TTZToQQ", "TTGJets"	// MC for comparison with data
  "TTGJets_ext","WW", "WW_ext",
  "TTJets_aMCatNLO", "DYJetsToLL_M10to50_aMCatNLO", 							// MC for control regions
  "DYJetsToLL_M10to50_aMCatNLO_ext", "TW", "TW_ext", "TbarW", "TbarW_ext",
  "WZTo3LNu", "WWTo2L2Nu", "ZZ", "ZZ_ext",
  "TTHonbb",																	// Signal samples
  "MuonEG", "DoubleMuon", "DoubleEG", "SingleElectron", "SingleMuon"			// Data samples
};
TString codepath 	= "/nfs/fanae/user/vrbouza/Documents/TFG/ttHAnalysis";
TString outputpath 	= "/nfs/fanae/user/vrbouza/Documents/TFG/Results";






//------------------------------------------------------------------------------
void ttHPlotter(TString samplename = "ZZ_ext") {
  TString codepath = "/nfs/fanae/user/vrbouza/Documents/TFG/ttHAnalysis";
  TString outputpath = "/nfs/fanae/user/vrbouza/Documents/TFG/Results";
  TString filename = "plot_" + samplename + ".pdf";
  TFile* f = TFile::Open(codepath + "/temp/" + "Tree_" + samplename + ".root");
  TH1F* h1;
  f->GetObject("fHDummy", h1);
  h1->SetTitle("Wololooooo");
  TCanvas* c = new TCanvas("c","c",800,600);
  h1->Draw("hist");
  c->Print(filename);
}

//void TTHPlotter(); // Execution of the function.
