#include <iomanip>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLine.h"
#include "TChain.h"
#include "TLatex.h"
#include "stdio.h"
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
#include <TGraph.h>
#include <TRandom3.h>
#include <fstream>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

TString codepath = "/nfs/fanae/user/vrbouza/Documents/TFG/TOP13TeV"
TString outputpath = "/nfs/fanae/user/vrbouza/Documents/TFG/TOP13TeV/TTHStuff/plots"

void TTHPlotter(TString samplename = "ZZ") {
  TFile* f;
  TFile f = TFile::Open(codepath + "/temp/" + "Tree_" + samplename + ".root");
  TH1F* h1;
  f->GetObject("fHDummy", h1);
  h1->SetTitle("Wololooooo");
  TCanvas* c = SetCanvas();
  

}

TTHPlotter(); // Execution of the function.
