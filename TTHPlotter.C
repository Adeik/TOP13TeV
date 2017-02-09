#include <iomanip>
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

void TTHPlotter(TString samplename = "ZZ") {
  TString codepath = "/nfs/fanae/user/vrbouza/Documents/TFG/TOP13TeV";
  TString outputpath = "/nfs/fanae/user/vrbouza/Documents/TFG/TOP13TeV/TTHStuff/plots";
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
