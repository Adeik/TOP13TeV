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

void ttHPlotter() {
	const UInt_t nmcSamples 	= 20;
	const UInt_t ndataSamples 	= 5;
	enum gCategories {
	    categories_begin,
	    twolSS = categories_begin,
	    threelSS,
	    Total,
	    gNCATEGORIES,
	};
	enum gChannel {
	    channels_begin,
	    MuMu = channels_begin,
	    ElEl,
	    ElMu,
	    All,
	    gNCHANNELS,
	};
	const TString gCatLabel	[gNCATEGORIES] 	= {"2lSS","3lSS","Total"};
	const TString gChanLabel[gNCHANNELS] 	= {"MuMu","ElEl","ElMu","All"};
	const TString mcsample	[nmcSamples] 	= {
		"TTWToLNu_ext2", "TTZToLLNuNu_ext", "TTZToLLNuNu_ext2", "TTZToQQ", "TTGJets"	// MC for comparison with data
	  	"TTGJets_ext","WW", "WW_ext",
	  	"TTJets_aMCatNLO", "DYJetsToLL_M10to50_aMCatNLO", 							// MC for control regions
	  	"DYJetsToLL_M10to50_aMCatNLO_ext", "TW", "TW_ext", "TbarW", "TbarW_ext",
	  	"WZTo3LNu", "WWTo2L2Nu", "ZZ", "ZZ_ext",
	  	"TTHonbb"																	// Signal samples
	};
	const TString datasample	[ndataSamples] 		= {
		"MuonEG", "DoubleMuon", "DoubleEG", "SingleElectron", "SingleMuon"			// Data samples
	};
	TString codepath 	= 	"/nfs/fanae/user/vrbouza/Documents/TFG/ttHAnalysis";
	TString outputpath 	= 	"/nfs/fanae/user/vrbouza/Documents/TFG/Results";
	TString filename	=	"Results.pdf";

	// Initializing THStacks of all the histograms.
	//------------------------------------------------------------------------------
	THStack*   fHSEvents    		[gNCATEGORIES][gNCHANNELS]; // Events
	THStack*   fHSTightLep			[gNCATEGORIES][gNCHANNELS]; // Yields
	THStack*   fHSFakeLep			[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSLooseLep			[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSTau				[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSJet				[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSMedBTagJet		[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSLosBTagJet		[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSPtLeading			[gNCATEGORIES][gNCHANNELS]; // Kinematic
	THStack*   fHSPtSubLeading		[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSPtSubSubLeading	[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSMET				[gNCATEGORIES][gNCHANNELS]; // MET
	THStack*   fHSMHT				[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSMETLD				[gNCATEGORIES][gNCHANNELS];
	THStack*   fHSChargeSum			[gNCATEGORIES][gNCHANNELS]; // Misc
	THStack*   fHSMass				[gNCATEGORIES][gNCHANNELS];

	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			fHSEvents    		[icat][ichan]	=	new THStack("HS_Events_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NEvents_"+gCatLabel[icat]+"_"+gChanLabel[ichan]); // Events
			fHSTightLep			[icat][ichan]	=	new THStack("HS_TightLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NTightLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan]); // Yields
			fHSFakeLep			[icat][ichan]	=	new THStack("HS_FakeLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NFakeLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSLooseLep			[icat][ichan]	=	new THStack("HS_LooseLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NLooseLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSTau				[icat][ichan]	=	new THStack("HS_Tau_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NTau_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSJet				[icat][ichan]	=	new THStack("HS_Jet_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSMedBTagJet		[icat][ichan]	=	new THStack("HS_MedBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NMedBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSLosBTagJet		[icat][ichan]	=	new THStack("HS_LosBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "NLosBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSPtLeading			[icat][ichan]	=	new THStack("HS_PtLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "PtLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan]); // Kinematic
			fHSPtSubLeading		[icat][ichan]	=	new THStack("HS_PtSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "PtSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSPtSubSubLeading	[icat][ichan]	=	new THStack("HS_PtSubSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "PtSubSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSMET				[icat][ichan]	=	new THStack("HS_MET_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "MET_"+gCatLabel[icat]+"_"+gChanLabel[ichan]); // MET
			fHSMHT				[icat][ichan]	=	new THStack("HS_MHT_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "MHT_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			fHSMETLD				[icat][ichan]	=	new THStack("HS_METLD_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "METLD_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
			if (icat == twolSS 		&& ichan != All) 	continue;
			fHSChargeSum			[icat][ichan]	=	new THStack("HS_ChargeSum_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "ChargeSum_"+gCatLabel[icat]+"_"+gChanLabel[ichan]); // Misc
			fHSMass				[icat][ichan]	=	new THStack("HS_Mass_"+gCatLabel[icat]+"_"+gChanLabel[ichan], "Mass_"+gCatLabel[icat]+"_"+gChanLabel[ichan]);
		}
	}

	TH1F* histEvents;
	TH1F* histTightLep;
	TH1F* histFakeLep;
	TH1F* histLooseLep;
	TH1F* histTau;
	TH1F* histJet;
	TH1F* histMedBTagJet;
	TH1F* histLosBTagJet;
	TH1F* histPtLeading;
	TH1F* histPtSubLeading;
	TH1F* histPtSubSubLeading;
	TH1F* histMET;
	TH1F* histMHT;
	TH1F* histMETLD;
	TH1F* histChargeSum;
	TH1F* histMass;

	for (UInt_t isample = 0; isample < nmcSamples; isample++) {
		TFile* f = TFile::Open(codepath + "/temp/" + "Tree_" + mcsample[isample] + ".root");
		for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
			for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
				if (icat == threelSS 	&& ichan != All) 	continue;
				if (icat == Total 		&& ichan != All) 	continue;
				f	->	GetObject("H_Events_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histEvents); // Events
				f	->	GetObject("H_TightLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histTightLep); // Yields
				f	->	GetObject("H_FakeLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histFakeLep);
				f	->	GetObject("H_LooseLep_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histLooseLep);
				f	->	GetObject("H_Tau_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histTau);
				f	->	GetObject("H_Jet_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histJet);
				f	->	GetObject("H_MedBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histMedBTagJet);
				f	->	GetObject("H_LosBTagJet_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histLosBTagJet);
				f	->	GetObject("H_PtLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histPtLeading); // Kinematic
				f	->	GetObject("H_PtSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histPtSubLeading);
				f	->	GetObject("H_PtSubSubLeading_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histPtSubSubLeading);
				f	->	GetObject("H_MET_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histMET); // MET
				f	->	GetObject("H_MHT_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histMHT);
				f	->	GetObject("H_METLD_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histMETLD);

				if (!(icat == twolSS 		&& ichan != All)) {
					f	->	GetObject("H_ChargeSum_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histChargeSum); // Misc
					f	->	GetObject("H_Mass_"+gCatLabel[icat]+"_"+gChanLabel[ichan],histMass);
				}

				fHSEvents    		[icat][ichan]	-> Add(histEvents); // Events
				fHSTightLep			[icat][ichan]	-> Add(histTightLep); // Yields
				fHSFakeLep			[icat][ichan]	-> Add(histFakeLep);
				fHSLooseLep			[icat][ichan]	-> Add(histLooseLep);
				fHSTau				[icat][ichan]	-> Add(histTau);
				fHSJet				[icat][ichan]	-> Add(histJet);
				fHSMedBTagJet		[icat][ichan]	-> Add(histMedBTagJet);
				fHSLosBTagJet		[icat][ichan]	-> Add(histLosBTagJet);
				fHSPtLeading		[icat][ichan]	-> Add(histPtLeading); // Kinematic
				fHSPtSubLeading		[icat][ichan]	-> Add(histPtSubLeading);
				fHSPtSubSubLeading	[icat][ichan]	-> Add(histPtSubSubLeading);
				fHSMET				[icat][ichan]	-> Add(histMET); // MET
				fHSMHT				[icat][ichan]	-> Add(histMHT);
				fHSMETLD			[icat][ichan]	-> Add(histMETLD);
				if (!(icat == twolSS 		&& ichan != All)) {
					fHSChargeSum		[icat][ichan]	-> Add(histChargeSum); // Misc
					fHSMass				[icat][ichan]	-> Add(histMass);
				}
				// histEvents				= 0;
				// histTightLep			= 0;
				// histFakeLep				= 0;
				// histLooseLep			= 0;
				// histTau					= 0;
				// histJet					= 0;
				// histMedBTagJet			= 0;
				// histLosBTagJet			= 0;
				// histPtLeading			= 0;
				// histPySubLeading		= 0;
				// histPtSubSubLeading		= 0;
				// histMET					= 0;
				// histMHT					= 0;
				// histMETLD				= 0;
				// histChargeSum			= 0;
				// histMass				= 0;
			}
		}
		// f						= 0;
	}

	// Drawing
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSEvents    		[icat][ichan]	-> Draw("hist"); // Events
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"Events"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"Events"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"Events"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSTightLep			[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"TightLep"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"TightLep"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"TightLep"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSFakeLep			[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"FakeLep"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"FakeLep"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"FakeLep"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSLooseLep			[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"LooseLep"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"LooseLep"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"LooseLep"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSTau				[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"Tau"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"Tau"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"Tau"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSJet				[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"Jet"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"Jet"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"Jet"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSMedBTagJet		[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"MedBTagJet"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"MedBTagJet"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"MedBTagJet"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSLosBTagJet		[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"LosBTagJet"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"LosBTagJet"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"LosBTagJet"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSPtLeading		[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"PtLeading"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"PtLeading"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"PtLeading"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSPtSubLeading		[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"PtSubLeading"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"PtSubLeading"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"PtSubLeading"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSPtSubSubLeading	[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"PtSubSubLeading"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"PtSubSubLeading"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"PtSubSubLeading"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSMET				[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"MET"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"MET"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"MET"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSMHT				[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"MHT"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"MHT"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"MHT"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSMETLD			[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"METLD"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"METLD"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"METLD"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			if (icat == twolSS 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSChargeSum		[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"ChargeSum"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"ChargeSum"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"ChargeSum"+".pdf");
			}
		}
	}
	for (UInt_t icat = 0; icat < gNCATEGORIES; icat++) {
		for (UInt_t ichan = 0; ichan < gNCHANNELS; ichan++) {
			if (icat == threelSS 	&& ichan != All) 	continue;
			if (icat == Total 		&& ichan != All) 	continue;
			if (icat == twolSS 		&& ichan != All) 	continue;
			TCanvas *c = new TCanvas("c", "c", 800, 600);
			fHSMass				[icat][ichan]	-> Draw("hist");
			if (icat == Total && ichan == All) {
				c->Print(outputpath+"/"+"Mass"+".pdf"+")");
			} else if (icat == categories_begin && ichan == channels_begin) {
				c->Print(outputpath+"/"+"Mass"+".pdf"+"(");
			}
			else {
				c->Print(outputpath+"/"+"Mass"+".pdf");
			}
		}
	}



}


//------------------------------------------------------------------------------
/*void ttHPlotter(TString samplename = "ZZ_ext") {
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
*/
