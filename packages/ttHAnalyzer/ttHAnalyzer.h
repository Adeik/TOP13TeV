/*==============================================================================

							ttHAnalyzer selector

==============================================================================*/
//------------------------------------------------------------------------------
//		Preprocessor directives
//------------------------------------------------------------------------------

#pragma once

// PAF inclusion
#include "PAFChainItemSelector.h"

// Analysis packages inclusion
//#include "GlobalVariables.h"
#include "mt2.h"
#include "PUWeight.h"
#include "BTagSFUtil.h"
#include "SusyLeptonSF.h"

// ROOT packages inclusion
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TRandom3.h"

// C++ packages inclusion
#include <vector>


//------------------------------------------------------------------------------
//		Enumerations, constants and other variable type declarations
//------------------------------------------------------------------------------
// Constants definitions
const Double_t pi = 3.1415926535897932384;
const Float_t Zm = 91.1876;
const TString gChanLabel[gNCHANNELS] = {"Muon","Elec","ElMu"};
const TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag","DYVeto","Exact1btag","Exact2btag"};

// New variable type definitions
enum gChannel {
    channels_begin,
    Muon = channels_begin,
    Elec,
    ElMu,
    gNCHANNELS,
};

//	Binning-related constants
const Int_t gNMuFPtBins = 6;
const Int_t gNMuPPtbins = 10;
const Int_t gNMuEtabins = 5;
const Int_t gNElFPtBins = 8;
const Int_t gNElPPtbins = 10;
const Int_t gNElEtabins = 5;
const Int_t gNElCMIdbins = 2;

const Double_t gMuFPtBins[gNMuFPtBins+1]	= {20., 25., 30., 35., 40., 50., 60.};						// Muon binning
const Double_t gMuPPtbins[gNMuPPtbins+1]	= {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.};
const Double_t gMuEtabins[gNMuEtabins+1]	= {0., 0.5, 1.0, 1.479, 2.0, 2.5};

const Double_t gElFPtBins[gNElFPtBins+1]  = {20., 25., 30., 40., 50., 60., 70., 80., 100.};			// Electron binning
const Double_t gElPPtbins[gNElPPtbins+1]  = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.};
const Double_t gElEtabins[gNElEtabins+1]  = {0., 0.5, 1.0, 1.479, 2.0, 2.5};

//------------------------------------------------------------------------------
//		Classes declarations
//------------------------------------------------------------------------------
// Particles class
class lepton {
    public:
        //lepton(){}
        //lepton(const lepton &l): p(l.p), charge(l.charge), type(l.type), index(l.index){ };
        lepton(TLorentzVector vec = TLorentzVector(0,0,0,0), Int_t ch = 0, Int_t ty = -1, Int_t ind = -1){
            p = vec;
            charge = ch;
            type = ty;
            index = ind;
        }
    	TLorentzVector p;
    	Int_t charge;
    	Int_t type; // -1(unknown), 0(tight muon), 1(tight electron); 2(fakeable muon), 3(fakeable electron); 4(loose muon), 5(loose electron)
    	Int_t index;
};

class jet {
	public:
    	jet(){};
    	jet(TLorentzVector vec, Bool_t btag, Int_t ind){
    		p = vec;
    		isbtag = btag;
    		index = ind;
    	};
    	TLorentzVector p;
    	Bool_t isbtag;
    	Int_t index;
};

// Analysis class
class ttHAnalyzer : public PAFChainItemSelector {
	public:
		////////////////////////////////////////////////////////////////////////
		//		Initial declarations
		////////////////////////////////////////////////////////////////////////
		// Constructor and destructor
		ttHAnalyzer();
		virtual ~ttHAnalyzer() {}

		// Core PAF-analysis methods
		virtual void Initialise();
		virtual void InsideLoop();
		virtual void Summary();

        // For printing
        void 	CoutEvent(ULong_t en = 0, TString t = " ");

		////////////////////////////////////////////////////////////////////////
		//		Trees-related declarations
		////////////////////////////////////////////////////////////////////////
		//	Methods
		//----------------------------------------------------------------------
		void 	GetTreeVariables();
		void 	GetParameters();

        //	Tree variables
		//----------------------------------------------------------------------
		Int_t   nLepGood;
        Int_t   nTauGood;
		Int_t   nJet;
        Long_t  evt;
		Float_t LepGood_px[30];
		Float_t LepGood_py[30];
		Float_t LepGood_pz[30];
		Float_t LepGood_energy[30];
		Float_t LepGood_pt[30];
		Float_t LepGood_eta[30];
		Float_t LepGood_dxy[30];
		Float_t LepGood_dz[30];
		Int_t   LepGood_charge[30];
		Int_t   LepGood_pdgId[30];
		Float_t LepGood_sip3d[30];				// NEW
		Float_t LepGood_miniRelIso[30];			// NEW
		Float_t LepGood_jetBTagCSV[30];			// NEW
		Int_t 	LepGood_mediumMuonId[30];		// NEW
		Float_t LepGood_mvaTTH[30];				// NEW
		Float_t LepGood_jetPtRatiov2[30];		// NEW
		Float_t LepGood_mvaIdSpring15[30];		// NEW
		Float_t LepGood_sigmaIEtaIEta[30];		// NEW
		Float_t LepGood_hadronicOverEm[30];		// NEW
		Float_t LepGood_dEtaScTrkIn[30];		// NEW
		Float_t LepGood_dPhiScTrkIn[30];		// NEW
		Float_t LepGood_eInvMinusPInv[30];		// NEW
		Float_t LepGood_convVeto[30];			// NEW
		Int_t 	LepGood_lostHits[30];			// NEW
		Int_t 	LepGood_tightCharge[30];		// NEW
		Float_t LepGood_jetDR[30];				// NEW
		Float_t Jet_px[50];
		Float_t Jet_py[50];
		Float_t Jet_pz[50];
		Float_t Jet_energy[50];
		Float_t Jet_eta[50];
		Float_t Jet_btagCSV[50];
		Int_t	TauGood_idDecayModeNewDMs[30];	// NEW
		Float_t	TauGood_pt[30];					// NEW
		Float_t	TauGood_eta[30];				// NEW
		Float_t	TauGood_phi[30];				// NEW
		Float_t	TauGood_mass[30];				// NEW
		Int_t	TauGood_idCI3hit[30];			// NEW

        ////////////////////////////////////////////////////////////////////////
		//		Histogram-related methods declarations
		////////////////////////////////////////////////////////////////////////
        // Initialising
		//----------------------------------------------------------------------
    	virtual void InitialiseYieldsHistos();

        //	Filling methods
		//----------------------------------------------------------------------
		void 	FillYieldsHistograms(gChannel);
		void 	FillYields();

    	// Saving
		//----------------------------------------------------------------------
    	void 	WriteHistos();                                                   // SIN DEFINIR
    	void 	WriteValidationsHistos();                                        // SIN DEFINIR

		////////////////////////////////////////////////////////////////////////
		//	   Leptons and jets selection
		////////////////////////////////////////////////////////////////////////
		Int_t  	getSelectedLeptons();
		std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);

        //  Muons
        //----------------------------------------------------------------------
		Bool_t	IsTightMuon(UInt_t, Float_t ptcut=20.);
		Bool_t	IsFakeableMuon(UInt_t, Float_t ptcut=20.);
		Bool_t	IsLooseMuon(UInt_t, Float_t ptcut=20.);

        //  Electrons
        //----------------------------------------------------------------------
		Bool_t	IsTightElectron(UInt_t,Float_t ptcut=20.,Int_t an=2);
		Bool_t	IsFakeableElectron(UInt_t,Float_t ptcut=20.);
		Bool_t	IsLooseElectron(UInt_t,Float_t ptcut=20.);

        //  Taus
        //----------------------------------------------------------------------
		Bool_t 	IsGoodTau(UInt_t iTau, Float_t ptcut);

        //  Jets
        //----------------------------------------------------------------------
		Int_t 	getSelectedJets();
		Bool_t	IsGoodJet(UInt_t, Float_t ptcut=25.);
		Bool_t 	IsGoodJetforprecuts(UInt_t, Float_t ptcut=25.);

		////////////////////////////////////////////////////////////////////////
		//	   Events selection
		////////////////////////////////////////////////////////////////////////
		Int_t  	IsDileptonEvent();
		Bool_t	IsMuMuEvent();                                                  // REDEF
		Bool_t	IsElMuEvent();                                                  // REDEF
		Bool_t	IsElElEvent();                                                  // REDEF
	    Bool_t 	IsSSEvent();					                                // NEw
	    Bool_t 	Is2lSSEvent();			                                        // NEw
	    Bool_t 	Is3lSSEvent();			                                        // NEw

		Bool_t 	PassesPreCuts();					                            // NEw

        ////////////////////////////////////////////////////////////////////////
		//	   Trigger methods
		////////////////////////////////////////////////////////////////////////
		Bool_t	triggermumuSS();
		Bool_t	triggereeSS();
		Bool_t	triggeremuSS();
		Bool_t	trigger3l4l();

        ////////////////////////////////////////////////////////////////////////
		//	   Set/reset methods
		////////////////////////////////////////////////////////////////////////
		void 	SetOriginalObjects();
		void 	SetEventObjects();
		void 	ResetOriginalObjects();
		void 	ResetHypLeptons();
        void    setMET(Float_t);

        ////////////////////////////////////////////////////////////////////////
		//	   Get methods
		////////////////////////////////////////////////////////////////////////
		Int_t   getNJets();
		Int_t   getNBTags();
		Float_t getMET();
		Float_t getHT();
		Float_t getMHT();					                                    // NEW
		Float_t getMETLD();
		Float_t getSF(gChannel);

	protected:
		////////////////////////////////////////////////////////////////////////
		//		Data members
		////////////////////////////////////////////////////////////////////////
		//	Input parameters
		//----------------------------------------------------------------------
	    TString gSampleName;
		Bool_t  gIsData;
		Float_t gWeight;
		Float_t gLumiForPU;
		Float_t gTotalLumi;
		Bool_t  gUseCSVM;

		//	PU and SF
		//----------------------------------------------------------------------
		PUWeight      *fPUWeight;      //The PU weight utility
		BTagSFUtil    *fBTagSFnom ;
		BTagSFUtil    *medfBTagSFnom ;				// NEW
		BTagSFUtil    *losfBTagSFnom ;
		SusyLeptonSF  *fLeptonSF;
		TRandom3      *fRand3;

		//	EventWeight
		//----------------------------------------------------------------------
		Float_t EventWeight;
		Float_t PUSF;

		//	Histograms and trees
		//----------------------------------------------------------------------
		TH1F*   fHDummy;
		TH1F*   hWeight;
		TH1F*   fHyields     [gNCHANNELS];
		TH1F*   fHSSyields   [gNCHANNELS];

		//	General
		//----------------------------------------------------------------------
		UInt_t 	nJets;
		UInt_t 	nBTags;
		UInt_t 	nTightMuon;					// NEW
		UInt_t 	nFakeableMuon;				// NEW
		UInt_t 	nLooseMuon;					// NEW
		UInt_t 	nTightElec;					// NEW
		UInt_t 	nFakeableElec;				// NEW
		UInt_t 	nLooseElec;					// NEW
		UInt_t 	nElec;
		UInt_t 	nLeptons;
		UInt_t 	njpt;						// NEW

		Float_t MET;
		Float_t MET_Phi;
		Float_t MHT;						// NEW

		lepton  fHypLepton1;
		lepton  fHypLepton2;
		std::vector<lepton> Lepton;
		std::vector<lepton> LooseLepton;	// NEW
		std::vector<lepton> FakeableLepton;	// NEW
		std::vector<lepton> TightLepton;	// NEW
		std::vector<jet>    Jet;

		std::vector<Float_t> JetEt;
		std::vector<Float_t> JetPt;
		std::vector<Float_t> JetPhi;
		std::vector<Float_t> MuPx;
		std::vector<Float_t> MuPy;
		std::vector<Float_t> MuPz;
		std::vector<Float_t> MuEnergy;
		std::vector<Float_t> ElPx;
		std::vector<Float_t> ElPy;
		std::vector<Float_t> ElPz;
		std::vector<Float_t> ElEnergy;

		ClassDef(ttHAnalyzer,0);  // ROOT definition as class
};
