/*==============================================================================

							ttHAnalyzer selector

==============================================================================*/
//------------------------------------------------------------------------------
//		Preprocessor directives
//------------------------------------------------------------------------------

#pragma once

// PAF inclusion
#include "PAFChainItemSelector.h"

// Packages inclusions
//#include "GlobalVariables.h"
#include "mt2.h"
#include "PUWeight.h"
#include "BTagSFUtil.h"
#include "SusyLeptonSF.h"

// ROOT inclusion
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
const Int_t nWeights = 248;
const Int_t nGenb = 0;
const Double_t pi = 3.1415926535897932384;
const Float_t Zm = 91.1876;

enum gChannel {
    channels_begin,
    Muon = channels_begin,
    Elec,
    ElMu,
    gNCHANNELS,
};

const TString gChanLabel[gNCHANNELS] = {"Muon","Elec","ElMu"};

enum gFPSwitch {
    SigSup,
    ZDecay,
    Sig
};

enum iCut {
    iDilepton,
    iZVeto,
    iMET,
    i2jets,
    i1btag,
    iDYVeto,
    iExact1btag,
    iExact2btag,
    iNCUTS
};

enum SR {
    AA, AB, AC, BA, BB, BC, CA, CB, CC,
    nSR
};

const TString SRlabel[nSR] = {
	"AA", "AB", "AC", "BA", "BB", "BC", "CA", "CB", "CC"
};

const TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag","DYVeto","Exact1btag","Exact2btag"};

enum gSystFlag {
    Norm,
    BtagUp,
    BtagDown,
    MisTagUp,
    MisTagDown,
    JESUp,
    JESDown,
    JER,
    LESUp,
    LESDown,
    /*  LepUp,
      LepDown,
      TrigUp,
      TrigDown,
    */
    PUUp,
    PUDown,
    TopPt,
    gNSYST
};

const TString SystName[gNSYST] = {
    "Normal",
    "BtagUp",
    "BtagDown",
    "MisTagUp",
    "MisTagDown",
    "JESUp",
    "JESDown",
    "JER",
    "LESUp",
    "LESDown",
    /*  "LepUp",
    "LepDown",
    "TrigUp",
    "TrigDown",*/
    //  "METUp",
    //  "METDown",
    "PUUp",
    "PUDown",
    "TopPt",
};
enum FakeSource {
    HF_mu,
    Other_mu,
    HF_el,
    Conv_el,
    Other_el,
    RightSign,
    WrongSign,
    gNFAKESOURCE
};

enum gNLOWeight {
    muR1muF1,
    muR1muF2,
    muR1muF05,
    muR2muF1,
    muR2muF2,
    muR2muF05,
    muR05muF1,
    muR05muF2,
    muR05muF05,
    gNWEIGHT
};

const TString WeiName[gNWEIGHT] = {
    "muR1muF1",
    "muR1muF2",
    "muR1muF05",
    "muR2muF1",
    "muR2muF2",
    "muR2muF05",
    "muR05muF1",
    "muR05muF2",
    "muR05muF05"
};

//	Binning
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

Int_t             getNFPtBins (gChannel chan);
const Double_t    *getFPtBins (gChannel chan);
Int_t             getNPPtBins (gChannel chan);
const Double_t    *getPPtBins (gChannel chan);
Int_t             getNEtaBins (gChannel chan);
const Double_t    *getEtaBins (gChannel chan);

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
		Int_t   ngenLep;
        Long_t  evt;
		Int_t   nJet;
        Int_t   nTauGood;
		Float_t genWeight;
		Float_t LepGood_px[30];
		Float_t LepGood_py[30];
		Float_t LepGood_pz[30];
		Float_t LepGood_energy[30];
		Float_t LepGood_pt[30];
		Float_t LepGood_etaSc[30];
		Float_t LepGood_eta[30];
		Float_t LepGood_dxy[30];
		Float_t LepGood_dz[30];
		Float_t LepGood_relIso03[30];
		Float_t LepGood_relIso04[30];
		Int_t   LepGood_charge[30];
		Int_t   LepGood_pdgId[30];
		Float_t LepGood_z[30];					// NEW
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
		Float_t genLep_pt[50];
		Float_t genLep_eta[50];
		Float_t genLep_phi[50];
		Float_t genLep_mass[50];
		Int_t   genLep_pdgId[50];
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
    	virtual void InitialiseGenHistos();
    	virtual void InitialiseDYHistos();
    	virtual void InitialiseYieldsHistos();
    	virtual void InitialiseKinematicHistos();
    	virtual void InitialiseSystematicHistos();

        //	Filling methods
		//----------------------------------------------------------------------
		void 	FillYieldsHistograms(gChannel, iCut, gSystFlag);
		void 	FillYields(gSystFlag sys=Norm);
		void 	FillDYHistograms();
		void 	FillKinematicHistos(gChannel,iCut);

    	// Saving
		//----------------------------------------------------------------------
    	void 	WriteHistos();
    	void 	WriteValidationsHistos(){};

		////////////////////////////////////////////////////////////////////////
		//	   Leptons, jets and MET selection
		////////////////////////////////////////////////////////////////////////
		void 	SelectedGenLepton();
		Int_t  	getSelectedLeptons();
		void 	ScaleLeptons(Int_t);
		std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);

        //  Muons
        //----------------------------------------------------------------------
		Bool_t	IsTightMuon(UInt_t, Float_t ptcut=20.);
		Bool_t	IsFakeableMuon(UInt_t, Float_t ptcut=20.);
		Bool_t	IsLooseMuon(UInt_t, Float_t ptcut=20.);
		Float_t getMuonIso(Int_t);

        //  Electrons
        //----------------------------------------------------------------------
		Bool_t	IsTightElectron(UInt_t,Float_t ptcut=20.,Int_t an=2);
		Bool_t	IsFakeableElectron(UInt_t,Float_t ptcut=20.);
		Bool_t	IsLooseElectron(UInt_t,Float_t ptcut=20.);
		Float_t getElecIso(Int_t);
		Float_t getEACorrection(Float_t);
		Bool_t	getMultiIso(UInt_t );

        //  Taus
        //----------------------------------------------------------------------
		Bool_t 	IsGoodTau(UInt_t iTau, Float_t ptcut);

        //  MET
        //----------------------------------------------------------------------
        Bool_t	METFilter();
        void 	propagateMET(TLorentzVector,TLorentzVector);
        void 	ScaleMET(Int_t);

        //  Jets
        //----------------------------------------------------------------------
		Int_t 	getSelectedJets();
		Bool_t	IsGoodJet(UInt_t, Float_t ptcut=25.);
		Bool_t 	IsGoodJetforprecuts(UInt_t, Float_t ptcut=25.);
		std::vector<Int_t> CleanedJetIndices(Float_t);
		void 	SmearJetPts(Int_t);

		////////////////////////////////////////////////////////////////////////
		//	   Events selection
		////////////////////////////////////////////////////////////////////////
		Int_t  	IsDileptonEvent();
		Bool_t	IsMuMuEvent();                                                  // REDEF
		Bool_t	IsElMuEvent();                                                  // REDEF
		Bool_t	IsElElEvent();                                                  // REDEF
	    Bool_t 	IsSSEvent();					                                        // NEw
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
		Int_t   getLeadingJetbTag();
		Float_t getDRClosestJet(TLorentzVector);
		Float_t getDPhiClosestJet(TLorentzVector);
		Float_t getMET();
		Float_t getMETPhi();
		Float_t getHT();
		Float_t getMHT();					                                    // NEW
		Float_t getMETLD();					                                    // NEW
		Float_t getJetPtIndex(UInt_t);
		Float_t getJetEtaIndex(UInt_t);
		Float_t getBtagJetPtIndex(UInt_t);
		Float_t getErrPt(Float_t,Float_t);
		Float_t getJERScaleUp(Int_t);
        Float_t getJERScale(Int_t);
		Float_t getJERScaleDown(Int_t);
		Float_t getSF(gChannel);
		Float_t getLeptonError(gChannel);
		Float_t getTriggerError(gChannel);
		Float_t getTopPtSF();
		Float_t getTopD();
		Float_t getDeltaPhillJet();
		Float_t weightNvtx(Int_t);
		Float_t getMT(gChannel);
		Float_t getMT2(TLorentzVector plep1, TLorentzVector plep2, TLorentzVector pmet, Float_t mass);
		Float_t getMT2ll(gChannel);
		Float_t getMT2b(gChannel);
		Float_t getMT2lb(gChannel);
		Float_t getMeff();
		Float_t getDPhiLepJet();
		Float_t getDelPhill();
		Float_t getDPhiJetMet();
		Float_t getDPhiLepMet();
		Float_t getDPhibMet();
		Float_t getMinDPhiMetJets();

		TLorentzVector getPtllb();

	protected:
		////////////////////////////////////////////////////////////////////////
		//		Data members
		////////////////////////////////////////////////////////////////////////
		//	Input parameters
		//----------------------------------------------------------------------
	    TString gSampleName;
		TString gfileSuffix;
		Float_t gWeight;
		Float_t gLumiForPU;
		Float_t gTotalLumi;
		Int_t   gSysSource;
		Int_t   gSysDirection;
		Bool_t  gDoSystStudies;
		Bool_t  gIsData;
		Bool_t  gUseCSVM;
		Bool_t  gIsMCatNLO;
		Bool_t  gIsT2tt;

		//	PU and SF
		//----------------------------------------------------------------------
		PUWeight *fPUWeight;      //The PU weight utility
		PUWeight *fPUWeightUp;    //The PU weight utility
		PUWeight *fPUWeightDown;  //The PU weight utility
		//BTagSFUtil *fBTagSF[5]; //The new BTag SF utility
		BTagSFUtil *fBTagSFnom ;
		BTagSFUtil *fBTagSFbUp ;
		BTagSFUtil *fBTagSFbDo ;
		BTagSFUtil *fBTagSFlUp ;
		BTagSFUtil *fBTagSFlDo ;
		BTagSFUtil *medfBTagSFnom ;				// NEW
		BTagSFUtil *medfBTagSFbUp ;
		BTagSFUtil *medfBTagSFbDo ;
		BTagSFUtil *medfBTagSFlUp ;
		BTagSFUtil *medfBTagSFlDo ;
		BTagSFUtil *losfBTagSFnom ;
		BTagSFUtil *losfBTagSFbUp ;
		BTagSFUtil *losfBTagSFbDo ;
		BTagSFUtil *losfBTagSFlUp ;
		BTagSFUtil *losfBTagSFlDo ;
		SusyLeptonSF *fLeptonSF;
		TRandom3 *fRand3;

		//	EventWeight
		//----------------------------------------------------------------------
		Float_t EventWeight;
		Float_t PUSF;
		Bool_t  fChargeSwitch;

		//	Histograms and trees
		//----------------------------------------------------------------------
		TH1F*   fHDummy;
		TH1F*   hWeight;
		TH1F*   fHyields     [gNCHANNELS][gNSYST];
		TH1F*   fHWeightyield[gNCHANNELS][gNWEIGHT];
		TH1F*   fHSSyields   [gNCHANNELS][gNSYST];
		TH1F*   fHTopPtWeight;

		//	Generation
		//----------------------------------------------------------------------
		std::vector<Double_t>       Gen_Muon_Charge;
		std::vector<Double_t>       Gen_Elec_Charge;
		std::vector<TLorentzVector> Gen_Muon;
		std::vector<TLorentzVector> Gen_Elec;

		std::vector<Int_t>          NGen_Jet;
		std::vector<Int_t>          NGen_b;

		std::vector<Double_t>       PtGen_Jet;
		std::vector<Double_t>       PtGen_b;

		Int_t 	nGenElec;
		Int_t	nGenMuon;
		Int_t	nGenTau;
		Int_t	nGenLepton;
		Int_t 	nTauElec;
		Int_t 	nTauMuon;
		Int_t 	nSGenMuon;
		Int_t 	nSGenElec;

		//	General
		//----------------------------------------------------------------------
		Int_t 	nGoodVertex;
		Float_t nVertex;
		Int_t 	nBtags;
		Int_t 	nJets;
		Int_t 	nTightMuon;					// NEW
		Int_t 	nFakeableMuon;				// NEW
		Int_t 	nLooseMuon;					// NEW
		Int_t 	nTightElec;					// NEW
		Int_t 	nFakeableElec;				// NEW
		Int_t 	nLooseElec;					// NEW
		Int_t 	nElec;
		Int_t 	nLeptons;
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
