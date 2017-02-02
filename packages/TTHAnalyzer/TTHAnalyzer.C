//========================================================
//  TTHAnalyzer selector
//========================================================

#include "TTHAnalyzer.h"
#include <iostream>
#include <math.h>

ClassImp(TTHAnalyzer);
const float gJetEtCut = 30.;

//#define DEBUG

//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void TTHAnalyzer::GetParameters(){
  gSampleName    = GetParam<TString>("sampleName");
  gIsData        = GetParam<bool>("IsData");
  gWeight        = GetParam<float>("weight"); // cross section / events in the sample
  gLumiForPU     = GetParam<float>("LumiForPU");
  gTotalLumi     = GetParam<float>("TotalLumi");
  gDoSystStudies = GetParam<bool>("DoSystStudies");
  gUseCSVM       = GetParam<bool>("UseCSVM");
  gStopMass      = GetParam<Int_t>("stopMass");
  gLspMass       = GetParam<Int_t>("lspMass");
  gIsT2tt        = false;
  if(gSampleName.BeginsWith("T2tt")) gIsT2tt = true;

  gIsMCatNLO     = GetParam<bool>("IsMCatNLO");
  gCreateTree    = GetParam<bool>("CreateTree");

  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gSampleName = %s",gSampleName.Data()));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gIsData = %d",gIsData ));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gWeight = %e", gWeight));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gLumiForPU = %f", gLumiForPU));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gTotalLumi = %f", gTotalLumi));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gDoSystStudies = %d", gDoSystStudies));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gUseCSVM = %d",gUseCSVM ));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gStopMass = %i", gStopMass));
  PAF_INFO("TTHAnalyzer::GetParameters()", Form("gLspMass = %i",gLspMass ));
}

//-----------------------------------------------------------------------------------
// GetTreeVariables
//-----------------------------------------------------------------------------------
void TTHAnalyzer::GetTreeVariables(){
  nLepGood             = Get<Int_t>("nLepGood");
	nJet                 = Get<Int_t>("nJet");
  evt                  = Get<ULong64_t>("evt");
	if(!gIsData){
		ngenLep              = Get<Int_t>("ngenLep");
		genWeight            = Get<Float_t>("genWeight");
	}
	for(int k = 0; k<nLepGood; k++){
    LepGood_px[k]      = Get<Float_t>("LepGood_px", k);
    LepGood_py[k]      = Get<Float_t>("LepGood_py", k);
    LepGood_pz[k]      = Get<Float_t>("LepGood_pz", k);
    LepGood_energy[k]  = Get<Float_t>("LepGood_energy", k);
    LepGood_pt[k]      = Get<Float_t>("LepGood_pt", k);
    LepGood_etaSc[k]   = Get<Float_t>("LepGood_etaSc", k);
    LepGood_eta[k]     = Get<Float_t>("LepGood_eta", k);
    LepGood_dxy[k]     = Get<Float_t>("LepGood_dxy", k);
    LepGood_dz[k]      = Get<Float_t>("LepGood_dz", k);
    LepGood_relIso03[k]= Get<Float_t>("LepGood_relIso03", k);
    LepGood_relIso04[k]= Get<Float_t>("LepGood_relIso04", k);
    LepGood_pdgId[k]   = Get<Int_t>("LepGood_pdgId", k);
    LepGood_charge[k]  = Get<Int_t>("LepGood_charge", k);
  }
  for(int k = 0; k<nJet; k++){
    Jet_px[k]          = Get<Float_t>("Jet_px", k);
    Jet_py[k]          = Get<Float_t>("Jet_py", k);
    Jet_pz[k]          = Get<Float_t>("Jet_pz", k);
    Jet_energy[k]      = Get<Float_t>("Jet_energy", k);
    Jet_eta[k]         = Get<Float_t>("Jet_eta", k);
    Jet_btagCSV[k]     = Get<Float_t>("Jet_btagCSV", k);
  }
	if(!gIsData){
		for(int k = 0; k<ngenLep; k++){
			genLep_pdgId[k]    = Get<Int_t>("genLep_pdgId", k);
			genLep_pt[k]       = Get<Float_t>("genLep_pt", k);
			genLep_eta[k]      = Get<Float_t>("genLep_eta", k);
			genLep_phi[k]      = Get<Float_t>("genLep_phi", k);
			genLep_mass[k]     = Get<Float_t>("genLep_mass", k);
		}
	}
}

int getNFPtBins(gChannel chan){ // fake ratios
  if(chan == Muon || chan == ElMu) return gNMuFPtBins;
  if(chan == Elec)                 return gNElFPtBins;
  else return -99;
};

const double *getFPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuFPtBins;
  else                             return gElFPtBins;
  //  if(chan == Elec)                 return gElFPtBins;
  //  else return *-999;
};

int getNPPtBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuPPtbins;
  if(chan == Elec)                 return gNElPPtbins;
  else return -99;
};

const double *getPPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuPPtbins;
  else                             return gElPPtbins;
  //  if(chan == Elec)                 return gElPPtbins;
  //  else return -99;
};

int getNEtaBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuEtabins;
  if(chan == Elec)                 return gNElEtabins;
  else return -99;
};

const double *getEtaBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuEtabins;
  else                             return gElEtabins;
  //  if(chan == Elec)                 return gElEtabins;
  //  else return *-99.;
};

//------------------------------------------------------------------------------------
// TTHAnalyzer class constructor (make sure the pointers are initialized to zero)
//------------------------------------------------------------------------------------
TTHAnalyzer::TTHAnalyzer() : PAFChainItemSelector() {
	fHDummy = 0;
	hWeight = 0;
	fHTopPtWeight = 0;
	// fHnGenEle = 0;
	// fHnGenMuo = 0;
	// fHGenElePt = 0;
	// fHGenMuoPt = 0;

	for (unsigned int ichan = 0; ichan < gNCHANNELS; ichan++) {
		for (unsigned int isyst = 0; isyst < gNSYST; isyst++) {
			fHyields     [ichan][isyst] = 0;
			fHSSyields   [ichan][isyst] = 0;
		}
		for (unsigned int iweight = 0; iweight < gNWEIGHT; iweight++) {
			fHWeightyield[ichan][iweight] = 0;
		}
	}
}

//-------------------------------------------------------------------
// Initialise
//-------------------------------------------------------------------
void TTHAnalyzer::Initialise() {
	PAF_INFO("TTHAnalyzer", "+ Initializing...");
	//PAF_INFO("TTHAnalyzer", "+ Initializing paramenters...");
	GetParameters();
	//PAF_INFO("TTHAnalyzer", "+ Sumw2 set for all histograms...");
	TH1::SetDefaultSumw2();
	fHDummy = CreateH1F("fHDummy","",1,0,1);
	//PAF_INFO("TTHAnalyzer", "+ Initialise Yield histograms...");
	InitialiseTree();
	InitialiseYieldsHistos();
	//PAF_INFO("TTHAnalyzer", "+ Initialise Kinematic histograms...");
	InitialiseKinematicHistos();
	if (!gIsData) {
		//PAF_INFO("TTHAnalyzer", "+ Initialise Gen histograms...");
		InitialiseGenHistos();
	}
	//PAF_INFO("TTHAnalyzer", "+ Initialise other histograms...");
	fHTopPtWeight  = CreateH1F("H_TopPtWeight" ,"TopPt Weight",100, 0, 2);


	if (gSampleName == "DoubleMuon"      ||
			gSampleName == "DoubleEG"        ||
			gSampleName == "SingleMu"        ||
			gSampleName == "SingleElectron"  ||
			gSampleName == "MuonEG"	       ||
			gSampleName.Contains("TTJets") || gSampleName.Contains("TTbar") ||
			gSampleName.Contains("DY") ||
			gSampleName.Contains("ZJets")){
		//	PAF_INFO("TTHAnalyzer", "+ Initialise Drell-Yan histograms...");
		InitialiseDYHistos();
	}
	PAF_INFO("TTHAnalyzer", "+ Initialise histograms for systematics studies...");
	InitialiseSystematicHistos();

	//	PU Reweight
	//--------------------------------------
	//PAF_INFO("TTHAnalyzer", "+ Initialise Pile-Up reweighting tool...");
  fPUWeight     = new PUWeight(gLumiForPU, Spring2016_25ns_poisson_OOTPU, "2016_ichep");
  if (!gIsData) {
    fPUWeightUp   = new PUWeight(18494.9,  Spring2016_25ns_poisson_OOTPU, "2016_ichep"); //  18494.9
    fPUWeightDown = new PUWeight(20441.7,  Spring2016_25ns_poisson_OOTPU, "2016_ichep"); //  20441.7
  }


	//if (gUseCSVM) fBTagSF   = new BTagSFUtil("CSVM","ABCD");//ReReco
	//else          fBTagSF   = new BTagSFUtil("CSVT","ABCD");//ReReco

	//PAF_INFO("TTHAnalyzer", "+ Initialise b-tag scale factors...");
	if (gUseCSVM){
		fBTagSFnom = new BTagSFUtil("mujets", "CSVv2", "Medium",  0);
		fBTagSFbUp = new BTagSFUtil("mujets", "CSVv2", "Medium",  1);
		fBTagSFbDo = new BTagSFUtil("mujets", "CSVv2", "Medium", -1);
		fBTagSFlUp = new BTagSFUtil("mujets", "CSVv2", "Medium",  3);
		fBTagSFlDo = new BTagSFUtil("mujets", "CSVv2", "Medium", -3);
	}
	else{
		fBTagSFnom = new BTagSFUtil("mujets", "CSVv2", "Tight",  0);
		fBTagSFbUp = new BTagSFUtil("mujets", "CSVv2", "Tight",  1);
		fBTagSFbDo = new BTagSFUtil("mujets", "CSVv2", "Tight", -1);
		fBTagSFlUp = new BTagSFUtil("mujets", "CSVv2", "Tight",  3);
		fBTagSFlDo = new BTagSFUtil("mujets", "CSVv2", "Tight", -3);
	}

	//PAF_INFO("TTHAnalyzer", "+ Initialise lepton scale factors...");
	fLeptonSF = new SusyLeptonSF();

	//PAF_INFO("TTHAnalyzer", "+ Initialise random 3...");
	fRand3 = new TRandom3(50);

	// No systematics activaded...
	gSysSource = Norm;
	PAF_INFO("TTHAnalyzer", "+ Initialisation DONE.");
}

void TTHAnalyzer::InitialiseTree(){
    fTree = CreateTree("sTopTree","Optimization tree");

    fTree->Branch("TWeight",      &TWeight,      "TWeight/F");
    fTree->Branch("TIsDoubleMuon",&TIsDoubleMuon,"TIsDoubleMuon/I");
    fTree->Branch("TIsDoubleElec",&TIsDoubleElec,"TIsDoubleElec/I");
    fTree->Branch("TIsElMu",      &TIsElMu,      "TIsElMu/I");
    fTree->Branch("TNJets",       &TNJets,       "TNJets/I");
    fTree->Branch("TNJetsBtag",   &TNJetsBtag,   "TNJetsBtag/I");

    fTree->Branch("TMET",         &TMET,         "TMET/F");
    fTree->Branch("TMT2ll",       &TMT2ll,       "TMT2ll/F");
    fTree->Branch("TMT2bb",       &TMT2bb,       "TMT2bb/F");
    fTree->Branch("TMT2lblb",     &TMT2lblb,     "TMT2lblb/F");
    fTree->Branch("TMll",         &TMll,         "TMll/F");
    fTree->Branch("TMeff",        &TMeff,        "TMeff/F");
    fTree->Branch("TPtllb",       &TPtllb,       "TPtllb/F");
    fTree->Branch("THT",          &THT,          "THT/F");

    fTree->Branch("TMET_Phi",     &TMET_Phi,     "TMET_Phi/F");
    fTree->Branch("TdPhiPtllbMET",&TdPhiPtllbMET,"TdPhiPtllbMET/F");
    fTree->Branch("TMinDPhiMetJets",&TMinDPhiMetJets,"TMinDPhiMetJets/F");
    fTree->Branch("TdPhiJetMet",  &TdPhiJetMet,  "TdPhiJetMet/F");
    fTree->Branch("TdPhiLepMet",  &TdPhiLepMet,  "TdPhiLepMet/F");
    fTree->Branch("TdPhiLepJet",  &TdPhiLepJet,  "TdPhiLepJet/F");
    fTree->Branch("TdPhill",      &TdPhill,      "TdPhill/F");

    fTree->Branch("TLep1_Px",     &TLep1_Px,     "TLep1_Px/F");
    fTree->Branch("TLep1_Py",     &TLep1_Py,     "TLep1_Py/F");
    fTree->Branch("TLep1_Pz",     &TLep1_Pz,     "TLep1_Pz/F");
    fTree->Branch("TLep1_E" ,     &TLep1_E ,     "TLep1_E/F");
    fTree->Branch("TLep1_Carge",  &TLep1_Charge, "TLep1_Charge/F");
    fTree->Branch("TLep2_Px",     &TLep2_Px,     "TLep2_Px/F");
    fTree->Branch("TLep2_Py",     &TLep2_Py,     "TLep2_Py/F");
    fTree->Branch("TLep2_Pz",     &TLep2_Pz,     "TLep2_Pz/F");
    fTree->Branch("TLep2_E" ,     &TLep2_E ,     "TLep2_E/F");
    fTree->Branch("TLep2_Carge",  &TLep2_Charge, "TLep2_Charge/F");

    fTree->Branch("TJet_isBJet",       TJet_isBJet,       "TJet_isBJet[TNJets]/I");
    fTree->Branch("TJet_Px",           TJet_Px,           "TJet_Px[TNJets]/F");
    fTree->Branch("TJet_Py",           TJet_Py,           "TJet_Py[TNJets]/F");
    fTree->Branch("TJet_Pz",           TJet_Pz,           "TJet_Pz[TNJets]/F");
    fTree->Branch("TJet_E",            TJet_E,            "TJet_E[TNJets]/F");
}

void TTHAnalyzer::InitialiseGenHistos(){

}

void TTHAnalyzer::InitialiseDYHistos(){

}

void TTHAnalyzer::InitialiseYieldsHistos(){
	hWeight = CreateH1F("hWeight","",200,0,1);
	//++ Yields histograms
		fHyields[Muon][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHyields[Elec][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHSSyields[Muon][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHSSyields[Elec][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHyields[ElMu][Norm]   = CreateH1F("H_Yields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHSSyields[ElMu][Norm] = CreateH1F("H_SSYields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);

	if (gDoSystStudies){
		for (size_t chan=0; chan<gNCHANNELS; chan++){
			for (size_t sys=1; sys<gNSYST; sys++){
				fHyields[chan][sys]   = CreateH1F("H_Yields_"+gChanLabel[chan]+"_"+SystName[sys],"",iNCUTS,-0.5,iNCUTS-0.5);
				fHSSyields[chan][sys] = CreateH1F("H_SSYields_"+gChanLabel[chan]+"_"+SystName[sys],"", iNCUTS, -0.5, iNCUTS-0.5);
			}
		}
	}
	for (size_t chan=0; chan<gNCHANNELS; chan++){
		for (int wei = 0; wei < gNWEIGHT; ++wei){
			fHWeightyield[chan][wei] = CreateH1F("H_Yields_Wei_"+gChanLabel[chan]+"_"+WeiName[wei],"",iNCUTS,-0.5,iNCUTS-0.5);
		}
	}
	for (size_t chan=0; chan<gNCHANNELS; chan++){
		for (size_t cut=0; cut<iNCUTS; cut++){
		}
	}
}

void TTHAnalyzer::InitialiseKinematicHistos(){
	//  PAF_DEBUG("TTHAnalyzer::InitialiseKinematicHistos()",Form("nWeights = %i", nWeights));
	//++ Kinematic histograms
	for (size_t ch=0; ch<gNCHANNELS; ch++){
		for (size_t cut=0; cut<iNCUTS; cut++){
			//PAF_DEBUG("TTHAnalyzer::InitialiseKinematicHistos()",Form("cut = %i", cut));
		}
	}
}

void TTHAnalyzer::InitialiseSystematicHistos(){

}

//---------------------------------------------------------------------------------------------------
// Set objets, to be called once per event, saving information in tmp vectors for systematic studies.
//---------------------------------------------------------------------------------------------------
void TTHAnalyzer::SetOriginalObjects(){
	ResetHypLeptons();
	gSysSource = Norm;

	// SAVING ORIGINAL VALUES FOR MET, JET, LEPTONS for SYST
	JetEt.clear();
	JetPt.clear();
	JetPhi.clear();
	MuPx.clear();
	MuPy.clear();
	MuPz.clear();
	MuEnergy.clear();
	ElPx.clear();
	ElPy.clear();
	ElPz.clear();
	ElEnergy.clear();
	MET  = 0.;
	MET_Phi = 0.;
	int k = 0;
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){
		j.SetPxPyPzE(Jet_px[i], Jet_py[i], Jet_pz[i], Jet_energy[i]);
		JetEt.push_back(j.Et());
		JetPt.push_back(j.Pt());
		JetPhi.push_back(j.Phi());
	}
	for (Int_t i=0; i < nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 11){
			ElPx.push_back(LepGood_px[i]);
			ElPy.push_back(LepGood_py[i]);
			ElPz.push_back(LepGood_pz[i]);
			ElEnergy.push_back(LepGood_energy[i]);
		}
	}
	for (Int_t i=0; i<nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 13){
			MuPx.push_back(LepGood_px[i]);
			MuPy.push_back(LepGood_py[i]);
			MuPz.push_back(LepGood_pz[i]);
			MuEnergy.push_back(LepGood_energy[i]);
		}
	}
	MET     = Get<Float_t>("met_pt"); //met
	MET_Phi = Get<Float_t>("met_phi"); //met
}

void TTHAnalyzer::SetEventObjects(){
	ResetHypLeptons();

	fChargeSwitch = false;
	EventWeight = 1.;

	// USEFUL COUNTERS
	nGenLepton = 0;
	nGenElec   = 0;
	nGenMuon   = 0;
	nGenTau    = 0;
	nTauElec   = 0;
	nTauMuon   = 0;
	nGoodVertex = 0;
	nVertex     = 0;
	nBtags      = 0;
	nJets       = 0;
	nMuon       = 0;
	nElec       = 0;
	nLeptons    = 0;

	//// READ AND SAVE OBJETS...
	Jet.clear();
	Lepton.clear();

	nLeptons = getSelectedLeptons();
	nJets    = getSelectedJets();
	nBtags   = getNBTags();
}

void TTHAnalyzer::ResetOriginalObjects(){
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){
		j.SetPxPyPzE(Jet_px[i],	Jet_py[i], Jet_pz[i],	Jet_energy[i]);
		JetEt[i]  = j.Et();
		JetPt[i]  = j.Pt();
		JetPhi[i] = j.Phi();
	}
	int k = 0;
	for (Int_t i=0; i<nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 11){
			ElPx[k] = LepGood_px[i];
			ElPy[k] = LepGood_py[i];
			ElPz[k] = LepGood_pz[i];
			ElEnergy[k] = LepGood_energy[i];
			k++;
		}
	}
	k = 0;
	for (Int_t i=0; i<nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 13){
			MuPx[k] = LepGood_px[i];
			MuPy[k] = LepGood_py[i];
			MuPz[k] = LepGood_pz[i];
			MuEnergy[k] = LepGood_energy[i];
			k++;
		}
	}
	setMET(Get<Float_t>("met_pt")); //met
}

void TTHAnalyzer::ResetHypLeptons(){
  TLorentzVector vec(0., 0., 0., 0.);
  fHypLepton1 = lepton(vec, 0, -1, -1);
  fHypLepton2 = lepton(vec, 0, -1, -1);
}
void TTHAnalyzer::CoutEvent(long unsigned int en, TString t){
  //if(en == 1000599168 || en == 1268707665 || en == 157395642 || en == 2726847580 || en == 42879335){
  //if(en == 1519610198 || en == 1559039433 || en == 998619292 || en == 1206329870 || en == 295644557 || en == 686746673 || en == 99957372 || en == 126485808 || en == 249297855){
  if(en == 1347253329 || en == 960559657){
    cout << t << endl;
  }
  else return;
}

void TTHAnalyzer::SetTreeVariables(gChannel chan){
  TWeight     = EventWeight;
  TNJets      = getNJets();
  TNJetsBtag  = getNBTags();

  TIsDoubleElec = 0; TIsDoubleMuon = 0; TIsElMu = 0;
  if(chan == Muon) TIsDoubleMuon = 1;
  if(chan == Elec) TIsDoubleElec = 1;
  if(chan == ElMu) TIsElMu       = 1;
  TMET          = getMET();
	TMT2ll        = getMT2ll(chan);
	TMT2bb        = getMT2b(chan);
	TMT2lblb      = getMT2lb(chan);
  TMll          = (fHypLepton1.p+fHypLepton2.p).M();
  TPtllb        = getPtllb().Pt();
  TMeff         = getMeff();
  THT           = getHT();
  TdPhiPtllbMET = getDPhibMet();
  TMinDPhiMetJets = getMinDPhiMetJets();
  TdPhiJetMet   = getDPhiJetMet();
  TdPhiLepMet   = getDPhiLepMet();
  TdPhiLepJet   = getDPhiLepJet();
  TdPhill       = getDelPhill();
  TMET_Phi      = getMETPhi();

  TLep1_Px      = fHypLepton1.p.Px();
  TLep1_Py      = fHypLepton1.p.Py();
  TLep1_Pz      = fHypLepton1.p.Pz();
  TLep1_E       = fHypLepton1.p.E();
  TLep1_Charge  = fHypLepton1.charge;
  TLep2_Px      = fHypLepton2.p.Px();
  TLep2_Py      = fHypLepton2.p.Py();
  TLep2_Pz      = fHypLepton2.p.Pz();
  TLep2_E       = fHypLepton2.p.E();
  TLep2_Charge  = fHypLepton2.charge;

  for(int k = 0; k<40; k++){
    if(k<TNJets){
      TJet_Px[k]           = Jet[k].p.Px();
      TJet_Py[k]           = Jet[k].p.Py();
      TJet_Pz[k]           = Jet[k].p.Pz();
      TJet_E[k]            = Jet[k].p.E();
      TJet_isBJet[k]       = Jet[k].isbtag;
    }
    else{
      TJet_Px[k]           = 0;
      TJet_Py[k]           = 0;
      TJet_Pz[k]           = 0;
      TJet_E[k]            = 0;
      TJet_isBJet[k]       = 0;
    }
  }
}


//-----------------------------------------------------------------------
// InsideLoop
//-----------------------------------------------------------------------
void TTHAnalyzer::InsideLoop() {
    if(gIsT2tt){
    //   cout << "GenSusyMStop = " << Get<Int_t>("GenSusyMStop") << ", gStopMass = " << gStopMass << ", GenSusyMNeutralino = " << Get<Int_t>("GenSusyMNeutralino") << ", gLspMass = " << gLspMass << endl;
    if((Get<Int_t>("GenSusyMStop") != gStopMass) || (Get<Int_t>("GenSusyMNeutralino") != gLspMass)) return;
    }
	fHDummy->Fill(0.5);
	if (!METFilter()) return;

    CoutEvent(evt, Form("Event number = %li", evt));
	// Calculate PU Weight
	if (!gIsData) PUSF = fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); //True       //nTruePU

	// Init data members ........................................................
	GetTreeVariables();
	SetOriginalObjects();
	SetEventObjects();

	// // Get number of generated leptons ........................................................

	// Fill Yields ...............................................................
#ifdef DEBUG
	cout << "N Leptons: " << Lepton.size() << endl;
	cout << "PassTriggerEMu/EE/MuMu= "
		<< triggermumuSS() << "/"
		<< triggermumuSS() << "/"
		<< triggermumuSS() << endl;
	cout << "Is ElMu/ElEl/MuMu Event= "
		<< IsElMuEvent() << "/"
		<< IsElElEvent() << "/"
		<< IsMuMuEvent() << endl;
#endif
	FillYields();

	// Get SS Yields...
	fChargeSwitch = true;
	FillYields();
	fChargeSwitch = false;

	// Fill DY DD histograms
	if (gSampleName == "DoubleMuon"      ||
			gSampleName == "DoubleEG"        ||
			gSampleName == "SingleMu"        ||
			gSampleName == "SingleElectron"  ||
			gSampleName == "MuonEG"	       ||
			gSampleName.Contains("ZJets")   ||
			gSampleName.Contains("DY"))  {
		FillDYHistograms();
	}

	// Fill Yields for syst. studies (only for MC) ..............................
  if (gIsData)         return;
  if (!gDoSystStudies) return;

  // B-tagging systematics .................................................................
  ResetOriginalObjects();
  gSysSource = BtagUp;
  SetEventObjects();
  FillYields(BtagUp);
  fChargeSwitch = true;
  FillYields(BtagUp); /* Get SS yields....*/
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = BtagDown;
  SetEventObjects();
  FillYields(BtagDown);
  fChargeSwitch = true;
  FillYields(BtagDown); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = MisTagUp;
  SetEventObjects();
  FillYields(MisTagUp);
  fChargeSwitch = true;
  FillYields(MisTagUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = MisTagDown;
  SetEventObjects();
  FillYields(MisTagDown);
  fChargeSwitch = true;
  FillYields(MisTagDown); /// Get SS yields....
  fChargeSwitch = false;

	// JES/JER sytematics ....................................................................
  ResetOriginalObjects();
  SmearJetPts(1);
  gSysSource = JESUp;
  SetEventObjects();
  FillYields(JESUp);
  fChargeSwitch = true;
  FillYields(JESUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  SmearJetPts(2);
  gSysSource = JESDown;
  SetEventObjects();
  FillYields(JESDown);
  fChargeSwitch = true;
  FillYields(JESDown); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  SmearJetPts(3);
  gSysSource = JER;
  SetEventObjects();
  FillYields(JER);
  fChargeSwitch = true;
  FillYields(JER); /// Get SS yields....
  fChargeSwitch = false;

  // Lepton Scale  sytematics ....................................................................
  ResetOriginalObjects();
  ScaleLeptons(1); //up
  gSysSource = LESUp;
  SetEventObjects();
  FillYields(LESUp);
  fChargeSwitch = true;
  FillYields(LESUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  ScaleLeptons(2); //down
  gSysSource = LESDown;
  SetEventObjects();
  FillYields(LESDown);
  fChargeSwitch = true;
  FillYields(LESDown); /// Get SS yields....
  fChargeSwitch = false;

  // Pile Up sytematics ....................................................................
  ResetOriginalObjects();
  if (!gIsData)
    PUSF = fPUWeightUp->GetWeight(Get<Float_t>("nTrueInt")); //nTruePU
  gSysSource = PUUp;
  SetEventObjects();
  FillYields(PUUp);

  ResetOriginalObjects();
  if (!gIsData)
    PUSF = fPUWeightDown->GetWeight(Get<Float_t>("nTrueInt")); //nTruePU
  gSysSource = PUDown;
  SetEventObjects();
  FillYields(PUDown);

  // Top PT ...............................................................................
  ResetOriginalObjects();
  gSysSource = TopPt;
  SetEventObjects();
  FillYields(TopPt);
  fChargeSwitch = true;
  FillYields(TopPt); /// Get SS yields....
  fChargeSwitch = false;
}

void TTHAnalyzer::Summary(){}

//------------------------------------------------------------------------------
// TRIGGER INFORMATION
//------------------------------------------------------------------------------
/*	The triggers for each region are defined exactly the same as in the CMS Draft Note of
03-03-2016 ("Search for ttH in multilepton final states at 13 TeV").
*/
bool TTHAnalyzer::triggermumuSS() {
    if(!gIsData) return true;
    //return true;
    Bool_t pass = false;
    if (gIsData){
        pass =  (Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") &&
                Get<Int_t>("HLT_BIT_IsoMu20_v") &&
                Get<Int_t>("HLT_BIT_HLT_IsoTkMu20_v"));
    }else{
        pass = false;
    }
    return pass;
}


bool TTHAnalyzer::triggereeSS(){
    if(!gIsData) return true;
    //return true;
    Bool_t pass = false;
    if (gIsData){
        pass =  (Get<Int_t>("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&
               Get<Int_t>("HLT_BIT_HLT_Ele25_WPTight_Gsf_v") &&
               Get<Int_t>("HLT_BIT_HLT_Ele45_WPLoose_Gsf_v"));
    }else{
        pass = false;
    }
    return pass;
}


bool TTHAnalyzer::triggeremuSS(){
    if(!gIsData) return true;
    //return true;
    Bool_t pass = false;
    if (gIsData){
        pass =  (Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVl_v") &&
                Get<Int_t>("HLT_BIT_IsoMu20_v") &&
                Get<Int_t>("HLT_BIT_HLT_IsoTkMu80_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele25_WPTight_Gsf_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele45_WPLoose_Gsf_v"));
    }else{
        pass = false;
    }
    return pass;
}

bool TTHAnalyzer::trigger3l4l(){
    if(!gIsData) return true;
    //return true;
    Bool_t pass = false;
    if (gIsData){
        pass =  (Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVl_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") &&
                Get<Int_t>("HLT_BIT_IsoMu20_v") &&
                Get<Int_t>("HLT_BIT_HLT_IsoTkMu20_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele25_WPTight_Gsf_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele45_WPLoose_Gsf_v") &&
                Get<Int_t>("HLT_BIT_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v") &&
                Get<Int_t>("HLT_BIT_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") &&
                Get<Int_t>("HLT_BIT_HLT_TripleMu_12_10_5_v") &&
                Get<Int_t>("HLT_BIT_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"));
    }else{
        pass = false;
    }
    return pass;
}


//------------------------------------------------------------------------------
// Get METHODS
//------------------------------------------------------------------------------
float TTHAnalyzer::getHT(){
  float ht(0);
  for (unsigned int i=0; i<Jet.size(); i++) ht+=Jet[i].p.Pt();
  return ht;
}
float TTHAnalyzer::getJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return Jet[ind].p.Pt();
}
float TTHAnalyzer::getJetEtaIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return TMath::Abs(Jet[ind].p.Eta());
}
float TTHAnalyzer::getBtagJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  Int_t btagInd = 0;
  if (ind==0) btagInd = getLeadingJetbTag();
  else  return -999.;
  return Jet[btagInd].p.Pt();
}

float TTHAnalyzer::getMT(gChannel chan){
  float ptl1 = fHypLepton1.p.Pt();
	float ptl2 = fHypLepton2.p.Pt();
	float dphi = getDelPhill();
	return TMath::Sqrt(2*ptl1*ptl2*(1-TMath::Cos(dphi)));
}

float TTHAnalyzer::getDPhiLepJet(){
	if (fHypLepton1.index == -1) return -999.; if (fHypLepton2.index == -1) return -999.;
	// Int_t ij = getLeadingJetbTag(); if (ij < 0) return -999.;
	if(Jet.size()<1) return -999.;
	TLorentzVector jet = Jet[0].p;
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(jet));
}

float TTHAnalyzer::getMinDPhiMetJets(){
	if (fHypLepton1.index == -1) return -999.; if (fHypLepton2.index == -1) return -999.;
	if(Jet.size()<2) return -999.;
	TLorentzVector jet1 = Jet[0].p;
  TLorentzVector jet2 = Jet[1].p;
  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
  float MinDelta = TMath::Min(TMath::Abs(pmet.DeltaPhi(jet1)), TMath::Abs(pmet.DeltaPhi(jet2)) );
  return MinDelta;
}

float TTHAnalyzer::getDelPhill(){ return fHypLepton1.p.DeltaPhi(fHypLepton2.p);}

float TTHAnalyzer::getDPhiJetMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return getDPhiClosestJet(pmet);
}

float TTHAnalyzer::getDPhiLepMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(pmet));
}

float TTHAnalyzer::getMT2(TLorentzVector plep1, TLorentzVector plep2, TLorentzVector pmet, float mass){
  double pa[3]; double pb[3]; double pmiss[3];
  pmiss[0] = 0.; // irrelevant
  pmiss[1] = pmet.Px(); pmiss[2] = pmet.Py();
  pa[0] = 0.; pa[1] = plep1.Px(); pa[2] = plep1.Py();
  pb[0] = 0.; pb[1] = plep2.Px(); pb[2] = plep2.Py();
  mt2 MT2bisect;
  MT2bisect.set_momenta(pa, pb, pmiss);
  MT2bisect.set_mn(mass); // testmass
  float MT2 = MT2bisect.get_mt2();
  return MT2;
}

float TTHAnalyzer::getMT2ll(gChannel chan){
  TLorentzVector plep1, plep2;
  if (chan == Muon) {
    plep1.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(1), MuPy.at(1), MuPz.at(1), MuEnergy.at(1));
  }
  if (chan == Elec) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(ElPx.at(1), ElPy.at(1), ElPz.at(1), ElEnergy.at(1));
  }
  if (chan == ElMu) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
  }
  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  return getMT2(plep1, plep2, pmet, 0.);
}

float TTHAnalyzer::getMT2b(gChannel chan){
  if (getNJets() < 2) return -1;
  TLorentzVector plep1, plep2;
  if (chan == Muon) {
    plep1.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(1), MuPy.at(1), MuPz.at(1), MuEnergy.at(1));
  }
  if (chan == Elec) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(ElPx.at(1), ElPy.at(1), ElPz.at(1), ElEnergy.at(1));
  }
  if (chan == ElMu) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
  }
  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  TLorentzVector pjet0 = Jet[0].p; TLorentzVector pjet1 = Jet[1].p;
  TLorentzVector lv;
  lv.SetPtEtaPhiM(TMath::Sqrt( pow(pmet.Px()+(plep1+plep2).Px(), 2) +  pow(pmet.Py()+(plep1+plep2).Py(), 2)), 0., TMath::ATan2(pmet.Py()+(plep1+plep2).Py(),pmet.Px()+(plep1+plep2).Px()), 0.);
  return getMT2(pjet0, pjet1, lv, 80.398);
}

float TTHAnalyzer::getMT2lb(gChannel chan){
  if (getNJets() < 2) return -1;
  TLorentzVector plep1, plep2;
  if (chan == Muon) {
    plep1.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(1), MuPy.at(1), MuPz.at(1), MuEnergy.at(1));
  }
  if (chan == Elec) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(ElPx.at(1), ElPy.at(1), ElPz.at(1), ElEnergy.at(1));
  }
  if (chan == ElMu) {
    plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
    plep2.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
  }
  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
  TLorentzVector pjet0 = Jet[0].p; TLorentzVector pjet1 = Jet[1].p;
  float MT2llbb00, MT2llbb01, MT2llbb;
  float METx = pmet.Px(); float METy = pmet.Py();
  float MET = pmet.Pt(); float MET_phi = pmet.Phi();
  TLorentzVector LepPlusBtagJet00 = plep1 + pjet0;
  TLorentzVector LepPlusBtagJet10 = plep2 + pjet0;
  TLorentzVector LepPlusBtagJet11 = plep2 + pjet1;
  TLorentzVector LepPlusBtagJet01 = plep1 + pjet1;
  if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && (LepPlusBtagJet10.M()>173 || LepPlusBtagJet01.M()>173))
    MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, pmet, 0.);
  else if ((LepPlusBtagJet11.M()>173 || LepPlusBtagJet00.M()>173) && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173)
    MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, pmet, 0.);
  else if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173) {
    if ( fabs(LepPlusBtagJet11.M()-LepPlusBtagJet00.M()) < fabs(LepPlusBtagJet10.M()-LepPlusBtagJet01.M()) )
      MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, pmet, 0.);
    else
      MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, pmet, 0.);
  }
  else
    MT2llbb=0;
  return MT2llbb;
}

float TTHAnalyzer::getMeff(){
	if(Jet.size()<2) return -999.;
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return fHypLepton1.p.Pt() + fHypLepton2.p.Pt() + Jet[0].p.Pt() + Jet[1].p.Pt() + pmet.Pt();
}

TLorentzVector TTHAnalyzer::getPtllb(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector pl1 = fHypLepton1.p; TLorentzVector pl2 = fHypLepton2.p;
	return pl1 + pl2 + pmet;
}

float TTHAnalyzer::getDPhibMet(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector Ptllb = getPtllb();
	return pmet.DeltaPhi(Ptllb);
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TTHAnalyzer::getJERScale(int jet){
	float eta = Jet_eta[ jet];
	// 8 TeV
	if(     TMath::Abs(eta) < 0.5) return 1.079;
	else if(TMath::Abs(eta) < 1.1) return 1.099;
	else if(TMath::Abs(eta) < 1.7) return 1.121;
	else if(TMath::Abs(eta) < 2.3) return 1.208;
	else if(TMath::Abs(eta) < 2.8) return 1.254;
	else if(TMath::Abs(eta) < 3.2) return 1.395;
	else                           return 1.056;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TTHAnalyzer::getJERScaleUp(int jet){
	float eta = Jet_eta[ jet];
	// up, 8 TeV
	if(	  TMath::Abs(eta) < 0.5) return 1.053;
	else if(TMath::Abs(eta) < 1.1) return 1.071;
	else if(TMath::Abs(eta) < 1.7) return 1.092;
	else if(TMath::Abs(eta) < 2.3) return 1.162;
	else if(TMath::Abs(eta) < 2.8) return 1.192;
	else if(TMath::Abs(eta) < 3.2) return 1.332;
	else  			 return 0.865;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TTHAnalyzer::getJERScaleDown(int jet){
	float eta = Jet_eta[ jet];
	// down, 8 TeV
	if(	  TMath::Abs(eta) < 0.5) return 1.105;
	else if(TMath::Abs(eta) < 1.1) return 1.127;
	else if(TMath::Abs(eta) < 1.7) return 1.150;
	else if(TMath::Abs(eta) < 2.3) return 1.254;
	else if(TMath::Abs(eta) < 2.8) return 1.316;
	else if(TMath::Abs(eta) < 3.2) return 1.458;
	else  			 return 1.247;
}

float TTHAnalyzer::getErrPt(float Pt, float Eta) {
	float InvPerr2;
	float N(0.), S(0.), C(0.), m(0.);

	if(TMath::Abs(Eta) < 0.5 ) {
		N = 3.96859;
		S = 0.18348;
		C = 0.;
		m = 0.62627;
	} else if( TMath::Abs(Eta) < 1.  ) {
		N = 3.55226;
		S = 0.24026;
		C = 0.;
		m = 0.52571;
	} else if( TMath::Abs(Eta) < 1.5  ) {
		N = 4.54826;
		S = 0.22652;
		C = 0.;
		m = 0.58963;
	} else if( TMath::Abs(Eta) < 2.  ) {
		N = 4.62622;
		S = 0.23664;
		C = 0.;
		m = 0.48738;
	} else if( TMath::Abs(Eta) < 2.5  ) {
		N = 2.53324;
		S = 0.34306;
		C = 0.;
		m = 0.28662;
	} else if( TMath::Abs(Eta) < 3.  ) {
		N = -3.33814;
		S = 0.73360;
		C = 0.;
		m = 0.08264;
	} else if( TMath::Abs(Eta) < 5.  ) {
		N = 2.95397;
		S = 0.11619;
		C = 0.;
		m = 0.96086;
	}
	// this is the absolute resolution (squared), not sigma(pt)/pt
	// so have to multiply by pt^2, thats why m+1 instead of m-1
	InvPerr2 =  (N * TMath::Abs(N) ) + (S * S) * pow(Pt, m+1) + (C * C) * Pt * Pt ;
	return sqrt(InvPerr2);
}

float TTHAnalyzer::getLeptonError(gChannel chan){
	float err1(0.), err2(0.);
	if (chan==Muon){
		err1 = fLeptonSF->GetTightMuonSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightMuonSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	if (chan==ElMu){
		err1 = fLeptonSF->GetTightMuonSF_err    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	if (chan==Elec){
		err1 = fLeptonSF->GetTightElectronSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	return TMath::Sqrt(err1*err1+err2*err2);
}

float TTHAnalyzer::getTriggerError(gChannel chan){
	float trig(0.);
	if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	if (chan==ElMu) trig = fLeptonSF->GetMuEGSF_err    (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
	if (chan==Elec) trig = fLeptonSF->GetDoubleElSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	return trig;
}

float TTHAnalyzer::getSF(gChannel chan) {
	if (gIsData)              return 1.; //Don't scale data
	float id1(1.),id2(1.), trig(1.);
	float err1(0.), err2(0.), err_trg(0.);
  float SF = 0; float FSSF = 0;
  float FS1 = 0; float FS2 = 0;
	if (chan == Muon){
		id1  = fLeptonSF->GetTightMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		id2  = fLeptonSF->GetTightMuonSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
		trig = fLeptonSF->GetDoubleMuSF (fHypLepton1.p.Pt(),fHypLepton1.p.Eta());
    FS1  = fLeptonSF->GetFastSimMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    FS2  = fLeptonSF->GetFastSimMuonSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  }
  else if (chan == Elec){
    id1  = fLeptonSF->GetTightElectronSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetDoubleElSF     (fHypLepton1.p.Pt(),fHypLepton1.p.Eta());
    FS1  = fLeptonSF->GetFastSimElectronSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    FS2  = fLeptonSF->GetFastSimElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  }
  else if (chan == ElMu){
    float leadingPt  = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Pt() : fHypLepton2.p.Pt();
    float leadingEta = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Eta() : fHypLepton2.p.Eta();
    id1  = fLeptonSF->GetTightMuonSF    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetMuEGSF         (leadingPt,leadingEta);
    FS1  = fLeptonSF->GetFastSimMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
    FS2  = fLeptonSF->GetFastSimElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
  SF = PUSF*id1*id2*trig;
  FSSF = FS1*FS2;
	if(gSampleName.BeginsWith("T2tt")) SF*= FSSF;
  return (SF);
}

float TTHAnalyzer::getTopPtSF(){
	// Return SF of the pt pt of the top
	// Only apply SF if the process is ttbar...
	if(!gSampleName.Contains("TTJets") && !gSampleName.Contains("TTbar_Powheg")) return 1.;
	if (gSysSource==TopPt) {
		TLorentzVector top;
		Float_t topSF = 0.;
		Float_t Weight = 1.;
		if (!gIsData) {
			Int_t nGenTop = Get<Int_t>("nGenTop");
			if (nGenTop != 2) return 1.;
			for (Int_t t=0; t<nGenTop; t++){
				top.SetPtEtaPhiM(Get<Float_t>("GenTop_pt",   t),
						Get<Float_t>("GenTop_eta",  t),
						Get<Float_t>("GenTop_phi",  t),
						Get<Float_t>("GenTop_mass", t));
				Float_t pt = TMath::Min(top.Pt(), 400.);
				topSF = TMath::Exp(0.156 - 0.00137 * pt);
				Weight *= topSF;
			}
			Weight = TMath::Sqrt(Weight);
		}
		return Weight;
	}
	return 1.;
}

//--------------------------------------------------------------------------
// Fill histograms
//------------------------------------------------------------------------
void TTHAnalyzer::FillDYHistograms(){
	float Mll = 0.;
	if (triggermumuSS()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu);
		if(gIsMCatNLO)  EventWeight = EventWeight * genWeight;
		Mll = (fHypLepton1.p+fHypLepton2.p).M();
	}
	ResetHypLeptons();
	if (triggermumuSS() && IsMuMuEvent()){
		EventWeight = gWeight * getSF(Muon);
		if(gIsMCatNLO)   EventWeight = EventWeight * genWeight;
		Mll = (fHypLepton1.p+fHypLepton2.p).M();

	}
	ResetHypLeptons();
	if (triggermumuSS()   && IsElElEvent()){
		EventWeight = gWeight * getSF(Elec);
		if(gIsMCatNLO)    EventWeight = EventWeight * genWeight;
		Mll = (fHypLepton1.p+fHypLepton2.p).M();

	}
	ResetHypLeptons();
}
void TTHAnalyzer::FillKinematicHistos(gChannel chan, iCut cut){
#ifdef DEBUG
	cout << "Filling KinematicHistos("<<chan<<","<<cut<<")... ";
	cout << fHypLepton1.index << " , " << fHypLepton2.index << endl;
#endif

	if (gSysSource != Norm)      return;  //only fill histograms for nominal distributions...
	if (fChargeSwitch == true  ) return;

	//++ jet info
	int njets = getNJets();
#ifdef DEBUG
	cout << " DONE!" << endl;
#endif

}

void TTHAnalyzer::FillYieldsHistograms(gChannel chan, iCut cut, gSystFlag sys){
#ifdef DEBUG
	cout << "FillYieldsHistograms("<<chan<<","<<cut<<","<<sys<<")...";
#endif
	if (fChargeSwitch){   fHSSyields[chan][sys]->Fill(cut, EventWeight);  }
	else {                fHyields[chan][sys]  ->Fill(cut, EventWeight);
	}

	/// FOR SYSTEMATIC STUDIES
	int njets  = 0; njets  = getNJets();
	int nbtags = 0; nbtags = getNBTags();


#ifdef DEBUG
	cout << " DONE! " << endl;
#endif
	return;
}
void TTHAnalyzer::FillYields(gSystFlag sys){
	ResetHypLeptons();

#ifdef DEBUG
	cout << "PassTriggerEMu= " << triggermumuSS() << endl;
	cout << "Is ElMu/ElEl/MuMu Event= "
		<< IsElMuEvent() << "/"
		<< IsElElEvent() << "/"
		<< IsMuMuEvent() << endl;
#endif

  CoutEvent(evt, Form(" PassTrigEmu: %i", triggermumuSS()));
	if (triggermumuSS()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu);// * getTopPtSF();
		hWeight -> Fill(EventWeight,1.);
#ifdef DEBUG
		cout << " pass trigger + emu, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115));

		if(
				(gCreateTree) && (sys==Norm)        && !(fChargeSwitch) &&
				PassesZVeto()      &&
				PassesMllVeto()    &&
				PassesNJetsCut()){
			SetTreeVariables(ElMu);
			fTree->Fill();
		}

#ifdef DEBUG
			cout << " pass mll, ";
#endif

		if (PassesMllVeto()){
			FillYieldsHistograms(ElMu, iDilepton, sys);
			if(sys==Norm) FillKinematicHistos(ElMu,iDilepton);

			FillYieldsHistograms(ElMu, iZVeto, sys);
			if(sys==Norm) FillKinematicHistos(ElMu, iZVeto);

			if(PassesMETCut()){
				FillYieldsHistograms(ElMu, iMET, sys);
				if(sys==Norm) FillKinematicHistos(ElMu,iMET);

				if (PassesNJetsCut()) {
#ifdef DEBUG
					cout << " pass njets with njets = "<<getNJets()<<", ";
#endif
					FillYieldsHistograms(ElMu, i2jets, sys);
					if(sys==Norm) FillKinematicHistos(ElMu,i2jets);
					if (PassesNBtagCut()) {
#ifdef DEBUG
						cout << " pass nbjets with nbtags = "<<getNBTags()<<", ";
#endif
						//if (sys == LESUp) cout << evt<< endl;  //LESup    8 //EventNumber
						FillYieldsHistograms(ElMu, i1btag, sys);
						if(sys==Norm) FillKinematicHistos(ElMu,i1btag);
						if(1){ // No DY Veto in emu
							FillYieldsHistograms(ElMu,iDYVeto, sys);
							if(sys==Norm) FillKinematicHistos(ElMu,iDYVeto);
						}
					}
				}
			}
			if (getNBTags() == 1){
#ifdef DEBUG
				cout << " pass nbjets=1";
#endif
				FillYieldsHistograms(ElMu, iExact1btag, sys);
				if(sys==Norm) FillKinematicHistos(ElMu,iExact1btag);
			}
			if (getNBTags() == 2){
#ifdef DEBUG
				cout << " pass nbjets=2";
#endif
				FillYieldsHistograms(ElMu, iExact2btag, sys);
				if(sys==Norm) FillKinematicHistos(ElMu,iExact2btag);
			}
		}
	}

	ResetHypLeptons();
	if (triggermumuSS() && IsMuMuEvent()){
		EventWeight = gWeight * getSF(Muon); //  * getTopPtSF();
		//EventWeight = 1.;
#ifdef DEBUG
		cout << " pass trigger + mumu, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115));

		if(
			(gCreateTree) &&	(sys==Norm) && !(fChargeSwitch) &&
				PassesZVeto()      &&
				PassesMllVeto()    &&
				PassesNJetsCut()){
			SetTreeVariables(Muon);
			fTree->Fill();
		}

		if (PassesMllVeto()){
#ifdef DEBUG
			cout << " pass mll, ";
#endif
			FillYieldsHistograms(Muon,iDilepton, sys);
			if(sys==Norm) FillKinematicHistos(Muon,iDilepton);
			if (PassesZVeto())    {
				FillYieldsHistograms(Muon,iZVeto, sys);
				if(sys==Norm) FillKinematicHistos(Muon,iZVeto);
				if (PassesMETCut())   {
					FillYieldsHistograms(Muon,iMET, sys);
					if(sys==Norm) FillKinematicHistos(Muon,iMET);

					if (getNBTags() == 1){
						FillYieldsHistograms(Muon, iExact1btag, sys);
						if(sys==Norm) FillKinematicHistos(Muon,iExact1btag);
					}
					if (getNBTags() == 2){
						FillYieldsHistograms(Muon, iExact2btag, sys);
						if(sys==Norm) FillKinematicHistos(Muon,iExact2btag);
					}
					if (PassesNJetsCut()) {
						FillYieldsHistograms(Muon,i2jets, sys);
						if(sys==Norm) FillKinematicHistos(Muon,i2jets);
						if (PassesNBtagCut()) {
							FillYieldsHistograms(Muon,i1btag, sys);
							if(sys==Norm) FillKinematicHistos(Muon,i1btag);
              if(PassesDYVetoCut()){
                FillYieldsHistograms(Muon,iDYVeto, sys);
                if(sys==Norm) FillKinematicHistos(Muon,iDYVeto);
              }
						}
					}
				}
            }
        }
    }

  ResetHypLeptons();
  if (triggermumuSS() && IsElElEvent()){
		EventWeight = gWeight * getSF(Elec);// * getTopPtSF();
		//EventWeight = 1.;
#ifdef DEBUG
		cout << " pass trigger + ee, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115));

		if(
			(gCreateTree) &&	(sys==Norm)        && !(fChargeSwitch) &&
				PassesZVeto()      &&
				PassesMllVeto()    &&
				PassesNJetsCut()){
			SetTreeVariables(Elec);
			fTree->Fill();
		}


}
  ResetHypLeptons();
#ifdef DEBUG
    cout << " DONE!"<<endl;
#endif
}

//----------------------------------------------------------------------
// Passes
//----------------------------------------------------------------------
bool TTHAnalyzer::PassesMuonEta2p1(gChannel chan){
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;

	if (chan == Muon){
		if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true;
		if (TMath::Abs(fHypLepton2.p.Eta()) < 2.1) return true;
	}
	else if (chan == ElMu){
		if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true;
	}
	else if (chan == Elec){
		return true;
	}
	return false;
}

bool TTHAnalyzer::Passes3rdLeptonVeto(){
	return true; // don't apply third lepton veto...
	// Return false if there are not 2 signal leptons
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;

	//  Int_t nvetoleptons = 0;
	for(Int_t i = 0; i < nLepGood; ++i){ //elec
		if (fHypLepton1.index == i && fHypLepton1.type == 0) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 0) continue;
		//    if (IsVetoMuon(i)) return false;
		//     nvetoleptons++;
	}

	for(Int_t i = 0; i < nLepGood; ++i){ // muon
		if (fHypLepton1.index == i && fHypLepton1.type == 1) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 1) continue;
		//     if (IsVetoElectron(i)) return false;
		//    nvetoleptons++;
	}
	//  if (nvetoleptons > 0) return false;
	return true;
}

bool TTHAnalyzer::PassesMllVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
	if (InvMass < 20.)            return false;
	return true;
}

bool TTHAnalyzer::PassesZVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
    if (InvMass > 76. && InvMass < 106.) return false;
    return true;
}

bool TTHAnalyzer::PassesTopDCut(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	if (getTopD() < 0.8) return false;
	return true;
}

bool TTHAnalyzer::PassesNJetsCut(){
	if (getNJets() <= 1) return false;
	return true;
}

bool TTHAnalyzer::PassesMETCut(){
	if (getMET() < 50.) return false;
	return true;
}

bool TTHAnalyzer::PassesNBtagCut(){
	if (getNBTags() < 1) return false;
	return true;
}

bool TTHAnalyzer::PassesDYVetoCut(){
    if(getMinDPhiMetJets() < 0.25) return false;
    if (getHT() == 0) return false;
    if(getMET()/TMath::Sqrt(getHT()) < 5.0) return false;
    return true;
}

bool TTHAnalyzer::IsElMuEvent(){
	if (fChargeSwitch){      return (IsDileptonEvent()  == 3);   }
	return (IsDileptonEvent() == -3);
}

bool TTHAnalyzer::IsMuMuEvent(){
	if (fChargeSwitch){  return (IsDileptonEvent()  == 1); }
	return (IsDileptonEvent() == -1);
}

bool TTHAnalyzer::IsElElEvent(){
    if (fChargeSwitch){    return (IsDileptonEvent()  == 2); }
    return (IsDileptonEvent() == -2);
}

int TTHAnalyzer::IsDileptonEvent(){
#ifdef DEBUG
	cout << "IsDileptonEvent(): NLeptons =" << Lepton.size()<< endl;
#endif
	if(Lepton.size() != 2) return 0; // xxx
    if(Lepton[0].p.Pt()<25) return 0;
	int select = Lepton[0].charge*Lepton[1].charge;
	int result = 0;
	if      (Lepton[0].type == 0 && Lepton[1].type == 0) result = 1; // mu/mu
	else if (Lepton[0].type == 1 && Lepton[1].type == 1) result = 2; // el/el
    else result = 3; // mu/el

	fHypLepton1 = lepton(Lepton[0]);
	fHypLepton2 = lepton(Lepton[1]);

	if(Lepton[0].type == 1 && Lepton[1].type == 0){
		fHypLepton1 = lepton(Lepton[1]);
		fHypLepton2 = lepton(Lepton[0]);
	}
	result *= select; // Add charge to result
#ifdef DEBUG
    cout << result;
	cout << " DONE!" << endl;
#endif
	return result;
}

//------------------------------------------------------------------------------
// LEPTON SELECTORS
//------------------------------------------------------------------------------
bool momentumComparator(lepton i, lepton j){ return (i.p.Pt()>j.p.Pt()); }

vector<lepton> TTHAnalyzer::SortLeptonsByPt(vector<lepton>& leptons){
    vector<lepton> theLep = leptons;
    sort (theLep.begin(), theLep.end(), momentumComparator);
    return theLep;
}

int TTHAnalyzer::getSelectedLeptons(){
    // Loops over the total number of muons and electrons and returns the number of leptons.
    if (Lepton.size() > 0) {
        cout << "[WARNING]: you have called this function previously... RESETTING..."<<endl;
        Lepton.clear();
    }
    vector<lepton> tmp_lepton;
    nMuon = 0;
    TLorentzVector lep;
    Int_t thetype = 0;
	for (Int_t i=0; i<nLepGood;i++){
		CoutEvent(evt, Form("---- Lepton %i", i));
		CoutEvent(evt, Form("   pdgId = %i ", LepGood_pdgId[i]));
		CoutEvent(evt, Form("   pT    = %f ", LepGood_pt[i]));
		CoutEvent(evt, Form("   Eta   = %f ", LepGood_eta[i]));
		CoutEvent(evt, Form("   Sip3D = %f ", Get<Float_t>("LepGood_sip3d", i)));
		CoutEvent(evt, Form("     dxy = %f ", LepGood_dxy[i]));
		CoutEvent(evt, Form("      dz = %f ", LepGood_dz[i]));
		CoutEvent(evt, Form("   LHits = %i ", Get<Int_t>("LepGood_lostHits", i)));
		CoutEvent(evt, Form("   mediumMuonId ---> %i", Get<Int_t>("LepGood_mediumMuonId",i)));
		CoutEvent(evt, Form("   tightElecId  ---> %i", Get<Int_t>("LepGood_tightId", i)));
		CoutEvent(evt, Form("   MULTIISO VT"));
		CoutEvent(evt, Form("   miniRelIso (<0.09) : %f", Get<Float_t>("LepGood_miniRelIso", i)));
		CoutEvent(evt, Form("      ptRatio (>0.84) : %f", Get<Float_t>("LepGood_jetPtRatiov2", i)));
		CoutEvent(evt, Form("        ptRel (>7.2 ) : %f", Get<Float_t>("LepGood_jetPtRelv2", i)));
        if(Get<Int_t>("LepGood_tightId", i) <3 && fabs(LepGood_pdgId[i]) == 11){
            CoutEvent(evt, Form("etaSC = %f", Get<Float_t>("LepGood_etaSc", i)));
            CoutEvent(evt, Form("sigmaIEtaIEta = %f", Get<Float_t>("LepGood_sigmaIEtaIEta", i)));
            CoutEvent(evt, Form("dEtaScTrkIn = %f", Get<Float_t>("LepGood_dEtaScTrkIn", i)));
            CoutEvent(evt, Form("dPhiScTrkIn = %f", Get<Float_t>("LepGood_dPhiScTrkIn", i)));
            CoutEvent(evt, Form("hadronicOverEm = %f", Get<Float_t>("LepGood_hadronicOverEm", i)));
            CoutEvent(evt, Form("eInvMinusPInv = %f", Get<Float_t>("LepGood_eInvMinusPInv", i)));
            CoutEvent(evt, Form("eInvMinusPInv_tkMom = %f", Get<Float_t>("LepGood_eInvMinusPInv_tkMom", i)));
	    }
		if(IsTightMuon(i)){
			CoutEvent(evt, "   Is Muon!!");
			thetype = 0;
			nMuon ++;
        }
        else if(IsTightElectron(i)){
			CoutEvent(evt, "   Is Electron!!");
            thetype = 1;
            nElec++;
        }
        else  continue;

        lep.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]);
        lepton tmpLepton(lep, LepGood_charge[i], thetype, i);
        tmp_lepton.push_back(tmpLepton);
    }
    Lepton = SortLeptonsByPt(tmp_lepton);
    CoutEvent(evt, Form("  ---> nselLeps = %i", Lepton.size()));
    if(Lepton.size() == 2) CoutEvent(evt, Form("  --->      Mll = %f", (Lepton[0].p+Lepton[1].p).M()));
    return Lepton.size();
}

bool TTHAnalyzer::METFilter(){
    if(gSampleName.BeginsWith("T2tt")) return true;
    if (Get<Int_t>("Flag_HBHENoiseFilter") &&
        Get<Int_t>("Flag_HBHENoiseIsoFilter") &&
        Get<Int_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
        Get<Int_t>("Flag_goodVertices") &&
        Get<Int_t>("Flag_eeBadScFilter")){

        //Get<Int_t>("Flag_BadMuon"))
        //Get<Int_t>("Flag_badChargedHadron"))
        //Get<Int_t>("Flag_globalTightHalo2016Filter"))
        CoutEvent(evt, "   Passes MET filters!!");
        return true;
    }
    else{
        CoutEvent(evt, "----->Does not pass MET filters!!");
        return false;
    }
}

//------------------------------------------------------------------------------
// Muon selectors
//------------------------------------------------------------------------------
bool TTHAnalyzer::IsTightMuon(unsigned int iMuon,float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return false;
	if (TMath::Abs(LepGood_eta[iMuon]) > 2.4) return false;
	if (LepGood_pt[iMuon] < 10) return false;
	if (TMath::Abs(LepGood_dxy[iMuon]) > 0.05) return false;
	if (TMath::Abs(LepGood_z[iMuon]) > 0.1) return false;
	if (LepGood_sip3d[iMuon] > 8) return false;
	if (LepGood_miniRelIso[iMuon] > 0.4) return false;
	// The condition of "loose muon" is not demanded as the trees which we are
	// doing the analysis already are filtered by it.
	if (LepGood_jetBTagCSV[iMuon] > 0.89) return false;
	if (LepGood_mediumMuonId[iMuon] != 1) return false;
	//if (LepGood_tightCharge[iMuon] != 2) return false; // WOLOLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
	if (LepGood_mvaTTH[iMuon] < 0.75) return false;
    return true;
}

bool TTHAnalyzer::IsFakeableMuon(unsigned int iMuon,float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return false;
	if (TMath::Abs(LepGood_eta[iMuon]) > 2.4) return false;
	if (LepGood_pt[iMuon] < 10) return false;
	if (TMath::Abs(LepGood_dxy[iMuon]) > 0.05) return false;
	if (TMath::Abs(LepGood_z[iMuon]) > 0.1) return false;
	if (LepGood_sip3d[iMuon] > 8) return false;
	if (LepGood_miniRelIso[iMuon] > 0.4) return false;
	// The condition of "loose muon" is not demanded as the trees which we are
	// doing the analysis already are filtered by it.
	if (lepGood_mvaTTH[iMuon] < 0.75) {
		if (LepGood_jetPtRatiov2[iMuon] < 0.3) return false;
		if (LepGood_jetBTagCSV[iMuon] > 0.605) return false;
	}
	else {
		if (LepGood_jetBTagCSV[iMuon] > 0.89) return false;
	}
    return true;
}

bool TTHAnalyzer::IsLooseMuon(unsigned int iMuon,float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return false;
	if (TMath::Abs(LepGood_eta[iMuon]) > 2.4) return false;
	if (LepGood_pt[iMuon] < 5) return false;
	if (TMath::Abs(LepGood_dxy[iMuon]) > 0.05) return false;
	if (TMath::Abs(LepGood_z[iMuon]) > 0.1) return false;
	if (LepGood_sip3d[iMuon] > 8) return false;
	if (LepGood_miniRelIso[iMuon] > 0.4) return false;
	// The condition of "loose muon" is not demanded as the trees which we are
	// doing the analysis already are filtered by it.
    return true;
}

float TTHAnalyzer::getMuonIso(int iMuon){
	if (iMuon < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return 9999.;
	return LepGood_relIso04[iMuon];
}

//------------------------------------------------------------------------------
// Electron selectors
//------------------------------------------------------------------------------
bool TTHAnalyzer::IsTightElectron(unsigned int iElec, float ptcut, int an){
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return false;
	if (TMath::Abs(LepGood_eta[iElec]) > 2.5) return false;
	if (an == 2) {
		if (LepGood_pt[iElec] < 15) return false;
	} 
	else {
		if (LepGood_pt[iElec] < 10) return false;
	}
	if (TMath::Abs(LepGood_dxy[iElec]) > 0.05) return false;
	if (TMath::Abs(LepGood_z[iElec]) > 0.1) return false;
	if (LepGood_sip3d[iElec] > 8) return false;
	if (LepGood_miniRelIso[iElec] > 0.4) return false;
	
	if (TMath::Abs(LepGood_eta[iElec]) < 0.8) {
		if (LepGood_mvaTTH[iElec] < -0.70) return false;
		if (LepGood_pt[iElec] > 30) {
			if (LepGood_sigmaIEtaIEta[iElec] > 0.011) return false;
			
		}
	} 
	else if ((TMath::Abs(LepGood_eta[iElec]) < 1.479) && (TMath::Abs(LepGood_eta[iElec]) >= 0.8)){
		if (LepGood_mvaTTH[iElec] < -0.83) return false;
		if (LepGood_pt[iElec] > 30) {
			if (LepGood_sigmaIEtaIEta[iElec] > 0.011) return false;
			
		}
	} 
	else if (TMath::Abs(LepGood_eta[iElec]) >= 1.479) {
		if (LepGood_mvaTTH[iElec] < -0.92) return false;
		if (LepGood_pt[iElec] > 30) {
			if (LepGood_sigmaIEtaIEta[iElec] > 0.030) return false;
			
		}
	}
	
		
		
		
		
		
		
	return true;
}

/*
bool TTHAnalyzer::IsFakeableElectron(unsigned int iElec, float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return false;
	if (LepGood_pt[iElec] < ptcut) return false;
	if (TMath::Abs(LepGood_eta[iElec]) > 2.4) return false;
    if (!getMultiIso(iElec)) return false;
    //if (getElecIso(iElec) > 0.0766) return false;
	if (TMath::Abs(LepGood_dxy[iElec]) >= 0.05) return false;
	if (TMath::Abs(LepGood_dz[ iElec]) >= 0.1 ) return false;
	if (Get<Int_t>("LepGood_lostHits", iElec)         > 0      ) return false;
	if (TMath::Abs(LepGood_etaSc[iElec]) > 1.4442 &&
			TMath::Abs(LepGood_etaSc[iElec]) < 1.566) return false;
	if (Get<Float_t>("LepGood_sip3d", iElec) > 4.0) return false;
	//if(TMath::Abs(Get<Int_t>("LepGood_tightId", iElec)) != 3) return false; // bin 3: tight ID electrons
    //'POG_SPRING15_25ns_v1_Tight' : [('dEtaIn', [0.00926, 0.00724]), ('dPhiIn', [0.0336, 0.0918]), ('sigmaIEtaIEta', [0.0101, 0.0279]), ('H/E', [0.0597, 0.0615]), ('1/E-1/p', [0.0120, 0.00999])],
    // Tight Electron Id:
    if(TMath::Abs(Get<Float_t>("LepGood_etaSc", iElec)) < 1.479){ // central
        if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.0101) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00926 ) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0336) return false;
        if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0597) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) > 0.012) return false;
	}
    else{ // forward
        if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.0279) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00724) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0918) return false;
        if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0615) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) >  	0.00999) return false;
    }
	return true;
}
bool TTHAnalyzer::IsLooseElectron(unsigned int iElec, float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return false;
	if (LepGood_pt[iElec] < ptcut) return false;
	if (TMath::Abs(LepGood_eta[iElec]) > 2.4) return false;
    if (!getMultiIso(iElec)) return false;
    //if (getElecIso(iElec) > 0.0766) return false;
	if (TMath::Abs(LepGood_dxy[iElec]) >= 0.05) return false;
	if (TMath::Abs(LepGood_dz[ iElec]) >= 0.1 ) return false;
	if (Get<Int_t>("LepGood_lostHits", iElec)         > 0      ) return false;
	if (TMath::Abs(LepGood_etaSc[iElec]) > 1.4442 &&
			TMath::Abs(LepGood_etaSc[iElec]) < 1.566) return false;
	if (Get<Float_t>("LepGood_sip3d", iElec) > 4.0) return false;
	//if(TMath::Abs(Get<Int_t>("LepGood_tightId", iElec)) != 3) return false; // bin 3: tight ID electrons
    //'POG_SPRING15_25ns_v1_Tight' : [('dEtaIn', [0.00926, 0.00724]), ('dPhiIn', [0.0336, 0.0918]), ('sigmaIEtaIEta', [0.0101, 0.0279]), ('H/E', [0.0597, 0.0615]), ('1/E-1/p', [0.0120, 0.00999])],
    // Tight Electron Id:
    if(TMath::Abs(Get<Float_t>("LepGood_etaSc", iElec)) < 1.479){ // central
        if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.0101) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00926 ) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0336) return false;
        if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0597) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) > 0.012) return false;
	}
    else{ // forward
        if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.0279) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00724) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0918) return false;
        if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0615) return false;
        if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) >  	0.00999) return false;
    }
	return true;
}
*/


float TTHAnalyzer::getElecIso(int iElec){
	if (iElec < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return 9999.;
	return LepGood_relIso03[iElec];
}

bool TTHAnalyzer::getMultiIso(unsigned int index){
    float max_mRelIso = 0.09;
    float min_ptRatio = 0.84; // Very tight
    float min_ptRel   = 7.2 ; // Tight
    //if      ((TMath::Abs(LepGood_pdgId[index])) == 11) max_mRelIso = 0.13;  // Tight for electrons
    //else if ((TMath::Abs(LepGood_pdgId[index])) == 13) max_mRelIso = 0.2;   // Medium for muons
    float mRelIso = Get<Float_t>("LepGood_miniRelIso", index);
    float ptRatio = Get<Float_t>("LepGood_jetPtRatiov2", index);
    float ptRel   = Get<Float_t>("LepGood_jetPtRelv2", index);
    // min_ptRatio = 0; min_ptRel = 0; // No MultiIso
    return (mRelIso < max_mRelIso && (ptRatio > min_ptRatio || ptRel > min_ptRel));
}

float TTHAnalyzer::getEACorrection(float eta){  // for a 0.3 CONE
	float abseta = TMath::Abs(eta);
	// numbers from https://indico.cern.ch/event/370494/contribution/2/material/slides/0.pdf
	if      (abseta < 0.8)                  return 0.1013;
	else if (abseta >= 0.8 && abseta < 1.3) return 0.0988;
	else if (abseta >= 1.3 && abseta < 2.0) return 0.0572;
	else if (abseta >= 2.0 && abseta < 2.2) return 0.0842;
	else if (abseta >= 2.2 && abseta < 5.0) return 0.1530;
	else                                    return 0.1530;

	cout << "[ERROR] getEACorrection(): No correction factor found!!" << endl;
	return -999999.;
}

void TTHAnalyzer::setMET(float newmet){ MET = newmet;}

float TTHAnalyzer::getMET(){ return MET; }

float TTHAnalyzer::getMETPhi(){ return MET_Phi;}

int TTHAnalyzer::getNJets(){ return nJets;}

float TTHAnalyzer::getDRClosestJet(TLorentzVector lep){
	float minDR = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDR > lep.DeltaR(Jet[i].p)) minDR = lep.DeltaR(Jet[i].p);
	}
	return minDR;
}

float TTHAnalyzer::getDPhiClosestJet(TLorentzVector lep){
	float minDphi = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDphi > TMath::Abs(lep.DeltaPhi(Jet[i].p))) minDphi = TMath::Abs(lep.DeltaPhi(Jet[i].p));
	}
	return minDphi;
}

int TTHAnalyzer::getLeadingJetbTag(){
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (Jet[i].isbtag) return i;
	}
	return  -1;
}

int TTHAnalyzer::getNBTags(){
	int ntags(0);
	for(UInt_t i = 0; i <Jet.size(); i++){
		if (Jet[i].isbtag) ntags++;
	}
	return ntags;
}

float TTHAnalyzer::getDeltaPhillJet(){
	if (fHypLepton1.index == -1) return -999.;
	if (fHypLepton2.index == -1) return -999.;
	Int_t ij = getLeadingJetbTag();
	if (ij < 0) return -999.;
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	TLorentzVector jet = Jet[ij].p;
	return TMath::Abs(dilep.DeltaPhi(jet));
}

float TTHAnalyzer::getTopD(){
	if (fHypLepton1.index == -1) return -999;
	if (fHypLepton2.index == -1) return -999;
	// Make Dilepton system
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	Float_t DeltaPhi(0.),TopD(0.);
	TLorentzVector jet;
	if (nJets == 0) return -999.;
	TopD     = 1 - (DeltaPhi/TMath::Pi()) * (1 - Jet_btagCSV[(Jet[0].index)]);
	return TopD;
}

int TTHAnalyzer::getSelectedJets(){
	int nj(0);
	if (Jet.size() > 0) {
		cout << "[WARNING]: you have called this function previously, RESETTING..."<<endl;
		Jet.clear();
	}
	//int btagSys = 0;
	//TLorentzVector jtDisc;
	//for (UInt_t i=0; i<nDiscJet; i++) {
	//  jtDisc.SetPtEtaPhiE(DiscJet_pt[i], DiscJet_eta[i], DiscJet_phi[i], DiscJet_energy[i]);
	//}
    TLorentzVector jt;
    for (Int_t i=0; i<nJet; i++) {
        if(!IsGoodJet(i,gJetEtCut)) continue;

        Float_t jetbtagi      = Jet_btagCSV[i];
        Float_t jetetai       = Jet_eta[i];
        Float_t jetenergyi    = Jet_energy[i];

        jt.SetPtEtaPhiE(JetPt.at(i), jetetai, JetPhi.at(i), jetenergyi);
        bool isbtag = false;
        if (gIsData) {
            isbtag = fBTagSFnom->IsTagged(Jet_btagCSV[i], -999999, JetPt.at(i), jetetai);
        }
        else {
			Int_t   jetmcflavouri = Get<Int_t>  ("Jet_mcFlavour", i);
			// official b-tag recommendation: use JetHadronFlavour instead of JetPartonFlavor
			/*
				 if(TMath::Abs(Jet_mcFlavour[i]) == 5 || TMath::Abs(Jet_mcFlavour[i]) == 4){
				 if (gSysSource == BtagUp)     btagSys =  1;
				 if (gSysSource == BtagDown)   btagSys = -1;
				 if (gSysSource == MisTagUp)   btagSys =  0;
				 if (gSysSource == MisTagDown) btagSys =  0;
				 }
				 else {
				 if (gSysSource == BtagUp)     btagSys =  0;
				 if (gSysSource == BtagDown)   btagSys =  0;
				 if (gSysSource == MisTagUp)   btagSys =  1;
				 if (gSysSource == MisTagDown) btagSys = -1;
				 }*/
            if      (gSysSource == BtagUp)      isbtag = fBTagSFbUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
            else if (gSysSource == BtagDown)    isbtag = fBTagSFbDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
            else if (gSysSource == MisTagUp)    isbtag = fBTagSFlUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
            else if (gSysSource == MisTagDown)  isbtag = fBTagSFlDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
            else                                isbtag = fBTagSFnom->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
            // Use this line only to get raw numbers for syncronization
            //if(Get<Float_t>("Jet_btagCSV",[i] > 0.89) isbtag=true;  // WP for 74
        }
        jet tmpjet(jt, isbtag, i);
        Jet.push_back(tmpjet);
        nj++;
    }
    return nj;
}

bool TTHAnalyzer::IsGoodJet(unsigned int ijet, float ptcut){
    float minDR = 0.4;
    TLorentzVector jet;
    jet.SetPtEtaPhiE(JetPt.at(ijet), Jet_eta[ijet], JetPhi.at(ijet), Jet_energy[ijet]);
    if (jet.Pt() < ptcut)     return false;
    if (abs(jet.Eta()) > 2.4) return false;
    if (Get<Int_t>("Jet_id",ijet) <= 0)    return false;
    //if (Jet_id[ijet] > 0) return true;
    // Remove jets close to all selected leptons...
    for(unsigned int i = 0; i < Lepton.size(); i++){
        if(jet.DeltaR(Lepton[i].p) < minDR) return false;
    }
    return true;
    //return false;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void TTHAnalyzer::SelectedGenLepton() {

}

void TTHAnalyzer::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
	TLorentzVector met;
	met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	// set the pfMET to the old MET minus original vector plus new vector
	setMET( (met+oVec-nVec).Pt() );
}

std::vector<int> TTHAnalyzer::CleanedJetIndices(float pt){
	std::vector<int> cleanJetsInd;
	for(Int_t i = 0; i <nJet; i++){
		if (IsGoodJet(i,pt)) cleanJetsInd.push_back(i);
	}
	return cleanJetsInd;
}

void TTHAnalyzer::SmearJetPts(int flag){
	// Modify the jet pt for systematics studies. Either shifted or smeared propagate to the MET!!
	if(gIsData)   return; // don't smear data
	if(flag == 0) return; // 0 makes no sense
	// select the jets you want to have...
	if (!gIsData) {
		std::vector<int> cleanJets = CleanedJetIndices(15.);
		TLorentzVector ojets, jets, tmp, genJet;                            // 4-vec of old jets, newjets and a tmp-vector

		std::vector<int>::const_iterator it = cleanJets.begin();

		for( it = cleanJets.begin(); it != cleanJets.end(); ++it) {
			tmp.SetPtEtaPhiE(JetPt.at(*it), Jet_eta[*it], JetPhi.at(*it), Jet_energy[*it]);         // set temp to the jet
			// Gen info for jets...
			genJet.SetPxPyPzE(Get<Float_t>("Jet_mcPx",*it),
					Get<Float_t>("Jet_mcPy",*it),
					Get<Float_t>("Jet_mcPz",*it),
					Get<Float_t>("Jet_mcEnergy",*it));
			ojets += tmp;
			if(flag == 1) JetPt.at(*it) *= Get<Float_t>("Jet_corr_JECUp",*it);   // vary up   for flag 1
			if(flag == 2) JetPt.at(*it) *= Get<Float_t>("Jet_corr_JECDown",*it); // vary down for flag 2;
			if(flag == 3){
				float jerScaleUp   = getJERScaleUp(*it);
				float jerScale     = getJERScale(*it);
				float factor1 = 0.;
				if (genJet.DeltaR(tmp) < 0.5) factor1 = max(0.0, genJet.Pt() + jerScale*(tmp.Pt() - genJet.Pt()) );
				else                          factor1 = tmp.Pt();
				float sigmaMC  = getErrPt(JetPt.at(*it), Jet_eta[*it]) / JetPt.at(*it);
				float factor   = fRand3->Gaus(1., sqrt(jerScale*jerScale -1.) * sigmaMC );
				JetPt.at(*it) = JetPt.at(*it) * factor;		// smear for flag 3
			}
			// set tmp to the scaled/smeared jet
			tmp.SetPtEtaPhiE(JetPt.at(*it), Jet_eta[*it], JetPhi.at(*it), Jet_energy[*it]);
			jets += tmp;  // add scaled/smeared jet to the new jets
		}
		propagateMET(jets, ojets);  // propagate this change to the MET
	}
}

void TTHAnalyzer::ScaleLeptons(int flag){
	// Shift the lepton pts for systematics studies
	if(gIsData) return; // don't smear data
	if(flag == 0) return;
	//float scale = 0.003; // 0.3% for muons
	float scale = 0.005;
	TLorentzVector oleps, leps, tmp;
	for(Int_t k = 0; k < nMuon; ++k){ //xxx Muon
		if(TMath::Abs(LepGood_pdgId[k]) != 13) continue;
		tmp.SetPxPyPzE(MuPx.at(k), MuPy.at(k), MuPz.at(k), MuEnergy.at(k));
		oleps += tmp;
		if(flag == 1) { MuPx.at(k) += scale*MuPx.at(k); MuPy.at(k) += scale*MuPy.at(k); }
		if(flag == 2) { MuPx.at(k) -= scale*MuPx.at(k); MuPy.at(k) -= scale*MuPy.at(k); }
		tmp.SetPxPyPzE(MuPx.at(k), MuPy.at(k), MuPz.at(k), MuEnergy.at(k));
		leps += tmp;
	}
	//scale = 0.0015; // 0.15% for electrons
	scale = 0.01;
	for(Int_t i = 0; i < nElec; ++i){ //xxx Elec
		if(TMath::Abs(LepGood_pdgId[i]) != 11) continue;
		tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), ElPz.at(i), ElEnergy.at(i));
		oleps += tmp;
		if(flag == 1) { ElPx.at(i) += scale*ElPx.at(i); ElPy.at(i) += scale*ElPy.at(i); }
		if(flag == 2) { ElPx.at(i) -= scale*ElPx.at(i); ElPy.at(i) -= scale*ElPy.at(i); }
		tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), ElPz.at(i), ElEnergy.at(i));
		leps += tmp;
	}
	propagateMET(leps, oleps);
	return;
}

float TTHAnalyzer::weightNvtx(int nvtx){
	float weight = 1.0;
	if(gIsData) return weight;
	// weights from single lepton region based on nVertex
	if(nvtx== 1)      weight = 3.55238;
	else if(nvtx== 2) weight = 3.55238;
	else if(nvtx== 3) weight = 3.55238;
	else if(nvtx== 4) weight = 3.55238;
	else if(nvtx== 5) weight = 3.55238;
	else if(nvtx== 6) weight = 1.14672;
	else if(nvtx== 7) weight = 1.23742;
	else if(nvtx== 8) weight = 1.18588;
	else if(nvtx== 9) weight = 1.18755;
	else if(nvtx==10) weight = 1.18737;
	else if(nvtx==11) weight = 1.17951;
	else if(nvtx==12) weight = 1.19148;
	else if(nvtx==13) weight = 1.23285;
	else if(nvtx==14) weight = 1.18385;
	else if(nvtx==15) weight = 1.19163;
	else if(nvtx==16) weight = 1.18142;
	else if(nvtx==17) weight = 1.12507;
	else if(nvtx==18) weight = 1.09837;
	else if(nvtx==19) weight = 1.0351;
	else if(nvtx==20) weight = 0.976103;
	else if(nvtx==21) weight = 0.876084;
	else if(nvtx==22) weight = 0.807204;
	else if(nvtx==23) weight = 0.71911;
	else if(nvtx==24) weight = 0.624918;
	else if(nvtx==25) weight = 0.512945;
	else if(nvtx==26) weight = 0.449579;
	else if(nvtx==27) weight = 0.356727;
	else if(nvtx==28) weight = 0.288227;
	else if(nvtx==29) weight = 0.237089;
	else if(nvtx==30) weight = 0.187185;
	else if(nvtx==31) weight = 0.155465;
	else if(nvtx==32) weight = 0.106072;
	else if(nvtx==33) weight = 0.0863709;
	else if(nvtx==34) weight = 0.0784446;
	else if(nvtx==35) weight = 0.0731524;
	else if(nvtx==36) weight = 0.0533728;
	else if(nvtx==37) weight = 0.0478281;
	else if(nvtx==38) weight = 0.0174096;
	else if(nvtx==39) weight = 0.0205018;
	else if(nvtx==40) weight = 0.0164379;
	else if(nvtx==42) weight = 0.0156376;
	else if(nvtx==43) weight = 0.0239829;
	return weight;
}

void TTHAnalyzer::ScaleMET(int flag){
	// first try on MET uncertainty
	if(gIsData) return; // don't scale data
	TLorentzVector umet, jets, leps, tmp;
	jets.SetPtEtaPhiM(0., 0., 0., 0.); // init
	leps.SetPtEtaPhiM(0., 0., 0., 0.); // init
	tmp.SetPtEtaPhiM(0., 0., 0., 0.);  // init
	umet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.); // add met
	// subtract uncleaned jets
	for (Int_t i=0; i<nJet; i++) {
		if (!IsGoodJet(i, 15.)) continue; // do this on all jets in the event, not only the good jets with pT > 40
		tmp.SetPxPyPzE(Jet_px[i], Jet_py[i], Jet_pz[i], Jet_energy[i]);
		umet += tmp;
		jets += tmp;
		tmp.SetPtEtaPhiE(0., 0., 0., 0.);
	}
	// subtract muons
	for (Int_t i=0; i<nLepGood; i++) {
		if (!IsTightMuon(i)) continue;
		tmp.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]);
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// subtract electrons
	for (Int_t i=0; i<nLepGood; i++) {
		if (!IsTightElectron(i)) continue;
		tmp.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]);
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// scale the unclustered energy by 10%
	if (flag == 0) tmp.SetPtEtaPhiE(1.1 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	if (flag == 1) tmp.SetPtEtaPhiE(0.9 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	// subtract the leptons and jets again
	tmp -= leps;
	tmp -= jets;
	setMET(tmp.Pt());
	return;
}
