// Load DatasetManager in ROOT 6
R__LOAD_LIBRARY(DatasetManager/DatasetManager.C+)
#include "DatasetManager/DatasetManager.h"

/*******************************************************************************
	* Main function
*******************************************************************************/
void RunttHAnalysis(TString		sampleName		=	"ZZ"	,
					Int_t		nSlots         	=  	8		,
					Bool_t  	DoSystStudies  	=  	false	,
					Long64_t 	nEvents        	= 	0		) {

	// Variables to be used as parameters
  	Float_t G_Total_Lumi    = 19664.225;
	Float_t G_Event_Weight  = 1.0;
	Bool_t  G_IsData        = false;
	Float_t G_LumiForPUData = 19468.3;	// luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
	Bool_t  G_IsMCatNLO     = false;

	// PAF mode
	//--------------------------------------------------------------------------
	cout << endl;
	PAFIExecutionEnvironment* pafmode = 0;
	if (nSlots <=1 ) {
    	PAF_INFO("RunttHAnalysis", "Sequential mode chosen");
    	pafmode = new PAFSequentialEnvironment();
  	}
  	else if (nSlots <= 8) {
    	PAF_INFO("RunttHAnalysis", "PROOF Lite mode chosen");
    	pafmode = new PAFPROOFLiteEnvironment(nSlots);
  	}
  	else {
    	PAF_INFO("RunttHAnalysis", "PoD mode chosen");
    	pafmode = new PAFPoDEnvironment(nSlots);
  	}

	// Create PAF Project whith that environment
	//--------------------------------------------------------------------------
	PAFProject* myProject = new PAFProject(pafmode);

	// Input data sample
	//--------------------------------------------------------------------------
	TString userhome = "/mnt_pool/fanae105/user/$USER/";
	DatasetManager* dm = DatasetManager::GetInstance();
	//dm->SetTab("DR80XasymptoticMiniAODv2");
	dm->SetTab("DR80XSummer16asymptoticMiniAODv2");

	// Deal with data samples
	if ((sampleName == "DoubleEG"   || 	sampleName == "DoubleMuon" ||
    sampleName == "MuonEG"	|| 	sampleName == "SingleEle"	||
    sampleName == "SingleMu")) {

    	cout << "   + Data..." << endl;
		TString datasuffix[] = { // 17.24
			"Run2016B_PromptReco_v2", // 5.86
			"Run2016C_PromptReco_v2", // 2.64
			"Run2016D_PromptReco_v2", // 4.35
			"Run2016G_PromptReco_v1", // 4.3c
			//"Run2015D_16Dec"
			//"Run2015C_05Oct",
			//"C_7016",
			//"D_7360"
    	};
	    const UInt_t nDataSamples = 4;
	    for (UInt_t i = 0; i < nDataSamples; i++) {
			TString asample = Form("Tree_%s_%s",sampleName.Data(), datasuffix[i].Data());
			cout << "   + Looking for " << asample << " trees..." << endl;
			myProject->AddDataFiles(dm->GetRealDataFiles(asample));
	    }
	    G_Event_Weight = 1.;
	    G_IsData = true;
  	}
  	else { // Deal with MC samples
	    G_IsData = false; //true;  // only for pseudodata
	    dm->LoadDataset(sampleName + "_ext");
	    if (sampleName != "TestHeppy" && !sampleName.Contains("T2tt")) myProject->AddDataFiles(dm->GetFiles());
    	if (sampleName == "TTWToLNu"  || sampleName == "TTWToQQ" || sampleName == "TTZToQQ" ||
	    sampleName == "WWZ" || sampleName == "WZZ" || sampleName == "ZZZ" ||
		sampleName.Contains("aMCatNLO") || sampleName.Contains("amcatnlo") ) {
			G_Event_Weight = dm->GetCrossSection() / dm->GetSumWeights();
			cout << endl;
			cout << " weightSum(MC@NLO) = " << dm->GetSumWeights()     << endl;
    	}
    	else if (sampleName.BeginsWith("T2tt")) {
      		TString lp = "/pool/ciencias/HeppyTreesSummer16/v1/";
      		myProject->AddDataFile(lp + "Tree_" + sampleName + "_ext_0.root");
    	}
    	else if (sampleName == "TestHeppy") {
			TString	localpath	=	"/pool/ciencias/users/user/palencia/";
			TString sample 		= 	"treeTtbar_jan19.root";
			myProject->AddDataFile(localpath + sample);
			G_Event_Weight = 1;
		}
		else {
			G_Event_Weight = dm->GetCrossSection() / dm->GetEventsInTheSample();
		}

    	if (nEvents == 0) nEvents = dm->GetEventsInTheSample();

	    cout << endl;
	    cout << " #===============================================" 	<< endl;
	    cout << " #     sampleName = " << sampleName                	<< endl;
		cout << " #	     x-section = " << dm->GetCrossSection()     	<< endl;
		cout << " #	       nevents = " << dm->GetEventsInTheSample()	<< endl;
    	cout << " # base file name = " << dm->GetBaseFileName()     	<< endl;
		cout << " #	        weight = " << G_Event_Weight	        	<< endl;
		cout << " #	        isData = " << G_IsData	                	<< endl;
    	cout << " #===============================================" 	<< endl;
		cout << endl;
  	}

	// Output file name
  	//--------------------------------------------------------------------------
	Bool_t G_Use_CSVM = true;
	TString outputDir = "./temp";
  	// if(sampleName.BeginsWith("T2tt")) outputDir += "/";

	gSystem->mkdir(outputDir, kTRUE);

	std::ostringstream oss;
	oss << G_Total_Lumi;

	TString LumiString = oss.str();
  	TString outputFile = outputDir;
  	outputFile += "/Tree_" + sampleName + "_ext.root";
  	if (outputFile.Contains("_ext2")) outputFile.ReplaceAll("_ext2","");
  	if (outputFile.Contains("_ext"))  outputFile.ReplaceAll("_ext","");

  	PAF_INFO("RunttHAnalysis", Form("Output file = %s", outputFile.Data()));
  	myProject->SetOutputFile(outputFile);

  	if (sampleName.Contains("aMCatNLO") || sampleName.Contains("amcatnlo") ||
	sampleName == "TTWToLNu"	|| sampleName == "TTWToQQ"	||
    sampleName == "TTZToQQ"		|| sampleName == "WWZ"      ||
    sampleName == "WZZ"			|| sampleName == "ZZZ" ){
    	PAF_INFO("RunttHAnalysis", "This is a MC@NLO sample!");
    	G_IsMCatNLO = true;
  	}

  	// Parameters for the analysis
  	//--------------------------------------------------------------------------
	myProject->SetInputParam("sampleName",    sampleName       );
	myProject->SetInputParam("IsData",        G_IsData         );
	myProject->SetInputParam("UseCSVM",       G_Use_CSVM       );
	myProject->SetInputParam("weight",        G_Event_Weight   );
	myProject->SetInputParam("LumiForPU",     G_LumiForPUData  );
	myProject->SetInputParam("TotalLumi",     G_Total_Lumi     );
	myProject->SetInputParam("DoSystStudies", DoSystStudies    );
	myProject->SetInputParam("IsMCatNLO"    , G_IsMCatNLO      );

	if(nEvents != 0) myProject->SetNEvents(nEvents);

	// Name of analysis class
	//--------------------------------------------------------------------------
	myProject->AddSelectorPackage("ttHAnalyzer");

	// Additional packages
	//--------------------------------------------------------------------------
	myProject->AddPackage("mt2");
	myProject->AddPackage("PUWeight");
	myProject->AddPackage("BTagSFUtil");
	myProject->AddPackage("SusyLeptonSF");

	// Let's rock!
	//--------------------------------------------------------------------------
	myProject->Run();
}
