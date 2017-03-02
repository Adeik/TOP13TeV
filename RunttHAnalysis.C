// Load DatasetManager in ROOT 6
R__LOAD_LIBRARY(DatasetManager/DatasetManager.C+)
#include "DatasetManager/DatasetManager.h"

/*******************************************************************************
	* Main function
*******************************************************************************/
void RunttHAnalysis(TString		sampleName		=	"ZZ_ext"	,
					Int_t		nSlots         	=  	1		,
					Long64_t 	nEvents        	= 	0		) {

	// Variables to be used as parameters
  	Float_t G_Total_Lumi    = 11300;  // (According with 03-03 Draft Note).
	Float_t G_Event_Weight  = 1.0;
	Bool_t  G_IsData        = false;

	// PAF mode choice and creation of project
	//--------------------------------------------------------------------------
	PAFIExecutionEnvironment* pafmode = 0;
    cout << endl;
	PAF_INFO("RunttHAnalysis", "======================================== Preprocess");
	if (nSlots <=1 ) {
    	PAF_INFO("RunttHAnalysis", "+ Sequential mode chosen");
    	pafmode = new PAFSequentialEnvironment();
  	}
  	else if (nSlots <= 64) {
    	PAF_INFO("RunttHAnalysis", "+ PROOF Lite mode chosen");
    	pafmode = new PAFPROOFLiteEnvironment(nSlots);
  	}
  	else {
    	PAF_INFO("RunttHAnalysis", "+ PoD mode chosen");
    	pafmode = new PAFPoDEnvironment(nSlots);
  	}

	PAFProject* myProject = new PAFProject(pafmode);

	// Input data sample using DatasetManager
	//--------------------------------------------------------------------------
	TString userhome = "/mnt_pool/fanae105/user/$USER/";
	PAF_INFO("RunttHAnalysis", "+ Configuring DatasetManager");
	DatasetManager* dm = DatasetManager::GetInstance();
	//dm->SetTab("DR80XasymptoticMiniAODv2");
	dm->SetTab("DR80XSummer16asymptoticMiniAODv2_2");

	cout << endl;
	PAF_INFO("RunttHAnalysis", "+ Obtaining sample info:");
	// Deal with data samples
	if (sampleName == "DoubleEG"   || 	sampleName == "DoubleMuon" ||
    sampleName == "MuonEG"	|| 	sampleName.BeginsWith("Single")) {

		PAF_INFO("RunttHAnalysis", "	> The sample is " + sampleName + " DATA");
		TString datasuffix[] = {
		"16B_03Feb2017",
		"16C_03Feb2017",
		"16D_03Feb2017",
		"16E_03Feb2017",
		"16F_03Feb2017",
		"16G_03Feb2017",
		"16H_03Feb2017_v2",
		"16H_03Feb2017_v3"
    	};
	    const UInt_t nDataSamples = 8;
	    for (UInt_t i = 0; i < nDataSamples; i++) {
			TString asample = Form("Tree_%s_%s",sampleName.Data(), datasuffix[i].Data());
			PAF_INFO("RunttHAnalysis", "	> Importing " + asample + " ...");
			myProject->AddDataFiles(dm->GetRealDataFiles(asample));
	    }
	    G_Event_Weight = 1.;
	    G_IsData = true;
  	}
  	else { // Deal with MC samples
		PAF_INFO("RunttHAnalysis", "	> The sample is MC SIMULATION ");
	    G_IsData = false;
	    dm->LoadDataset(sampleName);

	    if (sampleName != "TestHeppy") myProject->AddDataFiles(dm->GetFiles());

    	if (sampleName.Contains("TTWToLNu")	|| 	sampleName == "TTWToQQ" 		||
		sampleName.Contains("TTZToQQ")		||	sampleName == "WWZ" 			||
		sampleName == "WZZ" 				||	sampleName == "ZZZ" 			||
		sampleName.Contains("aMCatNLO") 	|| 	sampleName.Contains("amcatnlo") ||
		sampleName.Contains("TTZToLLNuNu")	||	sampleName.Contains("TTGJets") ) {
			G_Event_Weight 		= dm->GetCrossSection() / dm->GetSumWeights();
			cout << endl;
			cout << " weightSum(MC@NLO) = " << dm->GetSumWeights()     << endl;
    	}
    	else if (sampleName == "TestHeppy") {
			PAF_INFO("RunttHAnalysis", "	> This is a TEST SAMPLE");
			TString	localpath	=	"/pool/ciencias/users/user/palencia/";
			TString sample 		= 	"treeTtbar_jan19.root";
			myProject->AddDataFile(localpath + sample);
			G_Event_Weight 		= 1;
		}
		else {
			G_Event_Weight 		= dm->GetCrossSection() / dm->GetEventsInTheSample();
		}

    	if (nEvents ==0) nEvents= dm->GetEventsInTheSample();
	}
	
    cout << endl;
    cout << " #===============================================" 	<< endl;
    cout << " #          sampleName = " << sampleName               << endl;
	cout << " #	     x-section = " << dm->GetCrossSection()     	<< endl;
	cout << " #	       nevents = " << dm->GetEventsInTheSample()	<< endl;
	cout << " #      base file name = " << dm->GetBaseFileName()    << endl;
	cout << " #	        weight = " << G_Event_Weight	        	<< endl;
	cout << " #	        isData = " << G_IsData	                	<< endl;
	cout << " #===============================================" 	<< endl;
	cout << endl;

	// Output file name
  	//--------------------------------------------------------------------------
	TString outputDir = "./temp";

	gSystem->mkdir(outputDir, kTRUE);

	std::ostringstream oss;
	oss << G_Total_Lumi;

	TString LumiString = oss.str();
  	TString outputFile = outputDir;
  	outputFile += "/Tree_" + sampleName + ".root";

  	PAF_INFO("RunttHAnalysis", Form("+ Output file = %s", outputFile.Data()));
  	myProject->SetOutputFile(outputFile);

  	// Parameters for the analysis
  	//--------------------------------------------------------------------------
  	PAF_INFO("RunttHAnalysis", Form("+ Setting parameters"));
	myProject->SetInputParam("sampleName",    sampleName       );
	myProject->SetInputParam("IsData",        G_IsData         );
	myProject->SetInputParam("weight",        G_Event_Weight   );
	myProject->SetInputParam("TotalLumi",     G_Total_Lumi     );

	if(nEvents != 0) myProject->SetNEvents(nEvents);

	// Name of selector class
	PAF_INFO("RunttHAnalysis", "+ Loading selector and other packages");
	//--------------------------------------------------------------------------
	myProject->AddSelectorPackage("ttHAnalyzer");

	// Additional packages
	//--------------------------------------------------------------------------
	myProject->AddPackage("mt2");
	myProject->AddPackage("PUWeight");
	myProject->AddPackage("BTagSFUtil");
	myProject->AddPackage("SusyLeptonSF");

	// Start analysis
	//--------------------------------------------------------------------------
  	PAF_INFO("RunttHAnalysis", Form("+ Starting compilation"));
	myProject->Run();
}
