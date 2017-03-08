#===============================================================================
#
#                         Analysis of the ttH process
#
#===============================================================================

if [ "$1" == "an" ]; then
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ttH ANALYSIS EXECUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting up the environment"
    root6
    source /opt/PoD/PoD_env.sh
    source /opt/PAF/PAF_setup.sh
    resetpaf

    echo "%%%%%> DONE"

    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting analysis"
    if [ "$2" == "test" ]; then
        root -l -b -q "RunttHAnalysis.C(\"TestHeppy\"							, 1, 0)"
    else
    	echo "%%%%%> MC SAMPLES FOR COMPARISON WITH DATA"
        root -l -b -q "RunttHAnalysis.C(\"TTWToLNu_ext2\"						, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"TTZToLLNuNu_ext\"						, $2, 0)"
        #root -l -b -q "RunttHAnalysis.C(\"TTZToLLNuNu_ext2\"					, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"TTZToQQ\"								, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"WGToLNuG\"							, $2, 0)"	# ***************
        root -l -b -q "RunttHAnalysis.C(\"ZGTo2LG\"								, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"TGJets\"								, $2, 0)"	# ***************
        root -l -b -q "RunttHAnalysis.C(\"TGJets_ext\"							, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"TTGJets\"								, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"TTGJets_ext\"							, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"ZZZ\"									, $2, 0)"	# Rare # ***************
        root -l -b -q "RunttHAnalysis.C(\"WWZ\"									, $2, 0)"	# Rare # ***************
        #root -l -b -q "RunttHAnalysis.C(\"WW\"									, $2, 0)"	# Rare
        root -l -b -q "RunttHAnalysis.C(\"WW_ext\"								, $2, 0)"	# Rare
        root -l -b -q "RunttHAnalysis.C(\"tZq_ll\"								, $2, 0)"	# Rare # ***************
        root -l -b -q "RunttHAnalysis.C(\"TTTT\"								, $2, 0)"	# Rare # ***************

		echo "%%%%%> MC SAMPLES FOR CONTROL REGIONS"
        root -l -b -q "RunttHAnalysis.C(\"TTJets_aMCatNLO\"						, $2, 0)"
        #root -l -b -q "RunttHAnalysis.C(\"DYJetsToLL_M10to50_aMCatNLO\"			, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"DYJetsToLL_M10to50_aMCatNLO_ext\"		, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"WJetsToLNu_MLM\"						, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"WJetsToLNu_MLM\"						, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"TW\"									, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"TW_ext\"								, $2, 0)"
        #root -l -b -q "RunttHAnalysis.C(\"TbarW\"								, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"TbarW_ext\"							, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"T_tch\"								, $2, 0)"	# ***************
        root -l -b -q "RunttHAnalysis.C(\"Tbar_tch\"							, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"TToLeptons_sch_amcatnlo\"				, $2, 0)"	# ***************
        #root -l -b -q "RunttHAnalysis.C(\"WZTo3LNu\"							, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"WWTo2L2Nu\"							, $2, 0)"
        #root -l -b -q "RunttHAnalysis.C(\"ZZ\"									, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"ZZ_ext\"								, $2, 0)"

		echo "%%%%%> MC SIGNAL SAMPLES"
        root -l -b -q "RunttHAnalysis.C(\"TTHNonbb\"                       		, $2, 0)"

		echo "%%%%%> DATA SAMPLES"
        root -l -b -q "RunttHAnalysis.C(\"MuonEG\"                            	, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"DoubleMuon\"                        	, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"DoubleEG\"                          	, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"SingleElec\"                    		, $2, 0)"
        root -l -b -q "RunttHAnalysis.C(\"SingleMuon\"                        	, $2, 0)"
    fi
elif [ "$1" == "plot" ]; then
    echo "Ya lo haré"
else
    echo "ERROR - No valid arguments given"
    echo "Please, execute this script with a valid argument"
fi
