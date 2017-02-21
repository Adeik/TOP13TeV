#===============================================================================
#
#                         Analysis of the ttH process
#
#===============================================================================

if [ "$1" == "an" ]; then
    #------------------------------------ Setting up environment
    root6
    source /opt/PoD/PoD_env.sh
    source /opt/PAF/PAF_setup.sh
    resetpaf

    #------------------------------------ Initiating analysis
    if [ "$2" == "test" ]; then
        root -l -b -q 'RunttHAnalysis.C("TestHeppy", 1, 0)'
    else
		# MC
        root -l -b -q 'RunttHAnalysis.C("ZZ"                                , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("WW"                                , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("WZ"                                , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TW"                                , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TbarW"                             , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("DYJetsToLL_M50_aMCatNLO"           , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("DYJetsToLL_M10to50_aMCatNLO_ext"   , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("WJetsToLNu_aMCatNLO", 5, 0)'

        root -l -b -q 'RunttHAnalysis.C("TTWToLNu"	                        , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTZToQQ"	                        , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTGJets"	                        , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTWToQQ"	                        , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTZToLLNuNu"                       , $2, 0)'

		# Data
        root -l -b -q 'RunttHAnalysis.C("MuonEG"                            , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("DoubleMuon"                        , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("DoubleEG"                          , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("SingleElectron"                    , $2, 0)'
        root -l -b -q 'RunttHAnalysis.C("SingleMuon"                        , $2, 0)'
    fi
elif [ "$1" == "plot" ]; then
    echo "Ya lo haré"
else
    echo "ERROR - No valid arguments given"
    echo "Please, execute this script with a valid argument"
fi
