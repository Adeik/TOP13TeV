#===============================================================================
#
#                         Analysis of the ttH process
#
#===============================================================================
$ncores  =   1

if [ "$1" == "an" ]; then
    #------------------------------------ Setting up environment
    root6
    source /opt/PoD/PoD_env.sh
    source /opt/PAF/PAF_setup.sh
    resetpaf

    #------------------------------------ Initiating analysis
    if [ "$2" == "test" ]; then
        root -l -b -q 'RunTTHAnalysis.C("TestHeppy", 1, 0)'
    else
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_ext"                  , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_ext"                  , $ncores, 1)' #Semileptonicselection
        root -l -b -q 'RunTTHAnalysis.C("ZZ"                                , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("WW"                                , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("WZ"                                , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TW"                                , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TbarW"                             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("DYJetsToLL_M50_aMCatNLO"           , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("DYJetsToLL_M10to50_aMCatNLO_ext"   , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("WJetsToLNu_aMCatNLO", 5, 0)'

        root -l -b -q 'RunTTHAnalysis.C("TTWToLNu"	                        , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTZToQQ"	                        , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTGJets"	                        , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTWToQQ"	                        , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTZToLLNuNu"                       , $ncores, 0)'

        root -l -b -q 'RunTTHAnalysis.C("MuonEG"                            , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("DoubleMuon"                        , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("DoubleEG"                          , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("SingleElectron"                    , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("SingleMuon"                        , $ncores, 0)'

        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_scaleDown"            , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_scaleUp"              , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1665"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1695"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1715"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1735"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1755"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_mtop1785"             , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTJets_aMCatNLO"	                , $ncores, 0)'
        root -l -b -q 'RunTTHAnalysis.C("TTbar_Powheg_Herwig"               , $ncores, 0)'
    fi
elif [ "$1" == "plot" ]; then
    echo "Ya lo haré"

elif [ "$1" == "deepAnal" ]; then
    echo "Ya lo haré"

elif [ "$1" == "datacards" ]; then
    echo "Ya lo haré"
else
    echo "ERROR - No valid arguments given\n"
    echo "Please, execute this script with a valid argument"
fi
