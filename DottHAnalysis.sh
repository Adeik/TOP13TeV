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
        root -l -b -q 'RunttHAnalysis.C("TestHeppy", 1, 0)'
    else
        root -l -b -q 'RunttHAnalysis.C("TTbar_Powheg_ext"                  , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTbar_Powheg_ext"                  , $ncores, 1)' #Semileptonicselection
        root -l -b -q 'RunttHAnalysis.C("ZZ"                                , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("WW"                                , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("WZ"                                , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TW"                                , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TbarW"                             , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("DYJetsToLL_M50_aMCatNLO"           , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("DYJetsToLL_M10to50_aMCatNLO_ext"   , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("WJetsToLNu_aMCatNLO", 5, 0)'

        root -l -b -q 'RunttHAnalysis.C("TTWToLNu"	                        , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTZToQQ"	                        , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTGJets"	                        , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTWToQQ"	                        , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("TTZToLLNuNu"                       , $ncores, 0)'

        root -l -b -q 'RunttHAnalysis.C("MuonEG"                            , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("DoubleMuon"                        , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("DoubleEG"                          , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("SingleElectron"                    , $ncores, 0)'
        root -l -b -q 'RunttHAnalysis.C("SingleMuon"                        , $ncores, 0)'
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
