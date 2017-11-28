

OLDDIR=`pwd`/ntupler
cd $OLDDIR

#DAnalysis
DANALYSISPATH=$CMSSW_BASE/src/PhaseTwoAnalysis/delphesInterface/DAnalysis-v.1.1_rc2

#CMSSW and Delphes
cd $DANALYSISPATH/../
cd delphes
source DelphesEnv.sh
export DELPHES_PATH=`pwd`
cd $OLDDIR


export PYTHIA8=$CMSSW_RELEASE_BASE/../../../external/pythia8/212-ikhhed3
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
export DANALYSISPATH=$DANALYSISPATH
export LD_LIBRARY_PATH=$DANALYSISPATH:$LD_LIBRARY_PATH
export PATH=$PATH:$DANALYSISPATH
export LD_LIBRARY_PATH=$CMSSW_BASE/src/PhaseTwoAnalysis/delphesInterface/ntupler:$LD_LIBRARY_PATH

