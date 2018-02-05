cmsenv
source $CMSSW_BASE/src/PhaseTwoAnalysis/delphesInterface/env.sh
cd $CMSSW_BASE/src/PhaseTwoAnalysis/delphesInterface/ntupler/
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "${SCRIPT}")
export DANALYSISPATH=$BASEDIR
export LD_LIBRARY_PATH=$BASEDIR:$LD_LIBRARY_PATH
export PATH=$PATH:$BASEDIR
