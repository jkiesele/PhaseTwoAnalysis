Repository for collecting recipes for standard physics objects for analysis of simulated events with the CMS phase 2 detector.
=========================


Contains (in the future)

a) scripts to set up the correct CMSSW environment

b) a collection of latest recipes and a CMSSW config to apply them to get collections of standard objects suitable for physics analysis

c) an ntuple production configuration file that produces an ntuple that can be analysed in a similar manner als Delphes samples (with the DAnalysis framework)


Installation
--------------

```bash
cmsrel CMSSW_9_1_1_patch1
cd CMSSW_9_1_1_patch1/src
cmsenv
git cms-addpkg RecoEgamma/EgammaIsolationAlgos
cd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
cd ..
git clone git@github.com:jkiesele/PhaseTwoAnalysis.git
cp PhaseTwoAnalysis/RecoEgammaFix/* RecoEgamma/EgammaIsolationAlgos/plugins/
scram b -j8
```

How to run PAT on RECO datasets
----------------

The `PatProducer` folder contains a configuration file to produce miniAOD files from RECO files. Interactively, after updating the list of input files, one can run
```bash
cmsRun miniAOD-prod_PAT.py
```
A skeleton of crab configuration file is also provided in this folder. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

Plotting basic distributions from RECO collections
-----------------

A basic EDAnalyzer is available in the `BasicRecoDistrib` folder. Several private functions handle electron and forward muon ID. Lepton isolation is computed with a simple loop over neighbouring particles and there is no b-tagging information. More details are given in the `implementation` section of the `.cc` file.


Plotting basic distributions from PAT collections
-----------------

A basic EDAnalyzer is available in the `BasicPatDistrib` folder. Several private functions handle central electron and forward muon ID. A flag (`useDeepCSV`) can be set to true in the configuration file to use deepCSV rather than CSVv2 as b-tagging discriminant. More details are given in the `implementation` section of the `.cc` file.
After updating the list of input files, the analyzer can be run interactively from the `test` subfolder :
```bash
cmsRun ConFile_cfg.py
```
Plots in a pdf format can be obtained by running:
```bash
root -l 
.L plotIt.C++
plotIt()
```
