Repository for collecting recipes for standard physics objects for analysis of simulated events with the CMS phase 2 detector.
=========================


Contains (in the future)

a) scripts/recipes to set up the correct CMSSW environment

b) a collection of latest recipes and a CMSSW config to apply them to get collections of standard objects suitable for physics analysis as CMSSW producers (e.g. recommendedTightMuons = tightMuonProducer(slimmedMuons) etc.)

c) an ntuple production configuration file that produces an ntuple that can be analysed in a similar manner als Delphes samples (with the DAnalysis framework)


Installation
--------------

```bash
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git cms-init
git cms-merge-topic -u nsmith-:EgammaFromMultiCl_932v2
mkdir -p RecoEgamma && pushd RecoEgamma
git clone -b integrated https://github.com/nsmith-/Phase2InterimID.git
popd
git clone https://github.com/jkiesele/PhaseTwoAnalysis.git
scram b -j 8
```

As the global tag contains Run-2 JEC, you might want to download an SQLite file to rerun JEC. See the [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/Phase2HGCRecipes) for more details.

Producing flat ntuples (only recommended mode for the moment)
-----------------

Flat ntuples can be produced in the `NTupler` folder, either from PAT or RECO events, by running interactively:
```bash
cmsRun scripts/produceNtuples_cfg.py skim=False/True outFilename=MiniEvents.root inputFormat=RECO/PAT
```
If you want to rerun JEC, you can use the `updateJEC` argument with the path to the SQLite file.

The `skim` flag can be used to reduce the size of the output files. A histogram containing the number of events before the skim is then stored in the output files. By default, events are required to contain at least 1 lepton and 2 jets, but this can be easily modified ll.71-97 of `src/produceNtuples_cfg.py`.

The structure of the output tree can be seen/modified in `interface/MiniEvent.h` and `src/MiniEvent.cc`.

The main analyzers are:
   * `plugins/MiniFromPat.cc` -- to run over PAT events 
   * `plugins/MiniFromReco.cc` -- to run over RECO events (not recommended)

Details on the object definitions are given in the `implementation` section.

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

To adjust the input parameters of `scripts/produceNtuples_cfg.py`, the three following fields need to be modified consistently:
   * `config.JobType.pyCfgParams`
   * `config.JobType.inputFiles`
   * `config.JobType.outputFiles`

If you experience problems with multicrab, please consult the following link:
https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#Multiple_submission_fails_with_a


There is also an example multicrab script provided in the NTupler directory
