Repository for collecting recipes for standard physics objects for analysis of simulated events with the CMS phase 2 detector.
=========================


Contains (in the future)

a) scripts/recipes to set up the correct CMSSW environment

b) a collection of latest recipes and a CMSSW config to apply them to get collections of standard objects suitable for physics analysis as CMSSW producers (e.g. recommendedTightMuons = tightMuonProducer(slimmedMuons) etc.)

c) an ntuple production configuration file that produces an ntuple that can be analysed in a similar manner als Delphes samples (with the DAnalysis framework)


Installation
--------------

```bash
cmsrel CMSSW_9_1_1_patch3
cd CMSSW_9_1_1_patch3/src
cmsenv
git cms-addpkg RecoEgamma/EgammaIsolationAlgos
cd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
cd ..
git clone git@github.com:jkiesele/PhaseTwoAnalysis.git
cp PhaseTwoAnalysis/RecoEgammaFix/* RecoEgamma/EgammaIsolationAlgos/plugins/
scram b -j8
```

As the global tag contains Run-2 JEC, you might want to download an SQLite file to rerun JEC. See the [TWiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Jet_Energy_Corrections_JEC) for more details.

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

A basic EDAnalyzer is available in the `BasicRecoDistrib` folder. Several private functions handle electron and forward muon ID. Lepton isolation is computed with a simple loop over neighbouring particles and there is no b-tagging information. Normalization to luminosity is not handled. More details are given in the `implementation` section of the `.cc` file.
After updating the list of input files, the analyzer can be run interactively from the `test` subfolder :
```bash
cmsRun ConfFile_cfg.py
```

If you want to rerun JEC, you can use the `updateJEC` argument with the path to the SQLite file.
Befor the EDAnalyzer, PUPPI is run on the fly and jets are re-clustered. The MET is also recomputed but not exactly with the official recipe (that needs PAT collections).

Plots in a pdf format can be obtained by running:
```bash
root -l 
.L plotIt.C++
plotIt()
```

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```


Plotting basic distributions from PAT collections
-----------------

A basic EDAnalyzer is available in the `BasicPatDistrib` folder. Several private functions handle central electron and forward muon ID. A flag (`useDeepCSV`) can be set to true in the configuration file to use deepCSV rather than MVAv2 as b-tagging discriminant. Normalization to luminosity is not handled. More details are given in the `implementation` section of the `.cc` file.
After updating the list of input files, the analyzer can be run interactively from the `test` subfolder :
```bash
cmsRun ConfFile_cfg.py
```

If you want to rerun JEC, you can use the `updateJEC` argument with the path to the SQLite file.

Plots in a pdf format can be obtained by running:
```bash
root -l 
.L plotIt.C++
plotIt()
```

A skeleton of crab configuration file is also provided. The following fields need to be updated:
   * `config.General.requestName` 
   * `config.Data.inputDataset`
   * `config.Data.outLFNDirBase`

before running
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig.py
```

Producing flat ntuples
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
   * `plugins/MiniFromReco.cc` -- to run over RECO events 

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

Producing edm ntuples
-----------------

CMS edm ntuples can be produced in the `NTupler` folder too, either from PAT or RECO events, by running interactively:
```bash
cmsRun scripts/edmFilter_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.
If you want to rerun JEC, you can use the `updateJEC` argument with the path to the SQLite file.

To run over PAT events, the main producers are:
   * `../Electrons/plugins/PatElectronFilter.cc` 
   * `../Muons/plugins/PatMuonFilter.cc` 
   * `../Jets/plugins/PatJetFilter.cc`

while for RECO events, they are: 
   * `../Electrons/plugins/RecoElectronFilter.cc` 
   * `../Muons/plugins/RecoMuonFilter.cc` 
   * `../Jets/plugins/RecoJetFilter.cc` 

Details on the object definitions are given in the `implementation` section, but, please, already note that no JEC is applied and ak4 and PUPPI algorithms are considered everywhere.

When running over PAT events, the following collections are produced:
   * `doubles_electronfilter_LooseElectronRelIso_EDMFilter`
   * `doubles_electronfilter_MediumElectronRelIso_EDMFilter`
   * `doubles_electronfilter_TightElectronRelIso_EDMFilter`
   * `patElectrons_electronfilter_LooseElectrons_EDMFilter`
   * `patElectrons_electronfilter_MediumElectrons_EDMFilter`
   * `patElectrons_electronfilter_TightElectrons_EDMFilter`
   * `doubles_muonfilter_LooseMuonRelIso_EDMFilter`
   * `doubles_muonfilter_MediumMuonRelIso_EDMFilter`
   * `doubles_muonfilter_TightMuonRelIso_EDMFilter`
   * `patMuons_muonfilter_LooseMuons_EDMFilter`
   * `patMuons_muonfilter_MediumMuons_EDMFilter`
   * `patMuons_muonfilter_TightMuons_EDMFilter`
   * `patJets_jetfilter_Jets_EDMFilter`
   * `patJets_jetfilter_LooseMVAv2Jets_EDMFilter`
   * `patJets_jetfilter_MediumMVAv2Jets_EDMFilter`
   * `patJets_jetfilter_TightMVAv2Jets_EDMFilter`
   * `patJets_jetfilter_LooseDeepCSVJets_EDMFilter`
   * `patJets_jetfilter_MediumDeepCSVJets_EDMFilter`
   * `patJets_jetfilter_TightDeepCSVJets_EDMFilter`

while for RECO events, they are:
   * `doubles_electronfilter_LooseElectronRelIso_EDMFilter`
   * `doubles_electronfilter_MediumElectronRelIso_EDMFilter`
   * `doubles_electronfilter_TightElectronRelIso_EDMFilter`
   * `recoGsfElectrons_electronfilter_LooseElectrons_EDMFilter`
   * `recoGsfElectrons_electronfilter_MediumElectrons_EDMFilter`
   * `recoGsfElectrons_electronfilter_TightElectrons_EDMFilter`
   * `doubles_muonfilter_LooseMuonRelIso_EDMFilter`
   * `doubles_muonfilter_MediumMuonRelIso_EDMFilter`
   * `doubles_muonfilter_TightMuonRelIso_EDMFilter`
   * `recoMuons_muonfilter_LooseMuons_EDMFilter`
   * `recoMuons_muonfilter_MediumMuons_EDMFilter`
   * `recoMuons_muonfilter_TightMuons_EDMFilter`
   * `recoPFJets_jetfilter_Jets_EDMFilter`
   * `recoPFMETs_puppiMet__EDMFilter`

The initial vectors of electrons, muons, jets (and PFMETs) are dropped to avoid any confusion.
