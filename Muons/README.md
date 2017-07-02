Recipes for muon selection in analysis of simulated events with the CMS phase 2 detector.
=========================

Two vector of muons are added: one for muons passing the loose ID, one for muons passing the tight ID.
One should run:
```bash
cmsRun test/ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.

The main producers are:
   * `plugins/PatMuonFilter.cc` -- to run over PAT events 
   * `plugins/RecoMuonFilter.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

