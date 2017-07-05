Recipes for muon selection in analysis of simulated events with the CMS phase 2 detector.
=========================

One should run:
```bash
cmsRun ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.

The main producers are:
   * `plugins/PatMuonFilter.cc` -- to run over PAT events 
   * `plugins/RecoMuonFilter.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

For each ID quality, a vector of muons and a vector of double corresponding to the muon relative isolation are added: 
   * `doubles_muonfilter_LooseMuonRelIso_MuonFilter`
   * `doubles_muonfilter_MediumMuonRelIso_MuonFilter`
   * `doubles_muonfilter_TightMuonRelIso_MuonFilter`
   * `[reco|pat]Muons_muonfilter_LooseMuons_MuonFilter`
   * `[reco|pat]Muons_muonfilter_MediumMuons_MuonFilter`
   * `[reco|pat]Muons_muonfilter_TightMuons_MuonFilter`

The initial vector of muons is dropped to avoid any confusion.
