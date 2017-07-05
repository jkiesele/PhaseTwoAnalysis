Recipes for (b-tagged) jet selection in analysis of simulated events with the CMS phase 2 detector.
=========================

Please, note that the ak4 and PUPPI algorithms are considered here and no JEC is applied.

One should run:
```bash
cmsRun ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.

The main producers are:
   * `plugins/PatJetFilter.cc` -- to run over PAT events 
   * `plugins/RecoJetFilter.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

The following vectors of PF loose jets are added when running over PAT events: 
   * `patJets_jetfilter_Jets_JetFilter`
   * `patJets_jetfilter_LooseCSVv2Jets_JetFilter`
   * `patJets_jetfilter_MediumCSVv2Jets_JetFilter`
   * `patJets_jetfilter_TightCSVv2Jets_JetFilter`
   * `patJets_jetfilter_LooseDeepCSVJets_JetFilter`
   * `patJets_jetfilter_MediumDeepCSVJets_JetFilter`
   * `patJets_jetfilter_TightDeepCSVJets_JetFilter`

whereas only the following vector of jets is added when running over RECO events:
   * `recoJets_jetfilter_Jets_JetFilter`

The initial vector of jets is dropped to avoid any confusion.
