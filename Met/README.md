Recipes for met selection in analysis of simulated events with the CMS phase 2 detector.
=========================

Please, note that no JEC is applied.

For miniAOD files, one can simply use the `pat::MET` collection associated to the `slimmedMETsPuppi` input tag.

For RECO files, one should run:
```bash
cmsRun ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO
```

to recluster MET after applying PUPPI. One will get a `reco::PFMET` collection associated to the `puppiMet` input tag, while the initial vector of PF MET is dropped to avoid any confusion.
