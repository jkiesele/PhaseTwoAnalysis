Recipes for electron selection in analysis of simulated events with the CMS phase 2 detector.
=========================

One should run:
```bash
cmsRun ConfFile_cfg.py outFilename=FilteredEvents.root inputFormat=RECO/PAT
```

whether the input file format is RECO or miniAOD.

The main producers are:
   * `plugins/PatElectronFilter.cc` -- to run over PAT events 
   * `plugins/RecoElectronFilter.cc` -- to run over RECO events 

Details on the object definitions are given in the `implementation` section.

For each ID quality, a vector of muons and a vector of double corresponding to the muon relative isolation are added: 
   * `doubles_electronfilter_LooseElectronRelIso_ElectronFilter`
   * `doubles_electronfilter_MediumElectronRelIso_ElectronFilter`
   * `doubles_electronfilter_TightElectronRelIso_ElectronFilter`
   * `[recoGsf|pat]Electrons_electronfilter_LooseElectrons_ElectronFilter`
   * `[recoGsf|pat]Electrons_electronfilter_MediumElectrons_ElectronFilter`
   * `[recoGsf|pat]Electrons_electronfilter_TightElectrons_ElectronFilter`

The initial vector of electrons is dropped to avoid any confusion.
