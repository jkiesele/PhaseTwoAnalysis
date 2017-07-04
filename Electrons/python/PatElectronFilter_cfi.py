import FWCore.ParameterSet.Config as cms

electronfilter = cms.EDProducer('PatElectronFilter',
        electrons     = cms.InputTag("slimmedElectrons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
)
