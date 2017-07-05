import FWCore.ParameterSet.Config as cms

jetfilter = cms.EDProducer('PatJetFilter',
        electrons     = cms.InputTag("slimmedElectrons"),
        muons         = cms.InputTag("slimmedMuons"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
)
