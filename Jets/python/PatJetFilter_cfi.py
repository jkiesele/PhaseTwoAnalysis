import FWCore.ParameterSet.Config as cms

jetfilter = cms.EDProducer('PatJetFilter',
        pileup        = cms.uint32(200),
        electrons     = cms.InputTag("slimmedElectrons"),
        muons         = cms.InputTag("slimmedMuons"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
)
