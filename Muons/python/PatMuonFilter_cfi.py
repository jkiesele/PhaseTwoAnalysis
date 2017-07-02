import FWCore.ParameterSet.Config as cms

muonfilter = cms.EDProducer('PatMuonFilter',
        vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        muons         = cms.InputTag("slimmedMuons"),
)
