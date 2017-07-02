import FWCore.ParameterSet.Config as cms

muonfilter = cms.EDProducer('RecoMuonFilter',
        vertices      = cms.InputTag("offlinePrimaryVertices"),
        muons         = cms.InputTag("muons"),
)
