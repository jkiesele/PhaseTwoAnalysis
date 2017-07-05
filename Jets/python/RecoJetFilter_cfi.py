import FWCore.ParameterSet.Config as cms

jetfilter = cms.EDProducer('RecoJetFilter',
        electrons    = cms.InputTag("ecalDrivenGsfElectrons"),
        muons        = cms.InputTag("muons"),
        jets         = cms.InputTag("ak4PFJetsCHS"),
)
