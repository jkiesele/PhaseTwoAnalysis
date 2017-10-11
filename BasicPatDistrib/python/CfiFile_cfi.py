import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('BasicPatDistrib',
        pileup        = cms.uint32(200),
        vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        electrons     = cms.InputTag("slimmedElectrons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
#        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions"),
        muons         = cms.InputTag("slimmedMuons"),
#        jets          = cms.InputTag("slimmedJets"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
        useDeepCSV    = cms.bool(False),
#        mets          = cms.InputTag("slimmedMETs"),
        mets          = cms.InputTag("slimmedMETsPuppi"),
        genParts      = cms.InputTag("packedGenParticles"),
        genJets       = cms.InputTag("slimmedGenJets"),
)
