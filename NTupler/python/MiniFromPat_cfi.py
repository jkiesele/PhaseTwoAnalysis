import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer('MiniFromPat',
        vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        electrons     = cms.InputTag("slimmedElectrons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
        muons         = cms.InputTag("slimmedMuons"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
        mets          = cms.InputTag("slimmedMETsPuppi"),
        genParts      = cms.InputTag("packedGenParticles"),
        genJets       = cms.InputTag("slimmedGenJets"),
)
