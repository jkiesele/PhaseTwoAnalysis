import FWCore.ParameterSet.Config as cms

myana = cms.EDAnalyzer('BasicRecoDistrib',
        electrons    = cms.InputTag("ecalDrivenGsfElectrons"),
        beamspot     = cms.InputTag("offlineBeamSpot"),
        conversions  = cms.InputTag("particleFlowEGamma"),
        trackIsoValueMap = cms.InputTag("electronTrackIsolationLcone"),
        muons        = cms.InputTag("muons"),
        pfCands      = cms.InputTag("particleFlow"),
        pfCandsNoLep = cms.InputTag("particleFlow"),
        jets         = cms.InputTag("ak4PFJetsCHS"),
        met          = cms.InputTag("pfMet"),
        genParts     = cms.InputTag("genParticles"),
        genJets      = cms.InputTag("ak4GenJets"),
        vertices     = cms.InputTag("offlinePrimaryVertices"),
        HGCalIDToolConfig = cms.PSet(
            HGCBHInput = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
            HGCEEInput = cms.InputTag("HGCalRecHit","HGCEERecHits"),
            HGCFHInput = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
            HGCPFRecHits = cms.InputTag("particleFlowRecHitHGC::MyAna"),
            withPileup = cms.bool(True),
            debug = cms.bool(False),
        ),
        

)
