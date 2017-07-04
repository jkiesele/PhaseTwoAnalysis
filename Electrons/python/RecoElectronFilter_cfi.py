import FWCore.ParameterSet.Config as cms

electronfilter = cms.EDProducer('RecoElectronFilter',
        electrons    = cms.InputTag("ecalDrivenGsfElectrons"),
        beamspot     = cms.InputTag("offlineBeamSpot"),
        conversions  = cms.InputTag("particleFlowEGamma"),
        trackIsoValueMap = cms.InputTag("electronTrackIsolationLcone"),
        pfCandsNoLep = cms.InputTag("particleFlow"),
        genParts     = cms.InputTag("genParticles"),
        vertices     = cms.InputTag("offlinePrimaryVertices"),
        HGCalIDToolConfig = cms.PSet(
            HGCBHInput = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
            HGCEEInput = cms.InputTag("HGCalRecHit","HGCEERecHits"),
            HGCFHInput = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
            HGCPFRecHits = cms.InputTag("particleFlowRecHitHGC::ElectronFilter"),
            withPileup = cms.bool(True),
            debug = cms.bool(False),
        ),
)
