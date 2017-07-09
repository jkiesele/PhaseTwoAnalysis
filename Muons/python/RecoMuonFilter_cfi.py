import FWCore.ParameterSet.Config as cms

muonfilter = cms.EDProducer('RecoMuonFilter',
        vertices      = cms.InputTag("offlinePrimaryVertices"),
        muons         = cms.InputTag("muons"),
        puppiNoLepIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPINoLep","h+-DR040-ThresholdVeto000-ConeVeto000"),
        puppiNoLepIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPINoLep","h0-DR040-ThresholdVeto000-ConeVeto001"),
        puppiNoLepIsolationPhotons        = cms.InputTag("muonIsolationPUPPINoLep","gamma-DR040-ThresholdVeto000-ConeVeto001"),    
)

IsoConeDefinitions = cms.VPSet(
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.4),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.0001),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.4),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.01),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('MuonPFIsolationWithConeVeto'),
                  coneSize = cms.double(0.4),
                  VetoThreshold = cms.double(0.0),
                  VetoConeSize = cms.double(0.01),
                  isolateAgainst = cms.string('gamma'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),                  
)

muonIsolationPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
                srcToIsolate = cms.InputTag("muons"),
                srcForIsolationCone = cms.InputTag('packedPFCandidates'),
                puppiValueMap = cms.InputTag(''),
                usePUPPINoLepton = cms.bool(False),
                isolationConeDefinitions = IsoConeDefinitions
)
muonIsolationPUPPINoLep = muonIsolationPUPPI.clone(usePUPPINoLepton = cms.bool(True))
