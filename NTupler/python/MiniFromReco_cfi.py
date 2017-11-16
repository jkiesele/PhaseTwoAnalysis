import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer('MiniFromReco',
        # TODO: implement barrel and endcap collections
        electrons    = cms.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl"),
        beamspot     = cms.InputTag("offlineBeamSpot"),
        conversions  = cms.InputTag("particleFlowEGamma"),
        muons        = cms.InputTag("muons"),
        puppiNoLepIsolationChargedHadrons = cms.InputTag("muonIsolationPUPPINoLep","h+-DR040-ThresholdVeto000-ConeVeto000"),
        puppiNoLepIsolationNeutralHadrons = cms.InputTag("muonIsolationPUPPINoLep","h0-DR040-ThresholdVeto000-ConeVeto001"),
        puppiNoLepIsolationPhotons        = cms.InputTag("muonIsolationPUPPINoLep","gamma-DR040-ThresholdVeto000-ConeVeto001"),    
        pfCandsNoLep = cms.InputTag("particleFlow"),
        jets         = cms.InputTag("ak4PFJetsCHS"),
        met          = cms.InputTag("pfMet"),
        genParts     = cms.InputTag("genParticles"),
        genJets      = cms.InputTag("ak4GenJets"),
        vertices     = cms.InputTag("offlinePrimaryVertices"),
        photonsBarrel= cms.InputTag("gedPhotons"),
        phoBarrelMva = cms.InputTag("hgcPhotonMVAbarrel"),
        photonsEndcap= cms.InputTag("photonsFromMultiCl"),
        phoEndcapMva = cms.InputTag("hgcPhotonMVAendcap"),
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
