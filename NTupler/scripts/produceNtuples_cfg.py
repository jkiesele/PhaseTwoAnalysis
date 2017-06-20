import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFormat', 'PAT',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "format of the input files (PAT or RECO)"
                 )
options.register('skim', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "skim events with one lepton and 2 jets"
                 )
options.parseArguments()

process = cms.Process("MiniAnalysis")

# Geometry, GT, and other standard sequences
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v1', '')

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
)

# Input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*(
        '/store/user/ebouvier/RelValTTbar_14TeV/crab_UPG_CheckPat_miniAOD-prod_RelValTTbar/170612_140401/0000/miniAOD-prod_PAT_1.root',
    ))
)
if (options.inputFormat.lower() == "reco"):
    process.source.fileNames = cms.untracked.vstring(*(
        '/store/relval/CMSSW_9_1_1_patch1/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_91X_upgrade2023_realistic_v1_D17PU200r1-v1/10000/00052551-024E-E711-B071-0242AC130002.root',
    ))
process.source.inputCommands = cms.untracked.vstring("keep *")

# Pre-skim weight counter
process.weightCounter = cms.EDAnalyzer('WeightCounter')

# Skim filter
muonLabel = "slimmedMuons"
elecLabel = "slimmedElectrons"
jetLabel = "slimmedJets"
if (options.inputFormat.lower() == "reco"):
    muonLabel = "muons"
    elecLabel = "ecalDrivenGsfElectrons"
    jetLabel = "ak4PUPPIJets"
process.selectedMuons = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag(muonLabel),
                                     cut = cms.string("pt>10 && abs(eta)<3")
                                     )
process.selectedElectrons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag(elecLabel),
                                         cut = cms.string("pt>10 && abs(eta)<3")
                                         )
process.selectedJets = cms.EDFilter("CandPtrSelector",
                                    src = cms.InputTag(jetLabel),
                                    cut = cms.string("pt>20 && abs(eta)<5")
                                    )
process.allLeps = cms.EDProducer("CandViewMerger",
                                 src = cms.VInputTag(
                                                     cms.InputTag("selectedElectrons"),
                                                     cms.InputTag("selectedMuons")
                                                     )
                                 )
process.countLeps = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("allLeps"),
                                 minNumber = cms.uint32(1)
                                 )
process.countJets = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedJets"),
                                 minNumber = cms.uint32(2)
                                 )
process.preYieldFilter = cms.Sequence(process.selectedMuons+process.selectedElectrons+process.allLeps+process.countLeps+process.selectedJets+process.countJets)


# run Puppi 
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PhotonPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')
from CommonTools.PileupAlgos.PhotonPuppi_cff        import setupPuppiPhoton
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies(process)

# recluster jets
process.load('RecoJets/Configuration/RecoPFJets_cff')
process.ak4PUPPIJets  = process.ak4PFJets.clone(rParam=0.4, src = cms.InputTag('puppi'))

# recompute MET
process.load('RecoMET.METProducers.PFMET_cfi')
process.puppiMet = process.pfMet.clone()
process.puppiMet.src = cms.InputTag('puppi')

process.puSequence = cms.Sequence(process.pfNoLepPUPPI * process.puppi * process.puppiNoLep * process.ak4PUPPIJets * process.puppiMet)

# PF cluster producer for HFCal ID
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff")

# jurassic track isolation
# https://indico.cern.ch/event/27568/contributions/1618615/attachments/499629/690192/080421.Isolation.Update.RecHits.pdf
process.load("RecoEgamma.EgammaIsolationAlgos.electronTrackIsolationLcone_cfi")
process.electronTrackIsolationLcone.electronProducer = cms.InputTag("ecalDrivenGsfElectrons")
process.electronTrackIsolationLcone.intRadiusBarrel = 0.04
process.electronTrackIsolationLcone.intRadiusEndcap = 0.04

# analysis
moduleName = "MiniFromPat"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "MiniFromReco"
process.analysis = cms.EDAnalyzer(moduleName)
process.load("PhaseTwoAnalysis.NTupler."+moduleName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    process.analysis.jets = "ak4PUPPIJets"
    process.analysis.pfCands = "puppi"
    process.analysis.pfCandsNoLep = "puppiNoLep"
    process.analysis.met = "puppiMet"

# output
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )

# run
if options.skim:
    if (options.inputFormat.lower() == "reco"):
        process.p = cms.Path(process.weightCounter * process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.preYieldFilter * process.analysis)
    else:
        process.p = cms.Path(process.weightCounter*process.preYieldFilter*process.analysis)
else:
    if (options.inputFormat.lower() == "reco"):
        process.p = cms.Path(process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.analysis)
    else:
        process.p = cms.Path(process.analysis)
