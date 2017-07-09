import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'FilteredEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFormat', 'PAT',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "format of the input files (PAT or RECO)"
                 )
options.parseArguments()


process = cms.Process("MuonFilter")

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
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

# run Puppi 
process.load('CommonTools.PileupAlgos.Puppi_cff')
process.particleFlowNoLep = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
process.puppiNoLep = process.puppi.clone(candName = cms.InputTag('particleFlowNoLep'))

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")

# producer
moduleName = "PatMuonFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "RecoMuonFilter"
process.muonfilter = cms.EDProducer(moduleName)
process.load("PhaseTwoAnalysis.Muons."+moduleName+"_cfi")

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_*_*_*',
                                           'drop patMuons_slimmedMuons_*_*',
                                           'drop recoMuons_muons_*_*'),
    fileName = cms.untracked.string(options.outFilename)
)
  
if (options.inputFormat.lower() == "reco"):
    process.p = cms.Path(process.primaryVertexAssociation
                         +process.puppi
                         +process.particleFlowNoLep+process.puppiNoLep
                         +process.offlineSlimmedPrimaryVertices+process.packedPFCandidates
                         +process.muonIsolationPUPPI+process.muonIsolationPUPPINoLep * process.muonfilter)
else:
    process.p = cms.Path(process.muonfilter)

process.e = cms.EndPath(process.out)
