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
        '/store/mc/PhaseIITDRSpring17MiniAOD/TTToSemiLepton_TuneCUETP8M1_14TeV-powheg-pythia8/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v1/120000/008BFDF2-285E-E711-8055-001E674FC887.root',
    ))
)
if (options.inputFormat.lower() == "reco"):
    process.source.fileNames = cms.untracked.vstring(*(
        '/store/mc/PhaseIITDRSpring17DR/TTToSemiLepton_TuneCUETP8M1_14TeV-powheg-pythia8/AODSIM/PU200_91X_upgrade2023_realistic_v3-v1/120000/000CD008-7A58-E711-82DB-1CB72C0A3A61.root',
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
