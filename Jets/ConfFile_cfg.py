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


process = cms.Process("JetFilter")

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

process.puSequence = cms.Sequence(process.puppi * process.ak4PUPPIJets )

# producer
moduleName = "PatJetFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "RecoJetFilter"
process.jetfilter = cms.EDProducer(moduleName)
process.load("PhaseTwoAnalysis.Jets."+moduleName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    process.jetfilter.jets = "ak4PUPPIJets"    

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_*_*_*',
                                           'drop patJets_slimmedJetsPuppi_*_*',
                                           'drop reco*_ak4*Jets*_*_*'),
    fileName = cms.untracked.string(options.outFilename)
)
  
if (options.inputFormat.lower() == "reco"):
    process.p = cms.Path(process.puSequence * process.jetfilter)
else:
    process.p = cms.Path(process.jetfilter)

process.e = cms.EndPath(process.out)
