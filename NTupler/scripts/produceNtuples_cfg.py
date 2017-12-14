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
options.register('updateJEC', '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Name of the SQLite file (with path and extension) used to update the jet collection to the latest JEC and the era of the new JEC"
                )
options.register('noPU', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "noPU events: use ak4PFJets instead of Puppi"
                 )
options.parseArguments()

if len(options.updateJEC)==0:
    standardjec='/afs/cern.ch/work/v/vmilosev/HGCal/20171214/CMSSW_9_3_2/src/PhaseTwoAnalysis/NTupler/data/PhaseIIFall17_V3_MC.db'
    standardjec_tag='PhaseIIFall17_V3_MC'
    options.updateJEC=[standardjec,standardjec_tag]
    

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
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
)

# Input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) ) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/VBF_HToInvisible_M125_14TeV_powheg_pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/30000/145F80E4-A9BC-E711-9660-0242AC1C0502.root'
    ),
    secondaryFileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBF_HToInvisible_M125_14TeV_powheg_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/30000/6079EAF8-09B9-E711-87FA-0242AC110005.root'
        )
)

if (options.inputFormat.lower() == "reco"):
    process.source.fileNames = cms.untracked.vstring(*(
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBF_HToInvisible_M125_14TeV_powheg_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/30000/6079EAF8-09B9-E711-87FA-0242AC110005.root',
    ))

# HGCAL EGamma ID
if (options.inputFormat.lower() == "reco"):
    process.load("RecoEgamma.Phase2InterimID.phase2EgammaRECO_cff")
else:
    process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

# Get new JEC from an SQLite file rather than a GT
if options.updateJEC:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string('sqlite_file:'+options.updateJEC[0]),
                               toGet =  cms.VPSet(
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_"+options.updateJEC[1]+"_AK4PFPuppi"),
                     label = cms.untracked.string("AK4PFPuppi")),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_"+options.updateJEC[1]+"_AK4PF"),
                     label = cms.untracked.string("AK4PF"))
            )
                               )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

process.source.inputCommands = cms.untracked.vstring("keep *")

# Pre-skim weight counter
process.weightCounter = cms.EDAnalyzer('WeightCounter')

# Skim filter
muonLabel = "slimmedMuons"
elecLabel = "phase2Electrons"
photLabel = "phase2Photons"
if (options.inputFormat.lower() == "reco"):
    process.phoForYield = cms.EDProducer("CandViewMerger",
        src = cms.VInputTag(
            cms.InputTag("gedPhotons"),
            cms.InputTag("photonsFromMultiCl")
            )
        )
    photLabel = "phoForYield"
    process.eleForYield = cms.EDProducer("CandViewMerger",
        src = cms.VInputTag(
            cms.InputTag("gedGsfElectrons"),
            cms.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl")
            )
        )
    elecLabel = "eleForYield"
if options.updateJEC:
    if options.noPU:
        jetLabel = "updatedPatJetsUpdatedJECAK4PF"
    else:
        jetLabel = "updatedPatJetsUpdatedJECAK4PFPuppi"
else:    
    if options.noPU:
        jetLabel = "slimmedJets"
    else:
        jetLabel = "slimmedJetsPuppi"

if (options.inputFormat.lower() == "reco"):
    muonLabel = "muons"
    if options.updateJEC:
        if options.noPU:
            jetLabel = "ak4JetsL2L3"
        else:
            jetLabel = "ak4PUPPIJetsL1FastL2L3"
    else:
        if options.noPU:
            jetLabel = "ak4PFJets"
        else:
            jetLabel = "ak4PUPPIJets"


process.selectedMuons = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag(muonLabel),
                                     cut = cms.string("pt>10 && abs(eta)<3")
                                     )
process.selectedElectrons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag(elecLabel),
                                         cut = cms.string("pt>10 && abs(eta)<3")
                                         )
process.selectedPhotons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag(photLabel),
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
                                 minNumber = cms.uint32(0)
                                 )
process.countPhotons = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedPhotons"),
                                 minNumber = cms.uint32(0)
                                 )
process.countJets = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedJets"),
                                 minNumber = cms.uint32(2)
                                 )
process.preYieldFilter = cms.Sequence(process.selectedMuons+process.selectedElectrons+process.allLeps+process.countLeps+process.selectedPhotons+process.countPhotons+process.selectedJets+process.countJets)


# run Puppi 
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PhotonPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')
from CommonTools.PileupAlgos.PhotonPuppi_cff        import setupPuppiPhoton
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies(process)
process.particleFlowNoLep = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
process.puppiNoLep = process.puppi.clone(candName = cms.InputTag('particleFlowNoLep'))
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")

# recluster jets
process.load('RecoJets/Configuration/RecoPFJets_cff')
process.ak4PUPPIJets  = process.ak4PFJets.clone(rParam=0.4, src = cms.InputTag('puppi'))

# recompute MET
process.load('RecoMET.METProducers.PFMET_cfi')
process.puppiMet = process.pfMet.clone()
process.puppiMet.src = cms.InputTag('puppi')

# analysis
moduleName = "MiniFromPat"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "MiniFromReco"
process.ntuple = cms.EDAnalyzer(moduleName)
process.load("PhaseTwoAnalysis.NTupler."+moduleName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    if options.noPU:
        process.ntuple.pfCandsNoLep = "puppiNoLep"
        process.ntuple.met = "puppiMet"
    else:
        process.ntuple.pfCandsNoLep = "pfNoLep"
        process.ntuple.met = "pfMet"

    if options.updateJEC:
        # This will load several ESProducers and EDProducers which make the corrected jet collections
        # In this case the collection will be called ak4PUPPIJetsL1FastL2L3
        process.load('PhaseTwoAnalysis.Jets.JetCorrection_cff')
        if options.noPU:
            process.ntuple.jets = "ak4PFJetsL2L3"
        else:
            process.ntuple.jets = "ak4PUPPIJetsL1FastL2L3"
    else:
        # This simply switches the default AK4PFJetsCHS collection to the ak4PUPPIJets collection now that it has been produced
        if options.noPU:
            process.ntuple.jets = "ak4PFJets"

        else:
            process.ntuple.jets = "ak4PUPPIJets"

else:
    if options.noPU:
        process.ntuple.pileup = 0
        process.ntuple.jets = "slimmedJets"
        process.ntuple.mets = "slimmedMETs"
    else:
        process.ntuple.pileup = 200
        process.ntuple.jets = "slimmedJetsPuppi"
        process.ntuple.mets = "slimmedMETsPuppi"
    if options.updateJEC:
        # The updateJetCollection function will uncorred the jets from MiniAOD and then recorrect them using the current
        #  set of JEC in the event setup
        # The new name of the updated jet collection becomes updatedPatJetsUpdatedJECAK4PFPuppi
        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
        if options.noPU:
            updateJetCollection(process,
                                jetSource = cms.InputTag('slimmedJets'),
                                postfix = 'UpdatedJECAK4PF',
                                jetCorrections = ('AK4PF', ['L1FastJet','L2Relative','L3Absolute'], 'None')
                                )
            process.ntuple.jets = "updatedPatJetsUpdatedJECAK4PF"
        else:
            updateJetCollection(process,
                                jetSource = cms.InputTag('slimmedJetsPuppi'),
                                postfix = 'UpdatedJECAK4PFPuppi',
                                jetCorrections = ('AK4PFPuppi', ['L1FastJet','L2Relative','L3Absolute'], 'None')
                                )
            process.ntuple.jets = "updatedPatJetsUpdatedJECAK4PFPuppi"


# output
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )

# run
if (options.inputFormat.lower() == "reco"):
    if options.noPU:
        process.puSequence = cms.Sequence(process.primaryVertexAssociation * process.pfNoLep * process.offlineSlimmedPrimaryVertices * process.packedPFCandidates * process.muonIsolation * process.muonIsolationNoLep * process.ak4Jets * process.pfMet)
    else:
        process.puSequence = cms.Sequence(process.primaryVertexAssociation * process.pfNoLepPUPPI * process.puppi * process.particleFlowNoLep * process.puppiNoLep * process.offlineSlimmedPrimaryVertices * process.packedPFCandidates * process.muonIsolationPUPPI * process.muonIsolationPUPPINoLep * process.ak4PUPPIJets * process.puppiMet)

if options.skim:
    if (options.inputFormat.lower() == "reco"):
        if options.updateJEC:
            if options.noPU:
                process.p = cms.Path(process.weightCounter * process.phase2Egamma * process.puSequence * process.ak4PFL2L3CorrectorChain * process.ak4PFJetsL2L3 * process.preYieldFilter * process.ntuple)
            else:
                process.p = cms.Path(process.weightCounter * process.phase2Egamma * process.puSequence * process.ak4PFPuppiL1FastL2L3CorrectorChain * process.ak4PUPPIJetsL1FastL2L3 * process.preYieldFilter * process.ntuple)

        else:
            process.p = cms.Path(process.weightCounter * process.phase2Egamma * process.puSequence * process.preYieldFilter * process.ntuple)
    else:
        if options.updateJEC:
            if options.noPU:
                process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.patJetCorrFactorsUpdatedJECAK4PF * process.updatedPatJetsUpdatedJECAK4PF *process.preYieldFilter* process.ntuple)
            else:
                process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.patJetCorrFactorsUpdatedJECAK4PFPuppi * process.updatedPatJetsUpdatedJECAK4PFPuppi *process.preYieldFilter* process.ntuple)
        else:
            process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.preYieldFilter*process.ntuple)
else:
    if (options.inputFormat.lower() == "reco"):
        if options.updateJEC:
            if options.noPU:
                process.p = cms.Path(process.weightCounter * process.puSequence * process.ak4PFL2L3CorrectorChain * process.ak4JetsL2L3 * process.phase2Egamma * process.ntuple)
            else:
                process.p = cms.Path(process.weightCounter * process.puSequence * process.ak4PFPuppiL1FastL2L3CorrectorChain * process.ak4PUPPIJetsL1FastL2L3 * process.phase2Egamma * process.ntuple)

        else:
            process.p = cms.Path(process.weightCounter * process.puSequence * process.phase2Egamma * process.ntuple)
    else:
        if options.updateJEC:
            if options.noPU:
                process.p = cms.Path(process.weightCounter*process.patJetCorrFactorsUpdatedJECAK4PF * process.updatedPatJetsUpdatedJECAK4PF * process.phase2Egamma * process.ntuple)
            else:
                process.p = cms.Path(process.weightCounter*process.patJetCorrFactorsUpdatedJECAK4PFPuppi * process.updatedPatJetsUpdatedJECAK4PFPuppi * process.phase2Egamma * process.ntuple)
	else:    
            process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.ntuple)
