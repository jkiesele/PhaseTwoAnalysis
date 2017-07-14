import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *

ak4PFPuppiL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK4PFPuppi'),
    srcRho      = cms.InputTag('fixedGridRhoFastjetAll'),
    )
ak4PFPuppiL2Relative  = ak4PFL2Relative.clone( algorithm = 'AK4PFPuppi' )
ak4PFPuppiL3Absolute  = ak4PFL3Absolute.clone( algorithm = 'AK4PFPuppi' )
ak4PFPuppiL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak4PFPuppiL2Relative','ak4PFPuppiL3Absolute')
    )
ak4PFPuppiL1FastL2L3 = ak4PFPuppiL2L3.clone()
ak4PFPuppiL1FastL2L3.correctors.insert(0,'ak4PFPuppiL1Fastjet')

ak4PFPuppiL1FastjetCorrector = cms.EDProducer(
    'L1FastjetCorrectorProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK4PFPuppi'),
    srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )
ak4PFPuppiL2RelativeCorrector = cms.EDProducer(
    'LXXXCorrectorProducer',
    level     = cms.string('L2Relative'),
    algorithm = cms.string('AK4PFPuppi')
    )
ak4PFPuppiL3AbsoluteCorrector = cms.EDProducer(
    'LXXXCorrectorProducer',
    level     = cms.string('L3Absolute'),
    algorithm = cms.string('AK4PFPuppi')
    )
ak4PFPuppiL2L3Corrector = cms.EDProducer(
    'ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak4PFPuppiL2RelativeCorrector','ak4PFPuppiL3AbsoluteCorrector')
    )
ak4PFPuppiL1FastL2L3Corrector = ak4PFPuppiL2L3Corrector.clone()
ak4PFPuppiL1FastL2L3Corrector.correctors.insert(0,'ak4PFPuppiL1FastjetCorrector')
ak4PFPuppiL1FastL2L3CorrectorChain = cms.Sequence(
    ak4PFPuppiL1FastjetCorrector * ak4PFPuppiL2RelativeCorrector * ak4PFPuppiL3AbsoluteCorrector * ak4PFPuppiL1FastL2L3Corrector
    )

ak4PUPPIJetsL1 = cms.EDProducer(
    'CorrectedPFJetProducer',
    src         = cms.InputTag('ak4PUPPIJets'),
    correctors  = cms.VInputTag('ak4PUPPIL1FastjetCorrector')
    )
ak4PUPPIJetsL2L3 = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag('ak4PUPPIJets'),
    correctors  = cms.VInputTag('ak4PFPuppiL2L3Corrector')
    )
ak4PUPPIJetsL1FastL2L3  = ak4PUPPIJetsL2L3.clone(src = 'ak4PUPPIJets', correctors = ['ak4PFPuppiL1FastL2L3Corrector'])