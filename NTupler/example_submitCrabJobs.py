#
#
# example configuration for multicrab submission. Needs to be adapted to analysis needs
# (thanks to Sandhya for providing this)
#
#
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/afs/cern.ch/work/s/sandhya/Physics/Upgrade/VBF/JanFullSim/CMSSW_9_3_2/src/PhaseTwoAnalysis/NTupler/scripts/produceNtuples_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['MiniEvents.root']
#### if providing pyCfgParams, please replace submit(config) in multicrab by
##### p = Process(target=submit, args=(config,))
##### p.start()
##### p.join()
#config.JobType.inputFiles = ['PhaseIIFall17_V3_MC.db']
#config.JobType.pyCfgParams= ['updateJEC=PhaseIIFall17_V3_MC.db','updateJEC=PhaseIIFall17_V3_MC']

config.JobType.maxMemoryMB = 4000
config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/sandhya/HGCal/'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.secondaryInputDataset = ''

config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.storageSite = 'T2_CH_CERN'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

            
    config.Data.inputDataset = '/WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8_v2'
    config.Data.secondaryInputDataset = '/WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
     
    config.Data.inputDataset = '/WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WWJJ_SS_WToLNu_EWK_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v3/MINIAODSIM'
    config.General.requestName = 'TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8_v3'
    config.Data.secondaryInputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v2/MINIAODSIM'
    config.General.requestName = 'TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8_ext1'
    config.Data.secondaryInputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext2-v1/MINIAODSIM'
    config.General.requestName = 'TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8_ext2'
    config.Data.secondaryInputDataset = '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8_v2'
    config.Data.secondaryInputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
        
    config.Data.inputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17MiniAOD-PU200FEVT_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8_ext1'
    config.Data.secondaryInputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-PU200FEVT_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
######no extension/v2 for this dataset
    config.Data.inputDataset = '/DYToLL-M-50_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
    config.General.requestName = 'DYToLL-M-50_3J_14TeV-madgraphMLM-pythia8_v1'
    config.Data.secondaryInputDataset = '/DYToLL-M-50_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
        
    config.Data.inputDataset = '/WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
        
    config.Data.inputDataset = '/WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_ext1' 
    config.Data.secondaryInputDataset = '/WGToLNuG_PtG-40_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
############no extension for this dataset
    config.Data.inputDataset = '/WZJJTo3LNuJJ_MLL-4_EWK_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZJJTo3LNuJJ_MLL-4_EWK_14TeV-madgraph-pythia8_v2'
    config.Data.secondaryInputDataset = '/WZJJTo3LNuJJ_MLL-4_EWK_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
############no extension for this dataset
    config.Data.inputDataset = '/WZTo3LNu_0Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZTo3LNu_0Jets_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WZTo3LNu_0Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
        
    config.Data.inputDataset = '/WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8_v2' 
    config.Data.secondaryInputDataset = '/WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8_ext1' 
    config.Data.secondaryInputDataset = '/WZTo3LNu_1Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8_v2' 
    config.Data.secondaryInputDataset = '/WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8_ext1' 
    config.Data.secondaryInputDataset = '/WZTo3LNu_2Jets_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
#############no extension for this dataset
    config.Data.inputDataset = '/WZJJTo2L2QJJ_MLL-4_EWK_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZJJTo2L2QJJ_MLL-4_EWK_14TeV-madgraph-pythia8_v2'
    config.Data.secondaryInputDataset = '/WZJJTo2L2QJJ_MLL-4_EWK_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZZTo2Q2L_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ZZTo2Q2L_14TeV_powheg_pythia8_v2'
    config.Data.secondaryInputDataset = '/ZZTo2Q2L_14TeV_powheg_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZZTo2Q2L_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM' 
    config.General.requestName = 'ZZTo2Q2L_14TeV_powheg_pythia8_ext1'
    config.Data.secondaryInputDataset = '/ZZTo2Q2L_14TeV_powheg_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZZTo4L_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ZZTo4L_14TeV_powheg_pythia8_v2'
    config.Data.secondaryInputDataset = '/ZZTo4L_14TeV_powheg_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZZTo4L_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'ZZTo4L_14TeV_powheg_pythia8_ext1'
    config.Data.secondaryInputDataset = '/ZZTo4L_14TeV_powheg_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8_v2'
    config.Data.secondaryInputDataset = '/ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8_ext1'
    config.Data.secondaryInputDataset = '/ZGTo2LG_TuneCUETP8M1_14TeV-amcatnloFXFX-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8_v2' 
    config.Data.secondaryInputDataset = '/TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8_ext1' 
    config.Data.secondaryInputDataset = '/TTGamma_SingleLeptFromTbar_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8_v2'
    config.Data.secondaryInputDataset = '/TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8_ext1'
    config.Data.secondaryInputDataset = '/TTGamma_SingleLeptFromT_TuneCUETP8M1_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8_v2'
    config.Data.secondaryInputDataset = '/TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8_ext1' 
    config.Data.secondaryInputDataset = '/TTWJetsToLNu_TuneCUETP8M1_14TeV-amcatnloFXFX-madspin-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    

    config.Data.inputDataset = '/TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/TTZJets_TuneCUETP8M1_14TeV_madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
########no extension found for this dataset
    config.Data.inputDataset = '/WWW_4F_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WWW_4F_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/WWW_4F_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset =  '/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v3/MINIAODSIM'
    config.General.requestName = 'WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v3/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
#########no extension found for this dataset
    
    config.Data.inputDataset = '/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
#########no extension found for this dataset
    
    config.Data.inputDataset = '/WWG_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WWG_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/WWG_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
#########no extension found for this dataset
    
    config.Data.inputDataset = '/WZG_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WZG_TuneCUETP8M1_14TeV-amcatnlo-pythia8_v2'
    config.Data.secondaryInputDataset = '/WZG_TuneCUETP8M1_14TeV-amcatnlo-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ST_tW_DR_14TeV_top_incl-powheg-pythia8_v2'
    config.Data.secondaryInputDataset = '/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'ST_tW_DR_14TeV_top_incl-powheg-pythia8_ext1'
    config.Data.secondaryInputDataset = '/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'ST_tW_DR_14TeV_antitop_incl-powheg-pythia8_v2'
    config.Data.secondaryInputDataset = '/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'ST_tW_DR_14TeV_antitop_incl-powheg-pythia8_ext1'
    config.Data.secondaryInputDataset = '/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'
    submit(config)
    
    config.Data.inputDataset = '/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WToLNu_0J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    submit(config)

    config.Data.inputDataset = '/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WToLNu_0J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WToLNu_1J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WToLNu_1J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WToLNu_2J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WToLNu_2J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'WToLNu_3J_14TeV-madgraphMLM-pythia8_v2'
    config.Data.secondaryInputDataset = '/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'  
    submit(config)

    config.Data.inputDataset = '/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'
    config.General.requestName = 'WToLNu_3J_14TeV-madgraphMLM-pythia8_ext1'
    config.Data.secondaryInputDataset = '/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2_ext1-v1/GEN-SIM-RECO'  
    submit(config)
