from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'VBFHToInvisiblenoPU'
#config.General.requestName = 'VBFHToInvisible200PU'
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'scripts/produceNtuples_cfg.py'
#noPU
config.JobType.pyCfgParams= ['skim=False','outFilename=MiniEvents.root','inputFormat=PAT','noPU=True'] 
#PU
#config.JobType.pyCfgParams= ['skim=False','inputFormat=PAT','outFileName=MiniEvents.root']
# Uncomment the following line when running on PAT events
config.JobType.outputFiles = ['MiniEvents.root']

config.section_("Data")
#noPU
config.Data.inputDataset = '/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
#PU
#config.Data.inputDataset = '/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
# Uncomment to run on a fraction of the dataset
#config.Data.totalUnits = 5
config.Data.outLFNDirBase = '/store/user/vmilosev/'  
config.Data.publication = False

config.Data.useParent = True # need to run on GEN-SIM-RECO to apply photon ID

config.section_("Site")
config.Site.storageSite = 'T2_UK_London_IC'
config.Site.ignoreGlobalBlacklist = True
