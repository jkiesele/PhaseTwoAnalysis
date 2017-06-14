from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = <'requestName'>
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'
config.JobType.inputFiles = ['TMVAClassification_BDT.weights.xml']
config.JobType.outputFiles = ['histos.root']

config.section_("Data")
config.Data.inputDataset = <'inputDataset'>
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = <'dirBase'>  
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True
