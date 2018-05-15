from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = <'requestName'>
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'scripts/produceNtuples_cfg.py'
config.JobType.pyCfgParams= ['skim=False','inputFormat=PAT','outFileName=MiniEvents.root','pileup=200']
# Uncomment the following line when running on PAT events
config.JobType.outputFiles = ['MiniEvents.root']

config.section_("Data")
config.Data.inputDataset = <'inputDataset'>
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
# Uncomment to run on a fraction of the dataset
#config.Data.totalUnits = 5
config.Data.outLFNDirBase = <'dirBase'>  
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True
