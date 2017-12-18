from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'VBFHToInvisible_multicrab_test'
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'scripts/produceNtuples_cfg.py'
config.JobType.inputFiles  = ['data/PhaseIIFall17_V3_MC.db']
#noPU
#config.JobType.pyCfgParams= ['skim=False','outFilename=MiniEvents.root','inputFormat=PAT', 'updateJEC=PhaseIIFall17_V3_MC.db', 'updateJEC=PhaseIIFall17_V3_MC', 'noPU=True'] 
#PU
#config.JobType.pyCfgParams= ['skim=False','outFilename=MiniEvents.root','inputFormat=PAT', 'updateJEC=PhaseIIFall17_V3_MC.db', 'updateJEC=PhaseIIFall17_V3_MC', 'noPU=False']
config.JobType.maxMemoryMB = 4000
# Uncomment the following line when running on PAT events
config.JobType.outputFiles = ['MiniEvents.root']

config.section_("Data")
#noPU
#config.Data.inputDataset = '/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
#PU
#config.Data.inputDataset = '/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 4000
#TEST 10 EVENTS
config.Data.unitsPerJob = 10
# Uncomment to run on a fraction of the dataset
#comment out this for the full dataset
config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/vmilosev/'
#config.Data.outLFNDirBase = '/store/user/amagnan/'
config.Data.publication = False

config.Data.useParent = True # need to run on GEN-SIM-RECO to apply photon ID

config.section_("Site")
config.Site.storageSite = 'T2_UK_London_IC'
config.Site.ignoreGlobalBlacklist = True

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    
    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    
    def submit(inconfig):
        try:
            crabCommand('submit', config = inconfig)
        #            crabCommand('status')
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

    tasks=list()
    #Entries are in the format:
    #tasks.append((taskname,dataset name from das))
    
    # VBF_HToInvisible
    tasks.append(('DYJetsToLL_M-10to50_noPU','/DYJetsToLL_M-10to50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('DYJetsToLL_M-10to50_PU200ext1','/DYJetsToLL_M-10to50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('DYJetsToLL_M-10to50_PU200','/DYJetsToLL_M-10to50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_0J_noPU','/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_1J_PU200','/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_1J_noPU','/DYToLL-M-50_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_2J_PU200','/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_2J_PU200ext1','/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_2J_noPU','/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_3J_PU200','/DYToLL-M-50_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('DYToLL-M-50_3J_noPU','/DYToLL-M-50_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('EWKWMinus2Jets_WToLNu_PU200','/EWKWMinus2Jets_WToLNu_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('EWKWMinus2Jets_WToLNu_noPU','/EWKWMinus2Jets_WToLNu_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('EWKWPlus2Jets_WToLNu_M-50_PU200','/EWKWPlus2Jets_WToLNu_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('EWKWPlus2Jets_WToLNu_M-50_noPU','/EWKWPlus2Jets_WToLNu_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('EWKZ2Jets_ZToLL_M-50_PU200','/EWKZ2Jets_ZToLL_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('EWKZ2Jets_ZToLL_M-50_noPU','/EWKZ2Jets_ZToLL_M-50_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('EWKZ2Jets_ZToNuNu_PU200','/EWKZ2Jets_ZToNuNu_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('EWKZ2Jets_ZToNuNu_noPU','/EWKZ2Jets_ZToNuNu_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('QCD_Mdijet-1000toInf_PU200','/QCD_Mdijet-1000toInf_TuneCUETP8M1_14TeV-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('QCD_Mdijet-1000toInf_noPU','/QCD_Mdijet-1000toInf_TuneCUETP8M1_14TeV-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('TT_TuneCUETP8M2T4_PU200','/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v3/MINIAODSIM'))
    tasks.append(('TT_TuneCUETP8M2T4_noPU','/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('VBF_HToInvisible_M125_PU200','/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('VBF_HToInvisible_M125_noPU','/VBF_HToInvisible_M125_14TeV_powheg_pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('WToLNu_0J_PU200','/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('WToLNu_0J_noPU','/WToLNu_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('WToLNu_1J_PU200','/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('WToLNu_1J_PU200ext1','/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('WToLNu_1J_noPU','/WToLNu_1J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('WToLNu_2J_PU200','/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('WToLNu_2J_noPU','/WToLNu_2J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('WToLNu_3J_PU200','/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('WToLNu_3J_noPU','/WToLNu_3J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-100To200_PU200','/ZJetsToNuNu_HT-100To200_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-100To200_noPU','/ZJetsToNuNu_HT-100To200_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-1200To2500_PU200','/ZJetsToNuNu_HT-1200To2500_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-1200To2500_noPU','/ZJetsToNuNu_HT-1200To2500_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-200To400_PU200','/ZJetsToNuNu_HT-200To400_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-200To400_noPU','/ZJetsToNuNu_HT-200To400_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-400To600_PU200','/ZJetsToNuNu_HT-400To600_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-400To600_noPU','/ZJetsToNuNu_HT-400To600_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-600To800_PU200','/ZJetsToNuNu_HT-600To800_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-600To800_noPU','/ZJetsToNuNu_HT-600To800_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-800To1200_PU200','/ZJetsToNuNu_HT-800To1200_14TeV-madgraph/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ZJetsToNuNu_HT-800To1200_noPU','/ZJetsToNuNu_HT-800To1200_14TeV-madgraph/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_antitop_incl_noPU','/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_antitop_incl_PU200','/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_antitop_incl_PU200ext1','/ST_tW_DR_14TeV_antitop_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_top_incl_noPU','/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_top_incl_PU200','/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_top_incl_noPUext1','/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('ST_tW_DR_14TeV_top_incl_PU200ext1','/ST_tW_DR_14TeV_top_incl-powheg-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_antitop_incl_noPU','/ST_tch_14TeV_antitop_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_antitop_incl_PU200_old_v1','/ST_tch_14TeV_antitop_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_antitop_incl_PU200','/ST_tch_14TeV_antitop_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_top_incl_noPU','/ST_tch_14TeV_top_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_top_incl_PU200','/ST_tch_14TeV_top_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_top_incl_noPUext1','/ST_tch_14TeV_top_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))
    tasks.append(('ST_tch_14TeV_top_incl_PU200ext1','/ST_tch_14TeV_top_incl-powheg-pythia8-madspin/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2_ext1-v1/MINIAODSIM'))

    for task in tasks:
        print task[0]
        config.General.requestName = task[0]
        config.Data.inputDataset = task[1]
       
        if "PU200" in task[1]:
            config.JobType.pyCfgParams= ['skim=False','outFilename=MiniEvents.root','inputFormat=PAT', 'updateJEC=PhaseIIFall17_V3_MC.db', 'updateJEC=PhaseIIFall17_V3_MC', 'noPU=False']
        else :
            config.JobType.pyCfgParams= ['skim=False','outFilename=MiniEvents.root','inputFormat=PAT', 'updateJEC=PhaseIIFall17_V3_MC.db', 'updateJEC=PhaseIIFall17_V3_MC', 'noPU=True']


        submit(config)


