#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)

parser.add_option('-q', '--queue' ,  dest='queue'  , help='batch queue' , default='1nd')
parser.add_option('-d', '--outdir' ,  dest='outdir'  , help='output directory' , default='LooseVBFsel')
parser.add_option('-p', '--process' ,  dest='process'  , help='process short name' , default='')
parser.add_option('-j', '--dojes'   , action="store_true",  dest='dojes'     , help='do JES histograms')
(opt, args) = parser.parse_args()


workdir='/afs/cern.ch/work/a/amagnan/UPSGAna/'

processList=['DYToLL-M-50_0J','DYToLL-M-50_1J','DYToLL-M-50_2J','DYToLL-M-50_3J','EWKWMinus2Jets_WToLNu_M-50','EWKWPlus2Jets_WToLNu_M-50','EWKZ2Jets_ZToLL_M-50','EWKZ2Jets_ZToNuNu','QCD_Mdijet-1000toInf','ST_tW_DR_14TeV_top_ext1','ST_tW_DR_14TeV_top','ST_tW_DR','ST_tch_14TeV_antitop','ST_tch_14TeV_top_ext1','ST_tch_14TeV_top','TT_TuneCUETP8M2T4','VBFH','WToLNu_0J','WToLNu_1J','WToLNu_2J','WToLNu_3J','ZJetsToNuNu_HT-100To200','ZJetsToNuNu_HT-1200To2500','ZJetsToNuNu_HT-200To400','ZJetsToNuNu_HT-400To600','ZJetsToNuNu_HT-600To800']

if len(opt.process)>0:
    processList=[opt.process]

outDir='%s%s/'%(workdir,opt.outdir)
os.system('mkdir -p %s'%outDir)

for myproc in processList :

    #scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    #scriptFile.write('#!/bin/bash\n')
    #scriptFile.write('cd %s/../../../'%(os.getcwd()))
    #scriptFile.write('cmsenv')
    if opt.dojes :
        os.system('./bin/simplePlots %s %s 1'%(outDir,myproc))
    else:
        os.system('./bin/simplePlots %s %s'%(outDir,myproc))
    
    
