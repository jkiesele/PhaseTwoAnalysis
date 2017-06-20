#!/usr/bin/env python

import sys, os, string, shutil, datetime
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--skim", dest="skim", type="string", default=False, help="skim")
parser.add_option("-f", "--filelist", dest="filelist", type="string", default=False, help="file list directory")
parser.add_option("-v", "--version", dest="version", type="string", default=False, help="version of the day")
parser.add_option("-d", "--descr", dest="descr", type="string", default=False, help="description of the selection") 
(options, args) = parser.parse_args()

if not options.version:
    parser.error("you must specify the version")

if not options.filelist or not os.path.isdir(options.filelist):
    parser.error("you must specify a file list directory")

date = datetime.datetime.now().strftime("%y%b%d")

outRoot = date
if not os.path.isdir(outRoot):
    os.mkdir(outRoot)
outRoot = os.path.join(outRoot, options.version)   
if not os.path.isdir(outRoot):
    os.mkdir(outRoot)

try:    
    os.mkdir(outRoot)
except OSError:    
    raise OSError("this version for this channel already exists")

if not os.path.isdir("LogMyAna"):
    os.mkdir("LogMyAna/")
logFile = os.path.join("LogMyAna", date)    
logFile += "_"+options.version     
if options.skim:
    logFile += "_skim"
if options.descr:
    logFile += "_"+options.descr
logFile += ".log"
log = open(logFile, 'w')

fileListDir = options.filelist 

os.system("make clean")
os.system("make")

copy = "cp MyAna.cc " + outRoot + "/"
os.system(copy)

filelist = os.popen("ls  "+fileListDir+"/*.list").readlines()

for aFile in filelist:
   
    outRoo1 = outRoot+os.sep+string.split(string.strip(aFile),"/")[-1][:-5]+".root"
    outLog1 = outRoot+os.sep+string.split(string.strip(aFile),"/")[-1][:-5]+".log"
    
    option  = ""
    if options.skim
        option  = option+" -skim"

    if aFile.count("TTJets") == 1:
        option = option+" -sig"    
        
    cmd = "./runMyAna -filelist "+string.strip(aFile)+" -out "+outRoo1+option+" >& "+outLog1+"\n" 
    
    log.write(cmd)
    os.system(cmd)
    
log.close()
sys.exit()
