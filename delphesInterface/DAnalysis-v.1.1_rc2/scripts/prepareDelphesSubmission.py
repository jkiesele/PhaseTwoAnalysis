#! /usr/bin/python

import sys
import re
import os
import argparse

parser = argparse.ArgumentParser(description="Prepare necessary files for submission of Delphes simulations from a LHE file")
parser.add_argument("--inputLHE",help="Input LHE file name")
parser.add_argument("--outputDir",help="Name of the ouput directoy in eos where the files will be stored")
parser.add_argument("--submissionDir",default="submissionDir",help="Name of the submission area in the working directory where all submission files should be created.")
parser.add_argument("--nTotalEvents",default=100,help="Total number of events to simulate",type=int)
parser.add_argument("--nJobs",default=1,help="Number of jobs to be submitted to lxbatch",type=int)
parser.add_argument("--scenario",default='0PU',choices=['0PU', '200PU'],help="Scenario for the simulation. Choices: '0pu' or '200PU'")
args = parser.parse_args()

cmsswbase=os.environ.get('CMSSW_BASE')
delphespath=os.environ.get('DELPHES_PATH')
location=os.getcwd();
os.mkdir(args.submissionDir)
os.mkdir(args.outputDir)

######################################
# check number of events in LHE file #
######################################

fin = ""
try:
  fin = open(args.inputLHE)
except:
  print("Error: Input file: %s could not be opened, exiting." % args.inputFile)
  sys.exit(1)

if args.nTotalEvents%args.nJobs != 0:
  print("Error: %s is not an exact divider of %s, exiting." % (args.nJobs,args.nTotalEvents))
  sys.exit(1)

eventNum = 0
init = False
inFooter = False
footLines = []
for line in fin:
  if re.match(r"[^#]*</LesHouchesEvents>",line):
    inFooter = True
    footLines.append(line+"\n")
  elif inFooter:
    footLines.append(line+"\n")
  elif init:
    if re.match(r"[^#]*</event>",line):
      eventNum += 1
  elif re.match(r"[^#]*</init>",line):
    init = True

eventsTotal = eventNum
print "Total number of events in the input LHE file: %i" % eventsTotal

if args.nTotalEvents > eventsTotal:
  print("Error: There are less events in the LHE than you try to simulate, exiting.")
  sys.exit(1)

##################
# Split LHE file # 
##################

files = []
maxEventsFile = []
for i in range(args.nJobs):
  tmp = open(args.submissionDir+"/splitLHE"+str(i)+".lhe",'w')
  files.append(tmp)
  maxEventsFile.append(args.nTotalEvents/args.nJobs)
maxEventsFile[len(maxEventsFile)-1] += eventsTotal % args.nJobs

eventNum = 0
eventNumThisFile = 0
init = False
headLines = []
iFile = 0
fin.seek(0)
for line in fin:
  if init:
    files[iFile].write(line)
    if re.match(r"[^#]*</event>",line):
      eventNum += 1
      eventNumThisFile += 1
      if eventNumThisFile >= maxEventsFile[iFile]:
        files[iFile].writelines(footLines)
        iFile += 1
        eventNumThisFile = 0
        if iFile == args.nJobs:
          break
        files[iFile].writelines(headLines)
  elif re.match(r"[^#]*</init>",line):
    init = True
    headLines.append(line)
    files[iFile].writelines(headLines)
  else:
    headLines.append(line)

for f in files:
  f.close()

#####################
# Create cmnd files #
#####################

for i in range(args.nJobs):
  line = ""
  mycmnd = open(args.submissionDir+"/splitCmnd"+str(i)+".cmnd",'w')
  line = line + "! 1) Settings used in the main program." + "\n" 
  line = line + "Main:numberOfEvents = " + str(maxEventsFile[0]) + "! number of events to generate" + "\n"
  line = line + "Main:timesAllowErrors = 3          ! how many aborts before run stops" + "\n"  + "\n"
  line = line + "! 2) Settings related to output in init(), next() and stat()." + "\n"
  line = line + "Init:showChangedSettings = on      ! list changed settings" + "\n"
  line = line + "Init:showChangedParticleData = off ! list changed particle data" + "\n"
  line = line + "Next:numberCount = 100             ! print message every n events" + "\n"
  line = line + "Next:numberShowInfo = 1            ! print event information n times" + "\n"
  line = line + "Next:numberShowProcess = 1         ! print process record n times" + "\n"
  line = line + "Next:numberShowEvent = 0           ! print event record n times" + "\n" + "\n"
  line = line + "! 3) Set the input LHE file" + "\n"
  line = line + "Beams:frameType = 4" + "\n"
  line = line + "Beams:LHEF = "+args.submissionDir+"/splitLHE" + str(i) + ".lhe" + "\n"
  mycmnd.write(line)

#######################
# Create submit files #
#######################

myfinalfile = open("SubmitDelphes.sh","w")
lineDelphes = ""

for i in range(args.nJobs):
  line = ""
  mycsh = open(args.submissionDir+"/submitJob"+str(i)+".sh",'w')
  line = line + "#!/bin/sh" + "\n" + "\n"
  line = line + "cd " + cmsswbase + "\n"
  line = line + "eval `scramv1 runtime -sh`" + "\n"
  line = line + "cd " + delphespath +"\n" 
  line = line + "export PYTHONPATH=`pwd`/python:$PYTHONPATH\nexport LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH\n"
  line = line + "cd " + location +"\n"
  if args.scenario=="0PU":
     line = line + delphespath + "/DelphesPythia8 "+delphespath +"/cards/CMS_PhaseII/CMS_PhaseII_0PU.tcl " +args.submissionDir+"/splitCmnd" + str(i) + ".cmnd "  +args.outputDir+"/delphes_lhe_batch" + str(i) + ".root" + "\n"
  elif args.scenario=="200PU":
     line = line + delphespath + "/DelphesPythia8 "+delphespath +"/cards/CMS_PhaseII/CMS_PhaseII_200PU.tcl " +args.submissionDir+"/splitCmnd" + str(i) + ".cmnd " +args.outputDir+"/delphes_lhe_batch" + str(i) + ".root" + "\n"
  mycsh.write(line)
  lineDelphes = lineDelphes + "bsub -q 8nh -o /dev/null -e /dev/null < "+args.submissionDir+"/submitJob" + str(i) + ".sh" + "\n"

myfinalfile.write(lineDelphes)
os.system("chmod +x SubmitDelphes.sh")
