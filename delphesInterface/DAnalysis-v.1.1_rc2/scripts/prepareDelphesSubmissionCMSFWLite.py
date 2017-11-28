#! /usr/bin/python

import sys
import re
import os
import argparse

parser = argparse.ArgumentParser(description="Prepare necessary files for submission of Delphes simulations from a ROOT file")
parser.add_argument("--outputDir",help="Name of the ouput directoy in eos where the files will be stored")
parser.add_argument("--inputFile",help="Name of the input root file")
parser.add_argument("--scenario",default='0PU',choices=['0PU', '200PU'],help="Scenario for the simulation. Choices: '0pu' or '200PU'")
args = parser.parse_args()

location=os.getcwd();
os.mkdir(args.outputDir)

#######################
# Create submit files #
#######################

myfinalfile = open("SubmitDelphes.sh","w")
lineDelphes = ""

line = ""
mycsh = open("submitJob.csh",'w')
line = line + "#!/bin/csh" + "\n" + "\n"
line = line + "cd " + location + "\n"
line = line + "cmsenv" + "\n"
if args.scenario=="0PU":
   line = line + "./DelphesCMSFWLite cards/CMS_PhaseII/CMS_PhaseII_0PU.tcl " + args.outputDir + "/delphes_cmsfwlite_batch.root " +  args.inputFile +"\n"
elif args.scenario=="200PU":
   line = line + "./DelphesCMSFWLite cards/CMS_PhaseII/CMS_PhaseII_0PU.tcl " + args.outputDir + "/delphes_cmsfwlite_batch.root " + args.inputFile +"\n"
mycsh.write(line)
lineDelphes = lineDelphes + "bsub -q 8nh -o /dev/null -e /dev/null < submitJob.csh" + "\n"

myfinalfile.write(lineDelphes)

