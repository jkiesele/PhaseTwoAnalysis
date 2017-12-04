#!/bin/sh

ROOTFLAGS=`root-config --cflags`
ROOTLIBS=`root-config --libs`
LIBS="${ROOTLIBS} -L${DELPHES_PATH} -lDelphes -L${DANALYSISPATH} -lDAnalysis"

infile=$1

CPLUS_INCLUDE_PATH="-I${DELPHES_PATH} -I${DANALYSISPATH}"

g++ -g $ROOTFLAGS $CPLUS_INCLUDE_PATH -c -o $infile.o $infile
g++ -g -o $infile.exe  -Wall $LIBS $infile.o  
