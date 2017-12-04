#!/bin/bash
PPOLDDIR=`pwd`

echo $PPOLDDIR
wget https://github.com/delphes/delphes/archive/3.4.2pre07.tar.gz
tar xfz 3.4.2pre07.tar.gz
mv delphes-* delphes
cd delphes
pwd
source DelphesEnv.sh
./configure
sed -i -e 's/c++0x/c++1y/g' Makefile
make -j4  
cd $PPOLDDIR
source env.sh

cd $DANALYSISPATH
pwd
make -j3
cd $PPOLDDIR/ntupler
pwd
make -j3
cd -
