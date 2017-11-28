#!/bin/bash
wget https://github.com/delphes/delphes/archive/3.4.2pre07.tar.gz
tar xfz 3.4.2pre07.tar.gz
mv delphes-* delphes
cd delphes
pwd
source DelphesEnv.sh
./configure
sed -i -e 's/c++0x/c++1y/g' Makefile
make -j4  
cd -

source env.sh
cd ntupler
make -j3
cd -
