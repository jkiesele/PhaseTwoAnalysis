#!/bin/bash



ananame=$1
if [[ $2 ]]
then
    echo "USAGE: ${0} <analyser name>"
    exit -1
fi

if [[ $ananame == "" ]]
then
    echo "USAGE: ${0} <analyser name>"
    exit -1
fi
if [[ $ananame == *['!'@#\$%^\&*()_+-1234567890/=]* ]]
then
    echo "The name should not contain white spaces, numers or special characters"
    exit -1
fi

if [[ $DANALYSISPATH == "" ]]
then
	echo "the DAnalysis environment must be set. To do this, source the env.(c)sh script in the same directory this script is located in"
    exit -1 
fi

src=$DANALYSISPATH/templates/

sed -e 's;ANALYSER_TEMPL;'${ananame}';g' < $src/ANALYSER_TEMPL.h > ##workdir##/interface/${ananame}.h
sed -e 's;ANALYSER_TEMPL;'${ananame}';g' < $src/ANALYSER_TEMPL.cpp > ##workdir##/src/${ananame}.cpp
sed -e 's;ANALYSER_TEMPL;'${ananame}';g' < $src/ANALYSER_TEMPL_EXEC.cpp > ##workdir##/bin/${ananame}.cpp

echo "Skelton analyser/skimmer for the DAnalysis framework created."
echo "The source file can be found in the \"src\" directory. "
echo "It includes comments how to adapt it to the specific needs of the analysis."
echo "To compile the analyser, run \"make -j\""
echo "The analyser must be provided with a configuration file as input."
echo "Please have a look at the file \"config/testConfig.txt\" that gives indications about its structure"
echo "Other standard configuration files for the analyser can be found in the same directory."