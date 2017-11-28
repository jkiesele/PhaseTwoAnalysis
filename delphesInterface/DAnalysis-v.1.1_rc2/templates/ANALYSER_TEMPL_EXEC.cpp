/*
 * ANALYSER_TEMPL_EXEC.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */


#include <iostream>
#include "interface/ANALYSER_TEMPL.h"

int main(int argc, char* argv[]){

	if(argc!=2){
		std::cout << "ANALYSER_TEMPL: need exactly one input file" <<std::endl;
		exit (-1);
	}


	std::string inputfile=argv[1];

	ANALYSER_TEMPL analyser;


	analyser.readConfigFile(inputfile);

	analyser.start();

	return 1;
}
