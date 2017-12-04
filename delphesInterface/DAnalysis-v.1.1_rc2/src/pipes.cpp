/*
 * pipes.cpp
 *
 *  Created on: 1 Sep 2016
 *      Author: jkiesele
 */


#include "../interface/pipes.h"

#include <stdexcept>

size_t IPCPipeBase::openpipes=0;

IPCPipeBase::IPCPipeBase(){
	if(openpipes>511)
		throw std::runtime_error("IPCPipeBase: Too many open pipes!");
	openpipes++;
}

IPCPipeBase::~IPCPipeBase(){
	openpipes--;
}

