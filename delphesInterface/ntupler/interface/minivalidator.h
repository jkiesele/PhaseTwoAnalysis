/*
 * minivalidator.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef minivalidator_H_
#define minivalidator_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class minivalidator: public d_ana::basicAnalyzer{
public:
	minivalidator():d_ana::basicAnalyzer(){}
	~minivalidator(){}


private:
	void analyze(size_t id);

	void postProcess();
};





#endif /* minivalidator_H_ */
