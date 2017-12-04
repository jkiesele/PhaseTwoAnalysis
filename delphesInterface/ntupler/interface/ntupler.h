/*
 * ntupler.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef ntupler_H_
#define ntupler_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class ntupler: public d_ana::basicAnalyzer{
public:
	ntupler():d_ana::basicAnalyzer(){}
	~ntupler(){}


private:
	void analyze(size_t id);

	void postProcess();
};





#endif /* ntupler_H_ */
