/*
 * ANALYSER_TEMPL.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef ANALYSER_TEMPL_H_
#define ANALYSER_TEMPL_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class ANALYSER_TEMPL: public d_ana::basicAnalyzer{
public:
	ANALYSER_TEMPL():d_ana::basicAnalyzer(){}
	~ANALYSER_TEMPL(){}


private:
	void analyze(size_t id);

	void postProcess();
};





#endif /* ANALYSER_TEMPL_H_ */
