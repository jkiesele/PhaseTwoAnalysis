/*
 * helpers.h
 *
 *  Created on: 1 Jul 2016
 *      Author: jkiesele
 */

#ifndef CROSSSECTION_SUMMER16_HELPERS_H_
#define CROSSSECTION_SUMMER16_HELPERS_H_

#include <sstream>
#include "TString.h"

template<class t>
std::string toString(t in) {
	std::ostringstream s;
	s << in;
	std::string out = s.str();
	return out;
}


template<class t>
TString toTString(t in) {
	std::ostringstream s;
	s << in;
	TString out = s.str();
	return out;
}


#endif /* CROSSSECTION_SUMMER16_HELPERS_H_ */
