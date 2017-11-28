/*
 * sampleCollection.h
 *
 *  Created on: 8 Nov 2016
 *      Author: jkiesele
 */

#ifndef INTERFACE_SAMPLECOLLECTION_H_
#define INTERFACE_SAMPLECOLLECTION_H_

#include "histoCollection.h"
#include "TString.h"
#include <vector>

namespace d_ana{
class sampleCollection{
public:
	const histoCollection& getHistos(const TString&  legendname)const;

	void readFromFile(const TString& filename);

	std::vector<TString> listAllLegends()const;
private:

	std::vector<histoCollection> histosPerSample_;
};
}

#endif /* INTERFACE_SAMPLECOLLECTION_H_ */
