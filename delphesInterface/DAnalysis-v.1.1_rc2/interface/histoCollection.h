/*
 * histoCollection.h
 *
 *  Created on: 8 Nov 2016
 *      Author: jkiesele
 */

#ifndef INTERFACE_HISTOCOLLECTION_H_
#define INTERFACE_HISTOCOLLECTION_H_

#include "TH1.h"
#include "metaInfo.h"
#include "TString.h"
#include <vector>
#include "TH1.h"
#include <algorithm>

namespace d_ana{

class sampleCollection;

class histoCollection{
	friend class sampleCollection;
public:
	histoCollection(){}
	histoCollection(const histoCollection&);
	~histoCollection();
	histoCollection& operator=(const histoCollection&);

	const TH1* getHisto(const TString& histoname)const;
	TH1* cloneHisto(const TString& histoname)const;

	const int& getColor()const{return meta_.color;}
	const TString& getLegendName()const{return meta_.legendname;}
	const int& getLegendOrder()const{return meta_.legendorder;}

private:

	std::vector< std::pair<TString, TH1*> >histomap_;
	metaInfo meta_;

};


}


#endif /* INTERFACE_HISTOCOLLECTION_H_ */
