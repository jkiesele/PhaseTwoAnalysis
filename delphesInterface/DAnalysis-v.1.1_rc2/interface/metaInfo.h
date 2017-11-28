/*
 * metaInformation.h
 *
 *  Created on: 22 Aug 2016
 *      Author: jkiesele
 */

#ifndef INTERFACE_METAINFO_H_
#define INTERFACE_METAINFO_H_

#include "TString.h"
class TDirectory;

namespace d_ana{

class metaInfo{
public:

	metaInfo():
	color(0),
	legendorder(0),
	norm(0)
	{}
	~metaInfo(){}

	void Write()const;

	void extractFrom(const TString& str);

	void extractFrom( TDirectory* dir);

	TString legendname;
	int color;
	int legendorder;
	///only for book-keeping!
	double norm;

private:
	static const TString separator_;

};



}


#endif /* INTERFACE_METAINFO_H_ */
