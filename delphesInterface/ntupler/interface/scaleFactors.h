/*
 * scaleFactors.h
 *
 *  Created on: 5 Dec 2017
 *      Author: jkiesele
 */

#ifndef PHASETWOANALYSIS_DELPHESINTERFACE_NTUPLER_INTERFACE_SCALEFACTORS_H_
#define PHASETWOANALYSIS_DELPHESINTERFACE_NTUPLER_INTERFACE_SCALEFACTORS_H_

#include "TH2.h"
#include "TFile.h"
#include "TString.h"

class scaleFactors{
public:
	scaleFactors();
	void loadTH2D(const TH2D& in);
	void loadTH2D(const TString& filepath,const TString& histname);

	const float getSF( float eta_in,  float pt_in)const;

private:
	TH2D sfs_;
	float xmin,xmax,ymin,ymax;
	TString sourcepath;
};


#endif /* PHASETWOANALYSIS_DELPHESINTERFACE_NTUPLER_INTERFACE_SCALEFACTORS_H_ */
