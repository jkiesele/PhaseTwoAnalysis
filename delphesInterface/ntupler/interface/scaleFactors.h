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
#include <stdexcept>

class scaleFactors{
public:
	void loadTH2D(const TH2D& in){
		sfs_=in;
	}
	void loadTH2D(const TString& filepath,const TString& histname){
		TFile f(filepath,"READ");
		if(f.IsZombie())
			throw std::runtime_error("scaleFactors:loadTH2D: could not load histogram");
		TH2D* h=(TH2D*)f.Get(histname);
		if(h)
			sfs_ = *h;
		else
			throw std::runtime_error("scaleFactors:loadTH2D: could not load histogram");
	}

	inline const float getSF(const double& eta, const double& pt)const{
		if(sfs_.GetNbinsX()<1)return 1;
		int bin=sfs_.FindFixBin(eta,pt);
		int binx,biny,binz;
		sfs_.GetBinXYZ(bin,  binx,  biny,  binz);

		//if under/overflow go back to last usable bin
		if(binx==0)binx++;
		if(biny==0)biny++;
		if(binx>=sfs_.GetNbinsX()+1)binx--;
		if(biny>=sfs_.GetNbinsX()+1)biny--;

		return sfs_.GetBinContent(binx,biny);
	}

private:
	TH2D sfs_;

};


#endif /* PHASETWOANALYSIS_DELPHESINTERFACE_NTUPLER_INTERFACE_SCALEFACTORS_H_ */
