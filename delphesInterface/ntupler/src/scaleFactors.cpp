/*
 * scaleFactors.cpp
 *
 *  Created on: 22 Jan 2018
 *      Author: jkiesele
 */


#include "../interface/scaleFactors.h"
#include <iostream>
#include <stdexcept>
#include "TMath.h"
#include <cfloat>
#include <algorithm>


scaleFactors::scaleFactors():xmin(-1),xmax(1),ymin(-1),ymax(1){

}
void  scaleFactors::loadTH2D(const TH2D& in){
		sfs_=in;

		xmin=in.GetXaxis()->GetXmin()+FLT_EPSILON;
		xmax=in.GetXaxis()->GetXmax()-FLT_EPSILON;
		ymin=in.GetYaxis()->GetXmin()+FLT_EPSILON;
		ymax=in.GetYaxis()->GetXmax()-FLT_EPSILON;
}

void scaleFactors::loadTH2D(const TString& filepath,const TString& histname){
	TFile f(filepath,"READ");
	if(f.IsZombie())
		throw std::runtime_error("scaleFactors:loadTH2D: could not load histogram");
	TH2D* h=(TH2D*)f.Get(histname);
	if(h)
		sfs_ = *h;
	else
		throw std::runtime_error("scaleFactors:loadTH2D: could not load histogram");

	xmin=h->GetXaxis()->GetXmin()+FLT_EPSILON;
	xmax=h->GetXaxis()->GetXmax()-FLT_EPSILON;
	ymin=h->GetYaxis()->GetXmin()+FLT_EPSILON;
	ymax=h->GetYaxis()->GetXmax()-FLT_EPSILON;

	sourcepath=filepath;
}


const float scaleFactors::getSF( float eta_in,  float pt_in)const{
	if(sfs_.GetNbinsX()<2)return 1;

	double eta=std::max(std::min(eta_in,xmax),xmin);
	double pt =std::max(std::min(pt_in, ymax),ymin);

	int binx=sfs_.GetXaxis()->FindFixBin(eta);
	int biny=sfs_.GetYaxis()->FindFixBin(pt);

	float SF=sfs_.GetBinContent(binx,biny);

	return (float)sfs_.GetBinContent(binx,biny);
}
