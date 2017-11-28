/*
 * histoCollection.cpp
 *
 *  Created on: 8 Nov 2016
 *      Author: jkiesele
 */




#include "../interface/histoCollection.h"

namespace d_ana{


histoCollection::histoCollection(const histoCollection& rhs){
	for(size_t i=0;i<rhs.histomap_.size();i++){
		histomap_.push_back(std::pair<TString, TH1*> (rhs.histomap_.at(i).first, (TH1*)rhs.histomap_.at(i).second->Clone()));
	}
	meta_=rhs.meta_;
}

histoCollection::~histoCollection(){
	for(size_t i=0;i<histomap_.size();i++){
		if(histomap_.at(i).second)
			delete histomap_.at(i).second;
		histomap_.at(i).second=0;
	}
}

histoCollection& histoCollection::operator=(const histoCollection& rhs){
	if(&rhs==this) return *this;

	for(size_t i=0;i<histomap_.size();i++){
		if(histomap_.at(i).second)
			delete histomap_.at(i).second;
	}
	histomap_.clear();
	for(size_t i=0;i<rhs.histomap_.size();i++){
		histomap_.push_back(std::pair<TString, TH1*> (rhs.histomap_.at(i).first, (TH1*)rhs.histomap_.at(i).second->Clone()));
	}

	meta_=rhs.meta_;

	return *this;
}

const TH1* histoCollection::getHisto(const TString& histoname)const{

	//std::vector< std::pair<TString, TH1*> >histomap_;

	//by construction there are no different names!
	for(size_t i=0;i<histomap_.size();i++){
		if(histomap_.at(i).first == histoname)
			return histomap_.at(i).second;
	}
	std::string err="histoCollection::getHisto: ";
	err+=histoname.Data();
	err+=" not found";
	throw std::runtime_error(err);
}

TH1* histoCollection::cloneHisto(const TString& histoname)const{
	return (TH1*) getHisto(histoname)->Clone();
}

}
