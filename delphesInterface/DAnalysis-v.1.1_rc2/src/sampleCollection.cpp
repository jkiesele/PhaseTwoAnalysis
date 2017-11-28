/*
 * sampleCollection.cpp
 *
 *  Created on: 8 Nov 2016
 *      Author: jkiesele
 */

#include "../interface/sampleCollection.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TString.h"
#include "TKey.h"


namespace d_ana{

const histoCollection& sampleCollection::getHistos(const TString&  samplename)const{
	for(size_t i=0;i<histosPerSample_.size();i++){
		if(histosPerSample_.at(i).getLegendName() == samplename)
			return histosPerSample_.at(i);
	}
	std::string err="sampleCollection::getHistos: ";
	err+=samplename.Data();
	err+=" not found";
	throw std::runtime_error(err);
}

void sampleCollection::readFromFile(const TString& filename){
	TFile *fIn = new TFile(filename,"READ");
	if(!fIn || fIn->IsZombie())
		throw std::runtime_error("sampleCollection::readFromFile: could not read file " +(std::string)filename.Data());

	fIn->cd();


	//ownership by histoCollection
	TH1::AddDirectory(kFALSE);
	TDirectory::AddDirectory(kFALSE);

	TIter   dirIter(fIn->GetListOfKeys());
	TObject *cDirObj;
	TKey    *key;

	// iterate over directories and get all stacks
	while((key = (TKey *) dirIter())) {
		cDirObj=fIn->Get(key->GetName());
		if(!cDirObj->InheritsFrom(TDirectory::Class())) continue;

		TDirectory* cDir = (TDirectory*) cDirObj;
		metaInfo tMI;
		tMI.extractFrom(cDir);

		histoCollection tmphistcoll;
		tmphistcoll.meta_=tMI;

		TIter    histIter(cDir->GetListOfKeys());
		TObject* cHistObj;
		TKey*    cHistKey;

		while((cHistKey = (TKey*) histIter())) {
			cHistObj=cDir->Get(cHistKey->GetName());
			if(!cHistObj->InheritsFrom(TH1::Class())) continue;

			TH1* cHist = (TH1*) cHistObj->Clone();

			//takes ownership
			tmphistcoll.histomap_.push_back(std::pair<TString, TH1*> (cHist->GetName(), cHist));

		}

		histosPerSample_.push_back(tmphistcoll);
	}

	// cleanup
	fIn->Close();
	delete fIn;

}

std::vector<TString> sampleCollection::listAllLegends()const{
	std::vector<TString> out(histosPerSample_.size());
	for(size_t i=0;i<histosPerSample_.size();i++){
		out.at(i)=histosPerSample_.at(i).getLegendName();
	}
	return out;
}

}
