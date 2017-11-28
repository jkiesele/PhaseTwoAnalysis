/*
 * metaInfo.cpp
 *
 *  Created on: 22 Aug 2016
 *      Author: jkiesele
 */




#include "TObjString.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TIterator.h"
#include "TDirectory.h"
#include "../interface/metaInfo.h"
#include "../interface/textFormatter.h"

#include <iostream>

namespace d_ana{

const TString metaInfo::separator_="\n|%|\n";

void metaInfo::Write()const{

	TObjString str;

	TString information;
	information+=legendname+separator_;
	information+=color;
	information+=separator_;
	information+=legendorder;
	information+=separator_;
	information+=norm;


	str.SetString(information);
	str.Write();

}


void metaInfo::extractFrom(const TString& str){

	TObjArray *  arr=str.Tokenize(separator_);
	TIter iString(arr);
	TObjString* os=0;
	while ((os=(TObjString*)iString())) {
		legendname=  os->GetString();
		color = ((TObjString*)iString())->GetString().Atoi();
		legendorder= ((TObjString*)iString())->GetString().Atoi();
		norm= ((TObjString*)iString())->GetString().Atof();
	}

}

void metaInfo::extractFrom( TDirectory* dir){

	TIter next( dir->GetListOfKeys());
	TObject *obj;
	TKey *key;
	while ((key = (TKey *) next())) {
		obj = dir->Get(key->GetName());
		if(obj->InheritsFrom(TObjString::Class_Name())){
			TObjString * ostr= (TObjString *)obj;
			extractFrom(ostr->GetString());
			break;
		}
	}

}


}//ns
