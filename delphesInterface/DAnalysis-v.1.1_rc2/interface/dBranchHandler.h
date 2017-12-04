/*
 * dBranchHandler.h
 *
 *  Created on: 15 Aug 2016
 *      Author: jkiesele
 */

#ifndef INTERFACE_DBRANCHHANDLER_H_
#define INTERFACE_DBRANCHHANDLER_H_
/** \class d_ana::dBranchHandler
 *
 * Safe and easy way to read delphes classes from branches
 *
 * \original author: Jan Kieseler
 *
 * more docu
 *
 */

#include "tBranchHandler.h"
#include "TClonesArray.h"
#include "TBranchElement.h"
#include "TROOT.h"
#include <cxxabi.h>
#include "textFormatter.h"

namespace d_ana{

template<class T>
class dBranchHandler{
public:
	dBranchHandler(tTreeHandler * t, const TString& branchname):
		sizebranch_(t,branchname+"_size"),databranch_(t,branchname,1){

		if(databranch_.realcontent_->Capacity() ){//root doesn't want you to delete here
			databranch_.realcontent_->Clear();
		}
		databranch_.realcontent_=0;


		if(databranch_.branch_->IsA() == TBranchElement::Class())
		{
			TClonesArray * array=0;
			TBranchElement *element = static_cast<TBranchElement*>(databranch_.branch_);
			const char *className = element->GetClonesName();
			Int_t size = element->GetMaximum();
			TClass *cl = gROOT->GetClass(className);
			TString branchclass=className;
			//check for template
			T testobject;
			int status=0;
			TString templatename=abi::__cxa_demangle(typeid(testobject).name(),0,0,&status);

			if(cl)
			{
				array = new TClonesArray(cl, size);
				array->SetName(databranch_.getBranchName());
				databranch_.realcontent_=array;
				databranch_.tree()->tree()->SetBranchAddress(databranch_.getBranchName(),&databranch_.realcontent_,&databranch_.branch_);
			}
			else{
				throw std::runtime_error("dBranchVectorHandler: could not create TClonesArray");
			}
			//self-sonsitency check
			if(templatename != branchclass){
				//tBranchHandlerBase::debug=true;
				if(databranch_.realcontent_->Capacity() ){
					databranch_.realcontent_->Clear();
					delete databranch_.realcontent_;
				}
				databranch_.realcontent_=0;
				databranch_.pcontent_=0;
				throw std::runtime_error(("dBranchVectorHandler: trying to read branch of \"" +
						branchclass+ "\" with template argument \"" + templatename+"\"").Data());
			}
		}

	}

	~dBranchHandler(){
		if(databranch_.realcontent_->Capacity()){
			databranch_.realcontent_->Clear();
			delete databranch_.realcontent_;
		}
		databranch_.realcontent_=0;
		databranch_.pcontent_=0;
	}


	T  * at(const size_t & idx){
		size_t s=size();

		if(idx >= s)
			throw std::out_of_range(
					"dBranchVectorHandler::at(): index "+ textFormatter::toString(idx)+
					" out of range ("+ textFormatter::toString(s)+ ") in branch " +
					textFormatter::toString(databranch_.getBranchName()));

		TClonesArray * arr=(*databranch_.content());
		T* ret=(T*) arr->At(idx);
		return ret;
	}

	size_t size(){
		return *sizebranch_.content();
	} //cannot be const - read may be triggered

private:

	void getEntry(const Long64_t& entry){
		sizebranch_.getEntry(entry);
		tBranchHandler<TClonesArray*>::getEntry(entry);
	}
	tBranchHandler<int> sizebranch_;
	tBranchHandler<TClonesArray*> databranch_;


};
}






#endif /* INTERFACE_DBRANCHHANDLER_H_ */
