/*
 * tBranchHandler.cc
 *
 *  Created on: May 28, 2014
 *      Author: kiesej
 */


#include "../interface/tBranchHandler.h"

#include <algorithm>

namespace d_ana{



std::map< tTreeHandler* ,std::vector<TString> > tBranchHandlerBase::branchesfortree_;

bool tBranchHandlerBase::debug=false;

void tBranchHandlerBase::handleReturns(int ret,bool& missing,bool allowmissing)const{
	if(ret == -2 || ret == -1){
		std::cout << "tBranchHandler: Class type given for branch " << branchname_
				<< " does not match class type in tree. (root CheckBranchAddressType returned " << ret << ")" <<std::endl;
		throw std::runtime_error("tBranchHandler: Class type does not match class type in branch");
	}
	else if(ret == -4 || ret == -3){
		std::cout << "tBranchHandler: Internal error in branch " << branchname_
				<< " (root CheckBranchAddressType returned " << ret << ")" <<std::endl;
		throw std::runtime_error("tBranchHandler: Internal error in branch");
	}
	else if( ret == -5){
		if(allowmissing){
			std::cout << "tBranchHandler: allowed missing branch "<<branchname_<<std::endl;
			missing= true;
		}
		else{
			std::cout << "tBranchHandler: branch " << branchname_ << " does not exists!" << std::endl;
			throw std::runtime_error("tBranchHandler: branch does not exists!");
		}
	}

}

void tBranchHandlerBase::addTreeAndBranch(tTreeHandler * t, const TString& branchname){

	if(std::find(branchesfortree_[t].begin(),branchesfortree_[t].end(),branchname) != branchesfortree_[t].end()){
		std::string err="tBranchHandlerBase::addTreeAndBranch: Only one handler per branch allowed! (";
		err+=branchname.Data();
		err+=")";
		throw std::logic_error(err);
	}
	branchesfortree_[t].push_back(branchname);
	if(t_)
		t->associate(this);
}
void tBranchHandlerBase::removeTreeAndBranch( tTreeHandler * t, const TString& branchname){

	std::vector<TString> & branches=branchesfortree_[t];
	std::vector<TString>::iterator it=std::find(branches.begin(),branches.end(),branchname);
	if(it != branches.end()){
		branchesfortree_[t].erase(it);
	}
	if(t_)
		t->removeAsso(this);
	if(!missingbranch_)
		removeTree(t);
}


//all inlined

}

