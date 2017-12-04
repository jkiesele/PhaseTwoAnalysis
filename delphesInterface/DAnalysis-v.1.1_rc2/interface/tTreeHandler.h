#ifndef TTREEHANDLER_H_
#define TTREEHANDLER_H_
/** \class d_ana::tTreeHandler
 *
 * setting up a root tree in the tBranchHandler framework (more performant and safe)
 *
 * \original author: Jan Kieseler
 *
 * more docu
 *
 */

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <stdexcept>


namespace d_ana{

class tBranchHandlerBase;

class tTreeHandler{
public:
	tTreeHandler(const TString & filename, const TString &treename);
	tTreeHandler();
	~tTreeHandler();

	void load(const TString & filename, const TString &treename);
	void clear();

	const Long64_t& entries()const{return entries_;}
	void setEntry(const Long64_t& in);
	const Long64_t& currentEntry()const{return entry_;}
	void associate( tBranchHandlerBase*);
	void removeAsso( tBranchHandlerBase*);

	TTree * tree(){
		if(t_) return t_;
		throw std::logic_error("tTreeHandler::tree(): tree pointer 0");
	}

	/**
	 * Improves performance if many jobs access the same file
	 * Enables cashing of a few MB to RAM
	 */
	void setPreCache();

	void printStats()const;

	static bool debug;

protected:




private:

	tTreeHandler& operator = (const tTreeHandler& rhs){return *this;}
	tTreeHandler(const tTreeHandler&);

	TFile * file_;

	TTree * t_;

	Long64_t entry_;

	Long64_t entries_;

	std::vector< tBranchHandlerBase*> assobranches_;
};



}



#endif /* TTREEHANDLER_H_ */
