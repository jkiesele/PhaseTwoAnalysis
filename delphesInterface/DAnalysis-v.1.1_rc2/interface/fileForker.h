#ifndef TTZANALYSIS_TOOLS_INTERFACE_FILEFORKER_H_
#define TTZANALYSIS_TOOLS_INTERFACE_FILEFORKER_H_
/** \class d_ana::fileForker
 *
 * baseline implementation for forking a process
 *
 * \original author: Jan Kieseler
 *
 * more docu
 *
 */

/*
 *
 * think about omp tfile access again
 */
#define FF_PARENTINDEX 1e6

#include <string>
#include "pipes.h"
#include <vector>

namespace d_ana{
void segfault_callback_handler(int signum);
/*
 * Can do:
 * -> handle a list of input files
 * -> forks once for each of them
 * -> does something on them (virtual!)
 * -> handles an output file name
 * -> organises read/write to output file with locks
 *    (actual output function virtual)
 * -> handles errors/status messages from childs to main
 *
 *
 * NO ROOT/CMSSW DEPENDENCE HERE!
 *
 */
class fileForker{
	friend void segfault_callback_handler(int signum);
	friend void kill_callback_handler(int signum);
public:
	fileForker();
	virtual ~fileForker();

	enum fileforker_status{
		ff_status_child_hold=0,
		ff_status_child_busy,
		ff_status_child_writing,
		ff_status_child_success,
		ff_status_child_generror,
		ff_status_child_segfault,
		ff_status_child_exception,
		ff_status_child_aborted,

		ff_status_parent_childstospawn,
		ff_status_parent_success,
		ff_status_parent_busy,
		ff_status_parent_generror,
		ff_status_parent_exception,
		ff_status_parent_filewritten,

		ff_status_allowwrite,
		ff_status_askwrite
	};


	static bool debug;

	void setMaxChilds(size_t max){maxchilds_=max;}

	void setOutputFileName(std::string outname){outputfile_=outname;}
	const std::string& getOutputFileName()const{return outputfile_;}

protected:

	void setInputFiles(const std::vector<std::string>& in){inputfiles_=in;}

	//parent functions

	fileforker_status prepareSpawn();

	/*
	 * the following functions should remain in a loop until getStatus() returns a value different
	 * from ff_status_childstospawn or ff_status_busy
	 */
	fileforker_status spawnChildsAndUpdate();
	fileforker_status getStatus()const;


	fileforker_status getStatus(const size_t& childindex)const;

	//a number between 0 and 100 (%)
	int getBusyStatus(const size_t& childindex);


	std::string translateStatus(const fileforker_status&)const;

	/**
	 * implement the initial creation of an output file
	 * return true is successful
	 */
	virtual bool createOutFile()const=0;

	//child functions: use those!
	const std::string& getInputFileName()const{return inputfiles_.at(ownchildindex_);}
	void reportStatus(fileforker_status); //in case of an error or similar
	void reportBusyStatus(int bstat); //to report status of running job e.g. in %
	void processEndFunction(); // call at the end of the job

	const size_t ownChildIndex()const{return ownchildindex_;}

	/**
	 * implement what will be done with each input file (or string as input to the job)
	 */
	virtual void process()=0;

	/**
	 * implement the writing into the output file
	 * it will be protected by locks and only be allowed for one
	 * process at once
	 */
	virtual fileforker_status writeOutput()=0;

	const std::vector<pid_t>& getChildPids()const{return childPids_;}

	void abortChild(size_t idx, bool kill=true);

	/**
	 * Kills all childs
	 */
	void cleanUp();

	const size_t& getMaxChilds()const{return maxchilds_;}

private:
	int checkForWriteRequest();
	bool writeReady_block();
	void writeDone(fileforker_status);


	bool ischild_,spawnready_;
	std::vector<std::string> inputfiles_;
	std::string outputfile_;

	size_t getGlobalIndex(const size_t & runningidx)const;
	size_t getRunningIndex(const size_t & globalidx)const;
	void updateChildrenStatus();
	///communication pipes
	std::vector<int> runningidxs_;
	IPCPipes<int> p_idx;
	IPCPipes<int> p_busystatus;


	IPCPipes<int> p_askwrite;
	IPCPipes<int> p_allowwrite;
	IPCPipes<int> p_status;

	//parent members
	size_t maxchilds_,runningchilds_,donechilds_;
	int lastspawned_;
	std::vector<pid_t> childPids_;
	std::vector<int> busystatus_;
	std::vector<fileforker_status> childstatus_;

	//child members
	pid_t PID_;
	size_t ownchildindex_;
	size_t ownrunningindex_;
	bool processendcalled_;


};
} //ns
#endif /* TTZANALYSIS_TOOLS_INTERFACE_FILEFORKER_H_ */
