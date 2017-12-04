

#include "../interface/fileForker.h"
#include "../interface/helpers.h"
#include <sys/stat.h>
#include <iostream>
#include <sys/types.h>
#include <stdexcept>
#include "TString.h"
#include <sys/wait.h>
#include <signal.h>

#include <execinfo.h>
#include <algorithm>

namespace d_ana{



bool fileForker::debug=false;

fileForker * current;

void segfault_callback_handler(int signum){
	current->writeReady_block(); //pro-forma
	current->writeDone(fileForker::ff_status_child_segfault);
	const int lines=50;
	void *array[lines];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, lines);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", signum);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(signum);
}
void kill_callback_handler(int signum){
	current->cleanUp();
	exit(signum);
}


fileForker::fileForker():ischild_(false),spawnready_(false),
		maxchilds_(12),
		runningchilds_(0),
		donechilds_(0),
		lastspawned_(-1),
		PID_(0),
		ownchildindex_(FF_PARENTINDEX),
		ownrunningindex_(FF_PARENTINDEX),
		processendcalled_(false)
{

	current=this;

	signal(SIGHUP ,kill_callback_handler);
	signal(SIGINT ,kill_callback_handler);
}


fileForker::~fileForker(){
	cleanUp();
}


fileForker::fileforker_status fileForker::prepareSpawn(){
	if(!createOutFile())
		return ff_status_parent_generror;
	lastspawned_=-1;
	size_t filenumber=inputfiles_.size();
	if(!filenumber)
		throw std::logic_error("fileForker::prepareSpawn: no input");

	//open maxchilds pipes and communicate p_idx first
	p_idx.open(maxchilds_);
	p_allowwrite.open(maxchilds_);
	p_askwrite.open(maxchilds_);
	p_status.open(maxchilds_);
	p_busystatus.open(maxchilds_);
	runningidxs_.resize(maxchilds_,-100);//a bit of headroom

	childPids_.resize(filenumber,0);
	runningchilds_=0;
	busystatus_.resize(filenumber,0);
	childstatus_.resize(filenumber,ff_status_child_hold);
	ischild_=false;
	spawnready_=true;
	return ff_status_parent_success;
}

fileForker::fileforker_status fileForker::spawnChildsAndUpdate(){
	if(debug)
		std::cout << "fileForker::spawnChilds: " << inputfiles_.size()<<std::endl;
	using namespace std;
	using namespace d_ana;
	if(!spawnready_)
		throw std::logic_error("fileForker::spawnChilds: first prepare spawn");


	size_t filenumber=inputfiles_.size();
	//spawn
	if(filenumber>0){

		for(int i=lastspawned_+1;i<(int)filenumber;i++){

			if(runningchilds_>=maxchilds_)
				break;
			if(debug)
				std::cout << "last: " << lastspawned_<< " i: "<< i << std::endl;

			updateChildrenStatus();
			int newrunidx=-1;

			//search for free running index
			static size_t ri=0;
			while(1){
				if(runningidxs_.at(ri)<0){
					newrunidx=ri;
					runningidxs_.at(ri)=i;
					break;
				}
				ri++;
				if(ri>=runningidxs_.size())ri=0;
			}

			childPids_.at(i)=fork();
			if(childPids_.at(i)==0){//child process
				ischild_=true;
				current=this;
				signal(SIGSEGV, segfault_callback_handler);
				try{
					if(debug)
						std::cout << "fileForker::spawnChilds child "<<i << " running index: " << newrunidx << std::endl;
					ownrunningindex_=newrunidx;
					ownchildindex_=p_idx.get(ownrunningindex_)->pread(); //wait until parent got ok
					p_status.get(ownrunningindex_)->pwrite((int)ff_status_child_busy);
					reportBusyStatus(0);
					processendcalled_=false;
					process();
					if(!processendcalled_){
						throw std::logic_error("fileForker: implementation of process() must call processEndFunction() at its end!");
					}
					std::exit(EXIT_SUCCESS);
				}catch(std::exception& e){
					std::cout << "*************** Exception ***************" << std::endl;
					std::cout << e.what() << std::endl;
					std::cout << "in " << getInputFileName() << '\n'<<std::endl;
					writeReady_block(); //pro-forma
					writeDone(ff_status_child_exception);
					std::exit(EXIT_FAILURE);
				}
			}
			else{ //parent process
				if(debug)
					std::cout << "fileForker::spawnChilds parent "<< std::endl;
				ownchildindex_=FF_PARENTINDEX;//start processing
				p_idx.get(newrunidx)->pwrite(i); //send start signal
				//usleep(10);
				//updateBusyStatus();
				runningchilds_++;
				lastspawned_=i;

				//check for write requests
			}
		}
		//updateBusyStatus();

		int it_readytowrite=checkForWriteRequest();
		updateChildrenStatus();
		if(it_readytowrite >= 0){ // daugh ready to write

			size_t globalindex=getGlobalIndex(it_readytowrite);
			if(debug)
				std::cout << "fileForker::spawnChilds waiting for write " << globalindex<< " run: "<< it_readytowrite << std::endl;
			p_allowwrite.get(it_readytowrite)->pwrite(globalindex);
			childstatus_.at(globalindex)=ff_status_child_writing;

			childstatus_.at(globalindex)=(fileforker_status)p_status.get(it_readytowrite)->pread();
			//done
			updateChildrenStatus();
			runningidxs_.at(it_readytowrite)=-1;//free index again
			runningchilds_--;
			donechilds_++;
			//childPids_.at(globalindex)=0;
		}




		if(donechilds_ == filenumber)
			return ff_status_parent_success;

		if(it_readytowrite>=0)
			return ff_status_parent_filewritten;

	}
	else{ //only one process
		std::cout << "fileForker::spawnChilds: no input file, nothing to do"  << std::endl;
	}


	return ff_status_parent_busy;
}

fileForker::fileforker_status fileForker::getStatus(const size_t& childindex)const{
	return childstatus_.at(childindex);
}

fileForker::fileforker_status fileForker::getStatus()const{
	//update child status

	if(debug){
		std::cout << "fileForker::getStatus: "<< donechilds_<<"("<< runningchilds_<< ")/"<< inputfiles_.size()<<
				" last "<< lastspawned_<< std::endl;
		std::cout << "child stats: " <<std::endl;
		for(size_t c=0;c<inputfiles_.size();c++){
			std::cout << c << "\t"  << (int)childstatus_.at(c) <<"\t"<<inputfiles_.at(c) <<std::endl;
		}
	}
	if(! inputfiles_.size())
		return ff_status_parent_success;
	if(lastspawned_<(int)inputfiles_.size()-1)
		return ff_status_parent_childstospawn;
	if(donechilds_<inputfiles_.size())
		return ff_status_parent_busy;
	for(size_t i=0;i<childstatus_.size();i++){
		if(childstatus_.at(i) == ff_status_child_aborted
				|| childstatus_.at(i) == ff_status_child_exception
				|| childstatus_.at(i) == ff_status_child_generror)
			return ff_status_parent_generror;
	}

	return
			ff_status_parent_success;
}



int fileForker::getBusyStatus(const size_t& childindex){
	updateChildrenStatus();
	return busystatus_.at(childindex);
}

bool fileForker::writeReady_block(){
	if(debug)
		std::cout << "fileForker::writeReady_block: ischild: "<< ischild_ << std::endl;
	if(!ischild_)
		return true;

	p_askwrite.get(ownrunningindex_)->pwrite(ownchildindex_);
	p_allowwrite.get(ownrunningindex_)->pread();
	if(debug)
		std::cout << "fileForker::writeReady_block: done" << std::endl;
	return true;

}
void fileForker::writeDone(fileforker_status stat){
	if(debug)
		std::cout << "fileForker::writeDone" << std::endl;
	if(!ischild_)
		return;
	p_status.get(ownrunningindex_)->pwrite(stat);

}

size_t fileForker::getGlobalIndex(const size_t & runningidx)const{
	if(runningidx>=runningidxs_.size())
		throw std::out_of_range(("fileForker::getGlobalIndex "+toString(runningidx)+"/"+toString(runningidxs_.size())));
	return runningidxs_.at(runningidx);
}
size_t fileForker::getRunningIndex(const size_t & globalidx)const{
	size_t pos=std::find(runningidxs_.begin(),runningidxs_.end(),(int)globalidx)-runningidxs_.begin();
	if(pos<runningidxs_.size())
		return pos;
	throw std::out_of_range("fileForker::getRunningIndex");
}

void fileForker::updateChildrenStatus(){
	for(size_t i=0;i<runningidxs_.size();i++){
		if(runningidxs_.at(i)>=0){
			size_t globalidx=getGlobalIndex(i);
			while(p_busystatus.get(i)->preadready())
				busystatus_.at(globalidx)=p_busystatus.get(i)->pread();
			while(p_status.get(i)->preadready())
				childstatus_.at(globalidx)=(fileforker_status)p_status.get(i)->pread();
		}
	}
}
std::string fileForker::translateStatus(const fileforker_status& stat)const{
#define RETURNSTATUSSTRING(x) {if(stat==ff_status_ ## x) str= #x;}
	std::string str;
	RETURNSTATUSSTRING(child_hold)
	RETURNSTATUSSTRING(child_busy)
	RETURNSTATUSSTRING(child_writing)
	RETURNSTATUSSTRING(child_success)
	RETURNSTATUSSTRING(child_generror)
	RETURNSTATUSSTRING(child_segfault)
	RETURNSTATUSSTRING(child_exception)
	RETURNSTATUSSTRING(child_aborted)
	RETURNSTATUSSTRING(parent_childstospawn)
	RETURNSTATUSSTRING(parent_success)
	RETURNSTATUSSTRING(parent_busy)
	RETURNSTATUSSTRING(parent_generror)
	RETURNSTATUSSTRING(parent_exception)
	RETURNSTATUSSTRING(parent_filewritten)
	RETURNSTATUSSTRING(allowwrite)
	RETURNSTATUSSTRING(askwrite)
	return str;
#undef RETURNSTATUSSTRING
}

void fileForker::reportStatus(fileforker_status stat){
	p_status.get(ownrunningindex_)->pwrite(stat);
}
void fileForker::reportBusyStatus(int bstat){
	p_busystatus.get(ownrunningindex_)->pwrite(bstat);
}

void fileForker::processEndFunction(){
	if(debug)
		std::cout << "fileForker::processEndFunction" << std::endl;
	processendcalled_=true;
	writeReady_block() ;
	fileforker_status stat=writeOutput();
	writeDone(stat);
}

int fileForker::checkForWriteRequest(){
	for(size_t i=0;i<p_askwrite.size();i++){
		if(p_askwrite.get(i)->preadready()){
			p_askwrite.get(i)->pread(); //to free pipe again
			return i; //return first ready to write!
		}
	}
	return -1;
}

void fileForker::abortChild(size_t idx, bool dokill){
	if(idx>=childPids_.size())
		throw std::out_of_range("fileForker::abortChild: index");

	if(dokill)
		kill(childPids_.at(idx),-9);
	childstatus_.at(idx)=ff_status_child_aborted;
	runningidxs_.at(getRunningIndex(idx))=-1;
	runningchilds_--;
	donechilds_++;
}
void fileForker::cleanUp(){
	for(size_t i=0;i<childPids_.size();i++){
		if(childPids_.at(i)>0)
			kill(childPids_.at(i),15);
	}
	sleep(1);
	for(size_t i=0;i<childPids_.size();i++){
		if(childPids_.at(i)>0)
			kill(childPids_.at(i),9);
	}

}

}




