/*
 * Pipes.h
 *
 *  Created on: Jul 19, 2013
 *      Author: kiesej
 *
 *      Some template classes be used for inter program communication
 * IPCPipe and IPCPipes ("vector" of IPCPipe)
 */

#ifndef PIPES_H_
#define PIPES_H_


#include <sys/types.h>
#include <unistd.h>
#include <poll.h>
#include <vector>



class IPCPipeBase{
public:
	IPCPipeBase();
	~IPCPipeBase();
private:
	static size_t openpipes;
};

#define IPC_BUFFERSIZE 128
/**
 * class to enable communation between programs through pipes
 * it can only pass simple and basic data formats
 * the buffer size is limited to 128
 */
template<class T>
class IPCPipe: public IPCPipeBase{
public:
	/**
	 * opens a pipe with (p)read and (p)write end.
	 * additional function preadready checks if anything has been written to the buffer from write side
	 * buffersize is a define and set to a default of 256. In case it is to change, change IPC_BUFFERSIZE
	 * the preadready function polls for 2 ms until it returns
	 */
	IPCPipe():IPCPipeBase(){
		fds[0].events = POLLRDNORM | POLLIN;
		fds[1].events = POLLOUT | POLLWRBAND;;
		pipe(pfds);
	}
	~IPCPipe(){close(pfds[0]);close(pfds[1]);}

	T pwrite(T c){return write(pfds[1], &c , IPC_BUFFERSIZE);}
	T pread(){read(pfds[0], buf, IPC_BUFFERSIZE);return buf[0];}

	int preadready(){
		fds[0].fd=pfds[0];
		fds[1].fd=pfds[1];
		poll(fds, 2, 2);
		if (fds[0].revents & (POLLRDNORM | POLLIN))
			return 1;
		else
			return 0;
	}


private:
	int pfds[2];
	struct pollfd fds[2];
	T buf[IPC_BUFFERSIZE];

};

/**
 * class to enable communation between programs through a set of pipes
 * it can only pass simple and basic data formats
 * the buffer size of each pipe limited to 128
 */
template<class T>
class IPCPipes{
public:
	/**
	 * default constructor. Doesn't do anything
	 */
	IPCPipes(){size_=0;}
	/**
	 * constructor that opens <Size> pipes
	 */
	IPCPipes(size_t Size){size_=Size;open(size_);}
	~IPCPipes(){if(size_!=0)closePipes();}

	/**
	 * opens <numPipes> pipes (closes existing ones)
	 */
	void open(size_t numPipes){size_=numPipes;openPipes();}
	size_t size(){return size_;}

	IPCPipe<T> * get(size_t i){if(i<size_) return ipcs_[i]; else return 0;}
	// size_t pollAll(int timeout=100)

private:
	size_t size_;
	std::vector<IPCPipe<T> *> ipcs_;
	void openPipes(){
		closePipes();
		for(size_t i=0;i<size_;i++){
			ipcs_.push_back(new IPCPipe<T>);
		}
	}
	void closePipes(){
		for(size_t i=0;i<ipcs_.size();i++)
			if(ipcs_.at(i)) delete ipcs_.at(i);
		ipcs_.clear();
	}
};


#endif /* PIPES_H_ */
