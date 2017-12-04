/*
 * exampleanalyser.cpp
 *
 *  Created on: 5 Aug 2016
 *      Author: jkiesele
 */



#include "../interface/basicAnalyzer.h"
#include "../interface/tTreeHandler.h"
#include "../interface/tBranchHandler.h"

/*
 *
 * In reality there should be Delphes class includes
 *
 */

#include <iostream>

class exampleanalyser: public d_ana::basicAnalyzer{
public:
	exampleanalyser():d_ana::basicAnalyzer(){}
	~exampleanalyser(){}



private:
	/*
	 *
	 * the following have to be implemented here
	 *
	 */
	void analyze(size_t id){

		//open file, get some branches (commented since input files missing)
		//	d_ana::tTreeHandler tree(getSampleFile(),"treename");
		//	d_ana::tBranchHandler<std::vector<jet> > jets=d_ana::tBranchHandler<float>(&tree,"jets");




		/*
		 *
		 * The delphes trees have arrays instead of vectors, so there needs to be some
		 * king of explicit buffering. I did not realise this before. You can get ideas from
		 * https://github.com/jkiesele/TtZAnalysis/blob/master/Analysis/interface/wNTBaseInterface.h
		 *
		 * Something with dictionaries for the Delphes classes MIGHT be needed
		 * tTreeHandler does not work with TChains! This is due to a root limitation
		 * (at least from root 5). But it should not be a problem.
		 * Either hadd it or give it same legend name etc.
		 * For same legend name, output histogram files should just be added (taking into
		 * account the normalisation). This can be done in the writeOutput() function.
		 * This function is thread-safe. The file access is blocked while one process writes.
		 *
		 */

		std::cout << "event loop on " << getSampleFile()  <<std::endl;

		size_t nevents=1000;
		if(isTestMode())
			nevents=100;

		if(id>0)
			sleep(2); //just to hve different status % in this example
		for(size_t i=0;i<nevents;i++){

			//	tree.setEntry(i); //associate event entry

			reportStatus(i,nevents);
			usleep(4e4);

			//	size_t njets=jets.content()->size(); //how to access the branch content
		}

		processEndFunction(); //needs to be called
	}



};


int main(){

	exampleanalyser ana;


	ana.setDataSetDirectory("~/eos/cms/store/group/upgrade/delphes_framework");
	//adjust to what is needed. should aso be read in from file when we are done with testing

	ana.readConfigFile("exampleinput.txt");
	//	ana.setMaxChilds(1);

	//ana.setTestMode(true);//or not

	ana.start();

}

