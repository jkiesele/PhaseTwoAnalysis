/*
 * exampleanalyser.cpp
 *
 *  Created on: 5 Aug 2016
 *      Author: jkiesele
 */



#include "../interface/basicAnalyzer.h"
#include "../interface/tTreeHandler.h"
#include "../interface/tBranchHandler.h"
#include "../interface/dBranchHandler.h"


/*
 *
 * In reality there should be Delphes class includes
 *
 */
#include "classes/DelphesClasses.h"


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
		//std::cout << "analyse() on " << getSampleFile()  <<std::endl;

		//open file, get some branches (commented since input files missing)
		//	d_ana::tTreeHandler tree(getSampleFile(),"treename");
		//	d_ana::tBranchHandler<std::vector<jet> > jets=d_ana::tBranchHandler<float>(&tree,"jets");






		d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");

        //debug=true;

        // Add a histogram to the analysis
		TH1* histo=addPlot(new TH1D("histoname1","histotitle1",100,0,100) , "p_{T} [GeV]", "N_{e}");

        // Create tree in analysis (Default name is DAnalysis)
        TTree* anatree=addTree();

        // Add a branch to the tree
        Double_t elecPt=0;
        anatree->Branch("pte", &elecPt);

        // Add an object collection to the tree
        std::vector<Electron> skimmedelecs;
        anatree->Branch("Electrons",&skimmedelecs);


		size_t nevents=tree()->entries() ;
		if(isTestMode())
			nevents/=10;
		if(nevents<100)
			nevents=tree()->entries() ;

		for(size_t i=0;i<nevents;i++){

			tree()->setEntry(i); //associate event entry
			reportStatus(i,nevents);

			//std::cout << elecs.size() <<std::endl;
			if(elecs.size()>0) {
				histo->Fill(elecs.at(0)->PT);
                elecPt=elecs.at(0)->PT;
                skimmedelecs.clear();

                for(size_t elec=0;elec<elecs.size();elec++){
                	if(elecs.at(elec)->PT < 20 ) continue;
                	skimmedelecs.push_back(*elecs.at(elec));
                }

            }
			//usleep(4e4);

			//	size_t njets=jets.content()->size(); //how to access the branch content
            anatree->Fill();
		}

		processEndFunction(); //needs to be called
	}



};


int main(int argc, char* argv[]){

	exampleanalyser ana;
	if(argc!=2)exit(-1);

	ana.readConfigFile((std::string)argv[1]);
	//ana.readConfigFile("exampleinput.txt");
	ana.start();

}

