#ifndef TTZANALYSIS_ANALYSIS_INTERFACE_BASICANALYZER_H_
#define TTZANALYSIS_ANALYSIS_INTERFACE_BASICANALYZER_H_
/** \class d_ana::basicAnalyzer
 *
 * user-friendly parallel delphes analysis class
 *
 * \original author: Jan Kieseler
 *
 * more docu
 *
 */


#include "../interface/fileForker.h"
#include "../interface/textFormatter.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <vector>
#include <map>
//includes to make the classes already available, even though not needed here
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "../interface/basicAnalyzer.h"
#include "../interface/tTreeHandler.h"
#include "../interface/tBranchHandler.h"
#include "../interface/dBranchHandler.h"

/**
 * quick and dirty generic analysis interface
 */
namespace d_ana{
class basicAnalyzer : public fileForker{
public:
	basicAnalyzer();
	virtual ~basicAnalyzer();

	void readConfigFile(const std::string& );

	void setDataSetDirectory(const TString& dir){
		datasetdirectory_=dir;
		if(!datasetdirectory_.EndsWith("/"))
			datasetdirectory_+="/";
	}

	void setOutDir(const TString& dir);
	TString getOutDir()const{return outdir_;}

	TString getOutPath()const{return outdir_+getOutFileName();}
	TString getTreePath()const{
		return TString(textFormatter::stripFileExtension(getOutFileName().Data()))+"_ntuples.root";
	}


	//setters
	void setLumi(double Lumi){lumi_=Lumi;}
	void setSyst(const TString& syst){syst_=syst;}

	void setFilePostfixReplace(const TString& file,const TString& pf,bool clear=false);
	void setFilePostfixReplace(const std::vector<TString>& files,const std::vector<TString>& pf);

	void setTestMode(bool test){testmode_=test;}

	void setWriteTree(bool write=true){writeTree_=write;}

	//getters
	const TString& getSyst()const{return syst_;}


	virtual TString getOutFileName()const{
		if(syst_.Length()){
			return  (TString)getOutputFileName()+"_"+syst_;}
		else{
			return getOutputFileName();}
	}

	void start();


protected:

	/**
	 * Implements the event loop
	 */
	virtual void analyze(size_t )=0;

	/**
	 * Implements the extraction of parameters or similar from the output
	 */
	virtual void postProcess()=0;


	/**
	 * for child processes
	 * reports the Status (% of events already processed) to the main program
	 */
	void reportStatus(const Long64_t& entry,const Long64_t& nEntries);


	tTreeHandler* tree(){return tree_;}

	//adders
	TH1* addPlot(TH1* histo, const TString xaxis="", const TString yaxis="", const TString zaxis="");
	TTree* addTree(const TString& name="Delphes");

	const TString& getSampleFile()const{return thissample_.getInfile();}
	TString getSamplePath()const{return datasetdirectory_+thissample_.getInfile();}
	const TString& getLegendName()const{return thissample_.getLegend();}
	const int& getColor()const{return thissample_.getColor();}
	const double& getNorm()const{return thissample_.getNorm();}
	const double& getXsec()const{return thissample_.getXsec();}
        const int& getLegendOrder()const{return thissample_.getLegendorder();}
	const bool& getIsSignal()const{return thissample_.isSignal();}
	void setIsSignal(bool set){ thissample_.setSignal(set);}

	bool isTestMode()const{return testmode_;}

	void processEndFunction();

	void setTreeName(const TString& name){treename_=name;}
private:

	void process();

	void adjustNormalization(const tTreeHandler*);

	std::vector<TString> lsDirectory(const TString & dir, bool& hasmetadata, const TString mdtag="metaData",const TString sampleextension=".root")const;

	unsigned long getTotalEntries(const std::vector<TString>& infiles_indirectory)const;

	fileForker::fileforker_status  writeOutput();
	fileForker::fileforker_status  runParallels(int displaystatusinterval);

	/**
	 * returns false if already exists
	 */
	bool writeConfigHeader(std::fstream &file)const;
	void writeConfigFooter()const;

	TString makeNTupleFileName()const;



	bool createOutFile()const;

	class sampleDescriptor{
	public:
		sampleDescriptor(
				const TString& Infile,
				const TString& Legend,
				int Color,
				double Norm,
				long Direntries,
				int LegendOrder,
				bool Issignal,
				const TString& extraOpts
		);

		sampleDescriptor(const sampleDescriptor&);
		sampleDescriptor& operator=(const sampleDescriptor& rhs);

		const int& getColor() const {return color_;}
		const long& getDirentries() const {return direntries_;}
		void setDirentries(long entr){direntries_=entr;}
		const TString& getExtraopts() const {return extraopts_;}
		const TString& getInfile() const {return infile_;}
		const TString& getLegend() const {return legend_;}
		const int& getLegendorder() const {return legendorder_;}
		const double& getNorm() const {return norm_;}
		const double& getXsec()const{return xsec_;}
		void setNorm(double norm) {this->norm_ = norm;}
		const bool& isSignal() const {return signal_;}
		void setSignal(bool in) {signal_=in;}
	private:
		sampleDescriptor();
		TString infile_;
		TString legend_;
		int color_;
		double norm_;
		double xsec_;
		long direntries_;
		int legendorder_;
		bool signal_;
		TString extraopts_;
	};

	std::string configfile_;

	std::vector<sampleDescriptor> samples_;
	/*
	std::vector<TString> infiles_,legentries_;
	std::vector<int> colz_;
	std::vector<double> norms_;
	std::vector<long> direntriesv_;
	std::vector<size_t> legords_;
	std::vector<bool> issignal_;
	std::vector<TString> extraopts_;
	 */
	std::map<TString,TH1*> histos_;

	Bool_t rewriteoutfile_=true;
	Bool_t rewritentuple_=true;
	Bool_t writeTree_=true;
	TFile *ntuplefile_;
	TTree *ntuples_;

	sampleDescriptor thissample_;
	/*
	///child variables
	TString inputfile_;
	TString legendname_;
	int col_;
	double norm_;
	size_t legorder_;
	bool signal_;
	long direntries_;
	 */

	TString datasetdirectory_;

	TString syst_;
	double lumi_;
	std::vector<TString> fwithfix_,ftorepl_;
	int freplaced_;
	bool testmode_;

	bool isMC_;

	TString datalegend_;
	TString outdir_;

	TString treename_;
	tTreeHandler* tree_;

	bool runonoutputonly_;
};



}
#endif /* TTZANALYSIS_ANALYSIS_INTERFACE_BASICANALYZER_H_ */
