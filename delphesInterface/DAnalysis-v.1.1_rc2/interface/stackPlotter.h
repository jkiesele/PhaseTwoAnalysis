#ifndef TTZANALYSIS_ANALYSIS_INTERFACE_STACKPLOTTER_H_
#define TTZANALYSIS_ANALYSIS_INTERFACE_STACKPLOTTER_H_
/** \class d_ana::stackPlotter
 *
 * THStack script for formatted control plots  
 *
 * \original author: eacoleman 
 *
 * 
 *
 */


#include "../interface/metaInfo.h"
#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TKey.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "THStack.h"

/**
 * simple stack plotter for output files 
 */
namespace d_ana{
class stackPlotter {
public:
	stackPlotter();
	virtual ~stackPlotter();

    void plot();

    // runtime settings
    void rewriteOutfile(Bool_t rewrite=true){rewriteoutfile_=rewrite;} 
    void savePlots(Bool_t save=true){saveplots_=save;} 
    void saveCanvasRootFile(Bool_t save=true){savecanvases_=save;} 

	// setters
	void setInputFile(const TString& fIn){infile_=fIn;}
	void setOutDir(const TString& dir){outdir_=dir;}
	void setLumi(double Lumi){lumi_=Lumi;}
	void setTestMode(bool test){testmode_=test;}
    
protected:

	enum legendposition{lp_left,lp_top,lp_right};

	bool isTestMode()const{return testmode_;}

	virtual legendposition estimateBestLegendPosition(TH1* )const;

	virtual void applyStyleToAxis(THStack *, legendposition pos)const;
	virtual void applyStyleToCanvas(TVirtualPad* ,legendposition pos)const;
	virtual void applyStyleToTH1(TH1* ,legendposition pos)const;
	virtual void applyStyleToLegend(TLegend* ,legendposition )const;

	double yscaling=1.1;
	double yscaling_toplegend=1.4;

private:

    // utils for stacks
    std::map<TString, std::vector< std::pair<Int_t,TH1*> > > stacksLegEntries_;
    void moveDirHistsToStacks(TDirectory* tdir);
    void plotStack(const TString& key);

    // runtime settings
    Bool_t rewriteoutfile_=true;
    Bool_t saveplots_=true;
    Bool_t savecanvases_=false;
	Bool_t testmode_=false;
    Bool_t debug=false;
	Double_t lumi_=1;

    // input & output
	TString infile_;
	TString outdir_;
    TFile *outfile_;
};

}
#endif /* TTZANALYSIS_ANALYSIS_INTERFACE_STACKPLOTTER_H_ */
