/*
 * stackPlotter.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: eacoleman 
 */

#include "../interface/stackPlotter.h"
#include "TStyle.h"

namespace d_ana{

stackPlotter::stackPlotter()
{}

stackPlotter::~stackPlotter(){
	if(outfile_) {
		outfile_->Close();
		//delete outfile_;
	}
}


void stackPlotter::moveDirHistsToStacks(TDirectory* tdir){
	if(debug)
		std::cout << "stackPlotter::moveDirHistsToStacks" << std::endl;

	// get metainfo from directory, else exit TODO
	metaInfo tMI;
	tMI.extractFrom(tdir);

	if(debug) {
		std::cout << "stackPlotter::moveDirHistsToStacks || metaInfo color=" << tMI.color << std::endl;
		std::cout << "stackPlotter::moveDirHistsToStacks || metaInfo legendname=" << tMI.legendname<< std::endl;
		std::cout << "stackPlotter::moveDirHistsToStacks || metaInfo legendorder=" << tMI.legendorder << std::endl;
	}


	TIter    histIter(tdir->GetListOfKeys());
	TObject* cHistObj;
	TKey*    cHistKey;

	if(debug)
		std::cout << "stackPlotter::moveDirHistsToStacks || Iterating through histograms." << std::endl;

	// loop through keys in the directory
	while((cHistKey = (TKey*) histIter())) {
		cHistObj=tdir->Get(cHistKey->GetName());
		if(!cHistObj->InheritsFrom(TH1::Class())) continue;

		if(debug)
			std::cout << "stackPlotter::moveDirHistsToStacks || Found histogram "
			<< cHistKey->GetName() << std::endl;

		// prepare the histogram to be added to the stack
		TH1* cHist = (TH1*) cHistObj->Clone();
		cHist->SetDirectory(0);
		TString mapName = cHist->GetName();

		std::pair<Int_t,TH1*> newEntry(tMI.legendorder,cHist);

		// initialize the stack info if needed
		if(!stacksLegEntries_.count(mapName)) {
			std::vector<std::pair<Int_t,TH1*> > legInfo(0);
			legInfo.push_back(newEntry);
			stacksLegEntries_[mapName] = legInfo;
		}

		cHist->SetFillColor(tMI.color);
		cHist->SetFillStyle(1001);
		cHist->SetMarkerStyle(kNone);
		cHist->SetMarkerColor(kBlack);
		cHist->SetLineColor(kBlack);
		cHist->SetTitle(mapName);
		cHist->SetName(tMI.legendname);

		std::vector<std::pair<Int_t,TH1*> > legEntries = stacksLegEntries_[mapName];
		if(debug)
			std::cout << "stackPlotter::moveDirHistsToStacks || legEntries size is " << legEntries.size() << std::endl;
		for(size_t i=0; i < legEntries.size(); i++) {
			if(legEntries.at(i).second == cHist && legEntries.at(i).first == tMI.legendorder) break;

			if(legEntries.at(i).first >= tMI.legendorder) {
				if(debug)
					std::cout << "stackPlotter::moveDirHistsToStacks || i is " << i << std::endl;
				stacksLegEntries_[mapName].insert(stacksLegEntries_[mapName].begin()+i,newEntry);
				break;
			}

			if(i==legEntries.size()-1) {
				stacksLegEntries_[mapName].push_back(newEntry);
				break;
			}
		}

		if(debug)
			std::cout << "stackPlotter::moveDirHistsToStacks || legEntries size is " << legEntries.size() << std::endl;
	}

}


void stackPlotter::plotStack(const TString& key) {
	if(debug)
		std::cout << "stackPlotter::plotStack" << std::endl;

	std::vector<std::pair<Int_t,TH1*> > legEntries = stacksLegEntries_[key];
	if(legEntries.size()<1) return;

	TCanvas *c = new TCanvas(key,key,800,600);
	TLegend *leg = new TLegend(0.75,0.75,0.95,0.95); //will be overwritten by style!
	THStack *stack = new THStack();

	//determine legend position TBI (add all histos, look at mean and max values etc)
	//float max=stack->GetMaximum(); //unfortunately no mean function also not filled, yet.
	//needs another loop on histos
	TH1 * h=(TH1 *) legEntries.at(0).second->Clone();
	for(size_t i=1;i<legEntries.size();i++){
		h->Add(legEntries.at(i).second);
	}
//	int maxbin=h->GetMaximumBin();
	//int minbin=h->GetMinimumBin();
//	float ymax=h->GetBinContent(maxbin);
	//float ymin=h->GetBinContent(minbin);
	//pass this to the style functions
	legendposition legpos=estimateBestLegendPosition(h);
	delete h;

	//gStyle->SetOptTitle(0);//no title


	for(size_t i=0; i < legEntries.size(); i++) {
		TString legName = legEntries.at(i).second->GetName();
		if(legName == "") continue;

		applyStyleToTH1(legEntries.at(i).second,legpos);

		stack->Add(legEntries.at(i).second,"HIST");
		stack->SetTitle("");//legEntries.at(i).second->GetTitle());
		stack->SetName(legEntries.at(i).second->GetTitle());

		//mirror entry for legend
		leg->AddEntry(legEntries.at(legEntries.size()-1-i).second,
				legEntries.at(legEntries.size()-1-i).second->GetName(),
				"F");
	}

	// draw plot and legend

	c->cd();

	applyStyleToCanvas(c,legpos);

	stack->Draw();
	stack->GetXaxis()->SetTitle(legEntries.at(0).second->GetXaxis()->GetTitle());
	stack->GetYaxis()->SetTitle(legEntries.at(0).second->GetYaxis()->GetTitle());

	applyStyleToAxis(stack,legpos);
	applyStyleToLegend(leg,legpos);
	stack->Draw();
	leg->Draw();


	// save and exit
	if(saveplots_) {
		c->SaveAs(outdir_+"/"+key+".pdf");
		c->SaveAs(outdir_+"/"+key+".png");
	}

	if(savecanvases_ && outfile_) {
		outfile_->cd();
		c->Write();
	}

	delete stack;
	delete leg;
	delete c;

}


void stackPlotter::plot() {
	if(debug)
		std::cout << "stackPlotter::plot" << std::endl;

	TH1::AddDirectory(kFALSE);
	TDirectory::AddDirectory(kFALSE);
	gROOT->SetBatch(true);
	TFile *fIn = new TFile(infile_,"READ");
	fIn->cd();

	if(debug)
		std::cout << "stackPlotter::plot || input file '" << infile_ << "' is being read..." << std::endl;

	TIter   dirIter(fIn->GetListOfKeys());
	TObject *cDirObj;
	TKey    *key;

	// iterate over directories and get all stacks
	while((key = (TKey *) dirIter())) {
		cDirObj=fIn->Get(key->GetName());
		if(!cDirObj->InheritsFrom(TDirectory::Class())) continue;

		TDirectory* cDir = (TDirectory*) cDirObj;
		if(debug)
			std::cout << "stackPlotter::plot || Moving histograms from directory " << cDir->GetName()
			<< " to relevant maps." << std::endl;

		moveDirHistsToStacks(cDir);

	}

	if(debug)
		std::cout << "stackPlotter::plot || input file '" << infile_ << "' has been read. Closing..." << std::endl;

	// intermediate cleanup
	fIn->Close();
	delete fIn;

	if(debug)
		std::cout << "stackPlotter::plot || Closed. Saving output..." << std::endl;

	// create the outfile if need be
	if(savecanvases_) {
		if(debug)
			std::cout << "stackPlotter::plot || Opening output ROOT file" << std::endl;
		TString writeOption = rewriteoutfile_ ? "RECREATE" : "UPDATE";
		outfile_ = new TFile(outdir_+"/plotter.root",writeOption);
	}

	// plot all the stacks & save appropriately
	if(debug)
		std::cout << "stackPlotter::plot || Plotting all the canvases" << std::endl;
	for(const auto& it : stacksLegEntries_) {
		plotStack(it.first);
	}

	// close, save, and cleanup
	if(savecanvases_ && outfile_) {
		if(debug)
			std::cout << "stackPlotter::plot || Closing the outfile" << std::endl;
		outfile_->Close();
	}

	if(debug)
		std::cout << "stackPlotter::plot || Done!" << std::endl;
}


stackPlotter::legendposition stackPlotter::estimateBestLegendPosition(TH1* h)const{

//	int maxbin=h->GetMaximumBin();
//	int minbin=h->GetMinimumBin();
	//int nbins=h->GetNbinsX();

	//float ymax=h->GetBinContent(maxbin);
	//float ymin=h->GetBinContent(minbin);
	float mean=h->GetMean();

	float xmin=h->GetBinLowEdge(1);
	float xmax=h->GetBinLowEdge(h->GetNbinsX()+1);


	float range=xmax-xmin;
	float middle=range/2+xmin;

	float relative=(mean - middle)/fabs(range);

	if(fabs(relative) < 0.02)
		return lp_top;

	if(relative>0)
		return lp_left;

	if(relative<0)
		return lp_right;

	else
		return lp_right;
}


/*
 * These can be overwritten by inheriting classes
 *
 */


void stackPlotter::applyStyleToAxis(THStack * h, legendposition pos)const{

	h->GetYaxis()->SetTitleSize(0.07);
	h->GetYaxis()->SetLabelSize(0.06);
	h->GetYaxis()->SetTitleOffset(1.25);

	h->GetXaxis()->SetTitleSize(0.07);
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleOffset(0.95);

	double max=h->GetMaximum();
	double min=h->GetMinimum();
	if(min>0 && min/(max-min)<0.1)min=0;
	else min/=yscaling;

	if(pos==lp_top)
		max*=yscaling_toplegend;
	else
		max*=yscaling;
	//why like this?
	h->SetMaximum(max);
	h->SetMinimum(min);
}

void stackPlotter::applyStyleToCanvas(TVirtualPad* c,legendposition pos)const{
	c->SetBottomMargin(0.16);
	c->SetLeftMargin(0.17);
}
void stackPlotter::applyStyleToTH1(TH1* h,legendposition pos)const{


}

void stackPlotter::applyStyleToLegend(TLegend* leg ,legendposition pos)const{
	leg->SetColumnSeparation(0.1);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);

	if(pos==lp_right){
		leg->SetX1(0.62);
		leg->SetX2(0.86);
		leg->SetY1(0.32);
		leg->SetY2(0.88);
		leg->SetNColumns(1);
	}
	else if(pos==lp_top){
		leg->SetX1(0.2);
		leg->SetX2(0.9);
		leg->SetY1(0.95/yscaling_toplegend);
		leg->SetY2(0.89);
		leg->SetNColumns(3);
	}
	else if (pos==lp_left){
		leg->SetX1(0.22);
		leg->SetX2(0.46);
		leg->SetY1(0.32);
		leg->SetY2(0.88);
		leg->SetNColumns(1);
	}
}


}



int main(int argc, const char** argv){

	if(argc < 3) {
		std::cout << "***** stackPlotter ***************************************" << std::endl;
		std::cout << "**                                                       *" << std::endl;
		std::cout << "** Usage: ./stackPlotter <input file> <output directory> *" << std::endl;
		std::cout << "**                                                       *" << std::endl;
		std::cout << "**********************************************************\n\n" << std::endl;
		std::cout << "Incorrect usage: number of arguments is " << argc << std::endl;
		exit(EXIT_FAILURE);
	}

	system((TString("mkdir -p ")+TString(argv[2])).Data());
	d_ana::stackPlotter sPlots;

	sPlots.rewriteOutfile(true);
	sPlots.savePlots(true);
	sPlots.saveCanvasRootFile(true);

	sPlots.setInputFile(argv[1]);
	sPlots.setOutDir(argv[2]);
	sPlots.setLumi(1);

	sPlots.plot();
}
