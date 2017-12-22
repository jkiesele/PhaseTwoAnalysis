#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"

void fillROC(TGraphErrors *gr, TH1F **iso, TCanvas *myc){
  gr->SetTitle(";signal eff;(1-bkg eff)");
  //gr->SetMinimum(0);
  gr->SetMaximum(1);


  const unsigned nP = 100;
  const double step = 0.01;

  double totSig = iso[0]->Integral(0,iso[0]->GetNbinsX()+1);
  double totBkg = iso[1]->Integral(0,iso[1]->GetNbinsX()+1);

  for (unsigned iP(0); iP<nP; ++iP){
    //if (iP%10!=0) continue;
    double val = step*iP;//+step/2.;
    double sig,sigerr,bkg,bkgerr;
    int binS = iso[0]->GetXaxis()->FindBin(val+0.001);
    sig = iso[0]->Integral(0,binS)/totSig;
    bkg = iso[1]->Integral(0,binS)/totBkg;
    sigerr = sqrt(sig*(1-sig)/totSig);
    bkgerr = sqrt(bkg*(1-bkg)/totBkg);
    gr->SetPoint(iP,sig,1-bkg);
    gr->SetPointError(iP,sigerr,bkgerr);
    //gr->SetPointError(iP,0,0);
    if (iP%10==0){
      char buf[100];
      sprintf(buf,"iso=%3.1f",val);
      TLatex *latex = new TLatex(gr->GetX()[iP]+0.01, gr->GetY()[iP],buf);
      gr->GetListOfFunctions()->Add(latex);
    }

    if (fabs(sig-0.95)<0.002) std::cout << " SIG 95% sig eff = " << sig << " bkg = " << bkg << " isolation = " << val << std::endl;
    if (fabs(sig-0.98)<0.002) std::cout << " SIG 98% sig eff = " << sig << " bkg = " << bkg << " isolation = " << val << std::endl;
    if (fabs(bkg-0.3)<0.01) std::cout << " BKG 70% sig eff = " << sig << " bkg = " << bkg << " isolation = " << val << std::endl;

  }
  myc->cd();
  gr->Draw("AP");
  gr->GetXaxis()->SetNdivisions(520);
  myc->Update();

};

int rocLepton(){//main

  const std::string treename = "MuonTight";

  const unsigned nF = 2;
  TFile *fin[nF];
  fin[0] = TFile::Open("../rootfiles/DYToLL-M-50_0J_noPU.root");
  fin[1] = TFile::Open("../rootfiles/VBFH_noPU.root");
  //fin[2] = TFile::Open("../rootfiles/QCD.root");

  std::string label[nF] = {"DY","VBFH"};//,"QCD"};

  TCanvas *mycEB = new TCanvas("mycEB","mycEB",1);
  TCanvas *mycEE = new TCanvas("mycEE","mycEE",1);
  TCanvas *mycrocEB = new TCanvas("mycrocEB","mycrocEB",1);
  TCanvas *mycrocEE = new TCanvas("mycrocEE","mycrocEE",1);

  mycEB->cd();
  TH1F *isoEB[nF];
  TH1F *isoEE[nF];

  for (unsigned iF(0); iF<nF; ++iF){//loop on files
    fin[iF]->cd("ntuple");
    TTree *tree = (TTree*)gDirectory->Get(treename.c_str());
    isoEB[iF] = new TH1F(("isoEB_"+label[iF]).c_str(),";EB isolation;leptons",100,0,1);
    isoEB[iF]->Sumw2();
    tree->Draw(("IsolationVar>>isoEB_"+label[iF]).c_str(),iF==0?"Particle>=0 && TMath::Abs(Eta)<1.5 && PT>10":"TMath::Abs(Eta)<1.5 && PT>10");
    isoEB[iF]->SetLineColor(iF+1);
    isoEB[iF]->SetMarkerColor(iF+1);
    isoEB[iF]->SetMarkerStyle(iF+20);

    isoEE[iF] = new TH1F(("isoEE_"+label[iF]).c_str(),";EE isolation;leptons",100,0,1);
    isoEE[iF]->Sumw2();
    tree->Draw(("IsolationVar>>isoEE_"+label[iF]).c_str(),iF==0?"Particle>=0 && TMath::Abs(Eta)>1.5 && PT>10":"TMath::Abs(Eta)>1.5 && PT>10");
    isoEE[iF]->SetLineColor(iF+1);
    isoEE[iF]->SetMarkerColor(iF+1);
    isoEE[iF]->SetMarkerStyle(iF+20);
  }//loop on files


  TGraphErrors *grEB = new TGraphErrors();
  grEB->SetName("ROCEB");

  TGraphErrors *grEE = new TGraphErrors();
  grEE->SetName("ROCEE");

  std::cout << " -- EB -- " << std::endl;
  fillROC(grEB,isoEB,mycrocEB);
  mycrocEB->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  TLatex lat;
  lat.DrawLatexNDC(0.2,0.25,(treename+" EB").c_str());
  mycrocEB->Update();
  mycrocEB->Print(("ROC_EB_"+treename+".pdf").c_str());

  std::cout << " -- EE -- " << std::endl;
  fillROC(grEE,isoEE,mycrocEE);
  mycrocEE->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  lat.DrawLatexNDC(0.2,0.25,(treename+" EE").c_str());
  mycrocEE->Update();
  mycrocEE->Print(("ROC_EE_"+treename+".pdf").c_str());


  mycEB->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(0);
  isoEB[0]->Scale(1./isoEB[0]->GetEntries());
  isoEB[1]->Scale(1./isoEB[1]->GetEntries());
  isoEB[0]->Draw("PE");
  isoEB[1]->Draw("PEsame");
  lat.DrawLatexNDC(0.5,0.85,(treename+" EB").c_str());
  mycEB->Update();
  mycEB->Print(("Isolation_EB_"+treename+".pdf").c_str());

  mycEE->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(0);
  isoEE[0]->Scale(1./isoEE[0]->GetEntries());
  isoEE[1]->Scale(1./isoEE[1]->GetEntries());
  isoEE[0]->Draw("PE");
  isoEE[1]->Draw("PEsame");
  lat.DrawLatexNDC(0.5,0.85,(treename+" EE").c_str());
  mycEE->Update();
  mycEE->Print(("Isolation_EE_"+treename+".pdf").c_str());


  return 0;
};//main
