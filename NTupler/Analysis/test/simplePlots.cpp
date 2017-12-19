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

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};

unsigned readFileList(const std::string & datFile,
		      TChain* GenJet,TChain* Event,
		      TChain* Particle,TChain* GenPhoton,TChain* Vertex,
		      TChain* ElectronLoose,TChain* ElectronTight,
		      TChain* MuonLoose,TChain* MuonTight,
		      TChain* Jets,TChain* MissingET,
		      TChain* PhotonLoose,TChain* PhotonTight){
  std::ifstream inputdat;
  inputdat.open(datFile);
  if (!inputdat.is_open()) {
    std::cout << " ** ERROR, cannot open dat file " << datFile << std::endl;
    return 1;
  }
  //check number of files to read
  unsigned nFiles = 0;
  std::string lBuf;
  while (std::getline(inputdat, lBuf,'\n'))
    {
      if (lBuf.empty()) continue;
      nFiles++;
      
      TFile * inFile = 0;
      lBuf = "root://gfe02.grid.hep.ph.ic.ac.uk:1097/"+lBuf;
      if (!testInputFile(lBuf,inFile)) return 2;
      GenJet->AddFile(lBuf.c_str());
      Event->AddFile(lBuf.c_str());
      Particle->AddFile(lBuf.c_str());
      GenPhoton->AddFile(lBuf.c_str());
      Vertex->AddFile(lBuf.c_str());
      ElectronLoose->AddFile(lBuf.c_str());
      ElectronTight->AddFile(lBuf.c_str());
      MuonLoose->AddFile(lBuf.c_str());
      MuonTight->AddFile(lBuf.c_str());
      Jets->AddFile(lBuf.c_str());
      MissingET->AddFile(lBuf.c_str());
      PhotonLoose->AddFile(lBuf.c_str());
      PhotonTight->AddFile(lBuf.c_str());
    }
  std::cout << " -- Read " << nFiles << " input files. " << std::endl;
  inputdat.close();
  return 0;
};

bool checkEventSize(TChain* tree, const unsigned nEvts){
  if (tree->GetEntries()!=nEvts) 
    {
      std::cout << " -- Problem with tree " << tree->GetName() 
		<< " entries " <<tree->GetEntries() 
		<< " nevts = " << nEvts 
		<< std::endl;
      return false;
    }
  return true;
}

void plotVar(TCanvas *myc, TH1F *hist, 
	     const double & scale, 
	     const bool first,
	     const std::string & label,
	     const  std::string & file,
	     const std::string & plotDir){

  unsigned iP = first?0:1;
  myc->cd();
  gStyle->SetOptStat(0);
  hist->Sumw2();
  hist->Scale(scale);
  hist->SetLineColor(iP+1);
  if (first) hist->SetFillColor(5);
  hist->SetMarkerColor(iP+1);
  hist->SetMarkerStyle(iP+20);
  hist->Draw(first?"":"PEsame");
  myc->Update();
  if (!first) myc->Print((plotDir+label+"_"+file+".pdf").c_str());
  

}

void plotJES(TCanvas *myc, TH1F *hist, 
	     const double & scale, 
	     //const std::string & label,
	     //const  std::string & file,
	     //const std::string & plotDir,
	     double & mean,double & meanerr,
	     double & sigma, double & sigmaerr){

  myc->cd();
  gStyle->SetOptStat(0);//"eMRuo");
  gStyle->SetOptFit(0);//1111);
  hist->Sumw2();
  hist->Scale(scale);
  hist->SetLineColor(1);
  hist->SetMarkerColor(1);
  hist->SetMarkerStyle(2);
  hist->Draw("PE");
  hist->Fit("gaus","","",hist->GetMean()-2*hist->GetRMS(),hist->GetMean()+2*hist->GetRMS());
  TF1 *fit = (TF1*)hist->GetFunction("gaus");
  if (!fit) {
    mean = hist->GetMean();
    meanerr = hist->GetMeanError();
    sigma = hist->GetRMS();
    sigmaerr = hist->GetRMSError();
  } else {
    mean = fit->GetParameter(1);
    meanerr = fit->GetParError(1);
    sigma = fit->GetParameter(2);
    sigmaerr = fit->GetParError(2);
  }
  myc->Update();
  //myc->Print((plotDir+label+"_"+file+".pdf").c_str());
  

}


int makePlots(const std::string & aProcess, 
	      const bool doJES=false){//main


  const unsigned nPU = 2;

  std::string plotDir = "PLOTS_looseVBFsel/";

  std::string baseDir = "filelists";

  std::string pu[nPU] = {"noPU","200PU"};

  TFile *outFile = TFile::Open((plotDir+"HistosFile_"+aProcess+".root").c_str(),"RECREATE");

  const unsigned nC = 29;
  TCanvas *myc[nC];
  for (unsigned iC(0); iC<nC; ++iC){
    std::ostringstream label;
    label << "myc" << iC;
    myc[iC] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }
  outFile->cd();

  TH1F *hGenJet1_pt[nPU];
  TH1F *hGenJet2_pt[nPU];
  TH1F *hGenJet1_eta[nPU];
  TH1F *hGenJet2_eta[nPU];
  TH1F *hGen_Mjj[nPU];
  TH1F *hGen_detajj[nPU];
  TH1F *hGen_dphijj[nPU];

  TH1F *hnJets[nPU];

  TH1F *hJet1_pt[nPU];
  TH1F *hJet2_pt[nPU];
  TH1F *hJet1_eta[nPU];
  TH1F *hJet2_eta[nPU];
  TH1F *hJet1_ID[nPU];
  TH1F *hJet2_ID[nPU];
  TH1F *hJet1_genidx[nPU];
  TH1F *hJet2_genidx[nPU];
  TH1F *hJet1_parton[nPU];
  TH1F *hJet2_parton[nPU];
  TH1F *hJet1_deepcsv[nPU];
  TH1F *hJet2_deepcsv[nPU];
  TH1F *hMjj[nPU];
  TH1F *hdetajj[nPU];
  TH1F *hdphijj[nPU];

  TH1F *hmet[nPU];
  TH1F *hjetmetmindphi[nPU];

  //jes validation
  const unsigned neta = 10;
  const unsigned npt = 6;
  double etamin=0;
  double deta=0.5;
  double ptval[npt+1] = {10,30,50,80,120,200,14000};

  TH1F *hJet_checkDR[nPU];
  TH1F *hJet_JESall[nPU];
  TH1F *hJet_JES[nPU][neta][npt];

  TLatex lat;
  char buf[200];

  std::cout << " .. Processing file " << aProcess << std::endl;
  int nEvtsRef = 0;
  for (unsigned iP(0); iP<nPU; ++iP){
    std::cout << " ... Processing " <<  pu[iP] << std::endl;
    outFile->mkdir((aProcess+"/"+pu[iP]).c_str());
    outFile->cd((aProcess+"/"+pu[iP]).c_str());
    std::ostringstream hLabel;
    hLabel.str("");
    hLabel << "hGenJet1_pt_" << aProcess << "_" << pu[iP];
    hGenJet1_pt[iP] = new TH1F(hLabel.str().c_str(),";p_{T}^{genjet1} (GeV); Events",100,0,400);
    hLabel.str("");
    hLabel << "hGenJet1_eta_" << aProcess << "_" << pu[iP];
    hGenJet1_eta[iP] = new TH1F(hLabel.str().c_str(),";#eta^{genjet1}; Events",100,-5,5);
    hLabel.str("");
    hLabel << "hGenJet2_pt_" << aProcess << "_" << pu[iP];
    hGenJet2_pt[iP] = new TH1F(hLabel.str().c_str(),";p_{T}^{genjet2} (GeV); Events",100,0,400);
    hLabel.str("");
    hLabel << "hGenJet2_eta_" << aProcess << "_" << pu[iP];
    hGenJet2_eta[iP] = new TH1F(hLabel.str().c_str(),";#eta^{genjet2}; Events",100,-5,5);
    hLabel.str("");
    hLabel << "hGen_Mjj_" << aProcess << "_" << pu[iP];
    hGen_Mjj[iP] = new TH1F(hLabel.str().c_str(),";M_{jj}^{gen} (GeV); Events",100,0,5000);
    hLabel.str("");
    hLabel << "hGen_detajj_" << aProcess << "_" << pu[iP];
    hGen_detajj[iP] = new TH1F(hLabel.str().c_str(),";#Delta#eta_{jj}^{gen}; Events",100,0,10);
    hLabel.str("");
    hLabel << "hGen_dphijj_" << aProcess << "_" << pu[iP];
    hGen_dphijj[iP] = new TH1F(hLabel.str().c_str(),";#Delta#phi_{jj}^{gen}; Events",100,0,3.1416);

    hLabel.str("");
    hLabel << "hnJets_" << aProcess << "_" << pu[iP];
    hnJets[iP] = new TH1F(hLabel.str().c_str(),";n_{jets}; Events",50,0,50);
    hLabel.str("");
    hLabel << "hJet1_pt_" << aProcess << "_" << pu[iP];
    hJet1_pt[iP] = new TH1F(hLabel.str().c_str(),";p_{T}^{jet1} (GeV); Events",100,0,400);
    hLabel.str("");
    hLabel << "hJet1_eta_" << aProcess << "_" << pu[iP];
    hJet1_eta[iP] = new TH1F(hLabel.str().c_str(),";#eta^{jet1}; Events",100,-5,5);
    hLabel.str("");
    hLabel << "hJet1_ID_" << aProcess << "_" << pu[iP];
    hJet1_ID[iP] = new TH1F(hLabel.str().c_str(),";ID^{jet1}; Events",4,0,4);
    hLabel.str("");
    hLabel << "hJet1_genidx_" << aProcess << "_" << pu[iP];
    hJet1_genidx[iP] = new TH1F(hLabel.str().c_str(),";genIdx^{jet1}; Events",21,-1,20);
    hLabel.str("");
    hLabel << "hJet1_parton_" << aProcess << "_" << pu[iP];
    hJet1_parton[iP] = new TH1F(hLabel.str().c_str(),";jet 1 parton flavour; Events",28,-6,22);
    hLabel.str("");
    hLabel << "hJet1_deepcsv_" << aProcess << "_" << pu[iP];
    hJet1_deepcsv[iP] = new TH1F(hLabel.str().c_str(),";deepCSV^{jet1}; Events",10,0,10);

    hLabel.str("");
    hLabel << "hJet2_pt_" << aProcess << "_" << pu[iP];
    hJet2_pt[iP] = new TH1F(hLabel.str().c_str(),";p_{T}^{jet2} (GeV); Events",100,0,400);
    hLabel.str("");
    hLabel << "hJet2_eta_" << aProcess << "_" << pu[iP];
    hJet2_eta[iP] = new TH1F(hLabel.str().c_str(),";#eta^{jet2}; Events",100,-5,5);
    hLabel.str("");
    hLabel << "hJet2_ID_" << aProcess << "_" << pu[iP];
    hJet2_ID[iP] = new TH1F(hLabel.str().c_str(),";ID^{jet2}; Events",4,0,4);
    hLabel.str("");
    hLabel << "hJet2_genidx_" << aProcess << "_" << pu[iP];
    hJet2_genidx[iP] = new TH1F(hLabel.str().c_str(),";genIdx^{jet2}; Events",21,-1,20);
    hLabel.str("");
    hLabel << "hJet2_parton_" << aProcess << "_" << pu[iP];
    hJet2_parton[iP] = new TH1F(hLabel.str().c_str(),";jet 2 parton flavour; Events",28,-6,22);
    hLabel.str("");
    hLabel << "hJet2_deepcsv_" << aProcess << "_" << pu[iP];
    hJet2_deepcsv[iP] = new TH1F(hLabel.str().c_str(),";deepCSV^{jet2}; Events",10,0,10);
    hLabel.str("");
    hLabel << "hMjj_" << aProcess << "_" << pu[iP];
    hMjj[iP] = new TH1F(hLabel.str().c_str(),";M_{jj} (GeV); Events",100,0,5000);
    hLabel.str("");
    hLabel << "hdetajj_" << aProcess << "_" << pu[iP];
    hdetajj[iP] = new TH1F(hLabel.str().c_str(),";#Delta#eta_{jj}; Events",100,0,10);
    hLabel.str("");
    hLabel << "hdphijj_" << aProcess << "_" << pu[iP];
    hdphijj[iP] = new TH1F(hLabel.str().c_str(),";#Delta#phi_{jj}; Events",100,0,3.1416);


    hLabel.str("");
    hLabel << "hmet_" << aProcess << "_" << pu[iP];
    hmet[iP] = new TH1F(hLabel.str().c_str(),";MET (GeV); Events",100,0,500);
    hLabel.str("");
    hLabel << "hjetmetmindphi_" << aProcess << "_" << pu[iP];
    hjetmetmindphi[iP] = new TH1F(hLabel.str().c_str(),";min#Delta#phi(jet,MET); Events",100,0,3.1416);

    if (doJES){
      gDirectory->mkdir("JES");
      gDirectory->cd("JES");
      //outFile->mkdir((aProcess+"/JES").c_str());
      //outFile->cd((aProcess+"/JES").c_str());
	
      hLabel.str("");
      hLabel << "hJet_checkDR_" << aProcess << "_" << pu[iP];
      hJet_checkDR[iP] = new TH1F(hLabel.str().c_str(),";#DeltaR(genjet,PUPPiJet); Events",100,0,5);
      hJet_checkDR[iP]->SetStats(1);
      hLabel.str("");
      hLabel << "hJet_JESall_" << aProcess << "_" << pu[iP];
      hJet_JESall[iP] = new TH1F(hLabel.str().c_str(),";p_{T}^{reco}/p_{T}^{gen}; Events",150,0,3);
      hJet_JESall[iP]->SetStats(1);
      for (unsigned ieta(0); ieta<neta; ++ieta){
	for (unsigned ipt(0); ipt<npt; ++ipt){
	  hLabel.str("");
	  hLabel << "hJet_JES_" << aProcess << "_" << pu[iP] << "_" << ieta << "_" << ipt;
	  hJet_JES[iP][ieta][ipt] = new TH1F(hLabel.str().c_str(),";p_{T}^{reco}/p_{T}^{gen}; Events",300,0,6);
	  hJet_JES[iP][ieta][ipt]->SetStats(1);
	}
      }
    }

    std::ostringstream filePath;
    filePath << baseDir << "/" << aProcess << "_" << pu[iP] << ".dat";
    TChain* GenJet = new TChain("ntuple/GenJet");
    TChain* Event = new TChain("ntuple/Event");
    TChain* Particle = new TChain("ntuple/Particle");
    TChain* GenPhoton = new TChain("ntuple/GenPhoton");
    TChain* Vertex = new TChain("ntuple/Vertex");
    TChain* ElectronLoose = new TChain("ntuple/ElectronLoose");
    TChain* ElectronTight = new TChain("ntuple/ElectronTight");
    TChain* MuonLoose = new TChain("ntuple/MuonLoose");
    TChain* MuonTight = new TChain("ntuple/MuonTight");
    TChain* Jets = new TChain(iP==0?"ntuple/Jet":"ntuple/JetPUPPI");
    TChain* MissingET = new TChain(iP==0?"ntuple/MissingET":"ntuple/PuppiMissingET");
    TChain* PhotonLoose = new TChain("ntuple/PhotonLoose");
    TChain* PhotonTight = new TChain("ntuple/PhotonTight");

    if (readFileList(filePath.str(),GenJet,
		     Event,Particle,GenPhoton,
		     Vertex,ElectronLoose,ElectronTight,
		     MuonLoose,MuonTight,
		     Jets,MissingET,
		     PhotonLoose,PhotonTight)!=0) return 1;

    int nGenJets = 0;
    float genjet_pt[100];
    float genjet_eta[100];
    float genjet_phi[100];
    float genjet_mass[100];
    GenJet->SetBranchAddress("GenJet_size",&nGenJets);
    GenJet->SetBranchAddress("PT",&genjet_pt);
    GenJet->SetBranchAddress("Eta",&genjet_eta);
    GenJet->SetBranchAddress("Phi",&genjet_phi);
    GenJet->SetBranchAddress("Mass",&genjet_mass);

    int nJets = 0;
    float jet_pt[100];
    float jet_eta[100];
    float jet_phi[100];
    float jet_mass[100];
    int jet_ID[100];
    int jet_genidx[100];
    int jet_parton[100];
    int jet_DeepCSV[100];
    Jets->SetBranchAddress("JetPUPPI_size",&nJets);
    Jets->SetBranchAddress("PT",&jet_pt);
    Jets->SetBranchAddress("Eta",&jet_eta);
    Jets->SetBranchAddress("Phi",&jet_phi);
    Jets->SetBranchAddress("Mass",&jet_mass);
    Jets->SetBranchAddress("ID",&jet_ID);
    Jets->SetBranchAddress("GenJet",&jet_genidx);
    Jets->SetBranchAddress("DeepCSV",&jet_DeepCSV);
    Jets->SetBranchAddress("PartonFlavor",&jet_parton);

    float met = 0;
    float metphi = 0;
    MissingET->SetBranchAddress("MET",&met);
    MissingET->SetBranchAddress("Phi",&metphi);

    int nEvts = GenJet->GetEntries();
    if (iP==0) nEvtsRef = nEvts;
    if (!checkEventSize(GenJet,nEvts)) return 1;
    if (!checkEventSize(Particle,nEvts)) return 1;
    if (!checkEventSize(GenPhoton,nEvts)) return 1;
    if (!checkEventSize(Vertex,nEvts)) return 1;
    if (!checkEventSize(ElectronLoose,nEvts)) return 1;
    if (!checkEventSize(ElectronTight,nEvts)) return 1;
    if (!checkEventSize(MuonLoose,nEvts)) return 1;
    if (!checkEventSize(MuonTight,nEvts)) return 1;
    if (!checkEventSize(PhotonLoose,nEvts)) return 1;
    if (!checkEventSize(PhotonTight,nEvts)) return 1;
    if (!checkEventSize(Jets,nEvts)) return 1;
    if (!checkEventSize(MissingET,nEvts)) return 1;

    unsigned nDiffIdx = 0;

    //event loop
    for (int ievt(0); ievt<nEvts; ++ievt){
      if (ievt%1000==0) std::cout << ".... Processing entry " << ievt << std::endl;
      Event->GetEntry(ievt);
      GenJet->GetEntry(ievt);
      Jets->GetEntry(ievt);
      MissingET->GetEntry(ievt);

      //no selection for JES
      if (doJES){
	for (unsigned ij(0); ij<abs(nJets);++ij){
	  //redo jet-genjet matching
	  TLorentzVector rec;
	  rec.SetPtEtaPhiM(jet_pt[ij],jet_eta[ij],jet_phi[ij],jet_mass[ij]);
	  double mindr = 1000;
	  int saveGenIdx = jet_genidx[ij];
	  int newGenIdx = -1;
	  for (unsigned igj(0); igj<abs(nGenJets);++igj){
	    TLorentzVector gen;
	    gen.SetPtEtaPhiM(genjet_pt[igj],genjet_eta[igj],genjet_phi[igj],genjet_mass[igj]);
	    double dR = rec.DeltaR(gen);
	    if (dR<mindr) {
	      mindr=dR;
	      newGenIdx = igj;
	    }
	  }
	  if (newGenIdx!=saveGenIdx && saveGenIdx>=0){
	    //std::cout << " *** Found different index ! nGenJets=" << nGenJets << " nJets=" << nJets << " jet " << ij << " original index " << saveGenIdx << " new index " << newGenIdx << std::endl;
	    nDiffIdx++;
	    jet_genidx[ij] = newGenIdx;
	  }
	  if (jet_genidx[ij]<0) continue;
	  unsigned genidx = abs(jet_genidx[ij]);
	  hJet_checkDR[iP]->Fill(mindr);
	  if (mindr>0.4) continue;
	  hJet_JESall[iP]->Fill(jet_pt[ij]/genjet_pt[genidx]);
	  for (unsigned ieta(0); ieta<neta; ++ieta){
	    if (fabs(genjet_eta[genidx])<etamin+ieta*deta || 
		fabs(genjet_eta[genidx])>=etamin+(ieta+1)*deta) continue;
	    for (unsigned ipt(0); ipt<npt; ++ipt){
	      if (genjet_pt[genidx]<ptval[ipt] ||
		  genjet_pt[genidx]>=ptval[ipt+1]) continue;
	      hJet_JES[iP][ieta][ipt]->Fill(jet_pt[ij]/genjet_pt[genidx]);
	    }
	  }
	}
      }
      //apply some selection
      if (nGenJets <2) continue;
      if (genjet_pt[0]<30 || genjet_pt[1]<30 || fabs(genjet_eta[0])>5.0 || fabs(genjet_eta[1])>5.0) continue;

      TLorentzVector gen1;
      gen1.SetPtEtaPhiM(genjet_pt[0],genjet_eta[0],genjet_phi[0],genjet_mass[0]);
      TLorentzVector gen2;
      gen2.SetPtEtaPhiM(genjet_pt[1],genjet_eta[1],genjet_phi[1],genjet_mass[1]);
      TLorentzVector genPair = gen1 + gen2;

      double detajj = fabs(genjet_eta[0]-genjet_eta[1]);
      double dphijj = fabs(gen1.DeltaPhi(gen2));

      //apply loose VBF sel, gen level...
      if (genPair.M()<500 || detajj<1 || dphijj>2) continue;
      //apply tight VBF sel, gen level...
      //if (genPair.M()<1300 || detajj<4 || dphijj>2) continue;

      //gen jets info
      hGenJet1_pt[iP]->Fill(genjet_pt[0]);
      hGenJet1_eta[iP]->Fill(genjet_eta[0]);
      hGenJet2_pt[iP]->Fill(genjet_pt[1]);
      hGenJet2_eta[iP]->Fill(genjet_eta[1]);

      hGen_Mjj[iP]->Fill(genPair.M());
      hGen_detajj[iP]->Fill(detajj);
      hGen_dphijj[iP]->Fill(dphijj);

      //reco jets info
      if (nJets<2) continue;
      if (jet_pt[0]<40 || jet_pt[1]<40 || fabs(jet_eta[0])>4.7 || fabs(jet_eta[1])>4.7) continue;


      hJet1_pt[iP]->Fill(jet_pt[0]);
      hJet1_eta[iP]->Fill(jet_eta[0]);
      hJet1_ID[iP]->Fill(jet_ID[0]);
      hJet1_genidx[iP]->Fill(jet_genidx[0]);
      hJet1_parton[iP]->Fill(jet_parton[0]);
      hJet1_deepcsv[iP]->Fill(jet_DeepCSV[0]);
      hJet2_pt[iP]->Fill(jet_pt[1]);
      hJet2_eta[iP]->Fill(jet_eta[1]);
      hJet2_ID[iP]->Fill(jet_ID[1]);
      hJet2_genidx[iP]->Fill(jet_genidx[1]);
      hJet2_parton[iP]->Fill(jet_parton[1]);
      hJet2_deepcsv[iP]->Fill(jet_DeepCSV[1]);


      TLorentzVector rec1;
      rec1.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
      TLorentzVector rec2;
      rec2.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_mass[1]);
      TLorentzVector recPair = rec1 + rec2;
      hMjj[iP]->Fill(recPair.M());
      hdetajj[iP]->Fill(fabs(jet_eta[0]-jet_eta[1]));
      hdphijj[iP]->Fill(fabs(rec1.DeltaPhi(rec2)));

      //MET info
      hmet[iP]->Fill(met);
      TLorentzVector metvec;
      metvec.SetPtEtaPhiE(met,0,metphi,met);
      double mindphi=10;
      unsigned njets30 = 0;
      for (unsigned ij(0); ij<abs(nJets);++ij){
	//std::cout << ij << " pt " << jet_pt[ij] << std::endl;
	if (jet_pt[ij]<30) continue;
	njets30++;
	if (ij>4) continue;

	TLorentzVector recij;
	recij.SetPtEtaPhiM(jet_pt[ij],jet_eta[ij],jet_phi[ij],jet_mass[ij]);
	double dphi = fabs(recij.DeltaPhi(metvec));
	if (dphi<mindphi) mindphi = dphi;
      }
	
      hnJets[iP]->Fill(njets30);
      hjetmetmindphi[iP]->Fill(mindphi);
      
    }//event loop

    std::cout << " ** Number of jets with different gen index: " << nDiffIdx
	      << std::endl;


    double scale = 1.*nEvtsRef/nEvts;

    plotVar(myc[0],hGenJet1_pt[iP],scale,iP==0,"GenJet1_pt",aProcess,plotDir);
    plotVar(myc[1],hGenJet1_eta[iP],scale,iP==0,"GenJet1_eta",aProcess,plotDir);
    plotVar(myc[2],hGenJet2_pt[iP],scale,iP==0,"GenJet2_pt",aProcess,plotDir);
    plotVar(myc[3],hGenJet2_eta[iP],scale,iP==0,"GenJet2_eta",aProcess,plotDir);

    plotVar(myc[4],hGen_Mjj[iP],scale,iP==0,"Gen_Mjj",aProcess,plotDir);
    plotVar(myc[5],hGen_detajj[iP],scale,iP==0,"Gen_detajj",aProcess,plotDir);
    plotVar(myc[6],hGen_dphijj[iP],scale,iP==0,"Gen_dphijj",aProcess,plotDir);

    plotVar(myc[7],hnJets[iP],scale,iP==0,"nJets",aProcess,plotDir);

    plotVar(myc[8],hJet1_pt[iP],scale,iP==0,"Jet1_pt",aProcess,plotDir);
    plotVar(myc[9],hJet1_eta[iP],scale,iP==0,"Jet1_eta",aProcess,plotDir);
    plotVar(myc[10],hJet1_ID[iP],scale,iP==0,"Jet1_ID",aProcess,plotDir);
    plotVar(myc[11],hJet1_genidx[iP],scale,iP==0,"Jet1_genidx",aProcess,plotDir);
    plotVar(myc[12],hJet1_parton[iP],scale,iP==0,"Jet1_parton",aProcess,plotDir);
    plotVar(myc[13],hJet1_deepcsv[iP],scale,iP==0,"Jet1_deepcsv",aProcess,plotDir);
    plotVar(myc[14],hJet2_pt[iP],scale,iP==0,"Jet2_pt",aProcess,plotDir);
    plotVar(myc[15],hJet2_eta[iP],scale,iP==0,"Jet2_eta",aProcess,plotDir);
    plotVar(myc[16],hJet2_ID[iP],scale,iP==0,"Jet2_ID",aProcess,plotDir);
    plotVar(myc[17],hJet2_genidx[iP],scale,iP==0,"Jet2_genidx",aProcess,plotDir);
    plotVar(myc[18],hJet2_parton[iP],scale,iP==0,"Jet2_parton",aProcess,plotDir);
    plotVar(myc[19],hJet2_deepcsv[iP],scale,iP==0,"Jet2_deepcsv",aProcess,plotDir);


    plotVar(myc[20],hMjj[iP],scale,iP==0,"Mjj",aProcess,plotDir);
    plotVar(myc[21],hdetajj[iP],scale,iP==0,"detajj",aProcess,plotDir);
    plotVar(myc[22],hdphijj[iP],scale,iP==0,"dphijj",aProcess,plotDir);

    plotVar(myc[23],hmet[iP],scale,iP==0,"met",aProcess,plotDir);
    plotVar(myc[24],hjetmetmindphi[iP],scale,iP==0,"jetmetmindphi",aProcess,plotDir);

    if (doJES){
      double mean[npt][neta];
      double meanerr[npt][neta];
      double sigma[npt][neta];
      double sigmaerr[npt][neta];
      double res[npt][neta];
      double reserr[npt][neta];
      double etaval[neta];
      double etaerr[neta];
      myc[25+2*iP]->Print((plotDir+"JES_"+aProcess+pu[iP]+".pdf[").c_str());

      for (unsigned ieta(0); ieta<neta; ++ieta){
	etaval[ieta] = etamin +ieta*deta + deta/2.;
	etaerr[ieta] = deta/2.;
	for (unsigned ipt(0); ipt<npt; ++ipt){
	  //plotJES(myc[25+2*iP],hJet_JES[iP][ieta][ipt],1,"JES",(aProcess+pu[iP]).c_str(),plotDir,mean[ipt][ieta],meanerr[ipt][ieta],sigma[ipt][ieta],sigmaerr[ipt][ieta]);
	  plotJES(myc[25+2*iP],hJet_JES[iP][ieta][ipt],1,mean[ipt][ieta],meanerr[ipt][ieta],sigma[ipt][ieta],sigmaerr[ipt][ieta]);
	  myc[25+2*iP]->cd();
	  sprintf(buf,"%3.1f < |#eta| < %3.1f",etaval[ieta]-etaerr[ieta],etaval[ieta]+etaerr[ieta]);
	  lat.DrawLatexNDC(0.2,0.85,buf);
	  sprintf(buf,"%3.0f < p_{T} < %3.0f GeV",ptval[ipt],ptval[ipt+1]);
	  lat.DrawLatexNDC(0.2,0.8,buf);
	  myc[25+2*iP]->Update();
	  myc[25+2*iP]->Print((plotDir+"JES_"+aProcess+pu[iP]+".pdf").c_str());
	  res[ipt][ieta] = sigma[ipt][ieta]/mean[ipt][ieta];
	  reserr[ipt][ieta] = sigmaerr[ipt][ieta]/mean[ipt][ieta];

	}
      }
      myc[25+2*iP]->Print((plotDir+"JES_"+aProcess+pu[iP]+".pdf]").c_str());
      myc[25+2*iP]->Clear();
      std::ostringstream label;
      for (unsigned ipt(0); ipt<npt; ++ipt){
	myc[25+2*iP]->cd();
	TGraphErrors *grMean = new TGraphErrors(neta,etaval,mean[ipt],etaerr,meanerr[ipt]);
	label.str("");
	label << "MeanvsEta_" << aProcess << "_" << pu[iP] << "_pt" << ptval[ipt] << "_" << ptval[ipt+1];
	grMean->SetName(label.str().c_str());
	grMean->SetTitle(";|#eta|;<p_{T}^{rec}/p_{T}^{gen}>");
	grMean->SetMarkerColor(ipt+1);
	grMean->SetLineColor(ipt+1);
	grMean->SetMarkerStyle(ipt+20);
	grMean->SetMinimum(0);
	grMean->SetMaximum(2);
	grMean->Draw(ipt==0?"APL":"PLsame");

	myc[26+2*iP]->cd();
	TGraphErrors *grRes = new TGraphErrors(neta,etaval,res[ipt],etaerr,reserr[ipt]);
	label.str("");
	label << "ResovsEta_" << aProcess << "_" << pu[iP] << "_pt" << ptval[ipt] << "_" << ptval[ipt+1];
	grRes->SetName(label.str().c_str());
	grRes->SetTitle(";|#eta|;#sigma/mean");
	grRes->SetMarkerColor(ipt+1);
	grRes->SetLineColor(ipt+1);
	grRes->SetMarkerStyle(ipt+20);
	grRes->SetMinimum(0);
	grRes->SetMaximum(0.5);
	grRes->Draw(ipt==0?"APL":"PLsame");

	outFile->cd();
	grRes->Write();
      }

      myc[25+2*iP]->Update();
      label.str("");
      label << "MeanvsEta_" << aProcess << "_" << pu[iP];
      myc[25+2*iP]->Print((plotDir+label.str()+".pdf").c_str());

      myc[26+2*iP]->Update();
      label.str("");
      label << "ResovsEta_" << aProcess << "_" << pu[iP];
      myc[26+2*iP]->Print((plotDir+label.str()+".pdf").c_str());

    }//doJES
  }//PU loop


  outFile->Write();

  return 0;

}//process



int main(int argc, char** argv){//main

  if (argc!=2) {
    std::cout << " Usage: "
	      << argv[0] << " <process short name> "
	      << " (optional: <doJES>, default is false)"
	      << std::endl;
    return 1;
  }
  std::string lProcess = argv[1];
  bool doJES = false;
  if (argc>2) std::istringstream(argv[2])>>doJES;

  makePlots(lProcess,doJES);

  return 0;

}//main
