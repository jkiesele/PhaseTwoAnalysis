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

int makeTree(const std::string & plotDir,
	     const std::string & aProcess,
	     const std::string & pu){//main
  
  std::string baseDir = "filelists";
  
  std::string outname = plotDir+"/HistosFile_"+aProcess+"_"+pu+".root";
  TFile *outFile = TFile::Open(outname.c_str(),"RECREATE");

  if (!outFile) {
    std::cout << " ** ERROR, could not open output file " << outname << std::endl;
    return 1;
  }
  else {
    std::cout << " -- Output file " << outFile->GetName() << " successfully opened." << std::endl;
  }

  outFile->cd();
  TTree *outtree = new TTree("LightTree","Analysis light tree");
  double GenJet1_pt = 0;
  double GenJet2_pt = 0;
  double GenJet1_eta = 0;
  double GenJet2_eta = 0;
  double Gen_Mjj = 0;
  double Gen_detajj = 0;
  double Gen_dphijj = 0;

  unsigned njets = 0;
  unsigned nlooseEle = 0;
  unsigned ntightEle = 0;
  unsigned nlooseMu = 0;
  unsigned ntightMu = 0;
  unsigned nlooseGamma = 0;

  double Jet1_pt = 0;
  double Jet2_pt = 0;
  double Jet1_eta = 0;
  double Jet2_eta = 0;
  int Jet1_ID = 0;
  int Jet2_ID = 0;
  int Jet1_genidx = 0;
  int Jet2_genidx = 0;
  double Jet1_genjetDR = 0;
  double Jet1_JES = 0;
  double Jet2_genjetDR = 0;
  double Jet2_JES = 0;
  int Jet1_parton = 0;
  int Jet2_parton = 0;
  int Jet1_deepcsv = 0;
  int Jet2_deepcsv = 0;
  double Mjj = 0;
  double detajj = 0;
  double dphijj = 0;

  double met = 0;
  double jetmetmindphi = 0;

  outtree->Branch("GenJet1_pt",&GenJet1_pt);
  outtree->Branch("GenJet2_pt",&GenJet2_pt);
  outtree->Branch("GenJet1_eta",&GenJet1_eta);
  outtree->Branch("GenJet2_eta",&GenJet2_eta);
  outtree->Branch("Gen_Mjj",&Gen_Mjj);
  outtree->Branch("Gen_detajj",&Gen_detajj);
  outtree->Branch("Gen_dphijj",&Gen_dphijj);

  outtree->Branch("njets",&njets);
  outtree->Branch("nlooseEle",&nlooseEle);
  outtree->Branch("ntightEle",&ntightEle);
  outtree->Branch("nlooseMu",&nlooseMu);
  outtree->Branch("ntightMu",&ntightMu);
  outtree->Branch("nlooseGamma",&nlooseGamma);

  outtree->Branch("Jet1_pt",&Jet1_pt);
  outtree->Branch("Jet2_pt",&Jet2_pt);
  outtree->Branch("Jet1_eta",&Jet1_eta);
  outtree->Branch("Jet2_eta",&Jet2_eta);
  outtree->Branch("Jet1_ID",&Jet1_ID);
  outtree->Branch("Jet2_ID",&Jet2_ID);
  outtree->Branch("Jet1_genidx",&Jet1_genidx);
  outtree->Branch("Jet2_genidx",&Jet2_genidx);
  outtree->Branch("Jet1_genjetDR",&Jet1_genjetDR);
  outtree->Branch("Jet1_JES",&Jet1_JES);
  outtree->Branch("Jet2_genjetDR",&Jet2_genjetDR);
  outtree->Branch("Jet2_JES",&Jet2_JES);
  outtree->Branch("Jet1_parton",&Jet1_parton);
  outtree->Branch("Jet2_parton",&Jet2_parton);
  outtree->Branch("Jet1_deepcsv",&Jet1_deepcsv);
  outtree->Branch("Jet2_deepcsv",&Jet2_deepcsv);
  outtree->Branch("Mjj",&Mjj);
  outtree->Branch("detajj",&detajj);
  outtree->Branch("dphijj",&dphijj);

  outtree->Branch("met",&met);
  outtree->Branch("jetmetmindphi",&jetmetmindphi);

  std::cout << " .. Processing file " << aProcess << std::endl;
  std::cout << " ... Processing " <<  pu << std::endl;

  std::ostringstream filePath;
  filePath << baseDir << "/" << aProcess << "_" << pu << ".dat";
  TChain* GenJet = new TChain("ntuple/GenJet");
  TChain* Event = new TChain("ntuple/Event");
  TChain* Particle = new TChain("ntuple/Particle");
  TChain* GenPhoton = new TChain("ntuple/GenPhoton");
  TChain* Vertex = new TChain("ntuple/Vertex");
  TChain* ElectronLoose = new TChain("ntuple/ElectronLoose");
  TChain* ElectronTight = new TChain("ntuple/ElectronTight");
  TChain* MuonLoose = new TChain("ntuple/MuonLoose");
  TChain* MuonTight = new TChain("ntuple/MuonTight");
  TChain* Jets = new TChain((pu.find("noPU")!=pu.npos)?"ntuple/Jet":"ntuple/JetPUPPI");
  TChain* MissingET = new TChain((pu.find("noPU")!=pu.npos)?"ntuple/MissingET":"ntuple/PuppiMissingET");
  TChain* PhotonLoose = new TChain("ntuple/PhotonLoose");
  TChain* PhotonTight = new TChain("ntuple/PhotonTight");

  if (readFileList(filePath.str(),GenJet,
		   Event,Particle,GenPhoton,
		   Vertex,ElectronLoose,ElectronTight,
		   MuonLoose,MuonTight,
		   Jets,MissingET,
		   PhotonLoose,PhotonTight)!=0) return 1;
  
  int nLooseEle = 0;
  int nTightEle = 0;
  int nLooseMu = 0;
  int nTightMu = 0;
  int nLoosePhotons = 0;
  ElectronLoose->SetBranchAddress("ElectronLoose_size",&nLooseEle);
  ElectronTight->SetBranchAddress("ElectronTight_size",&nTightEle);
  MuonLoose->SetBranchAddress("MuonLoose_size",&nLooseMu);
  MuonTight->SetBranchAddress("MuonTight_size",&nTightMu);
  PhotonLoose->SetBranchAddress("PhotonLoose_size",&nLoosePhotons);
  
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
  
  float metpf = 0;
  float metphi = 0;
  MissingET->SetBranchAddress("MET",&metpf);
  MissingET->SetBranchAddress("Phi",&metphi);
  
  int nEvts = GenJet->GetEntries();
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
    Jets->GetEntry(ievt);
    //reco jets info
    //basic selection
    if (nJets<2) continue;
    if (jet_pt[0]<80 || jet_pt[1]<40 || fabs(jet_eta[0])>4.7 || fabs(jet_eta[1])>4.7) continue;


    TLorentzVector rec1;
    rec1.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
    TLorentzVector rec2;
    rec2.SetPtEtaPhiM(jet_pt[1],jet_eta[1],jet_phi[1],jet_mass[1]);
    TLorentzVector recPair = rec1 + rec2;
    Mjj = recPair.M();
    detajj = fabs(jet_eta[0]-jet_eta[1]);

    if (Mjj<500 || detajj<1) continue;

    Event->GetEntry(ievt);
    GenJet->GetEntry(ievt);
    MissingET->GetEntry(ievt);
    ElectronLoose->GetEntry(ievt);
    ElectronTight->GetEntry(ievt);
    MuonLoose->GetEntry(ievt);
    MuonTight->GetEntry(ievt);
    PhotonLoose->GetEntry(ievt);

    dphijj = fabs(rec1.DeltaPhi(rec2));
    
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE(metpf,0,metphi,metpf);
    double mindphi=10;
    unsigned njets30 = 0;
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
      if (ij==0) {
	Jet1_genjetDR = mindr;
	Jet1_JES = jet_pt[ij]/genjet_pt[genidx];
      }
      else if (ij==1){
	Jet2_genjetDR = mindr;
	Jet2_JES = jet_pt[ij]/genjet_pt[genidx];
      }
      
      if (jet_pt[ij]<30) continue;
      njets30++;
      if (ij>4) continue;
      
      double dphi = fabs(rec.DeltaPhi(metvec));
      if (dphi<mindphi) mindphi = dphi;
      
    }//loop on jets

    TLorentzVector gen1;
    gen1.SetPtEtaPhiM(genjet_pt[0],genjet_eta[0],genjet_phi[0],genjet_mass[0]);
    TLorentzVector gen2;
    gen2.SetPtEtaPhiM(genjet_pt[1],genjet_eta[1],genjet_phi[1],genjet_mass[1]);
    TLorentzVector genPair = gen1 + gen2;
    
    //gen jets info
    GenJet1_pt = genjet_pt[0];
    GenJet1_eta = genjet_eta[0];
    GenJet2_pt = genjet_pt[1];
    GenJet2_eta = genjet_eta[1];
    
    Gen_Mjj = genPair.M();
    Gen_detajj = fabs(genjet_eta[0]-genjet_eta[1]);
    Gen_dphijj = fabs(gen1.DeltaPhi(gen2));
    
    //recojets info
    Jet1_pt = jet_pt[0];
    Jet1_eta = jet_eta[0];
    Jet1_ID = jet_ID[0];
    Jet1_genidx = jet_genidx[0];
    Jet1_parton = jet_parton[0];
    Jet1_deepcsv = jet_DeepCSV[0];
    Jet2_pt = jet_pt[1];
    Jet2_eta = jet_eta[1];
    Jet2_ID = jet_ID[1];
    Jet2_genidx = jet_genidx[1];
    Jet2_parton = jet_parton[1];
    Jet2_deepcsv = jet_DeepCSV[1];
    
    
    //MET info
    met = metpf;
    
    njets = njets30;
    nlooseEle = nLooseEle;
    nlooseMu = nLooseMu;
    nlooseGamma = nLoosePhotons;
    ntightEle = nTightEle;
    ntightMu = nTightMu;
    jetmetmindphi = mindphi;
    
    outtree->Fill();   
  }//event loop
  
  std::cout << " ** Number of jets with different gen index: " << nDiffIdx
	    << std::endl;
  
  
  outFile->cd();
  outtree->Write();
  outFile->Write();
  
  return 0;
  
}//process



int main(int argc, char** argv){//main

  if (argc!=3) {
    std::cout << " Usage: "
	      << argv[0] << " <outputDirName> <process short name> "
	      << std::endl;
    return 1;
  }
  std::string lPlotDir = argv[1];
  std::string lProcess = argv[2];

  makeTree(lPlotDir,lProcess,"noPU");
  makeTree(lPlotDir,lProcess,"PU200");


  return 0;

}//main
