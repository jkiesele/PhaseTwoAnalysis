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

double isoCutnoPU(const unsigned id, const double & eta){
  //0,1,2 are loose, medium, tight electrons
  //3,4 are loose, tight muons
  bool isEB = false;
  if (fabs(eta)<1.5) isEB = true;
  if (isEB){
    //95% signal eff wrt ID
    if (id==0) return 0.34;
    else if (id==1) return 0.28;
    else if (id==2) return 0.23;
    else if (id==3) return 0.6;
    else if (id==4) return 0.6;
  }
  else {
    //70% bkg rej wrt ID.
    if (id==0) return 0.13;
    else if (id==1) return 0.11;
    else if (id==2) return 0.09;
    else if (id==3) return 0.35;
    else if (id==4) return 0.35;
  }
  return 0;
};

double isoCutPU200(const unsigned id, const double & eta){
  //0,1,2 are loose, medium, tight electrons
  //3,4 are loose, tight muons
  bool isEB = false;
  if (fabs(eta)<1.5) isEB = true;
  if (isEB){
    //95% signal eff wrt ID
    if (id==0) return 0.12;
    else if (id==1) return 0.18;
    else if (id==2) return 0.11;
    else if (id==3) return 0.26;
    else if (id==4) return 0.49;
  }
  else {
    //70% bkg rej wrt ID.
    if (id==0) return 0.46;
    else if (id==1) return 0.42;
    else if (id==2) return 0.36;
    else if (id==3) return 0.35;
    else if (id==4) return 0.28;
  }
  return 0;
};

double isoCut(const unsigned id, const double & eta,
	      const bool pu200){
  if (!pu200) return isoCutnoPU(id,eta);
  else return isoCutPU200(id,eta);
};



bool checkEventSize(TTree* tree, const unsigned nEvts){
  if (tree->GetEntries()!=nEvts || tree->GetEntries()==0) 
    {
      std::cout << " ---- Problem with tree " << tree->GetName() 
		<< " entries " <<tree->GetEntries() 
		<< " nevts = " << nEvts 
		<< ". Skipping..."
		<< std::endl;
      return false;
    }
  return true;
};

bool checkEventSize(TFile* file, const std::string & pu){
  if (!file->cd("ntuple")) return false;
  TTree* myGenJet = (TTree*)gDirectory->Get("GenJet");
  TTree* myEvent = (TTree*)gDirectory->Get("Event");
  TTree* myParticle = (TTree*)gDirectory->Get("Particle");
  TTree* myGenPhoton = (TTree*)gDirectory->Get("GenPhoton");
  TTree* myVertex = (TTree*)gDirectory->Get("Vertex");
  TTree* myElectronLoose = (TTree*)gDirectory->Get("ElectronLoose");
  TTree* myElectronMedium = (TTree*)gDirectory->Get("ElectronMedium");
  TTree* myElectronTight = (TTree*)gDirectory->Get("ElectronTight");
  TTree* myMuonLoose = (TTree*)gDirectory->Get("MuonLoose");
  TTree* myMuonTight = (TTree*)gDirectory->Get("MuonTight");
  TTree* myJets = (TTree*)gDirectory->Get((pu.find("noPU")!=pu.npos)?"Jet":"JetPUPPI");
  TTree* myMissingET = (TTree*)gDirectory->Get((pu.find("noPU")!=pu.npos)?"MissingET":"PuppiMissingET");
  TTree* myPhotonLoose = (TTree*)gDirectory->Get("PhotonLoose");
  TTree* myPhotonTight = (TTree*)gDirectory->Get("PhotonTight");

  int nEvts = myEvent->GetEntries();
  if (!checkEventSize(myGenJet,nEvts)) return false;
  if (!checkEventSize(myParticle,nEvts)) return false;
  if (!checkEventSize(myGenPhoton,nEvts)) return false;
  if (!checkEventSize(myVertex,nEvts)) return false;
  if (!checkEventSize(myElectronLoose,nEvts)) return false;
  if (!checkEventSize(myElectronMedium,nEvts)) return false;
  if (!checkEventSize(myElectronTight,nEvts)) return false;
  if (!checkEventSize(myMuonLoose,nEvts)) return false;
  if (!checkEventSize(myMuonTight,nEvts)) return false;
  if (!checkEventSize(myPhotonLoose,nEvts)) return false;
  if (!checkEventSize(myPhotonTight,nEvts)) return false;
  if (!checkEventSize(myJets,nEvts)) return false;
  if (!checkEventSize(myMissingET,nEvts)) return false;
  return true;
}

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else {
    std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
    if (!checkEventSize(file,input)) return false;
  }
  return true;
};

unsigned readFileList(const std::string & datFile,
		      TChain* GenJet,TChain* Event,
		      TChain* Particle,TChain* GenPhoton,TChain* Vertex,
		      TChain* ElectronLoose,TChain* ElectronMedium,TChain* ElectronTight,
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
  unsigned nFail = 0;
  std::string lBuf;
  while (std::getline(inputdat, lBuf,'\n'))
    {
      if (lBuf.empty()) continue;
      
      TFile * inFile = 0;
      lBuf = "root://gfe02.grid.hep.ph.ic.ac.uk:1097/"+lBuf;
      if (!testInputFile(lBuf,inFile)) {
	nFail++;
	continue;
      }
      GenJet->AddFile(lBuf.c_str());
      Event->AddFile(lBuf.c_str());
      Particle->AddFile(lBuf.c_str());
      GenPhoton->AddFile(lBuf.c_str());
      Vertex->AddFile(lBuf.c_str());
      ElectronLoose->AddFile(lBuf.c_str());
      ElectronMedium->AddFile(lBuf.c_str());
      ElectronTight->AddFile(lBuf.c_str());
      MuonLoose->AddFile(lBuf.c_str());
      MuonTight->AddFile(lBuf.c_str());
      Jets->AddFile(lBuf.c_str());
      MissingET->AddFile(lBuf.c_str());
      PhotonLoose->AddFile(lBuf.c_str());
      PhotonTight->AddFile(lBuf.c_str());
      nFiles++;
    }
  std::cout << " -- Read " << nFiles << " input files. Failed: " << nFail << std::endl;
  inputdat.close();
  return 0;
};

bool checkEventSize(TChain* tree, const unsigned nEvts){
  if (tree->GetEntries()!=nEvts) 
    {
      std::cout << " ---- Problem with tree " << tree->GetName() 
		<< " entries " <<tree->GetEntries() 
		<< " nevts = " << nEvts 
		<< ". Skipping..."
		<< std::endl;
      return false;
    }
  return true;
};


int makeTree(const std::string & plotDir,
	     const std::string & aProcess,
	     const std::string & pu){//main
  
  bool isPU = pu.find("no")==pu.npos;

  std::string baseDir = "filelists_20171221/";
  
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

  unsigned njets30 = 0;
  double ht30 = 0;

  unsigned nlooseEle = 0;
  unsigned nmediumEle = 0;
  unsigned ntightEle = 0;
  unsigned nlooseMu = 0;
  unsigned ntightMu = 0;
  unsigned nlooseGamma = 0;
  unsigned ntightGamma = 0;

  double ele_mt = 0;
  double mu_mt = 0;
  double Mee = 0;
  double Mmumu = 0;

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
  double metnolep = 0;
  double metphi = 0;
  double jetmetmindphi = 0;
  double jetmetnolepmindphi = 0;

  outtree->Branch("GenJet1_pt",&GenJet1_pt);
  outtree->Branch("GenJet2_pt",&GenJet2_pt);
  outtree->Branch("GenJet1_eta",&GenJet1_eta);
  outtree->Branch("GenJet2_eta",&GenJet2_eta);
  outtree->Branch("Gen_Mjj",&Gen_Mjj);
  outtree->Branch("Gen_detajj",&Gen_detajj);
  outtree->Branch("Gen_dphijj",&Gen_dphijj);

  outtree->Branch("njets30",&njets30);
  outtree->Branch("ht30",&ht30);

  outtree->Branch("nlooseEle",&nlooseEle);
  outtree->Branch("nmediumEle",&nmediumEle);
  outtree->Branch("ntightEle",&ntightEle);
  outtree->Branch("nlooseMu",&nlooseMu);
  outtree->Branch("ntightMu",&ntightMu);
  outtree->Branch("nlooseGamma",&nlooseGamma);
  outtree->Branch("ntightGamma",&ntightGamma);

  outtree->Branch("ele_mt",&ele_mt);
  outtree->Branch("mu_mt",&mu_mt);
  outtree->Branch("Mee",&Mee);
  outtree->Branch("Mmumu",&Mmumu);

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
  outtree->Branch("metnolep",&metnolep);
  outtree->Branch("metphi",&metphi);
  outtree->Branch("jetmetmindphi",&jetmetmindphi);
  outtree->Branch("jetmetnolepmindphi",&jetmetnolepmindphi);

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
  TChain* ElectronMedium = new TChain("ntuple/ElectronMedium");
  TChain* ElectronTight = new TChain("ntuple/ElectronTight");
  TChain* MuonLoose = new TChain("ntuple/MuonLoose");
  TChain* MuonTight = new TChain("ntuple/MuonTight");
  TChain* Jets = new TChain((pu.find("noPU")!=pu.npos)?"ntuple/Jet":"ntuple/JetPUPPI");
  TChain* MissingET = new TChain((pu.find("noPU")!=pu.npos)?"ntuple/MissingET":"ntuple/PuppiMissingET");
  TChain* PhotonLoose = new TChain("ntuple/PhotonLoose");
  TChain* PhotonTight = new TChain("ntuple/PhotonTight");

  if (readFileList(filePath.str(),GenJet,
		   Event,Particle,GenPhoton,
		   Vertex,ElectronLoose,ElectronMedium,ElectronTight,
		   MuonLoose,MuonTight,
		   Jets,MissingET,
		   PhotonLoose,PhotonTight)!=0) return 1;
  
  int nLooseEle = 0;
  int nMediumEle = 0;
  int nTightEle = 0;
  int nLooseMu = 0;
  int nTightMu = 0;
  int nLoosePhotons = 0;
  int nTightPhotons = 0;

  float looseEle_pt[100];
  float mediumEle_pt[100];
  float tightEle_pt[100];
  float looseEle_eta[100];
  float mediumEle_eta[100];
  float tightEle_eta[100];
  float looseEle_phi[100];
  float mediumEle_phi[100];
  float tightEle_phi[100];
  float looseEle_iso[100];
  float mediumEle_iso[100];
  float tightEle_iso[100];
  float looseEle_tkiso[100];
  float mediumEle_tkiso[100];
  float tightEle_tkiso[100];
  float looseEle_mass[100];
  float mediumEle_mass[100];
  float tightEle_mass[100];

  float looseMu_pt[100];
  float tightMu_pt[100];
  float looseMu_eta[100];
  float tightMu_eta[100];
  float looseMu_phi[100];
  float tightMu_phi[100];
  float looseMu_iso[100];
  float tightMu_iso[100];
  float looseMu_mass[100];
  float tightMu_mass[100];

  float loosePhoton_pt[100];
  float loosePhoton_pt_multi[100];
  float loosePhoton_eta[100];
  float loosePhoton_eta_multi[100];
  float loosePhoton_phi[100];
  float loosePhoton_iso[100];
  float loosePhoton_mass[100];
  int   loosePhoton_isEB[100];

  float tightPhoton_pt[100];

  ElectronLoose->SetBranchAddress("ElectronLoose_size",&nLooseEle);
  ElectronMedium->SetBranchAddress("ElectronMedium_size",&nMediumEle);
  ElectronTight->SetBranchAddress("ElectronTight_size",&nTightEle);
  ElectronLoose->SetBranchAddress("PT",&looseEle_pt);
  ElectronLoose->SetBranchAddress("Eta",&looseEle_eta);
  ElectronLoose->SetBranchAddress("Phi",&looseEle_phi);
  ElectronLoose->SetBranchAddress("Mass",&looseEle_mass);
  ElectronLoose->SetBranchAddress("IsolationVar",&looseEle_iso);
  ElectronLoose->SetBranchAddress("TrackIso",&looseEle_tkiso);
  ElectronMedium->SetBranchAddress("PT",&mediumEle_pt);
  ElectronMedium->SetBranchAddress("Eta",&mediumEle_eta);
  ElectronMedium->SetBranchAddress("Phi",&mediumEle_phi);
  ElectronMedium->SetBranchAddress("Mass",&mediumEle_mass);
  ElectronMedium->SetBranchAddress("IsolationVar",&mediumEle_iso);
  ElectronMedium->SetBranchAddress("TrackIso",&mediumEle_tkiso);
  ElectronTight->SetBranchAddress("PT",&tightEle_pt);
  ElectronTight->SetBranchAddress("Eta",&tightEle_eta);
  ElectronTight->SetBranchAddress("Phi",&tightEle_phi);
  ElectronTight->SetBranchAddress("Mass",&tightEle_mass);
  ElectronTight->SetBranchAddress("IsolationVar",&tightEle_iso);
  ElectronTight->SetBranchAddress("TrackIso",&tightEle_tkiso);

  MuonLoose->SetBranchAddress("MuonLoose_size",&nLooseMu);
  MuonLoose->SetBranchAddress("PT",&looseMu_pt);
  MuonLoose->SetBranchAddress("Eta",&looseMu_eta);
  MuonLoose->SetBranchAddress("Phi",&looseMu_phi);
  MuonLoose->SetBranchAddress("Mass",&looseMu_mass);
  MuonLoose->SetBranchAddress("IsolationVar",&looseMu_iso);

  MuonTight->SetBranchAddress("MuonTight_size",&nTightMu);
  MuonTight->SetBranchAddress("PT",&tightMu_pt);
  MuonTight->SetBranchAddress("Eta",&tightMu_eta);
  MuonTight->SetBranchAddress("Phi",&tightMu_phi);
  MuonTight->SetBranchAddress("Mass",&tightMu_mass);
  MuonTight->SetBranchAddress("IsolationVar",&tightMu_iso);

  PhotonLoose->SetBranchAddress("PhotonLoose_size",&nLoosePhotons);
  PhotonLoose->SetBranchAddress("PT",&loosePhoton_pt);
  PhotonLoose->SetBranchAddress("PT_multi",&loosePhoton_pt_multi);
  PhotonLoose->SetBranchAddress("Eta",&loosePhoton_eta);
  PhotonLoose->SetBranchAddress("Eta_multi",&loosePhoton_eta_multi);
  PhotonLoose->SetBranchAddress("Phi",&loosePhoton_phi);
  PhotonLoose->SetBranchAddress("Mass",&loosePhoton_mass);
  PhotonLoose->SetBranchAddress("IsolationVar",&loosePhoton_iso);
  PhotonLoose->SetBranchAddress("IsEB",&loosePhoton_isEB);

  PhotonTight->SetBranchAddress("PhotonTight_size",&nTightPhotons);
  PhotonTight->SetBranchAddress("PT",&tightPhoton_pt);

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
  float metphipf = 0;
  MissingET->SetBranchAddress("MET",&metpf);
  MissingET->SetBranchAddress("Phi",&metphipf);
  
  int nEvts = GenJet->GetEntries();
  if (!checkEventSize(GenJet,nEvts)) return 1;
  if (!checkEventSize(Particle,nEvts)) return 1;
  if (!checkEventSize(GenPhoton,nEvts)) return 1;
  if (!checkEventSize(Vertex,nEvts)) return 1;
  if (!checkEventSize(ElectronLoose,nEvts)) return 1;
  if (!checkEventSize(ElectronMedium,nEvts)) return 1;
  if (!checkEventSize(ElectronTight,nEvts)) return 1;
  if (!checkEventSize(MuonLoose,nEvts)) return 1;
  if (!checkEventSize(MuonTight,nEvts)) return 1;
  if (!checkEventSize(PhotonLoose,nEvts)) return 1;
  if (!checkEventSize(PhotonTight,nEvts)) return 1;
  if (!checkEventSize(Jets,nEvts)) return 1;
  if (!checkEventSize(MissingET,nEvts)) return 1;
  
  unsigned nDiffIdx = 0;
  
  std::cout << " -- Processing " << nEvts << " events." << std::endl;
  //event loop
  for (int ievt(0); ievt<nEvts; ++ievt){
    if (ievt%1000==0) std::cout << ".... Processing entry " << ievt << std::endl;

    GenJet1_pt = 0;
    GenJet2_pt = 0;
    GenJet1_eta = 0;
    GenJet2_eta = 0;
    Gen_Mjj = 0;
    Gen_detajj = 0;
    Gen_dphijj = 0;
    
    njets30 = 0;
    ht30 = 0;
    
    nlooseEle = 0;
    nmediumEle = 0;
    ntightEle = 0;
    nlooseMu = 0;
    ntightMu = 0;
    nlooseGamma = 0;
    ntightGamma = 0;
    
    ele_mt = 0;
    mu_mt = 0;
    Mee = 0;
    Mmumu = 0;
    
    Jet1_pt = 0;
    Jet2_pt = 0;
    Jet1_eta = 0;
    Jet2_eta = 0;
    Jet1_ID = 0;
    Jet2_ID = 0;
    Jet1_genidx = 0;
    Jet2_genidx = 0;
    Jet1_genjetDR = 0;
    Jet1_JES = 0;
    Jet2_genjetDR = 0;
    Jet2_JES = 0;
    Jet1_parton = 0;
    Jet2_parton = 0;
    Jet1_deepcsv = 0;
    Jet2_deepcsv = 0;
    Mjj = 0;
    detajj = 0;
    dphijj = 0;
    
    met = 0;
    metnolep = 0;
    metphi = 0;
    jetmetmindphi = 0;
    jetmetnolepmindphi = 0;
    
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
    dphijj = fabs(rec1.DeltaPhi(rec2));

    if (Mjj<500 || detajj<1) continue;

    Event->GetEntry(ievt);
    GenJet->GetEntry(ievt);
    MissingET->GetEntry(ievt);
    ElectronLoose->GetEntry(ievt);
    ElectronMedium->GetEntry(ievt);
    ElectronTight->GetEntry(ievt);
    MuonLoose->GetEntry(ievt);
    MuonTight->GetEntry(ievt);
    PhotonLoose->GetEntry(ievt);
    PhotonTight->GetEntry(ievt);

    
    //loose leptons, get px and py component to add to the met
    //apply loose selections on them
    double pxSum = 0;
    double pySum = 0;
    nlooseEle = 0;
    bool first = true;
    bool second = true;
    unsigned ele1 = 0;
    unsigned ele2 = 0;
    for (unsigned il(0); il<abs(nLooseEle);++il){
      if (looseEle_iso[il] >= isoCut(0,looseEle_eta[il],isPU)) continue;
      if (looseEle_pt[il]<10 || fabs(looseEle_eta[il])>2.8) continue;
      if (fabs(looseEle_eta[il])>1.444 && fabs(looseEle_eta[il])<1.566) continue;
      double px = looseEle_pt[il]*cos(looseEle_phi[il]);
      double py = looseEle_pt[il]*sin(looseEle_phi[il]);
      pxSum += px;
      pySum += py;
      if (first){
	ele1 = il;
	first = false;
      }
      else if (second){
	ele2 = il;
	second = false;
      }
      nlooseEle++;
    }
    nlooseMu = 0;
    first = true;
    second = true;
    unsigned mu1 = 0;
    unsigned mu2 = 0;
    for (unsigned il(0); il<abs(nLooseMu);++il){
      if (looseMu_iso[il] >= isoCut(3,looseMu_eta[il],isPU)) continue;
      if (looseMu_pt[il]<10 || fabs(looseMu_eta[il])>3.0) continue;
      double px = looseMu_pt[il]*cos(looseMu_phi[il]);
      double py = looseMu_pt[il]*sin(looseMu_phi[il]);
      pxSum += px;
      pySum += py;
      if (first){
	mu1 = il;
	first = false;
      }
      else if (second){
	mu2 = il;
	second = false;
      }
      nlooseMu++;
    }
    //apply tight selections on tight leptons
    nmediumEle = 0;
    for (unsigned il(0); il<abs(nMediumEle);++il){
      if (mediumEle_iso[il] >= isoCut(1,mediumEle_eta[il],isPU)) continue;
      if (mediumEle_pt[il]<10 || fabs(mediumEle_eta[il])>2.8) continue;
      if (fabs(mediumEle_eta[il])>1.444 && fabs(mediumEle_eta[il])<1.566) continue;
      nmediumEle++;
    }
    ntightEle = 0;
    for (unsigned il(0); il<abs(nTightEle);++il){
      if (tightEle_iso[il] >= isoCut(2,tightEle_eta[il],isPU)) continue;
      if (tightEle_pt[il]<20 || fabs(tightEle_eta[il])>2.8) continue;
      if (fabs(tightEle_eta[il])>1.444 && fabs(tightEle_eta[il])<1.566) continue;
      ntightEle++;
    }

    ntightMu = 0;
    for (unsigned il(0); il<abs(nTightMu);++il){
      if (tightMu_iso[il] >= isoCut(4,tightMu_eta[il],isPU)) continue;
      if (tightMu_pt[il]<20 || fabs(tightMu_eta[il])>3.0) continue;
      double px = tightMu_pt[il]*cos(tightMu_phi[il]);
      double py = tightMu_pt[il]*sin(tightMu_phi[il]);
      pxSum += px;
      pySum += py;
      ntightMu++;
    }

    //photons 
    nlooseGamma = 0;
    for (unsigned il(0); il<abs(nLoosePhotons);++il){
      //if (loosePhoton_iso[il] >= isoCut(4,loosePhoton_eta[il])) continue;
      if (loosePhoton_pt[il]<15) continue;
      nlooseGamma++;
    }
    ntightGamma = 0;
    for (unsigned il(0); il<abs(nTightPhotons);++il){
      //if (tightPhoton_iso[il] >= isoCut(4,tightPhoton_eta[il])) continue;
      if (tightPhoton_pt[il]<15) continue;
      ntightGamma++;
    }

    //met and jets
    TLorentzVector metvec;
    metvec.SetPtEtaPhiE(metpf,0,metphipf,metpf);
    if (nlooseEle==1) ele_mt = sqrt(2*looseEle_pt[ele1]*metpf*(1-cos(looseEle_phi[ele1]-metphipf)));
    if (nlooseEle==2){
      TLorentzVector ele1vec ;
      ele1vec.SetPtEtaPhiM(looseEle_pt[ele1],looseEle_eta[ele1],looseEle_phi[ele1],looseEle_mass[ele1]);
      TLorentzVector ele2vec ;
      ele2vec.SetPtEtaPhiM(looseEle_pt[ele2],looseEle_eta[ele2],looseEle_phi[ele2],looseEle_mass[ele2]);
      Mee = (ele1vec+ele2vec).M();
    }
    if (nlooseMu==1) mu_mt = sqrt(2*looseMu_pt[mu1]*metpf*(1-cos(looseMu_phi[mu1]-metphipf)));
    if (nlooseMu==2){
      TLorentzVector mu1vec ;
      mu1vec.SetPtEtaPhiM(looseMu_pt[mu1],looseMu_eta[mu1],looseMu_phi[mu1],looseMu_mass[mu1]);
      TLorentzVector mu2vec ;
      mu2vec.SetPtEtaPhiM(looseMu_pt[mu2],looseMu_eta[mu2],looseMu_phi[mu2],looseMu_mass[mu2]);
      Mmumu = (mu1vec+mu2vec).M();
    }

    //calculate met no lepton
    TLorentzVector metnolepvec;
    pxSum += metpf*cos(metphipf);
    pySum += metpf*sin(metphipf);
    metnolep = sqrt(pxSum*pxSum+pySum*pySum);
    double metnolepphi = acos(pxSum/metnolep);
    metnolepvec.SetPtEtaPhiE(metnolep,0,metnolepphi,metnolep);

    jetmetmindphi=10;
    jetmetnolepmindphi=10;
    njets30 = 0;
    ht30 = 0;
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
      ht30 += jet_pt[ij];
      if (ij>4) continue;
      
      double dphi = fabs(rec.DeltaPhi(metvec));
      if (dphi<jetmetmindphi) jetmetmindphi = dphi;
      double dphinolep = fabs(rec.DeltaPhi(metnolepvec));
      if (dphinolep<jetmetnolepmindphi) jetmetnolepmindphi = dphinolep;
    
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
    metphi = metphipf;
    
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

  if (argc!=4) {
    std::cout << " Usage: "
	      << argv[0] << " <outputDirName> <process short name> <pu: noPU or PU200>"
	      << std::endl;
    return 1;
  }
  std::string lPlotDir = argv[1];
  std::string lProcess = argv[2];
  std::string lPu = argv[3];

  makeTree(lPlotDir,lProcess,lPu);

  return 0;

}//main
