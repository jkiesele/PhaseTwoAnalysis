// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/BasicRecoDistrib
// Class:      BasicRecoDistrib
// 
/**\class BasicRecoDistrib BasicRecoDistrib.cc PhaseTwoAnalysis/BasicRecoDistrib/plugins/BasicRecoDistrib.cc

Description: produces histograms of basic quantities from RECO collections

Implementation:
   - muon isolation comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_isolation 
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_identification
   - electron isolation needs to be refined
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
   - no jet ID nor JEC are applied
   - b-tagging is not available 

*/
//
// Original Author:  Elvire Bouvier
//         Created:  Wed, 14 Jun 2017 14:16:16 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMVA/Reader.h" 
#include "TMVA/MethodBDT.h" 

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BasicRecoDistrib : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit BasicRecoDistrib(const edm::ParameterSet&);
    ~BasicRecoDistrib();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  

  private:
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endJob() override;

    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
    bool isME0MuonSelNew(reco::Muon, double, double, double);    
    bool isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    int matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles);
    void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
//    float evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

//    std::unique_ptr<HGCalIDTool> hgcEmId_; 
//    TMVA::Reader tmvaReader_;
//    float hgcId_startPosition, hgcId_lengthCompatibility, hgcId_sigmaietaieta, hgcId_deltaEtaStartPosition, hgcId_deltaPhiStartPosition, hOverE_hgcalSafe, hgcId_cosTrackShowerAngle, trackIsoR04jurassic_D_pt, ooEmooP, d0, dz, pt, etaSC, phiSC, nPV, expectedMissingInnerHits, passConversionVeto, isTrue;

    unsigned int pileup_;
    edm::EDGetTokenT<std::vector<reco::GsfElectron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<edm::ValueMap<double>> trackIsoValueMapToken_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_charged_hadrons_;
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_neutral_hadrons_;
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPINoLeptonsIsolation_photons_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsToken_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsNoLepToken_;
    edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken_;
    edm::EDGetTokenT<std::vector<reco::PFMET>> metToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    const ME0Geometry* ME0Geometry_; 
    double muThres_;

    // Electrons
    TH1D* h_allElecs_n_;
    TH1D* h_allElecs_pt_;
    TH1D* h_allElecs_eta_;
    TH1D* h_allElecs_phi_;
    TH1D* h_allElecs_iso_;
    TH1D* h_allElecs_id_;
    // after cut ID
    TH1D* h_elecs_n_;
    TH1D* h_elecs_pt_;
    TH1D* h_elecs_eta_;
    TH1D* h_elecs_phi_;
    TH1D* h_elecs_iso_;
    TH1D* h_PFElecs_n_;
    TH1D* h_PFElecs_pt_;
    TH1D* h_PFElecs_eta_;
    TH1D* h_PFElecs_phi_;
    TH1D* h_PFElecs_iso_;
    // ... that are isolated
    TH1D* h_goodElecs_n_;
    TH1D* h_goodElecs_pt_;
    TH1D* h_goodElecs_eta_;
    TH1D* h_goodElecs_phi_;
    TH1D* h_goodElecs_iso_;
    TH1D* h_goodPFElecs_n_;
    TH1D* h_goodPFElecs_pt_;
    TH1D* h_goodPFElecs_eta_;
    TH1D* h_goodPFElecs_phi_;
    TH1D* h_goodPFElecs_iso_;

    // Muons
    TH1D* h_allMuons_n_;
    TH1D* h_allMuons_pt_;
    TH1D* h_allMuons_eta_;
    TH1D* h_allMuons_phi_;
    TH1D* h_allMuons_iso_;
    TH1D* h_allMuons_id_;
    // ... that are tight
    TH1D* h_muons_n_;
    TH1D* h_muons_pt_;
    TH1D* h_muons_eta_;
    TH1D* h_muons_phi_;
    TH1D* h_muons_iso_;
    TH1D* h_PFMuons_n_;
    TH1D* h_PFMuons_pt_;
    TH1D* h_PFMuons_eta_;
    TH1D* h_PFMuons_phi_;
    TH1D* h_PFMuons_iso_;
    // ... that are isolated
    TH1D* h_goodMuons_n_;
    TH1D* h_goodMuons_pt_;
    TH1D* h_goodMuons_eta_;
    TH1D* h_goodMuons_phi_;
    TH1D* h_goodMuons_iso_;
    TH1D* h_goodPFMuons_n_;
    TH1D* h_goodPFMuons_pt_;
    TH1D* h_goodPFMuons_eta_;
    TH1D* h_goodPFMuons_phi_;
    TH1D* h_goodPFMuons_iso_;

    // Jets
    // ... with p_T > 20 GeV
    TH1D* h_jet20_n_;
    TH1D* h_jet20_pt_;
    TH1D* h_jet20_eta_;
    TH1D* h_jet20_phi_;
    // ... with p_T > 30 GeV
    TH1D* h_jet30_n_;
    TH1D* h_jet30_pt_;
    TH1D* h_jet30_eta_;
    TH1D* h_jet30_phi_;
    // ... with p_T > 40 GeV
    TH1D* h_jet40_n_;
    TH1D* h_jet40_pt_;
    TH1D* h_jet40_eta_;
    TH1D* h_jet40_phi_;
    // ... with p_T > 50 GeV
    TH1D* h_jet50_n_;
    TH1D* h_jet50_pt_;
    TH1D* h_jet50_eta_;
    TH1D* h_jet50_phi_;

    // MET
    TH1D* h_met_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BasicRecoDistrib::BasicRecoDistrib(const edm::ParameterSet& iConfig): 
  pileup_(iConfig.getParameter<unsigned int>("pileup")),
  elecsToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
  trackIsoValueMapToken_(consumes<edm::ValueMap<double>>(iConfig.getParameter<edm::InputTag>("trackIsoValueMap"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  pfCandsToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  pfCandsNoLepToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandsNoLep"))),
  jetsToken_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("met"))),
  genPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices")))
{
  //now do what ever initialization is needed
  PUPPINoLeptonsIsolation_charged_hadrons_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNoLepIsolationChargedHadrons"));
  PUPPINoLeptonsIsolation_neutral_hadrons_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNoLepIsolationNeutralHadrons"));
  PUPPINoLeptonsIsolation_photons_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiNoLepIsolationPhotons"));

  if (pileup_ == 0) {
    muThres_ = 0.152;
  } else if (pileup_ == 140) {
    muThres_ = 0.204;
  } else if (pileup_ == 200) {
    muThres_ = 0.212;
  } else 
    muThres_ = 0.;

  usesResource("TFileService");

//  const edm::ParameterSet& hgcIdCfg = iConfig.getParameterSet("HGCalIDToolConfig");
//  auto cc = consumesCollector();
//  hgcEmId_.reset( new HGCalIDTool(hgcIdCfg, cc) );

/*  tmvaReader_.SetOptions("!Color:Silent:!Error");
  tmvaReader_.AddVariable("hgcId_startPosition", &hgcId_startPosition);
  tmvaReader_.AddVariable("hgcId_lengthCompatibility", &hgcId_lengthCompatibility);
  tmvaReader_.AddVariable("hgcId_sigmaietaieta", &hgcId_sigmaietaieta);
  tmvaReader_.AddVariable("abs(hgcId_deltaEtaStartPosition)", &hgcId_deltaEtaStartPosition);
  tmvaReader_.AddVariable("abs(hgcId_deltaPhiStartPosition)", &hgcId_deltaPhiStartPosition);
  tmvaReader_.AddVariable("hOverE_hgcalSafe", &hOverE_hgcalSafe);
  tmvaReader_.AddVariable("hgcId_cosTrackShowerAngle", &hgcId_cosTrackShowerAngle);
  tmvaReader_.AddVariable("trackIsoR04jurassic_D_pt := trackIsoR04jurassic/pt", &trackIsoR04jurassic_D_pt);
  tmvaReader_.AddVariable("abs(ooEmooP)", &ooEmooP);
  tmvaReader_.AddVariable("abs(d0)", &d0);
  tmvaReader_.AddVariable("abs(dz)", &dz);
  tmvaReader_.AddVariable("expectedMissingInnerHits", &expectedMissingInnerHits);
  tmvaReader_.AddSpectator("pt",  &pt);
  tmvaReader_.AddSpectator("nPV",  &nPV);
  tmvaReader_.AddSpectator("etaSC",  &etaSC);
  tmvaReader_.AddSpectator("phiSC",  &phiSC);
  tmvaReader_.AddSpectator("isTrue",  &isTrue);
  tmvaReader_.AddSpectator("passConversionVeto", &passConversionVeto);

  tmvaReader_.BookMVA("PhaseIIEndcapHGCal","TMVAClassification_BDT.weights.xml");
*/
  // Electrons
  h_allElecs_n_ = fs_->make<TH1D>("AllElecsN",";Number of electrons;Events / 1",10,0.,10.);
  h_allElecs_pt_ = fs_->make<TH1D>("AllElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)",50,0.,250.);
  h_allElecs_eta_ = fs_->make<TH1D>("AllElecsEta",";#eta(e);Events / 0.2",30,-3.,3.);
  h_allElecs_phi_ = fs_->make<TH1D>("AllElecsPhi",";#phi(e);Events / 0.2",30,-3.,3.);
  h_allElecs_iso_ = fs_->make<TH1D>("AllElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 40, 0., 0.4);
  h_allElecs_id_ = fs_->make<TH1D>("AllElecsID",";;Electrons / 1", 4, 0., 4.);
  h_allElecs_id_->SetOption("bar");
  h_allElecs_id_->SetBarWidth(0.75);
  h_allElecs_id_->SetBarOffset(0.125);
  h_allElecs_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allElecs_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allElecs_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allElecs_id_->GetXaxis()->SetBinLabel(4,"Tight");
  
  // after cut ID
  h_elecs_n_ = fs_->make<TH1D>("ElecsN",";Number of electrons;Events / 1",4,0.,4.);
  h_elecs_pt_ = fs_->make<TH1D>("ElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)",50,0.,250.);
  h_elecs_eta_ = fs_->make<TH1D>("ElecsEta",";#eta(e);Events / 0.2",30,-3.,3.);
  h_elecs_phi_ = fs_->make<TH1D>("ElecsPhi",";#phi(e);Events / 0.2",30,-3.,3.);
  h_elecs_iso_ = fs_->make<TH1D>("ElecsIso",";I_{rel}^{PUPPI}(iso e);Events / 0.02",200,0.,4.);
  h_PFElecs_n_ = fs_->make<TH1D>("PFElecsN",";Number of electrons;Events / 1",4,0.,4.);
  h_PFElecs_pt_ = fs_->make<TH1D>("PFElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)",50,0.,250.);
  h_PFElecs_eta_ = fs_->make<TH1D>("PFElecsEta",";#eta(e);Events / 0.2",30,-3.,3.);
  h_PFElecs_phi_ = fs_->make<TH1D>("PFElecsPhi",";#phi(e);Events / 0.2",30,-3.,3.);
  h_PFElecs_iso_ = fs_->make<TH1D>("PFElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.02",200,0.,4.);
  //... that are isolated
  h_goodElecs_n_ = fs_->make<TH1D>("GoodElecsN",";Number of isolated electrons;Events / 1",4,0.,4.);
  h_goodElecs_pt_ = fs_->make<TH1D>("GoodElecsPt",";p_{T}(iso e) (GeV);Events / (5 GeV)",50,0.,250.);
  h_goodElecs_eta_ = fs_->make<TH1D>("GoodElecsEta",";#eta(iso e);Events / 0.2",30,-3.,3.);
  h_goodElecs_phi_ = fs_->make<TH1D>("GoodElecsPhi",";#phi(iso e);Events / 0.2",30,-3.,3.);
  h_goodElecs_iso_ = fs_->make<TH1D>("GoodElecsIso",";I_{rel}^{PUPPI}(iso e);Events / 0.01",20,0.,0.2);
  h_goodPFElecs_n_ = fs_->make<TH1D>("GoodPFElecsN",";Number of isolated electrons;Events / 1",4,0.,4.);
  h_goodPFElecs_pt_ = fs_->make<TH1D>("GoodPFElecsPt",";p_{T}(iso e) (GeV);Events / (5 GeV)",50,0.,250.);
  h_goodPFElecs_eta_ = fs_->make<TH1D>("GoodPFElecsEta",";#eta(iso e);Events / 0.2",30,-3.,3.);
  h_goodPFElecs_phi_ = fs_->make<TH1D>("GoodPFElecsPhi",";#phi(iso e);Events / 0.2",30,-3.,3.);
  h_goodPFElecs_iso_ = fs_->make<TH1D>("GoodPFElecsIso",";I_{rel}^{PUPPI}(iso e);Events / 0.01",20,0.,0.2);

  // Muons
  h_allMuons_n_ = fs_->make<TH1D>("AllMuonsN",";Number of muons;Events / 1",10,0.,10.);
  h_allMuons_pt_ = fs_->make<TH1D>("AllMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)",50,0.,250.);
  h_allMuons_eta_ = fs_->make<TH1D>("AllMuonsEta",";#eta(#mu);Events / 0.2",30,-3.,3.);
  h_allMuons_phi_ = fs_->make<TH1D>("AllMuonsPhi",";#phi(#mu);Events / 0.2",30,-3.,3.);
  h_allMuons_iso_ = fs_->make<TH1D>("AllMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 40, 0., 0.4);
  h_allMuons_id_ = fs_->make<TH1D>("AllMuonsID",";;Muons / 1", 4, 0., 4.);
  h_allMuons_id_->SetOption("bar");
  h_allMuons_id_->SetBarWidth(0.75);
  h_allMuons_id_->SetBarOffset(0.125);
  h_allMuons_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allMuons_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allMuons_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allMuons_id_->GetXaxis()->SetBinLabel(4,"Tight");
  // after tight ID
  h_muons_n_ = fs_->make<TH1D>("MuonsN",";Number of muons;Events / 1",4,0.,4.);
  h_muons_pt_ = fs_->make<TH1D>("MuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)",50,0.,250.);
  h_muons_eta_ = fs_->make<TH1D>("MuonsEta",";#eta(#mu);Events / 0.2",30,-3.,3.);
  h_muons_phi_ = fs_->make<TH1D>("MuonsPhi",";#phi(#mu);Events / 0.2",30,-3.,3.);
  h_muons_iso_ = fs_->make<TH1D>("MuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.02",200,0.,4.);
  h_PFMuons_n_ = fs_->make<TH1D>("PFMuonsN",";Number of muons;Events / 1",4,0.,4.);
  h_PFMuons_pt_ = fs_->make<TH1D>("PFMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)",50,0.,250.);
  h_PFMuons_eta_ = fs_->make<TH1D>("PFMuonsEta",";#eta(#mu);Events / 0.2",30,-3.,3.);
  h_PFMuons_phi_ = fs_->make<TH1D>("PFMuonsPhi",";#phi(#mu);Events / 0.2",30,-3.,3.);
  h_PFMuons_iso_ = fs_->make<TH1D>("PFMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.02",200,0.,4.);
  //... that are isolated
  h_goodMuons_n_ = fs_->make<TH1D>("GoodMuonsN",";Number of isolated muons;Events / 1",4,0.,4.);
  h_goodMuons_pt_ = fs_->make<TH1D>("GoodMuonsPt",";p_{T}(iso #mu) (GeV);Events / (5 GeV)",50,0.,250.);
  h_goodMuons_eta_ = fs_->make<TH1D>("GoodMuonsEta",";#eta(iso #mu);Events / 0.2",30,-3.,3.);
  h_goodMuons_phi_ = fs_->make<TH1D>("GoodMuonsPhi",";#phi(iso #mu);Events / 0.2",30,-3.,3.);
  h_goodMuons_iso_ = fs_->make<TH1D>("GoodMuonsIso",";I_{rel}^{PUPPI}(iso #mu);Events / 0.01",20,0.,0.2);
  h_goodPFMuons_n_ = fs_->make<TH1D>("GoodPFMuonsN",";Number of isolated muons;Events / 1",4,0.,4.);
  h_goodPFMuons_pt_ = fs_->make<TH1D>("GoodPFMuonsPt",";p_{T}(iso #mu) (GeV);Events / (5 GeV)",50,0.,250.);
  h_goodPFMuons_eta_ = fs_->make<TH1D>("GoodPFMuonsEta",";#eta(iso #mu);Events / 0.2",30,-3.,3.);
  h_goodPFMuons_phi_ = fs_->make<TH1D>("GoodPFMuonsPhi",";#phi(iso #mu);Events / 0.2",30,-3.,3.);
  h_goodPFMuons_iso_ = fs_->make<TH1D>("GoodPFMuonsIso",";I_{rel}^{PUPPI}(iso #mu);Events / 0.01",20,0.,0.2);

  // Jets
  // ... with p_T > 20 GeV
  h_jet20_n_ = fs_->make<TH1D>("Jets20N",";Number of jets;Events / 1",14,0.,14.);
  h_jet20_pt_ = fs_->make<TH1D>("Jets20Pt",";p_{T}(jet) (GeV);Events / (5 GeV)",46,20.,250.);
  h_jet20_eta_ = fs_->make<TH1D>("Jets20Eta",";#eta(jet);Events / 0.2",40,-4.,4.);
  h_jet20_phi_ = fs_->make<TH1D>("Jets20Phi",";#phi(jet);Events / 0.2",30,-3.,3.);
  // ... with p_T > 30 GeV
  h_jet30_n_ = fs_->make<TH1D>("Jets30N",";Number of jets;Events / 1",14,0.,14.);
  h_jet30_pt_ = fs_->make<TH1D>("Jets30Pt",";p_{T}(jet) (GeV);Events / (5 GeV)",44,30.,250.);
  h_jet30_eta_ = fs_->make<TH1D>("Jets30Eta",";#eta(jet);Events / 0.2",40,-4.,4.);
  h_jet30_phi_ = fs_->make<TH1D>("Jets30Phi",";#phi(jet);Events / 0.2",30,-3.,3.);
  // ... with p_T > 40 GeV
  h_jet40_n_ = fs_->make<TH1D>("Jets40N",";Number of jets;Events / 1",12,0.,12);
  h_jet40_pt_ = fs_->make<TH1D>("Jets40Pt",";p_{T}(jet) (GeV);Events / (5 GeV)",42,40.,250.);
  h_jet40_eta_ = fs_->make<TH1D>("Jets40Eta",";#eta(jet);Events / 0.2",40,-4.,4.);
  h_jet40_phi_ = fs_->make<TH1D>("Jets40Phi",";#phi(jet);Events / 0.2",30,-3.,3.);
  // ... with p_T > 50 GeV
  h_jet50_n_ = fs_->make<TH1D>("Jets50N",";Number of jets;Events / 1",12,0.,12);
  h_jet50_pt_ = fs_->make<TH1D>("Jets50Pt",";p_{T}(jet) (GeV);Events / (5 GeV)",40,50.,250);
  h_jet50_eta_ = fs_->make<TH1D>("Jets50Eta",";#eta(jet);Events / 0.2",40,-4.,4.);
  h_jet50_phi_ = fs_->make<TH1D>("Jets50Phi",";#phi(jet);Events / 0.2",30,-3.,3.);

  // MET
  h_met_ = fs_->make<TH1D>("METPt",";p_{T}(MET) (GeV); Events / (5 GeV)",60,0.,300.);

}


BasicRecoDistrib::~BasicRecoDistrib()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
BasicRecoDistrib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

//  hgcEmId_->getEventSetup(iSetup);
//  hgcEmId_->getEvent(iEvent);

  Handle<std::vector<reco::GsfElectron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
//  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<ValueMap<double>> trackIsoValueMap;
  iEvent.getByToken(trackIsoValueMapToken_, trackIsoValueMap);

  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_charged_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_neutral_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPINoLeptonsIsolation_photons;
  iEvent.getByToken(PUPPINoLeptonsIsolation_charged_hadrons_, PUPPINoLeptonsIsolation_charged_hadrons);
  iEvent.getByToken(PUPPINoLeptonsIsolation_neutral_hadrons_, PUPPINoLeptonsIsolation_neutral_hadrons);
  iEvent.getByToken(PUPPINoLeptonsIsolation_photons_, PUPPINoLeptonsIsolation_photons);  

  Handle<std::vector<reco::PFCandidate>> pfCands;
  iEvent.getByToken(pfCandsToken_, pfCands);

  Handle<std::vector<reco::PFCandidate>> pfCandsNoLep;
  iEvent.getByToken(pfCandsNoLepToken_, pfCandsNoLep);

  Handle<std::vector<reco::PFJet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<reco::PFMET>> met;
  iEvent.getByToken(metToken_, met);

  Handle<std::vector<reco::GenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);


  int prVtx = -1;
  for(size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4.) continue;
    if (prVtx < 0) prVtx = i;
  }
  if (prVtx < 0.) return;


  // Electrons
  int nElec = 0;
  int nGoodElec = 0;
  h_allElecs_n_->Fill(elecs->size());
  for(size_t i = 0; i < elecs->size(); i++) { 
    h_allElecs_pt_->Fill(elecs->at(i).pt());
    h_allElecs_eta_->Fill(elecs->at(i).eta());
    h_allElecs_phi_->Fill(elecs->at(i).phi());
    double isoEl = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(elecs->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.4) continue;
      isoEl += pfCandsNoLep->at(k).pt();
    }
    if (elecs->at(i).pt() > 0.) isoEl = isoEl / elecs->at(i).pt(); 
    else isoEl = -1.;
    h_allElecs_iso_->Fill(isoEl);
    Ptr<const reco::GsfElectron> el4iso(elecs,i);
//    double eljurassicIso = (*trackIsoValueMap)[el4iso];
//    double elpt = elecs->at(i).pt();
//    double elMVAVal = -1.;
//    if (hgcEmId_->setElectronPtr(&(elecs->at(i)))) 
//      elMVAVal = (double)evalMVAElec(elecs->at(i),vertices->at(prVtx),conversions,beamspot,genParts,eljurassicIso/elpt,vertices->size());
    h_allElecs_id_->Fill(0.);
//    if (isLooseElec(elecs->at(i),conversions,beamspot,elMVAVal)) h_allElecs_id_->Fill(1.);    
//    if (isMediumElec(elecs->at(i),conversions,beamspot,elMVAVal)) h_allElecs_id_->Fill(2.);    
//    if (isTightElec(elecs->at(i),conversions,beamspot,elMVAVal)) h_allElecs_id_->Fill(3.);    

//    if (!isTightElec(elecs->at(i),conversions,beamspot,elMVAVal)) continue;
    if (fabs(elecs->at(i).eta()) > 2.8) continue;
    if (elecs->at(i).pt() < 20.) continue;
    h_elecs_pt_->Fill(elecs->at(i).pt());
    h_elecs_eta_->Fill(elecs->at(i).eta());
    h_elecs_phi_->Fill(elecs->at(i).phi());
    h_elecs_iso_->Fill(isoEl);
    ++nElec;

    if (elecs->at(i).pt() < 30. || isoEl > 0.15 ||
        (fabs(elecs->at(i).eta()) > 1.479 && fabs(elecs->at(i).eta()) < 1.5660)) continue;
    h_goodElecs_pt_->Fill(elecs->at(i).pt());
    h_goodElecs_eta_->Fill(elecs->at(i).eta());
    h_goodElecs_phi_->Fill(elecs->at(i).phi());
    h_goodElecs_iso_->Fill(isoEl);
    ++nGoodElec;
  }
  h_elecs_n_->Fill(nElec);
  h_goodElecs_n_->Fill(nGoodElec);
  //PF elecs
  int nPFElec =0;
  int nGoodPFElec = 0;
  for (size_t i = 0; i < pfCands->size(); i++) {
    if (abs(pfCands->at(i).pdgId()) != 11) continue;
    double isoPFEl = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(pfCands->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.4) continue;
      isoPFEl += pfCandsNoLep->at(k).pt();
    }
    isoPFEl = isoPFEl / pfCands->at(i).pt();
    if (fabs(pfCands->at(i).eta()) > 2.8) continue;
    if (pfCands->at(i).pt() < 20.) continue;
    h_PFElecs_pt_->Fill(pfCands->at(i).pt());
    h_PFElecs_eta_->Fill(pfCands->at(i).eta());
    h_PFElecs_phi_->Fill(pfCands->at(i).phi());
    h_PFElecs_iso_->Fill(isoPFEl);
    ++nPFElec;

    if (pfCands->at(i).pt() < 30. || isoPFEl > 0.15 ||
        (fabs(pfCands->at(i).eta()) > 1.479 && fabs(pfCands->at(i).eta()) < 1.5660)) continue;
    h_goodPFElecs_pt_->Fill(pfCands->at(i).pt());
    h_goodPFElecs_eta_->Fill(pfCands->at(i).eta());
    h_goodPFElecs_phi_->Fill(pfCands->at(i).phi());
    h_goodPFElecs_iso_->Fill(isoPFEl);
    ++nGoodPFElec;
  }
  h_PFElecs_n_->Fill(nPFElec);
  h_goodPFElecs_n_->Fill(nGoodPFElec);

  // Muons
  int nMuon =0;
  int nGoodMuon = 0;
  h_allMuons_n_->Fill(muons->size());
  for(size_t i = 0; i < muons->size(); i++){
    h_allMuons_pt_->Fill(muons->at(i).pt());
    h_allMuons_eta_->Fill(muons->at(i).eta());
    h_allMuons_phi_->Fill(muons->at(i).phi());
    Ptr<const reco::Muon> muref(muons,i);
    double muon_puppiIsoNoLep_ChargedHadron = (*PUPPINoLeptonsIsolation_charged_hadrons)[muref];
    double muon_puppiIsoNoLep_NeutralHadron = (*PUPPINoLeptonsIsolation_neutral_hadrons)[muref];
    double muon_puppiIsoNoLep_Photon = (*PUPPINoLeptonsIsolation_photons)[muref];
    double isoMu = (muon_puppiIsoNoLep_ChargedHadron+muon_puppiIsoNoLep_NeutralHadron+muon_puppiIsoNoLep_Photon)/muons->at(i).pt();
    h_allMuons_iso_->Fill(isoMu);
    
    // Loose ID
    double dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.056);
    double dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0096);    
    bool isLooseMuon = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut));

    // Medium ID -- needs to be updated
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    if (muons->at(i).innerTrack().isNonnull()){
    	ipxy = std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
    	ipz = std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
    	validPxlHit = muons->at(i).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    	highPurity = muons->at(i).innerTrack()->quality(reco::Track::highPurity);
    }    
    bool isMediumMuon = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    // Tight ID
    dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.032);
    dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0041);
    bool isTightMuon = (fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.048, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    h_allMuons_id_->Fill(0.);
    if (isLooseMuon) h_allMuons_id_->Fill(1.);
    if (isMediumMuon) h_allMuons_id_->Fill(2.);
    if (isTightMuon) h_allMuons_id_->Fill(3.);

    if (!isTightMuon) continue;
    if (fabs(muons->at(i).eta()) > 2.8) continue;
    if (muons->at(i).pt() < 10.) continue;
    h_muons_pt_->Fill(muons->at(i).pt());
    h_muons_eta_->Fill(muons->at(i).eta());
    h_muons_phi_->Fill(muons->at(i).phi());
    h_muons_iso_->Fill(isoMu);
    ++nMuon;

    if (muons->at(i).pt() < 26. || isoMu > muThres_) continue;
    h_goodMuons_pt_->Fill(muons->at(i).pt());
    h_goodMuons_eta_->Fill(muons->at(i).eta());
    h_goodMuons_phi_->Fill(muons->at(i).phi());
    h_goodMuons_iso_->Fill(isoMu);
    ++nGoodMuon;
  }
  h_muons_n_->Fill(nMuon);
  h_goodMuons_n_->Fill(nGoodMuon);
  //PF muons
  int nPFMuon =0;
  int nGoodPFMuon = 0;
  for (size_t i = 0; i < pfCands->size(); i++) {
    if (abs(pfCands->at(i).pdgId()) != 13) continue;
    double isoPFMu = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(pfCands->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.4) continue;
      isoPFMu += pfCandsNoLep->at(k).pt();
    }
    isoPFMu = isoPFMu / pfCands->at(i).pt();
    if (fabs(pfCands->at(i).eta()) > 2.8) continue;
    if (pfCands->at(i).pt() < 10.) continue;
    h_PFMuons_pt_->Fill(pfCands->at(i).pt());
    h_PFMuons_eta_->Fill(pfCands->at(i).eta());
    h_PFMuons_phi_->Fill(pfCands->at(i).phi());
    h_PFMuons_iso_->Fill(isoPFMu);
    ++nPFMuon;

    if (pfCands->at(i).pt() < 26. || isoPFMu > muThres_) continue;
    h_goodPFMuons_pt_->Fill(pfCands->at(i).pt());
    h_goodPFMuons_eta_->Fill(pfCands->at(i).eta());
    h_goodPFMuons_phi_->Fill(pfCands->at(i).phi());
    h_goodPFMuons_iso_->Fill(isoPFMu);
    ++nGoodPFMuon;
  }
  h_PFMuons_n_->Fill(nPFMuon);
  h_goodPFMuons_n_->Fill(nGoodPFMuon);

  // Jets
  int nJet20 = 0;
  int nJet30 = 0;
  int nJet40 = 0;
  int nJet50 = 0;
  for(size_t i = 0; i < jets->size(); i++){
    bool overlaps = false;
    for (size_t j = 0; j < elecs->size(); j++) {
      if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    for (size_t j = 0; j < muons->size(); j++) {
      if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    if (fabs(jets->at(i).eta()) > 4.7) continue;

    if (jets->at(i).pt() < 20.) continue;
    h_jet20_pt_->Fill(jets->at(i).pt());
    h_jet20_eta_->Fill(jets->at(i).eta());
    h_jet20_phi_->Fill(jets->at(i).phi());
    ++nJet20;

    if (jets->at(i).pt() < 30.) continue;
    h_jet30_pt_->Fill(jets->at(i).pt());
    h_jet30_eta_->Fill(jets->at(i).eta());
    h_jet30_phi_->Fill(jets->at(i).phi());
    ++nJet30;

    if (jets->at(i).pt() < 40.) continue;
    h_jet40_pt_->Fill(jets->at(i).pt());
    h_jet40_eta_->Fill(jets->at(i).eta());
    h_jet40_phi_->Fill(jets->at(i).phi());
    ++nJet40;

    if (jets->at(i).pt() < 50.) continue;
    h_jet50_pt_->Fill(jets->at(i).pt());
    h_jet50_eta_->Fill(jets->at(i).eta());
    h_jet50_phi_->Fill(jets->at(i).phi());
    ++nJet50;
  }
  h_jet20_n_->Fill(nJet20);
  h_jet30_n_->Fill(nJet30);
  h_jet40_n_->Fill(nJet40);
  h_jet50_n_->Fill(nJet50);

  // MET 
  h_met_->Fill(met->at(0).pt());

}

// ------------ method to improve ME0 muon ID ----------------
  bool 
BasicRecoDistrib::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = std::abs(chamber->x - segment->x);
          deltaY   = std::abs(chamber->y - segment->y);
          pullX    = std::abs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = std::abs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = std::abs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

bool 
BasicRecoDistrib::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      if (chamber->detector() == 5){

        for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

          LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
          LocalPoint seg_loc_coord(segment->x, segment->y, 0);
          LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
          LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

          const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);

          GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
          GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);

          //double segDPhi = segment->me0SegmentRef->deltaPhi();
          // need to check if this works
          double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
          double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);

          deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
          deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
          deltaPhiBend = std::abs(segDPhi - trackDPhi);

          if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;

        }
      }
    }

  }

  return result;

}

// ------------ loose elec ID -----------
bool 
BasicRecoDistrib::isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  bool isLoose = false;
  double Ooemoop = 999.;
  if (recoEl.ecalEnergy()==0) Ooemoop = 999.;
  else if (!std::isfinite(recoEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1.0/recoEl.ecalEnergy() - recoEl.eSuperClusterOverP()/recoEl.ecalEnergy());
  if (fabs(recoEl.superCluster()->eta()) < 1.479
      && recoEl.full5x5_sigmaIetaIeta() < 0.02992
      && fabs(recoEl.deltaEtaSuperClusterTrackAtVtx()) < 0.004119
      && fabs(recoEl.deltaPhiSuperClusterTrackAtVtx()) < 0.05176
      && recoEl.hcalOverEcal() < 6.741
      && Ooemoop < 73.76
      && recoEl.pfIsolationVariables().sumChargedHadronPt / recoEl.pt() < 2.5
      && !ConversionTools::hasMatchedConversion(recoEl, conversions, beamspot.position())) 
    isLoose = true;

  return (fabs(recoEl.superCluster()->eta()) < 1.556 ? isLoose : (MVAVal > -0.01)); 
}

// ------------ medium elec ID -----------
bool 
BasicRecoDistrib::isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  bool isMedium = false;
  double Ooemoop = 999.;
  if (recoEl.ecalEnergy()==0) Ooemoop = 999.;
  else if (!std::isfinite(recoEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1.0/recoEl.ecalEnergy() - recoEl.eSuperClusterOverP()/recoEl.ecalEnergy());
  if (fabs(recoEl.superCluster()->eta()) < 1.479
      && recoEl.full5x5_sigmaIetaIeta() < 0.01609
      && fabs(recoEl.deltaEtaSuperClusterTrackAtVtx()) < 0.001766
      && fabs(recoEl.deltaPhiSuperClusterTrackAtVtx()) < 0.03130
      && recoEl.hcalOverEcal() < 7.371
      && Ooemoop < 22.6
      && recoEl.pfIsolationVariables().sumChargedHadronPt / recoEl.pt() < 1.325
      && !ConversionTools::hasMatchedConversion(recoEl, conversions, beamspot.position())) 
    isMedium = true;

  return (fabs(recoEl.superCluster()->eta()) < 1.556 ? isMedium : (MVAVal > 0.03));  
}

// ------------ tight elec ID -----------
bool 
BasicRecoDistrib::isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  bool isTight = false;
  double Ooemoop = 999.;
  if (recoEl.ecalEnergy()==0) Ooemoop = 999.;
  else if (!std::isfinite(recoEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1.0/recoEl.ecalEnergy() - recoEl.eSuperClusterOverP()/recoEl.ecalEnergy());
  if (fabs(recoEl.superCluster()->eta()) < 1.479
      && recoEl.full5x5_sigmaIetaIeta() < 0.01614
      && fabs(recoEl.deltaEtaSuperClusterTrackAtVtx()) < 0.001322
      && fabs(recoEl.deltaPhiSuperClusterTrackAtVtx()) < 0.06129
      && recoEl.hcalOverEcal() < 4.492
      && Ooemoop < 18.26
      && recoEl.pfIsolationVariables().sumChargedHadronPt / recoEl.pt() < 1.255
      && !ConversionTools::hasMatchedConversion(recoEl, conversions, beamspot.position())) 
    isTight = true;

  return (fabs(recoEl.superCluster()->eta()) < 1.556 ? isTight : (MVAVal > 0.1)); 
}

// ------------ match reco elec to gen elec ------------
int 
BasicRecoDistrib::matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles) {
  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    // Drop everything that is not electron or not status 1
    if(abs(particle->pdgId()) != 11 || particle->status() != 1)
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( recoEl.p4(), particle->p4() );
    if(dRtmp < dR){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if(!(closestElectron != 0 && dR < 0.1)) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if(ancestorPID == -999 && ancestorStatus == -999){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("SimpleElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }

  if(abs(ancestorPID) > 50 && ancestorStatus == 2)
    return TRUE_NON_PROMPT_ELECTRON;

  if(abs(ancestorPID) == 15 && ancestorStatus == 2)
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}
void 
BasicRecoDistrib::findFirstNonElectronMother(const reco::Candidate *particle,
    int &ancestorPID, int &ancestorStatus) {

  if(particle == 0){
    printf("SimpleElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if(abs(particle->pdgId()) == 11){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

// ------------ tight HGCal electron ID --------------
/*float 
BasicRecoDistrib::evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize) {

  if (fabs(recoEl.superCluster()->eta()) < 1.556) return -1.;

  bool isHGCal = hgcEmId_->setElectronPtr(&recoEl);
  if (isHGCal)  {
    hgcId_startPosition = std::abs(hgcEmId_->getClusterStartPosition().z());
    hgcId_lengthCompatibility = hgcEmId_->getClusterLengthCompatibility();
    hgcId_sigmaietaieta = hgcEmId_->getClusterSigmaEtaEta();
    hgcId_deltaEtaStartPosition = recoEl.trackPositionAtCalo().eta() - hgcEmId_->getClusterStartPosition().eta();
    hgcId_deltaPhiStartPosition = reco::deltaPhi(recoEl.trackPositionAtCalo().phi(), hgcEmId_->getClusterStartPosition().phi());
    hOverE_hgcalSafe = hgcEmId_->getClusterHadronFraction();
    hgcId_cosTrackShowerAngle = recoEl.trackMomentumOut().Unit().Dot(hgcEmId_->getClusterShowerAxis().Unit());
  } else {
    hgcId_startPosition = -1.;
    hgcId_lengthCompatibility = -1.;
    hgcId_sigmaietaieta = recoEl.full5x5_sigmaIetaIeta();
    hgcId_deltaEtaStartPosition = -1.;
    hgcId_deltaPhiStartPosition = -1.;
    hOverE_hgcalSafe = recoEl.hcalOverEcal();
    hgcId_cosTrackShowerAngle = -1.;
  }
  trackIsoR04jurassic_D_pt = (float)isoEl;
  ooEmooP = 1e30;
  if (recoEl.ecalEnergy() == 0) ooEmooP = 1e30;
  else if (!std::isfinite(recoEl.ecalEnergy())) ooEmooP = 1e30;
  else ooEmooP = fabs(1.0/recoEl.ecalEnergy() - recoEl.eSuperClusterOverP()/recoEl.ecalEnergy());
  d0 = recoEl.gsfTrack()->dxy(recoVtx.position());
  dz = recoEl.gsfTrack()->dz(recoVtx.position());
  pt = recoEl.pt();
  etaSC = recoEl.superCluster()->eta();
  phiSC = recoEl.superCluster()->phi();

  expectedMissingInnerHits = (float)recoEl.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  isTrue = (float)matchToTruth(recoEl, genParticles);
  nPV = (float)vertexSize;
  if (!ConversionTools::hasMatchedConversion(recoEl, conversions, beamspot.position())) passConversionVeto = 1.;
  else passConversionVeto = 0.;

  return (isHGCal ? tmvaReader_.EvaluateMVA("PhaseIIEndcapHGCal") : -1.);
}*/


// ------------ method called once each job just before starting event loop  ------------
  void 
BasicRecoDistrib::beginJob()
{
}

// ------------ method called once each run ----------------
  void
BasicRecoDistrib::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
  void
BasicRecoDistrib::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
BasicRecoDistrib::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BasicRecoDistrib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicRecoDistrib);
