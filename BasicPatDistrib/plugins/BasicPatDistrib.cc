// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/BasicPatDistrib
// Class:      BasicPatDistrib
// 
/**\class BasicPatDistrib BasicPatDistrib.cc PhaseTwoAnalysis/BasicPatDistrib/plugins/BasicPatDistrib.cc

Description: produces histograms of basic quantities from PAT collections

Implementation:
   - muon isolation comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_isolation
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_identification
   - electron isolation might need to be refined
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
      /!\ no ID is implemented for forward electrons as:
      - PFClusterProducer does not run on miniAOD
      - jurassic isolation needs tracks
   - PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
   - no JEC applied
   - b-tagging WPs come from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#B_tagging 
*/

//
// Original Author:  Elvire Bouvier
//         Created:  Wed, 14 Jun 2017 14:16:22 GMT
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BasicPatDistrib : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit BasicPatDistrib(const edm::ParameterSet&);
    ~BasicPatDistrib();

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

    bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
    bool isME0MuonSelNew(reco::Muon, double, double, double);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    unsigned int pileup_;
    bool useDeepCSV_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    PFJetIDSelectionFunctor jetIDTight_;
    edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
    edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    const ME0Geometry* ME0Geometry_; 
    double mvaThres_;
    double deepThres_;
    double muThres_;

    // MC truth in fiducial phase space
    TH1D* h_genMuons_n_;
    TH1D* h_genMuons_pt_;
    TH1D* h_genMuons_phi_;
    TH1D* h_genMuons_eta_;
    TH1D* h_genMuons_iso_;  
    TH1D* h_genElecs_n_;
    TH1D* h_genElecs_pt_;
    TH1D* h_genElecs_phi_;
    TH1D* h_genElecs_eta_;
    TH1D* h_genElecs_iso_; 
    TH1D* h_genJets_n_;
    TH1D* h_genJets_pt_;
    TH1D* h_genJets_phi_;
    TH1D* h_genJets_eta_;

    // Vertices
    TH1D* h_allVertices_n_;
    // ... that pass ID
    TH1D* h_goodVertices_n_;

    // Muons
    TH1D* h_allMuons_n_;
    TH1D* h_allMuons_pt_;
    TH1D* h_allMuons_phi_;
    TH1D* h_allMuons_eta_;
    TH1D* h_allMuons_iso_;
    TH1D* h_allMuons_id_;
    // ... that pass kin cuts, tight ID, and are isolated
    TH1D* h_goodMuons_n_;
    TH1D* h_goodMuons_pt_;
    TH1D* h_goodMuons_phi_;
    TH1D* h_goodMuons_eta_;
    TH1D* h_goodMuons_iso_;

    // Elecs
    TH1D* h_allElecs_n_;
    TH1D* h_allElecs_pt_;
    TH1D* h_allElecs_phi_;
    TH1D* h_allElecs_eta_;
    TH1D* h_allElecs_iso_;
    TH1D* h_allElecs_id_;
    // ... that pass kin cuts, tight ID, and are isolated
    TH1D* h_goodElecs_n_;
    TH1D* h_goodElecs_pt_;
    TH1D* h_goodElecs_phi_;
    TH1D* h_goodElecs_eta_;
    TH1D* h_goodElecs_iso_;

    // Jets
    TH1D* h_allJets_n_;
    TH1D* h_allJets_pt_;
    TH1D* h_allJets_phi_;
    TH1D* h_allJets_eta_;
    TH1D* h_allJets_disc_;
    TH1D* h_allJets_id_;
    // ... that pass kin cuts, loose ID
    TH1D* h_goodJets_n_;
    TH1D* h_goodJets_nb_;
    TH1D* h_goodJets_pt_;
    TH1D* h_goodJets_phi_;
    TH1D* h_goodJets_eta_;
    TH1D* h_goodJets_disc_;
    TH1D* h_goodLJets_n_;
    TH1D* h_goodLJets_nb_;
    TH1D* h_goodLJets_pt_;
    TH1D* h_goodLJets_phi_;
    TH1D* h_goodLJets_eta_;
    TH1D* h_goodLJets_disc_;
    TH1D* h_goodBJets_n_;
    TH1D* h_goodBJets_nb_;
    TH1D* h_goodBJets_pt_;
    TH1D* h_goodBJets_phi_;
    TH1D* h_goodBJets_eta_;
    TH1D* h_goodBJets_disc_;

    TH1D* h_goodMET_pt_;
    TH1D* h_goodMET_phi_;      
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
BasicPatDistrib::BasicPatDistrib(const edm::ParameterSet& iConfig):
  pileup_(iConfig.getParameter<unsigned int>("pileup")),
  useDeepCSV_(iConfig.getParameter<bool>("useDeepCSV")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),  
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
  jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
  metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
  genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets")))
{
  //now do what ever initialization is needed
  if (pileup_ == 0) {
     mvaThres_ = 0.822;
     deepThres_ = 0.741;
     muThres_ = 0.152;
  } else if (pileup_ == 140) {
    mvaThres_ = 0.864;
    deepThres_ = 0.799;
    muThres_ = 0.204;
  } else if (pileup_ == 200) {
    mvaThres_ = 0.878;
    deepThres_ = 0.821;
    muThres_ = 0.212;
  } else {
    mvaThres_ = -1.;
    deepThres_ = 0.;
    muThres_ = 0.;
  }  

  usesResource("TFileService");

  // MC truth in fiducial phase space
  h_genMuons_n_ = fs_->make<TH1D>("GenMuonsN",";Muon multiplicity;Events / 1", 4, 0., 4.);
  h_genMuons_pt_ = fs_->make<TH1D>("GenMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_genMuons_phi_ = fs_->make<TH1D>("GenMuonsPhi",";#phi(#mu);Events / 0.2", 30, -3., 3.);
  h_genMuons_eta_ = fs_->make<TH1D>("GenMuonsEta",";#eta(#mu);Events / 0.2", 30, -3., 3.);
  h_genMuons_iso_ = fs_->make<TH1D>("GenMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 22, 0., 0.22); 
  h_genElecs_n_ = fs_->make<TH1D>("GenElecsN",";Electron multiplicity;Events / 1", 4, 0., 4.);
  h_genElecs_pt_ = fs_->make<TH1D>("GenElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_genElecs_phi_ = fs_->make<TH1D>("GenElecsPhi",";#phi(e);Events / 0.2", 30, -3., 3.);
  h_genElecs_eta_ = fs_->make<TH1D>("GenElecsEta",";#eta(e);Events / 0.2", 30, -3., 3.);
  h_genElecs_iso_ = fs_->make<TH1D>("GenElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 20, 0., 0.2); 
  h_genJets_n_ = fs_->make<TH1D>("GenJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_genJets_pt_ = fs_->make<TH1D>("GenJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_genJets_phi_ = fs_->make<TH1D>("GenJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_genJets_eta_ = fs_->make<TH1D>("GenJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);

  // Vertices
  if (pileup_ == 200)  h_allVertices_n_ = fs_->make<TH1D>("AllVertices",";Vertex multiplicity;Events / 1", 20, 0., 200.);
  else h_allVertices_n_ = fs_->make<TH1D>("AllVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);
  // ... that pass ID
  if (pileup_ == 200)  h_goodVertices_n_ = fs_->make<TH1D>("GoodVertices",";Vertex multiplicity;Events / 1", 20, 0., 200.);
  else h_goodVertices_n_ = fs_->make<TH1D>("GoodVertices",";Vertex multiplicity;Events / 1", 7, 0., 7.);

  // Muons
  h_allMuons_n_ = fs_->make<TH1D>("AllMuonsN",";Muon multiplicity;Events / 1", 6, 0., 6.);
  h_allMuons_pt_ = fs_->make<TH1D>("AllMuonsPt",";p_{T}(#mu) (GeV);Events / (2 GeV)", 75, 0., 150.);
  h_allMuons_phi_ = fs_->make<TH1D>("AllMuonsPhi",";#phi(#mu);Events / 0.1", 60, -3., 3.);
  h_allMuons_eta_ = fs_->make<TH1D>("AllMuonsEta",";#eta(#mu);Events / 0.1", 60, -3., 3.);
  h_allMuons_iso_ = fs_->make<TH1D>("AllMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 40, 0., 0.4);
  h_allMuons_id_ = fs_->make<TH1D>("AllMuonsID",";;Muons / 1", 4, 0., 4.);
  h_allMuons_id_->SetOption("bar");
  h_allMuons_id_->SetBarWidth(0.75);
  h_allMuons_id_->SetBarOffset(0.125);
  h_allMuons_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allMuons_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allMuons_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allMuons_id_->GetXaxis()->SetBinLabel(4,"Tight");
  // ... that pass kin cuts, tight ID, and are isolated
  h_goodMuons_n_ = fs_->make<TH1D>("GoodMuonsN",";Muon multiplicity;Events / 1", 4, 0., 4.);
  h_goodMuons_pt_ = fs_->make<TH1D>("GoodMuonsPt",";p_{T}(#mu) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_goodMuons_phi_ = fs_->make<TH1D>("GoodMuonsPhi",";#phi(#mu);Events / 0.2", 30, -3., 3.);
  h_goodMuons_eta_ = fs_->make<TH1D>("GoodMuonsEta",";#eta(#mu);Events / 0.2", 30, -3., 3.);
  h_goodMuons_iso_ = fs_->make<TH1D>("GoodMuonsIso",";I_{rel}^{PUPPI}(#mu);Events / 0.01", 22, 0., 0.22);

  // Elecs
  h_allElecs_n_ = fs_->make<TH1D>("AllElecsN",";Electron multiplicity;Events / 1", 6, 0., 6.);
  h_allElecs_pt_ = fs_->make<TH1D>("AllElecsPt",";p_{T}(e) (GeV);Events / (2 GeV)", 75, 0., 150.);
  h_allElecs_phi_ = fs_->make<TH1D>("AllElecsPhi",";#phi(e);Events / 0.1", 60, -3., 3.);
  h_allElecs_eta_ = fs_->make<TH1D>("AllElecsEta",";#eta(e);Events / 0.1", 60, -3., 3.);
  h_allElecs_iso_ = fs_->make<TH1D>("AllElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 40, 0., 0.4);
  h_allElecs_id_ = fs_->make<TH1D>("AllElecsID",";;Electrons / 1", 4, 0., 4.);
  h_allElecs_id_->SetOption("bar");
  h_allElecs_id_->SetBarWidth(0.75);
  h_allElecs_id_->SetBarOffset(0.125);
  h_allElecs_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allElecs_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allElecs_id_->GetXaxis()->SetBinLabel(3,"Medium");
  h_allElecs_id_->GetXaxis()->SetBinLabel(4,"Tight");
  // ... that pass kin cuts, tight ID, and are isolated
  h_goodElecs_n_ = fs_->make<TH1D>("GoodElecsN",";Electron multiplicity;Events / 1", 4, 0., 4.);
  h_goodElecs_pt_ = fs_->make<TH1D>("GoodElecsPt",";p_{T}(e) (GeV);Events / (5 GeV)", 26, 20., 150.);
  h_goodElecs_phi_ = fs_->make<TH1D>("GoodElecsPhi",";#phi(e);Events / 0.2", 30, -3., 3.);
  h_goodElecs_eta_ = fs_->make<TH1D>("GoodElecsEta",";#eta(e);Events / 0.2", 30, -3., 3.);
  h_goodElecs_iso_ = fs_->make<TH1D>("GoodElecsIso",";I_{rel}^{PUPPI}(e);Events / 0.01", 20, 0., 0.2);

  // Jets
  h_allJets_n_ = fs_->make<TH1D>("AllJetsN",";Jet multiplicity;Events / 1", 15, 0., 15.);
  h_allJets_pt_ = fs_->make<TH1D>("AllJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 100, 0., 200.);
  h_allJets_phi_ = fs_->make<TH1D>("AllJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_allJets_eta_ = fs_->make<TH1D>("AllJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_allJets_disc_ = fs_->make<TH1D>("AllJetsDisc",";b-tagging discriminant;Events / 0.02", 50, 0., 1.);
  h_allJets_id_ = fs_->make<TH1D>("AllJetsID",";;Jets / 1", 3, 0., 3.);
  h_allJets_id_->SetOption("bar");
  h_allJets_id_->SetBarWidth(0.75);
  h_allJets_id_->SetBarOffset(0.125);
  h_allJets_id_->GetXaxis()->SetBinLabel(1,"All");
  h_allJets_id_->GetXaxis()->SetBinLabel(2,"Loose");
  h_allJets_id_->GetXaxis()->SetBinLabel(3,"Tight");
  // ... that pass kin cuts, loose ID
  h_goodJets_n_ = fs_->make<TH1D>("GoodJetsN",";Jet multiplicity;Events / 1", 14, 0., 14.);
  h_goodJets_nb_ = fs_->make<TH1D>("GoodJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodJets_pt_ = fs_->make<TH1D>("GoodJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodJets_phi_ = fs_->make<TH1D>("GoodJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodJets_eta_ = fs_->make<TH1D>("GoodJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodJets_disc_ = fs_->make<TH1D>("GoodJetsDisc",";b-tagging discriminant;Events / 0.02", 50, 0., 1.);
  h_goodLJets_n_ = fs_->make<TH1D>("GoodLightJetsN",";Jet multiplicity;Events / 1", 12, 0., 12.);
  h_goodLJets_nb_ = fs_->make<TH1D>("GoodLightJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodLJets_pt_ = fs_->make<TH1D>("GoodLightJetsPt",";p_{T}(jet) (GeV);Events / (2 GeV)", 90, 20., 200.);
  h_goodLJets_phi_ = fs_->make<TH1D>("GoodLightJetsPhi",";#phi(jet);Events / 0.1", 60, -3., 3.);
  h_goodLJets_eta_ = fs_->make<TH1D>("GoodLightJetsEta",";#eta(jet);Events / 0.1", 100, -5., 5.);
  h_goodLJets_disc_ = fs_->make<TH1D>("GoodLightJetsDisc",";b-tagging discriminant;Events / 0.02", 50, 0., 1.);
  h_goodBJets_n_ = fs_->make<TH1D>("GoodBtaggedJetsN",";Jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_nb_ = fs_->make<TH1D>("GoodBtaggedJetsNb",";b jet multiplicity;Events / 1", 5, 0., 5.);
  h_goodBJets_pt_ = fs_->make<TH1D>("GoodBtaggedJetsPt",";p_{T}(jet) (GeV);Events / (5 GeV)", 36, 20., 200.);
  h_goodBJets_phi_ = fs_->make<TH1D>("GoodBtaggedJetsPhi",";#phi(jet);Events / 0.2", 30, -3., 3.);
  h_goodBJets_eta_ = fs_->make<TH1D>("GoodBtaggedJetsEta",";#eta(jet);Events / 0.2", 50, -5., 5.);
  h_goodBJets_disc_ = fs_->make<TH1D>("GoodBtaggedJetsDisc",";b-tagging discriminant;Events / 0.01", 20, 0.8, 1.);

  // MET
  h_goodMET_pt_ = fs_->make<TH1D>("GoodMETPt",";p_{T}(MET) (GeV);Events / (5 GeV)", 60, 0., 300.);
  h_goodMET_phi_ = fs_->make<TH1D>("GoodMETPhi",";#phi(MET);Events / 0.2", 30, -3., 3.);


}


BasicPatDistrib::~BasicPatDistrib()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
BasicPatDistrib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);

  Handle<std::vector<pat::Electron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();  

  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<std::vector<pat::MET>> mets;
  iEvent.getByToken(metsToken_, mets);

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<pat::PackedGenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  // Vertices
  int prVtx = -1;
  size_t nVtx = 0;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
    ++ nVtx;
  }
  if (prVtx < 0) return;
  h_goodVertices_n_->Fill(nVtx);
  h_allVertices_n_->Fill(vertices->size());
   
  // MC truth in fiducial phase space
  std::vector<size_t> jGenJets;
  size_t nGenJets = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < genParts->size(); j++) {
      if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    jGenJets.push_back(i);

    if (genJets->at(i).pt() < 30.) continue;
    if (fabs(genJets->at(i).eta()) > 4.7) continue;
    h_genJets_pt_->Fill(genJets->at(i).pt());
    h_genJets_phi_->Fill(genJets->at(i).phi());
    h_genJets_eta_->Fill(genJets->at(i).eta());
    ++nGenJets;
  }
  h_genJets_n_->Fill(nGenJets);

  size_t nGenMuons = 0;
  size_t nGenElecs = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (fabs(genParts->at(i).eta()) > 2.8) continue;
    double genIso = 0.;
    for (size_t j = 0; j < jGenJets.size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
      std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
      for (size_t k = 0; k < jconst.size(); k++) {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
        if (deltaR < 0.01 || deltaR > 0.4) continue;
        genIso = genIso + jconst[k]->pt();
      }
    }
    genIso = genIso / genParts->at(i).pt();
    if (abs(genParts->at(i).pdgId()) == 13) {
      if (genIso > muThres_) continue;
      if (genParts->at(i).pt() < 26.) continue;
      h_genMuons_pt_->Fill(genParts->at(i).pt());
      h_genMuons_phi_->Fill(genParts->at(i).phi());
      h_genMuons_eta_->Fill(genParts->at(i).eta());
      h_genMuons_iso_->Fill(genIso); 
      ++nGenMuons;
    }
    if (abs(genParts->at(i).pdgId()) == 11) {
      if (genIso > 0.15) continue;
      if (genParts->at(i).pt() < 30.) continue;
      h_genElecs_pt_->Fill(genParts->at(i).pt());
      h_genElecs_phi_->Fill(genParts->at(i).phi());
      h_genElecs_eta_->Fill(genParts->at(i).eta());
      h_genElecs_iso_->Fill(genIso); 
      ++nGenElecs;
    }
  }
  h_genMuons_n_->Fill(nGenMuons);
  h_genElecs_n_->Fill(nGenElecs);

  // Muons
  size_t nGoodMuons = 0;
  for (size_t i = 0; i < muons->size(); i++) {
    h_allMuons_pt_->Fill(muons->at(i).pt());
    h_allMuons_phi_->Fill(muons->at(i).phi());
    h_allMuons_eta_->Fill(muons->at(i).eta());
    h_allMuons_iso_->Fill((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt());

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

    if (muons->at(i).pt() < 26.) continue;
    if (fabs(muons->at(i).eta()) > 2.8) continue;
    if ((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt() > muThres_) continue;
    if (!isTightMuon) continue;
    h_goodMuons_pt_->Fill(muons->at(i).pt());
    h_goodMuons_phi_->Fill(muons->at(i).phi());
    h_goodMuons_eta_->Fill(muons->at(i).eta());
    h_goodMuons_iso_->Fill((muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt());
    ++nGoodMuons;
  }
  h_goodMuons_n_->Fill(nGoodMuons);
  h_allMuons_n_->Fill(muons->size());
 
  // Electrons
  size_t nGoodElecs = 0;
  for (size_t i = 0; i < elecs->size(); i++) {
    h_allElecs_pt_->Fill(elecs->at(i).pt());
    h_allElecs_phi_->Fill(elecs->at(i).phi());
    h_allElecs_eta_->Fill(elecs->at(i).eta());
    h_allElecs_iso_->Fill((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt());
    h_allElecs_id_->Fill(0.);
    if (isLooseElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(1.);    
    if (isMediumElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(2.);    
    if (isTightElec(elecs->at(i),conversions,beamspot)) h_allElecs_id_->Fill(3.);    

    if (elecs->at(i).pt() < 30.) continue;
    if (fabs(elecs->at(i).eta()) > 2.8) continue;
    if ((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt() > 0.15) continue;
    if (!isTightElec(elecs->at(i),conversions,beamspot)) continue;    
    h_goodElecs_pt_->Fill(elecs->at(i).pt());
    h_goodElecs_phi_->Fill(elecs->at(i).phi());
    h_goodElecs_eta_->Fill(elecs->at(i).eta());
    h_goodElecs_iso_->Fill((elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt());
    ++nGoodElecs;
  }
  h_goodElecs_n_->Fill(nGoodElecs);
  h_allElecs_n_->Fill(elecs->size());
  
  // Jets
  size_t nGoodJets = 0;
  size_t nbGoodJets = 0;
  size_t nGoodLightJets = 0;
  size_t nbGoodLightJets = 0;
  size_t nGoodBtaggedJets = 0;
  size_t nbGoodBtaggedJets = 0;
  for (size_t i =0; i < jets->size(); i++) {
    bool overlaps = false;
    for (size_t j = 0; j < elecs->size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    for (size_t j = 0; j < muons->size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;

    double btagDisc = -1.;
    if (useDeepCSV_)
        btagDisc = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    else
        btagDisc = jets->at(i).bDiscriminator("pfCombinedMVAV2BJetTags");
    
    h_allJets_pt_->Fill(jets->at(i).pt());
    h_allJets_phi_->Fill(jets->at(i).phi());
    h_allJets_eta_->Fill(jets->at(i).eta());
    h_allJets_disc_->Fill(btagDisc); 
    h_allJets_id_->Fill(0.);
    pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    retLoose.set(false);
    if (jetIDLoose_(jets->at(i), retLoose)) h_allJets_id_->Fill(1.);
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    if (jetIDTight_(jets->at(i), retTight)) h_allJets_id_->Fill(2.);

    if (jets->at(i).pt() < 30.) continue;
    if (fabs(jets->at(i).eta()) > 4.7) continue;
    if (!jetIDLoose_(jets->at(i), retLoose)) continue;
    h_goodJets_pt_->Fill(jets->at(i).pt());
    h_goodJets_phi_->Fill(jets->at(i).phi());
    h_goodJets_eta_->Fill(jets->at(i).eta());
    h_goodJets_disc_->Fill(btagDisc); 
    ++nGoodJets;
    if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodJets;
    if ((useDeepCSV_ && btagDisc > deepThres_)
            || (!useDeepCSV_ && btagDisc > mvaThres_)){  
      h_goodBJets_pt_->Fill(jets->at(i).pt());
      h_goodBJets_phi_->Fill(jets->at(i).phi());
      h_goodBJets_eta_->Fill(jets->at(i).eta());
      h_goodBJets_disc_->Fill(btagDisc);
      ++nGoodBtaggedJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodBtaggedJets;
    } else {
      h_goodLJets_pt_->Fill(jets->at(i).pt());
      h_goodLJets_phi_->Fill(jets->at(i).phi());
      h_goodLJets_eta_->Fill(jets->at(i).eta());
      h_goodLJets_disc_->Fill(btagDisc); 
      ++nGoodLightJets;
      if (jets->at(i).genParton() && fabs(jets->at(i).genParton()->pdgId()) == 5) ++nbGoodLightJets;
    }
  }
  h_goodLJets_n_->Fill(nGoodLightJets);
  h_goodLJets_nb_->Fill(nbGoodLightJets);
  h_goodBJets_n_->Fill(nGoodBtaggedJets);
  h_goodBJets_nb_->Fill(nbGoodBtaggedJets);
  h_goodJets_n_->Fill(nGoodJets);
  h_goodJets_nb_->Fill(nbGoodJets);
  h_allJets_n_->Fill(jets->size());
  
  // MET
  if (mets->size() > 0) {
    h_goodMET_pt_->Fill(mets->at(0).pt());
    h_goodMET_phi_->Fill(mets->at(0).phi());
  }

}


// ------------ method check that an e passes loose ID ----------------------------------
  bool
BasicPatDistrib::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
  if (patEl.hcalOverEcal() > 6.741) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 73.76) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
  bool
BasicPatDistrib::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
  if (patEl.hcalOverEcal() > 7.371) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 22.6) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
  bool
BasicPatDistrib::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
  if (patEl.hcalOverEcal() > 4.492) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 18.26) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method to improve ME0 muon ID ----------------
  bool 
BasicPatDistrib::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
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
BasicPatDistrib::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
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

// ------------ method called once each job just before starting event loop  ------------
  void 
BasicPatDistrib::beginJob()
{
}

// ------------ method called once each run ----------------
void
BasicPatDistrib::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
  void
BasicPatDistrib::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
BasicPatDistrib::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BasicPatDistrib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicPatDistrib);
