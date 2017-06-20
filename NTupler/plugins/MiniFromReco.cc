// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/NTupler
// Class:      MiniFromReco
// 
/**\class MiniFromReco MiniFromReco.cc PhaseTwoAnalysis/NTupler/plugins/MiniFromReco.cc

Description: produces flat ntuples from RECO collections
   - storing gen, reco, and pf leptons with pT > 10 GeV and |eta| < 3
   - storing gen and reco jets with pT > 20 GeV and |eta| < 5

Implementation:
   - lepton isolation needs to be refined
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/UPGTrackerTDRStudies#Muon_identification
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
   - no jet ID is stored
   - b-tagging is not available 


*/
//
// Original Author:  Elvire Bouvier
//         Created:  Tue, 20 Jun 2017 11:27:06 GMT
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
#include "DataFormats/Math/interface/deltaR.h"

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

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

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

class MiniFromReco : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit MiniFromReco(const edm::ParameterSet&);
    ~MiniFromReco();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  

  private:
    virtual void beginJob() override;
    void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
    bool isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    int matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles);
    void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
    float evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    std::unique_ptr<HGCalIDTool> hgcEmId_; 
    TMVA::Reader tmvaReader_;
    float hgcId_startPosition, hgcId_lengthCompatibility, hgcId_sigmaietaieta, hgcId_deltaEtaStartPosition, hgcId_deltaPhiStartPosition, hOverE_hgcalSafe, hgcId_cosTrackShowerAngle, trackIsoR04jurassic_D_pt, ooEmooP, d0, dz, pt, etaSC, phiSC, nPV, expectedMissingInnerHits, passConversionVeto, isTrue;

    edm::EDGetTokenT<std::vector<reco::GsfElectron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<edm::ValueMap<double>> trackIsoValueMapToken_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsToken_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsNoLepToken_;
    edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken_;
    edm::EDGetTokenT<std::vector<reco::PFMET>> metToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;

    TTree *tree_;
    MiniEvent_t ev_;

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
MiniFromReco::MiniFromReco(const edm::ParameterSet& iConfig): 
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
  usesResource("TFileService");

  const edm::ParameterSet& hgcIdCfg = iConfig.getParameterSet("HGCalIDToolConfig");
  auto cc = consumesCollector();
  hgcEmId_.reset( new HGCalIDTool(hgcIdCfg, cc) );

  tmvaReader_.SetOptions("!Color:Silent:!Error");
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

  tree_ = fs_->make<TTree>("events","events");
  createMiniEventTree(tree_,ev_);
}


MiniFromReco::~MiniFromReco()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method to fill gen level event -------------
  void
MiniFromReco::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<reco::GenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  // Jets
  std::vector<size_t> jGenJets;
  ev_.ngj = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
    if (genJets->at(i).pt() < 25.) continue;
    if (fabs(genJets->at(i).eta()) > 5) continue;

    bool overlaps = false;
    for (size_t j = 0; j < genParts->size(); j++) {
      if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
      if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    jGenJets.push_back(i);

    ev_.gj_pt[ev_.ngj]   = genJets->at(i).pt();
    ev_.gj_phi[ev_.ngj]  = genJets->at(i).phi();
    ev_.gj_eta[ev_.ngj]  = genJets->at(i).eta();
    ev_.gj_mass[ev_.ngj] = genJets->at(i).mass();
    ev_.gj_pid[ev_.ngj]  = genJets->at(i).pdgId();
    ev_.ngj++;
  }

  // Leptons
  ev_.ngl = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (genParts->at(i).pt() < 20.) continue;
    if (fabs(genParts->at(i).eta()) > 3.) continue;
    double genIso = 0.;
    for (size_t j = 0; j < jGenJets.size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
      std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
      for (size_t k = 0; k < jconst.size(); k++) {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
        if (deltaR < 0.01) continue;
        if (abs(genParts->at(i).pdgId()) == 13 && deltaR > 0.4) continue;
        if (abs(genParts->at(i).pdgId()) == 11 && deltaR > 0.3) continue;
        genIso = genIso + jconst[k]->pt();
      }
    }
    genIso = genIso / genParts->at(i).pt();
    ev_.gl_pid[ev_.ngl]    = genParts->at(i).pdgId();
    ev_.gl_pt[ev_.ngl]     = genParts->at(i).pt();
    ev_.gl_phi[ev_.ngl]    = genParts->at(i).phi();
    ev_.gl_eta[ev_.ngl]    = genParts->at(i).eta();
    ev_.gl_mass[ev_.ngl]   = genParts->at(i).mass();
    ev_.gl_relIso[ev_.ngl] = genIso; 
    ev_.ngl++;
  }

}

// ------------ method to fill reco level pat -------------
  void
MiniFromReco::recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  hgcEmId_->getEventSetup(iSetup);
  hgcEmId_->getEvent(iEvent);

  Handle<std::vector<reco::GsfElectron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<ValueMap<double>> trackIsoValueMap;
  iEvent.getByToken(trackIsoValueMapToken_, trackIsoValueMap);

  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);

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

  // Leptons
  ev_.nl = 0;

  for(size_t i = 0; i < muons->size(); i++){
    if (muons->at(i).pt() < 10.) continue;
    if (fabs(muons->at(i).eta()) > 3.) continue;

    double isoMu = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(muons->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.3) continue;
      isoMu += pfCandsNoLep->at(k).pt();
    }
    if (muons->at(i).pt() > 0.) isoMu = isoMu / muons->at(i).pt();
    else isoMu = -1.;

    bool isLoose  = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.5));
    bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.3));
    bool isTight  = (fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSel(muons->at(i), 3, 4, 3, 4, 0.1));

    ev_.l_id[ev_.nl]     = (isTight | (isMedium<<1) | (isLoose<<2));
    ev_.l_pid[ev_.nl]    = muons->at(i).charge()*13;
    ev_.l_pt[ev_.nl]     = muons->at(i).pt();
    ev_.l_phi[ev_.nl]    = muons->at(i).phi();
    ev_.l_eta[ev_.nl]    = muons->at(i).eta();
    ev_.l_relIso[ev_.nl] = isoMu;
    ev_.l_g[ev_.nl] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != abs(ev_.l_pid[ev_.nl])) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.l_eta[ev_.nl],ev_.l_phi[ev_.nl]) > 0.4) continue;
      ev_.l_g[ev_.nl]    = ig;
    }
    ev_.nl++;
  }

  for(size_t i = 0; i < elecs->size(); i++) { 
    if (elecs->at(i).pt() < 10.) continue;
    if (fabs(elecs->at(i).eta()) > 3.) continue;

    double isoEl = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(elecs->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.4) continue;
      isoEl += pfCandsNoLep->at(k).pt();
    }
    if (elecs->at(i).pt() > 0.) isoEl = isoEl / elecs->at(i).pt(); 
    else isoEl = -1.;

    Ptr<const reco::GsfElectron> el4iso(elecs,i);
    double eljurassicIso = (*trackIsoValueMap)[el4iso];
    double elpt = elecs->at(i).pt();
    double elMVAVal = -1.;
    if (hgcEmId_->setElectronPtr(&(elecs->at(i)))) 
      elMVAVal = (double)evalMVAElec(elecs->at(i),vertices->at(prVtx),conversions,beamspot,genParts,eljurassicIso/elpt,vertices->size());
    bool isLoose  = isLooseElec(elecs->at(i),conversions,beamspot,elMVAVal);    
    bool isMedium = isMediumElec(elecs->at(i),conversions,beamspot,elMVAVal);    
    bool isTight  = isTightElec(elecs->at(i),conversions,beamspot,elMVAVal);    

    ev_.l_id[ev_.nl]     = (isTight | (isMedium<<1) | (isLoose<<2));
    ev_.l_pid[ev_.nl]    = elecs->at(i).charge()*11;
    ev_.l_pt[ev_.nl]     = elecs->at(i).pt();
    ev_.l_phi[ev_.nl]    = elecs->at(i).phi();
    ev_.l_eta[ev_.nl]    = elecs->at(i).eta();
    ev_.l_relIso[ev_.nl] = isoEl;
    ev_.l_g[ev_.nl] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != abs(ev_.l_pid[ev_.nl])) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.l_eta[ev_.nl],ev_.l_phi[ev_.nl]) > 0.4) continue;
      ev_.l_g[ev_.nl]    = ig;
    }
    ev_.nl++;

  }

  // Jets
  ev_.nj = 0;
  for(size_t i = 0; i < jets->size(); i++){
    if (jets->at(i).pt() < 20.) continue;
    if (fabs(jets->at(i).eta()) > 5) continue;

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

    ev_.j_id[ev_.nl]      = -1;
    ev_.j_pt[ev_.nj]      = jets->at(i).pt();
    ev_.j_phi[ev_.nj]     = jets->at(i).phi();
    ev_.j_eta[ev_.nj]     = jets->at(i).eta();
    ev_.j_mass[ev_.nj]    = jets->at(i).mass();
    ev_.j_csvv2[ev_.nj]   = -1; 
    ev_.j_deepcsv[ev_.nj] = -1;
    ev_.j_flav[ev_.nj]    = -1;
    ev_.j_hadflav[ev_.nj] = -1;
    ev_.j_pid[ev_.nj]     = jets->at(i).pdgId();
    ev_.j_g[ev_.nj] = -1;
    for (int ig = 0; ig < ev_.ngj; ig++) {
      if (reco::deltaR(ev_.gj_eta[ig],ev_.gj_phi[ig],ev_.j_eta[ev_.nj],ev_.j_phi[ev_.nj]) > 0.4) continue;
      ev_.j_g[ev_.nj]     = ig;
      break;
    }	
    ev_.nj++;

  }

  // MET 
  ev_.nmet = 0;
  if (met->size() > 0) {
    ev_.met_pt[ev_.nmet]  = met->at(0).pt();
    ev_.met_phi[ev_.nmet] = met->at(0).phi();
    ev_.nmet++;
  }

  // PF leptons
  ev_.npf = 0;
  for (size_t i = 0; i < pfCands->size(); i++) {
    if (abs(pfCands->at(i).pdgId()) != 11 and abs(pfCands->at(i).pdgId()) != 13) continue;
    if (pfCands->at(i).pt() < 10.) continue;
    if (fabs(pfCands->at(i).eta()) > 3.) continue;

    double isoPF = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(pfCands->at(i).p4(),pfCandsNoLep->at(k).p4()) > (abs(pfCands->at(i).pdgId()) == 11 ? 0.4 : 0.3)) continue;
      isoPF += pfCandsNoLep->at(k).pt();
    }
    isoPF = isoPF / pfCands->at(i).pt();

    ev_.pf_pid[ev_.npf]    = pfCands->at(i).pdgId();
    ev_.pf_pt[ev_.npf]     = pfCands->at(i).pt();
    ev_.pf_eta[ev_.npf]    = pfCands->at(i).eta();
    ev_.pf_phi[ev_.npf]    = pfCands->at(i).phi();
    ev_.pf_mass[ev_.npf]   = pfCands->at(i).mass();
    ev_.pf_relIso[ev_.npf] = isoPF;
    ev_.pf_hp[ev_.npf]     = (pfCands->at(i).trackRef().isNonnull() && pfCands->at(i).trackRef()->quality(reco::Track::highPurity));
    ev_.npf++;

  }
  
}

// ------------ method called for each event  ------------
  void
MiniFromReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //analyze the event
  if(!iEvent.isRealData()) genAnalysis(iEvent, iSetup);
  recoAnalysis(iEvent, iSetup);
  
  //save event if at least one lepton at gen or reco level
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  tree_->Fill();

}

// ------------ method to improve ME0 muon ID ----------------
  bool 
MiniFromReco::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
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

          deltaX   = fabs(chamber->x - segment->x);
          deltaY   = fabs(chamber->y - segment->y);
          pullX    = fabs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = fabs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = fabs(atan(chamber->dXdZ) - atan(segment->dXdZ));

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

// ------------ loose elec ID -----------
bool 
MiniFromReco::isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
MiniFromReco::isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
MiniFromReco::isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
MiniFromReco::matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles) {
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
MiniFromReco::findFirstNonElectronMother(const reco::Candidate *particle,
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
float 
MiniFromReco::evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize) {

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
}


// ------------ method called once each job just before starting event loop  ------------
  void 
MiniFromReco::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MiniFromReco::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniFromReco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniFromReco);
