// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/NTupler
// Class:      MiniFromPat
// 
/**\class MiniFromPat MiniFromPat.cc PhaseTwoAnalysis/NTupler/plugins/MiniFromPat.cc

Description: produces flat ntuples from PAT collections
   - storing gen, reco, and pf leptons with pT > 10 GeV and |eta| < 3
   - storing gen and reco jets with pT > 20 GeV and |eta| < 5

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
//         Created:  Tue, 20 Jun 2017 11:27:12 GMT
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
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
#include "DataFormats/Math/interface/deltaR.h"

#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"


static inline float matchClosestDr(
		float* geneta,float* genphi,float* genpt,
		int ngen,
		float* recoeta,float* recophi,float* recopt,
		int recoidx,
		int& matched,
		float maxDr,
		float maxDeltaRelPt=0.5,
		int * genpdgid=0,
		int reqabspdgid=-1){

	matched=-1;
	float mindr=maxDr;
	for (int ig = 0; ig < ngen; ig++) {
		if(genpdgid){
			if(std::abs(genpdgid[ig])!=reqabspdgid)continue;
		}
		float thisdr=reco::deltaR(geneta[ig],genphi[ig],recoeta[recoidx],recophi[recoidx]);
		if (thisdr > mindr) continue;
		if(fabs(genpt[ig]-recopt[recoidx]) / genpt[ig] > maxDeltaRelPt)continue;
		mindr=thisdr;
		matched     = ig;
	}
	return mindr;
}

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniFromPat : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit MiniFromPat(const edm::ParameterSet&);
    ~MiniFromPat();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  


  private:
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endJob() override;

    bool isME0MuonSelNew(const reco::Muon&, double, double, double, edm::EventSetup const& );

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    unsigned int pileup_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Tau>> tausToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    PFJetIDSelectionFunctor jetIDTight_;
    edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
    edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
    edm::EDGetTokenT<std::vector<pat::Photon>> photonsToken_;
    edm::EDGetTokenT<EcalRecHitCollection> ecalRecHitsToken_;
    const ME0Geometry* ME0Geometry_; 
    double mvaThres_[3];
    double deepThres_[3];
    TH2F * photonEcorr_{nullptr};

    TTree *t_event_, *t_genParts_, *t_vertices_, *t_genJets_, *t_genPhotons_, *t_looseElecs_, *t_mediumElecs_, *t_tightElecs_, *t_looseMuons_, *t_tightMuons_, *t_allTaus_,*t_puppiJets_, *t_puppiMET_, *t_loosePhotons_, *t_tightPhotons_;

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
MiniFromPat::MiniFromPat(const edm::ParameterSet& iConfig):
  pileup_(iConfig.getParameter<unsigned int>("pileup")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),  
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  tausToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
  jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
  metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  photonsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  ecalRecHitsToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalRecHits")))
{
	ME0Geometry_=0;
  //now do what ever initialization is needed
  if (pileup_ == 0) {
    mvaThres_[0] = -0.694;
    mvaThres_[1] = 0.128;
    mvaThres_[2] = 0.822;
    deepThres_[0] = 0.131;
    deepThres_[1] = 0.432;
    deepThres_[2] = 0.741;
  } else if (pileup_ == 140) {
    mvaThres_[0] = -0.654;
    mvaThres_[1] = 0.214;
    mvaThres_[2] = 0.864;
    deepThres_[0] = 0.159;
    deepThres_[1] = 0.507;
    deepThres_[2] = 0.799;
  } else if (pileup_ == 200) {
    mvaThres_[0] = -0.642;
    mvaThres_[1] = 0.236;
    mvaThres_[2] = 0.878;
    deepThres_[0] = 0.170;
    deepThres_[1] = 0.527;
    deepThres_[2] = 0.821;
  } else {
    mvaThres_[0] = -1.;
    mvaThres_[1] = -1.;
    mvaThres_[2] = -1.;
    deepThres_[0] = 0.;
    deepThres_[1] = 0.;
    deepThres_[2] = 0.;
  }  

  // Load photon correction
  if ( iConfig.existsAs<edm::FileInPath>("photonEcorr") ) {
    auto photonEcorrFile = iConfig.getParameter<edm::FileInPath>("photonEcorr");
    TFile * photonEcorr = TFile::Open(photonEcorrFile.fullPath().c_str());
    photonEcorr_ = (TH2F*) photonEcorr->Get("combinedECorrection");
    photonEcorr_->SetDirectory(0);  // don't delete
    delete photonEcorr;
  }

  usesResource("TFileService");

  t_event_      = fs_->make<TTree>("Event","Event");
  t_genParts_   = fs_->make<TTree>("Particle","Particle");
  t_genPhotons_ = fs_->make<TTree>("GenPhoton","GenPhoton"); 
  t_vertices_   = fs_->make<TTree>("Vertex","Vertex");
  t_genJets_    = fs_->make<TTree>("GenJet","GenJet");
  t_looseElecs_ = fs_->make<TTree>("ElectronLoose","ElectronLoose");
  t_mediumElecs_ = fs_->make<TTree>("ElectronMedium","ElectronMedium");
  t_tightElecs_ = fs_->make<TTree>("ElectronTight","ElectronTight");
  t_looseMuons_ = fs_->make<TTree>("MuonLoose","MuonLoose");
  t_tightMuons_ = fs_->make<TTree>("MuonTight","MuonTight");
  t_allTaus_ = fs_->make<TTree>("TauAll","TauAll");
  if (pileup_!=0)  t_puppiJets_  = fs_->make<TTree>("JetPUPPI","JetPUPPI");
  else   t_puppiJets_  = fs_->make<TTree>("Jet","Jet");
  if (pileup_!=0)  t_puppiMET_   = fs_->make<TTree>("PuppiMissingET","PuppiMissingET");
  else   t_puppiMET_   = fs_->make<TTree>("MissingET","MissingET");
  t_loosePhotons_ = fs_->make<TTree>("PhotonLoose","PhotonLoose");
  t_tightPhotons_ = fs_->make<TTree>("PhotonTight","PhotonTight");
  createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_genPhotons_, t_looseElecs_,
		  t_mediumElecs_,t_tightElecs_, t_looseMuons_, t_tightMuons_, t_allTaus_,t_puppiJets_, t_puppiMET_,
		  t_loosePhotons_, t_tightPhotons_, ev_);

}


MiniFromPat::~MiniFromPat()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method to fill gen level pat -------------
  void
MiniFromPat::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<pat::PackedGenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  // Jets
  std::vector<size_t> jGenJets;
  ev_.ngj = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
	if (ev_.ngj>=MiniEvent_t::maxjets) break;
    if (genJets->at(i).pt() < 20.) continue;
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
    ev_.ngj++;
  }

  // Leptons
  ev_.ngl = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
	if (ev_.ngl>=MiniEvent_t::maxpart) break;
	if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (genParts->at(i).pt() < 10.) continue;
    if (fabs(genParts->at(i).eta()) > 5.) continue;
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
    ev_.gl_pid[ev_.ngl]    = genParts->at(i).pdgId();
    ev_.gl_ch[ev_.ngl]     = genParts->at(i).charge();
    ev_.gl_st[ev_.ngl]     = genParts->at(i).status();
    ev_.gl_p[ev_.ngl]      = genParts->at(i).p();
    ev_.gl_px[ev_.ngl]     = genParts->at(i).px();
    ev_.gl_py[ev_.ngl]     = genParts->at(i).py();
    ev_.gl_pz[ev_.ngl]     = genParts->at(i).pz();
    ev_.gl_nrj[ev_.ngl]    = genParts->at(i).energy();
    ev_.gl_pt[ev_.ngl]     = genParts->at(i).pt();
    ev_.gl_phi[ev_.ngl]    = genParts->at(i).phi();
    ev_.gl_eta[ev_.ngl]    = genParts->at(i).eta();
    ev_.gl_mass[ev_.ngl]   = genParts->at(i).mass();
    ev_.gl_relIso[ev_.ngl] = genIso; 
    ev_.ngl++;
  }

  // Photons
  ev_.ngp = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
	if (ev_.ngp>=MiniEvent_t::maxpart) break;
    if (abs(genParts->at(i).pdgId()) != 22) continue;
    if (genParts->at(i).pt() < 10.) continue;
    if (fabs(genParts->at(i).eta()) > 3.) continue;

    ev_.gp_st[ev_.ngp]     = genParts->at(i).status();
    ev_.gp_p[ev_.ngp]      = genParts->at(i).p();
    ev_.gp_px[ev_.ngp]     = genParts->at(i).px();
    ev_.gp_py[ev_.ngp]     = genParts->at(i).py();
    ev_.gp_pz[ev_.ngp]     = genParts->at(i).pz();
    ev_.gp_nrj[ev_.ngp]    = genParts->at(i).energy();
    ev_.gp_pt[ev_.ngp]     = genParts->at(i).pt();
    ev_.gp_phi[ev_.ngp]    = genParts->at(i).phi();
    ev_.gp_eta[ev_.ngp]    = genParts->at(i).eta();
    ev_.ngp++;
  }
  
  // Generator weights
  ev_.g_nw = 0; ev_.g_w[0] = 1.0;
  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);
  if(evt.isValid()) {
    ev_.g_w[0] = evt->weight();
    ev_.g_nw++;
    for (unsigned int i = 1; i<evt->weights().size(); i++) {
      if (ev_.g_nw>=MiniEvent_t::maxweights) break;
      ev_.g_w[ev_.g_nw]=evt->weights()[i];
      ev_.g_nw++;
    }
  }
  // LHE weights
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);
  if(evet.isValid()) {
    double asdd=evet->originalXWGTUP();
    for(unsigned int i=0  ; i<evet->weights().size();i++) {
      if (ev_.g_nw>=MiniEvent_t::maxweights) break;
      double asdde=evet->weights()[i].wgt;
      ev_.g_w[ev_.g_nw]=ev_.g_w[0]*asdde/asdd;
      ev_.g_nw++;
    }
  }
}

// ------------ method to fill reco level pat -------------
  void
MiniFromPat::recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  Handle<std::vector<pat::Tau>> taus;
  iEvent.getByToken(tausToken_, taus);

  Handle<std::vector<pat::MET>> mets;
  iEvent.getByToken(metsToken_, mets);

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<pat::Photon>> photons;
  iEvent.getByToken(photonsToken_, photons);

  Handle<EcalRecHitCollection> ecalRecHits;
  iEvent.getByToken(ecalRecHitsToken_, ecalRecHits);
   
  // Vertices
  int prVtx = -1;
  ev_.nvtx = 0;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
    ev_.v_pt2[ev_.nvtx] = vertices->at(i).p4().pt();
    ev_.nvtx++;
  }
  if (prVtx < 0) return;
  auto primaryVertex=vertices->at(prVtx);

  // Taus 
  ev_.ntau=0;

  for (size_t i = 0; i < taus->size(); i++) {
    if (taus->at(i).pt()<15.) continue; 
    if (fabs(taus->at(i).eta()) > 3.0) continue;
    if (taus->at(i).tauID("decayModeFinding")<0) continue;     
 
    if (ev_.ntau<MiniEvent_t::maxpart){

       ev_.tau_ch[ev_.ntau]     = taus->at(i).charge();
       ev_.tau_pt[ev_.ntau]     = taus->at(i).pt();
       ev_.tau_phi[ev_.ntau]    = taus->at(i).phi();
       ev_.tau_eta[ev_.ntau]    = taus->at(i).eta();
       ev_.tau_mass[ev_.ntau]   = taus->at(i).mass();
       ev_.tau_dm[ev_.ntau]   = taus->at(i).decayMode();
       ev_.tau_chargedIso[ev_.ntau] = taus->at(i).tauID("chargedIsoPtSum");
       ev_.tau_g[ev_.ntau] = -1;
       for (int ig = 0; ig < ev_.ngl; ig++) {
         // I need to fix it matching to genTau
         if (abs(ev_.gl_pid[ig]) != 15) continue;
         if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.tau_eta[ev_.ntau],ev_.tau_phi[ev_.ntau]) > 0.4) continue;
         ev_.tau_g[ev_.ntau]    = ig;
       }
       ev_.ntau++;
    }
  }

  // Muons
  ev_.nlm = 0;
  ev_.ntm = 0;

  for (size_t i = 0; i < muons->size(); i++) {
    if (muons->at(i).pt() < 2.) continue;
    if (fabs(muons->at(i).eta()) > 2.8) continue;

    // Loose ID
    double dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.056);
    double dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0096);    
    bool isLoose = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i)))
    		|| (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut,iSetup));

    // Medium ID -- needs to be updated
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    float dz=0, dxy=0;
    if (muons->at(i).innerTrack().isNonnull()){
    	dxy=std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position()));
    	dz= std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position()));
    	ipxy = dxy < 0.2;
    	ipz = dz < 0.5;
    	validPxlHit = muons->at(i).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    	highPurity = muons->at(i).innerTrack()->quality(reco::Track::highPurity);
    }    
    // bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    // Tight ID
    dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.032);
    dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0041);
    bool isTight = (fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx)))
    		|| (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.048, dPhiCut, dPhiBendCut,iSetup) && ipxy && ipz && validPxlHit && highPurity);

    float trackIso03=muons->at(i).trackIso()/muons->at(i).pt();

    if (!isLoose) continue;
    if (ev_.nlm<MiniEvent_t::maxpart){

       ev_.lm_ch[ev_.nlm]     = muons->at(i).charge();
       ev_.lm_pt[ev_.nlm]     = muons->at(i).pt();
       ev_.lm_phi[ev_.nlm]    = muons->at(i).phi();
       ev_.lm_eta[ev_.nlm]    = muons->at(i).eta();
       ev_.lm_mass[ev_.nlm]   = muons->at(i).mass();
       ev_.lm_dxy[ev_.nlm]   = dxy;
       ev_.lm_dz[ev_.nlm]   = dz;
       ev_.lm_relIso[ev_.nlm] = trackIso03;//(muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
       ev_.lm_g[ev_.nlm] = -1;
       matchClosestDr(ev_.gl_eta,
           		ev_.gl_phi,
       			ev_.gl_pt,
       			ev_.ngl,
       			ev_.lm_eta,
       			ev_.lm_phi,
       			ev_.lm_pt,
       			ev_.nlm,
       			ev_.lm_g[ev_.nlm],
       			0.1,
       			1, ev_.gl_pid, 13);

       ev_.nlm++;
    }

    if (!isTight) continue;

    if (ev_.ntm>=MiniEvent_t::maxpart) break;

    ev_.tm_ch[ev_.ntm]     = muons->at(i).charge();
    ev_.tm_pt[ev_.ntm]     = muons->at(i).pt();
    ev_.tm_phi[ev_.ntm]    = muons->at(i).phi();
    ev_.tm_eta[ev_.ntm]    = muons->at(i).eta();
    ev_.tm_mass[ev_.ntm]   = muons->at(i).mass();
    ev_.tm_dxy[ev_.ntm]   = dxy;
    ev_.tm_dz[ev_.ntm]   = dz;
    ev_.tm_relIso[ev_.ntm] = trackIso03;//(muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
    ev_.tm_g[ev_.ntm] = -1;
    matchClosestDr(ev_.gl_eta,
        		ev_.gl_phi,
    			ev_.gl_pt,
    			ev_.ngl,
    			ev_.tm_eta,
    			ev_.tm_phi,
    			ev_.tm_pt,
    			ev_.ntm,
    			ev_.tm_g[ev_.ntm],
    			0.1,
    			1, ev_.gl_pid, 13);
    ev_.ntm++;
  }

  // Electrons

  ev_.nle = 0;
  ev_.nme = 0;
  ev_.nte = 0;

  for (size_t i = 0; i < elecs->size(); i++) {
    if (elecs->at(i).pt() < 10.) continue;
    if (fabs(elecs->at(i).eta()) > 3.) continue;

    float mvaValue = elecs->at(i).userFloat("mvaValue");
    bool isEB = elecs->at(i).isEB();
     
    bool isLoose = 0;
    bool isMedium = 0;
    bool isTight = 0;

    if( isEB ) {
      if (elecs->at(i).pt() < 20.) {
        isLoose = (mvaValue > -0.661);
        isMedium = mvaValue > 0.855;
        isTight = (mvaValue > 0.986);
      }
      else {
        isLoose = (mvaValue > -0.797);
        isMedium = mvaValue > 0.723;
        isTight = (mvaValue > 0.988);
      }
    }
    else {
      if (not (elecs->at(i).userFloat("hgcElectronID:ecEnergy") > 0)) continue;
      if (not (elecs->at(i).userFloat("hgcElectronID:sigmaUU") > 0)) continue;
      if (not (elecs->at(i).fbrem() > -1)) continue;
      if (not (elecs->at(i).userFloat("hgcElectronID:measuredDepth") < 40)) continue;
      if (not (elecs->at(i).userFloat("hgcElectronID:nLayers") > 20)) continue;
      if (elecs->at(i).pt() < 20.) {
        isLoose = (mvaValue > -0.320);
        isMedium = mvaValue > 0.777;
        isTight = (mvaValue > 0.969);
      }
      else {
        isLoose = (mvaValue > -0.919);
        isMedium = mvaValue > 0.591;
        isTight = (mvaValue > 0.983);
      }
    }
    float dxy=0;
    float dz=0;
    if(elecs->at(i).gsfTrack().isNonnull()){
    	dxy=std::abs(elecs->at(i).gsfTrack()->dxy(primaryVertex.position()));
    	dz=std::abs(elecs->at(i).gsfTrack()->dz(primaryVertex.position()));
    }

    if (!isLoose) continue;

    if (ev_.nle<MiniEvent_t::maxpart){
       ev_.le_ch[ev_.nle]     = elecs->at(i).charge();
       ev_.le_pt[ev_.nle]     = elecs->at(i).pt();
       ev_.le_phi[ev_.nle]    = elecs->at(i).phi();
       ev_.le_eta[ev_.nle]    = elecs->at(i).eta();
       ev_.le_mass[ev_.nle]   = elecs->at(i).mass();
       ev_.le_bdt[ev_.nle]   = mvaValue;
       ev_.le_dz[ev_.nle]    = dz;
       ev_.le_dxy[ev_.nle]   = dxy;
       if( isEB )
         ev_.le_relIso[ev_.nle] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
       else
         ev_.le_relIso[ev_.nle] = (elecs->at(i).userFloat("hgcElectronID:caloIsoRing1") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing2") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing3") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing4")) / elecs->at(i).energy();
       ev_.le_relTkIso[ev_.nle] = elecs->at(i).dr03TkSumPt()/elecs->at(i).pt();

       matchClosestDr(ev_.gl_eta,
           		ev_.gl_phi,
       			ev_.gl_pt,
       			ev_.ngl,
       			ev_.le_eta,
       			ev_.le_phi,
       			ev_.le_pt,
       			ev_.nle,
       			ev_.le_g[ev_.nle],
       			0.1,
       			1, ev_.gl_pid, 11);
       ev_.nle++;
    }

    if (!isMedium) continue;
    if (ev_.nme<MiniEvent_t::maxpart){

    	ev_.me_ch[ev_.nme]     = elecs->at(i).charge();
    	ev_.me_pt[ev_.nme]     = elecs->at(i).pt();
    	ev_.me_phi[ev_.nme]    = elecs->at(i).phi();
    	ev_.me_eta[ev_.nme]    = elecs->at(i).eta();
    	ev_.me_mass[ev_.nme]   = elecs->at(i).mass();
        ev_.me_bdt[ev_.nme]   = mvaValue;
        ev_.me_dz[ev_.nme]    = dz;
        ev_.me_dxy[ev_.nme]   = dxy;
    	if( isEB )
    		ev_.me_relIso[ev_.nme] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    	else
    		ev_.me_relIso[ev_.nme] = (elecs->at(i).userFloat("hgcElectronID:caloIsoRing1") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing2") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing3") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing4")) / elecs->at(i).energy();
    	ev_.me_g[ev_.nme] = -1;
    	ev_.me_relTkIso[ev_.nme] = elecs->at(i).dr03TkSumPt()/elecs->at(i).pt();
    	matchClosestDr(ev_.gl_eta,
    	           		ev_.gl_phi,
    	       			ev_.gl_pt,
    	       			ev_.ngl,
    	       			ev_.me_eta,
    	       			ev_.me_phi,
    	       			ev_.me_pt,
    	       			ev_.nme,
    	       			ev_.me_g[ev_.nme],
    	       			0.1,
    	       			1, ev_.gl_pid, 11);
    	ev_.nme++;
    }
    if (!isTight) continue;
    if (ev_.nte>=MiniEvent_t::maxpart) break;

    ev_.te_ch[ev_.nte]     = elecs->at(i).charge();
    ev_.te_pt[ev_.nte]     = elecs->at(i).pt();
    ev_.te_phi[ev_.nte]    = elecs->at(i).phi();
    ev_.te_eta[ev_.nte]    = elecs->at(i).eta();
    ev_.te_mass[ev_.nte]   = elecs->at(i).mass();
    ev_.te_bdt[ev_.nte]   = mvaValue;
    ev_.te_dz[ev_.nte]    = dz;
    ev_.te_dxy[ev_.nte]   = dxy;
    if( isEB )
      ev_.te_relIso[ev_.nte] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    else
      ev_.te_relIso[ev_.nte] = (elecs->at(i).userFloat("hgcElectronID:caloIsoRing1") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing2") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing3") + elecs->at(i).userFloat("hgcElectronID:caloIsoRing4")) / elecs->at(i).energy();
    ev_.te_g[ev_.nte] = -1;
    ev_.te_relTkIso[ev_.nte] = elecs->at(i).dr03TkSumPt()/elecs->at(i).pt();
    matchClosestDr(ev_.gl_eta,
    		ev_.gl_phi,
			ev_.gl_pt,
			ev_.ngl,
			ev_.te_eta,
			ev_.te_phi,
			ev_.te_pt,
			ev_.nte,
			ev_.te_g[ev_.nte],
			0.1,
			1, ev_.gl_pid, 11);
    ev_.nte++;
  }

  // Jets
  ev_.nj = 0;
  for (size_t i =0; i < jets->size(); i++) {
	if (ev_.nj>=MiniEvent_t::maxjets) break;
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

    pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    retLoose.set(false);
    bool isLoose = jetIDLoose_(jets->at(i), retLoose);
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    bool isTight = jetIDTight_(jets->at(i), retTight);

    double mvav2   = jets->at(i).bDiscriminator("pfCombinedMVAV2BJetTags"); 
    bool isLooseMVAv2  = mvav2 > mvaThres_[0];
    bool isMediumMVAv2 = mvav2 > mvaThres_[1];
    bool isTightMVAv2  = mvav2 > mvaThres_[2];
    double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
                            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    bool isLooseDeepCSV  = deepcsv > deepThres_[0];
    bool isMediumDeepCSV = deepcsv > deepThres_[1];
    bool isTightDeepCSV  = deepcsv > deepThres_[2];

    ev_.j_id[ev_.nj]      = (isTight | (isLoose<<1));
    ev_.j_pt[ev_.nj]      = jets->at(i).pt();
    ev_.j_phi[ev_.nj]     = jets->at(i).phi();
    ev_.j_eta[ev_.nj]     = jets->at(i).eta();
    ev_.j_mass[ev_.nj]    = jets->at(i).mass();
    ev_.j_jecf[ev_.nj]    = jets->at(i).pt()/jets->at(i).correctedJet("Uncorrected").pt();
    ev_.j_mvav2[ev_.nj]   = (isTightMVAv2 | (isMediumMVAv2<<1) | (isLooseMVAv2<<2)); 
    ev_.j_deepcsv[ev_.nj] = (isTightDeepCSV | (isMediumDeepCSV<<1) | (isLooseDeepCSV<<2));
    ev_.j_flav[ev_.nj]    = jets->at(i).partonFlavour();
    ev_.j_hadflav[ev_.nj] = jets->at(i).hadronFlavour();
    ev_.j_pid[ev_.nj]     = (jets->at(i).genParton() ? jets->at(i).genParton()->pdgId() : 0);
    ev_.j_g[ev_.nj] = -1;

    matchClosestDr(ev_.gj_eta,
    		ev_.gj_phi,
			ev_.gj_pt,
			ev_.ngj,
			ev_.j_eta,
			ev_.j_phi,
			ev_.j_pt,
			ev_.nj,
			ev_.j_g[ev_.nj],
			0.2,
			1);
    ev_.nj++;

  }
  
  // MET
  ev_.nmet = 0;
  if (mets->size() > 0) {
    ev_.met_pt[ev_.nmet]  = mets->at(0).pt();
    ev_.met_eta[ev_.nmet] = mets->at(0).eta();
    ev_.met_phi[ev_.nmet] = mets->at(0).phi();
    ev_.nmet++;
  }

  // Photons

  ev_.nlp = 0;
  ev_.ntp = 0;

  for (size_t i = 0; i < photons->size(); i++) {
    if (photons->at(i).pt() < 10.) continue;
    if (fabs(photons->at(i).eta()) > 3.) continue;

    float mvaValue = photons->at(i).userFloat("mvaValue");
    bool isEB = photons->at(i).isEB();

    float photonEnergy = -1.;
    // Change this if the default primary vertex choice is not satisfactory
    float photonEta = photons->at(i).eta();
    if ( isEB ) {
      // nMax15 algorithm, see https://indico.cern.ch/event/659834/contributions/2691448/attachments/1508099/2357330/170814_upsg_ledovskoy.pdf
      // barrel energy algorithm to compensate for high pileup
      // sort crystals in cluster, then add top 15 to find energy

      auto& pho = photons->at(i);

      // Build map of unique crystals to their total energy used in supercluster
      std::map<DetId, double> unique_cells;
      double nCellsEffective{0.};
      for( auto&& detid_frac : pho.superCluster()->hitsAndFractions() ) {
        auto hit = ecalRecHits->find(detid_frac.first);
        if ( hit == ecalRecHits->end() and detid_frac.first.subdetId() == DetId::Ecal ) {
          std::cout << "uh oh, missing a hit" << std::endl;
          continue;
        }
        else if ( hit == ecalRecHits->end() ) {
          continue;
        }
        double frac = detid_frac.second;
        unique_cells[detid_frac.first] += hit->energy() * frac;
        nCellsEffective += frac;
      }

      // Copy map to vector and sort
      std::vector<double> cells;
      for( auto&& det_e : unique_cells ) {
        cells.push_back(det_e.second);
      }
      std::sort(cells.begin(), cells.end());

      // Sum top N crystals
      double energy_nmax15{0.};
      size_t nCells{0u};
      for(auto it=cells.rbegin(); it!=cells.rend(); ++it) {
        if ( nCells < 15 ) energy_nmax15 += *it;
        if ( nCells == 15 ) break;
        nCells++;
      }

      photonEnergy = energy_nmax15;

    } else { // isEB

      photonEnergy = photons->at(i).superCluster()->seed()->energy();
      if ( photonEcorr_ != nullptr ) {
        float eCorr = photonEcorr_->GetBinContent(photonEcorr_->FindBin(std::abs(photons->at(i).eta()), photonEnergy));
        if ( eCorr != 0. ) photonEnergy *= eCorr;
      }
    }
     
    bool isLoose = 0;
    bool isTight = 0;

    if( isEB )
      {
	 isLoose = (mvaValue > 0.00);
	 isTight = (mvaValue > 0.56);
      }     
     else
      {
	 isLoose = (mvaValue > 0.20);
	 isTight = (mvaValue > 0.68);
      }          

    if (!isLoose) continue;
    if (ev_.nlp<MiniEvent_t::maxpart){

       ev_.lp_isEB[ev_.nlp]   = isEB ? 1 : 0;
       ev_.lp_pt[ev_.nlp]     = photonEnergy / std::cosh(photonEta);
       ev_.lp_phi[ev_.nlp]    = photons->at(i).phi();
       ev_.lp_eta[ev_.nlp]    = photonEta;
       ev_.lp_nrj[ev_.nlp]    = photonEnergy;
       ev_.lp_g[ev_.nlp] = -1;
       ev_.lp_bdt[ev_.nlp] = mvaValue;
       // add multicluster quantities too
       // for endcap it's seed multi, otherwise just use supercluster
       if (isEB) {
         ev_.lp_pt_multi[ev_.nlp]     = photons->at(i).superCluster()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
         ev_.lp_phi_multi[ev_.nlp]    = photons->at(i).superCluster()->phi();
         ev_.lp_eta_multi[ev_.nlp]    = photons->at(i).superCluster()->eta();
         ev_.lp_nrj_multi[ev_.nlp]    = photons->at(i).superCluster()->energy();
       }
       else {
         ev_.lp_pt_multi[ev_.nlp]     = photons->at(i).superCluster()->seed()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
         ev_.lp_phi_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->phi();
         ev_.lp_eta_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->eta();
         ev_.lp_nrj_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->energy();
       }
       matchClosestDr(ev_.gp_eta,
           		ev_.gp_phi,
       			ev_.gp_pt,
       			ev_.ngp,
       			ev_.lp_eta,
       			ev_.lp_phi,
       			ev_.lp_pt,
       			ev_.nlp,
       			ev_.lp_g[ev_.nlp],
       			0.1,
       			1);
       ev_.nlp++;
    }

    if (!isTight) continue;
    if (ev_.ntp>=MiniEvent_t::maxpart) break;

    ev_.tp_isEB[ev_.ntp]   = isEB ? 1 : 0;
    ev_.tp_pt[ev_.ntp]     = photonEnergy / std::cosh(photonEta);
    ev_.tp_phi[ev_.ntp]    = photons->at(i).phi();
    ev_.tp_eta[ev_.ntp]    = photonEta;
    ev_.tp_nrj[ev_.ntp]    = photonEnergy;
    ev_.tp_g[ev_.ntp] = -1;
    ev_.tp_bdt[ev_.ntp] = mvaValue;
    // add multicluster quantities too
    // for endcap it's seed multi, otherwise just use supercluster
    if (isEB) {
      ev_.tp_pt_multi[ev_.ntp]     = photons->at(i).superCluster()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
      ev_.tp_phi_multi[ev_.ntp]    = photons->at(i).superCluster()->phi();
      ev_.tp_eta_multi[ev_.ntp]    = photons->at(i).superCluster()->eta();
      ev_.tp_nrj_multi[ev_.ntp]    = photons->at(i).superCluster()->energy();
    }
    else {
      ev_.tp_pt_multi[ev_.ntp]     = photons->at(i).superCluster()->seed()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
      ev_.tp_phi_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->phi();
      ev_.tp_eta_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->eta();
      ev_.tp_nrj_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->energy();
    }

    matchClosestDr(ev_.gp_eta,
        		ev_.gp_phi,
    			ev_.gp_pt,
    			ev_.ngp,
    			ev_.tp_eta,
    			ev_.tp_phi,
    			ev_.tp_pt,
    			ev_.ntp,
    			ev_.tp_g[ev_.nlp],
    			0.1,
    			1);

    ev_.ntp++;
  }
   
}

// ------------ method called for each event  ------------
  void
MiniFromPat::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //analyze the event
  if(!iEvent.isRealData()) genAnalysis(iEvent, iSetup);
  recoAnalysis(iEvent, iSetup);
  
  //save event if at least one lepton at gen or reco level
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  t_event_->Fill();
  t_genParts_->Fill();
  t_genPhotons_->Fill();
  t_vertices_->Fill();
  t_genJets_->Fill();
  t_looseElecs_->Fill();
  t_mediumElecs_->Fill();
  t_tightElecs_->Fill();
  t_looseMuons_->Fill();
  t_tightMuons_->Fill();
  t_allTaus_->Fill();
  t_puppiJets_->Fill();
  t_puppiMET_->Fill();
  t_loosePhotons_->Fill();
  t_tightPhotons_->Fill();

}



bool 
MiniFromPat::isME0MuonSelNew(const reco::Muon& muon, double dEtaCut, double dPhiCut, double dPhiBendCut, edm::EventSetup const& iSetup)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;

    if(!ME0Geometry_){
    	edm::ESHandle<ME0Geometry> hGeom;
    	iSetup.get<MuonGeometryRecord>().get(hGeom);
    	ME0Geometry_ =( &*hGeom);
    	if(!ME0Geometry_)
    		return false;
    }

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      if (chamber->detector() == 5){

        for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

          LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
          LocalPoint seg_loc_coord(segment->x, segment->y, 0);
          LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
          LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

          const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);
          if(!me0chamber)continue;

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
MiniFromPat::beginJob()
{
}

// ------------ method called once each run ----------------
void
MiniFromPat::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
  void
MiniFromPat::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MiniFromPat::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniFromPat::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniFromPat);
