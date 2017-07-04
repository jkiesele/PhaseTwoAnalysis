// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/RecoElectronFilter
// Class:      RecoElectronFilter
// 
/**\class RecoElectronFilter RecoElectronFilter.cc PhaseTwoAnalysis/RecoElectronFilter/plugins/RecoElectronFilter.cc

Description: adds a vector of reco electrons

Implementation:
- lepton isolation needs to be refined
- electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
*/
//
// Original Author:  Elvire Bouvier
//         Created:  Sun, 02 Jul 2017 20:54:15 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include "Math/GenVector/VectorUtil.h"
#include "TMVA/Reader.h" 
#include "TMVA/MethodBDT.h" 

//
// class declaration
//

class RecoElectronFilter : public edm::stream::EDProducer<> {
  public:
    explicit RecoElectronFilter(const edm::ParameterSet&);
    ~RecoElectronFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};      

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    bool isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    int matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles);
    void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
    float evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize);

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    std::unique_ptr<HGCalIDTool> hgcEmId_; 
    TMVA::Reader tmvaReader_;
    float hgcId_startPosition, hgcId_lengthCompatibility, hgcId_sigmaietaieta, hgcId_deltaEtaStartPosition, hgcId_deltaPhiStartPosition, hOverE_hgcalSafe, hgcId_cosTrackShowerAngle, trackIsoR04jurassic_D_pt, ooEmooP, d0, dz, pt, etaSC, phiSC, nPV, expectedMissingInnerHits, passConversionVeto, isTrue;

    edm::EDGetTokenT<std::vector<reco::GsfElectron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<edm::ValueMap<double>> trackIsoValueMapToken_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsNoLepToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;    

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
RecoElectronFilter::RecoElectronFilter(const edm::ParameterSet& iConfig):
  elecsToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
  trackIsoValueMapToken_(consumes<edm::ValueMap<double>>(iConfig.getParameter<edm::InputTag>("trackIsoValueMap"))),
  pfCandsNoLepToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandsNoLep"))),
  genPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices")))  
{
  produces<std::vector<reco::GsfElectron>>("LooseElectrons");
  produces<std::vector<double>>("LooseElectronRelIso");
  produces<std::vector<reco::GsfElectron>>("MediumElectrons");
  produces<std::vector<double>>("MediumElectronRelIso");
  produces<std::vector<reco::GsfElectron>>("TightElectrons");
  produces<std::vector<double>>("TightElectronRelIso");

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

}


RecoElectronFilter::~RecoElectronFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
  void
RecoElectronFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  hgcEmId_->getEventSetup(iSetup);
  hgcEmId_->getEvent(iEvent);  

  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);
  // Vertices
  int prVtx = -1;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
  }

  Handle<std::vector<reco::GsfElectron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<ValueMap<double>> trackIsoValueMap;
  iEvent.getByToken(trackIsoValueMapToken_, trackIsoValueMap);
  Handle<std::vector<reco::PFCandidate>> pfCandsNoLep;
  iEvent.getByToken(pfCandsNoLepToken_, pfCandsNoLep);  
  Handle<std::vector<reco::GenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);
  std::unique_ptr<std::vector<reco::GsfElectron>> filteredLooseElectrons;
  std::unique_ptr<std::vector<double>> filteredLooseElectronRelIso;
  std::unique_ptr<std::vector<reco::GsfElectron>> filteredMediumElectrons;
  std::unique_ptr<std::vector<double>> filteredMediumElectronRelIso;
  std::unique_ptr<std::vector<reco::GsfElectron>> filteredTightElectrons;
  std::unique_ptr<std::vector<double>> filteredTightElectronRelIso;
  std::vector<reco::GsfElectron> looseVec;
  std::vector<double> looseIsoVec;
  std::vector<reco::GsfElectron> mediumVec;
  std::vector<double> mediumIsoVec;
  std::vector<reco::GsfElectron> tightVec;
  std::vector<double> tightIsoVec;

  for(size_t i = 0; i < elecs->size(); i++) { 
    if (elecs->at(i).pt() < 10.) continue;
    if (fabs(elecs->at(i).eta()) > 3.) continue;

    double relIso = 0.;
    for (size_t k = 0; k < pfCandsNoLep->size(); k++) {
      if (ROOT::Math::VectorUtil::DeltaR(elecs->at(i).p4(),pfCandsNoLep->at(k).p4()) > 0.4) continue;
      relIso += pfCandsNoLep->at(k).pt();
    }
    if (elecs->at(i).pt() > 0.) relIso = relIso / elecs->at(i).pt(); 
    else relIso = -1.;

    Ptr<const reco::GsfElectron> el4iso(elecs,i);
    double eljurassicIso = (*trackIsoValueMap)[el4iso];
    double elpt = elecs->at(i).pt();
    double elMVAVal = -1.;
    if (prVtx > -0.5 && hgcEmId_->setElectronPtr(&(elecs->at(i)))) 
      elMVAVal = (double)evalMVAElec(elecs->at(i),vertices->at(prVtx),conversions,beamspot,genParts,eljurassicIso/elpt,vertices->size());
    bool isLoose  = isLooseElec(elecs->at(i),conversions,beamspot,elMVAVal);    
    bool isMedium = isMediumElec(elecs->at(i),conversions,beamspot,elMVAVal);    
    bool isTight  = isTightElec(elecs->at(i),conversions,beamspot,elMVAVal);    

    if (!isLoose) continue;
    looseVec.push_back(elecs->at(i));
    looseIsoVec.push_back(relIso);

    if (!isMedium) continue;
    mediumVec.push_back(elecs->at(i));
    mediumIsoVec.push_back(relIso);

    if (!isTight) continue;
    tightVec.push_back(elecs->at(i));
    tightIsoVec.push_back(relIso);

  }

  filteredLooseElectrons.reset(new std::vector<reco::GsfElectron>(looseVec));
  filteredLooseElectronRelIso.reset(new std::vector<double>(looseIsoVec));
  filteredMediumElectrons.reset(new std::vector<reco::GsfElectron>(mediumVec));
  filteredMediumElectronRelIso.reset(new std::vector<double>(mediumIsoVec));
  filteredTightElectrons.reset(new std::vector<reco::GsfElectron>(tightVec));
  filteredTightElectronRelIso.reset(new std::vector<double>(tightIsoVec));

  iEvent.put(std::move(filteredLooseElectrons), "LooseElectrons");
  iEvent.put(std::move(filteredLooseElectronRelIso), "LooseElectronRelIso");
  iEvent.put(std::move(filteredMediumElectrons), "MediumElectrons");
  iEvent.put(std::move(filteredMediumElectronRelIso), "MediumElectronRelIso");
  iEvent.put(std::move(filteredTightElectrons), "TightElectrons");
  iEvent.put(std::move(filteredTightElectronRelIso), "TightElectronRelIso");

  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
RecoElectronFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
RecoElectronFilter::endStream() {
}


// ------------ loose elec ID -----------
bool 
RecoElectronFilter::isLooseElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
RecoElectronFilter::isMediumElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
RecoElectronFilter::isTightElec(const reco::GsfElectron & recoEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
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
RecoElectronFilter::matchToTruth(const reco::GsfElectron & recoEl, const edm::Handle<std::vector<reco::GenParticle>> & genParticles) {
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
RecoElectronFilter::findFirstNonElectronMother(const reco::Candidate *particle,
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
RecoElectronFilter::evalMVAElec(const reco::GsfElectron & recoEl, const reco::Vertex & recoVtx, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, const edm::Handle<std::vector<reco::GenParticle>> & genParticles, double isoEl, int vertexSize) {

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

// ------------ method called when starting to processes a run  ------------
/*
   void
   RecoElectronFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   RecoElectronFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   RecoElectronFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   RecoElectronFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoElectronFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoElectronFilter);
