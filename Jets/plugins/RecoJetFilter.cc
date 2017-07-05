// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/RecoJetFilter
// Class:      RecoJetFilter
// 
/**\class RecoJetFilter RecoJetFilter.cc PhaseTwoAnalysis/RecoJetFilter/plugins/RecoJetFilter.cc

Description: adds a vector of reco ak4 PUPPI jets

Implementation:
- no jet ID is applied
- b-tagging is not available 
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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <vector>
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

class RecoJetFilter : public edm::stream::EDProducer<> {
  public:
    explicit RecoJetFilter(const edm::ParameterSet&);
    ~RecoJetFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<reco::GsfElectron>> elecsToken_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken_;

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
RecoJetFilter::RecoJetFilter(const edm::ParameterSet& iConfig):
  elecsToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets")))
{
  produces<std::vector<reco::PFJet>>("Jets");

}


RecoJetFilter::~RecoJetFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
  void
RecoJetFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  Handle<std::vector<reco::GsfElectron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<std::vector<reco::PFJet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  std::unique_ptr<std::vector<reco::PFJet>> filteredJets;
  std::vector<reco::PFJet> Vec;
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
    Vec.push_back(jets->at(i));

  }

  filteredJets.reset(new std::vector<reco::PFJet>(Vec));

  iEvent.put(std::move(filteredJets), "Jets");

  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
RecoJetFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
RecoJetFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   RecoJetFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   RecoJetFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   RecoJetFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   RecoJetFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoJetFilter);
