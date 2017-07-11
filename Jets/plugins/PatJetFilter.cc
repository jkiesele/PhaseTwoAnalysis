// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/PatJetFilter
// Class:      PatJetFilter
// 
/**\class PatJetFilter PatJetFilter.cc PhaseTwoAnalysis/PatJetFilter/plugins/PatJetFilter.cc

Description: adds vectors of pat ak4 PUPPI loose PF jets

Implementation:
- PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
- b-tagging WPs come from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#B_tagging 
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <vector>

//
// class declaration
//

class PatJetFilter : public edm::stream::EDProducer<> {
  public:
    explicit PatJetFilter(const edm::ParameterSet&);
    ~PatJetFilter();

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
    unsigned int pileup_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    double mvaThres_[3];
    double deepThres_[3];
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
PatJetFilter::PatJetFilter(const edm::ParameterSet& iConfig):
  pileup_(iConfig.getParameter<unsigned int>("pileup")),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE) 
{
  produces<std::vector<pat::Jet>>("Jets");
  produces<std::vector<pat::Jet>>("LooseMVAv2Jets");
  produces<std::vector<pat::Jet>>("MediumMVAv2Jets");
  produces<std::vector<pat::Jet>>("TightMVAv2Jets");
  produces<std::vector<pat::Jet>>("LooseDeepCSVJets");
  produces<std::vector<pat::Jet>>("MediumDeepCSVJets");
  produces<std::vector<pat::Jet>>("TightDeepCSVJets");

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
}


PatJetFilter::~PatJetFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
  void
PatJetFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  Handle<std::vector<pat::Electron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  std::unique_ptr<std::vector<pat::Jet>> filteredJets;
  std::unique_ptr<std::vector<pat::Jet>> filteredLooseMVAv2Jets;
  std::unique_ptr<std::vector<pat::Jet>> filteredMediumMVAv2Jets;
  std::unique_ptr<std::vector<pat::Jet>> filteredTightMVAv2Jets;
  std::unique_ptr<std::vector<pat::Jet>> filteredLooseDeepCSVJets;
  std::unique_ptr<std::vector<pat::Jet>> filteredMediumDeepCSVJets;
  std::unique_ptr<std::vector<pat::Jet>> filteredTightDeepCSVJets;
  std::vector<pat::Jet> Vec;
  std::vector<pat::Jet> looseMVAv2Vec;
  std::vector<pat::Jet> mediumMVAv2Vec;
  std::vector<pat::Jet> tightMVAv2Vec;
  std::vector<pat::Jet> looseDeepCSVVec;
  std::vector<pat::Jet> mediumDeepCSVVec;
  std::vector<pat::Jet> tightDeepCSVVec;

  for (size_t i = 0; i < jets->size(); i++) {
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

    if (!isLoose) continue;
    Vec.push_back(jets->at(i));

    double mvav2   = jets->at(i).bDiscriminator("pfCombinedMVAV2BJetTags"); 
    if (mvav2 > mvaThres_[0]) looseMVAv2Vec.push_back(jets->at(i));
    if (mvav2 > mvaThres_[1]) mediumMVAv2Vec.push_back(jets->at(i));
    if (mvav2 > mvaThres_[2]) tightMVAv2Vec.push_back(jets->at(i));

    double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
      jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    if (deepcsv > deepThres_[0]) looseDeepCSVVec.push_back(jets->at(i));
    if (deepcsv > deepThres_[1]) mediumDeepCSVVec.push_back(jets->at(i));
    if (deepcsv > deepThres_[2]) tightDeepCSVVec.push_back(jets->at(i));

  }

  filteredJets.reset(new std::vector<pat::Jet>(Vec));
  filteredLooseMVAv2Jets.reset(new std::vector<pat::Jet>(looseMVAv2Vec));
  filteredMediumMVAv2Jets.reset(new std::vector<pat::Jet>(mediumMVAv2Vec));
  filteredTightMVAv2Jets.reset(new std::vector<pat::Jet>(tightMVAv2Vec));
  filteredLooseDeepCSVJets.reset(new std::vector<pat::Jet>(looseDeepCSVVec));
  filteredMediumDeepCSVJets.reset(new std::vector<pat::Jet>(mediumDeepCSVVec));
  filteredTightDeepCSVJets.reset(new std::vector<pat::Jet>(tightDeepCSVVec));

  iEvent.put(std::move(filteredJets), "Jets");
  iEvent.put(std::move(filteredLooseMVAv2Jets), "LooseMVAv2Jets");
  iEvent.put(std::move(filteredMediumMVAv2Jets), "MediumMVAv2Jets");
  iEvent.put(std::move(filteredTightMVAv2Jets), "TightMVAv2Jets");
  iEvent.put(std::move(filteredLooseDeepCSVJets), "LooseDeepCSVJets");
  iEvent.put(std::move(filteredMediumDeepCSVJets), "MediumDeepCSVJets");
  iEvent.put(std::move(filteredTightDeepCSVJets), "TightDeepCSVJets");

  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
PatJetFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PatJetFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   PatJetFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   PatJetFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   PatJetFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   PatJetFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatJetFilter);
