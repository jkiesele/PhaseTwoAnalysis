// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/PatElectronFilter
// Class:      PatElectronFilter
// 
/**\class PatElectronFilter PatElectronFilter.cc PhaseTwoAnalysis/PatElectronFilter/plugins/PatElectronFilter.cc

Description: adds a vector of pat electrons

Implementation:
- electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
   /!\ no ID is implemented for forward electrons as:
   - PFClusterProducer does not run on miniAOD
   - jurassic isolation needs tracks
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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include <vector>

//
// class declaration
//

class PatElectronFilter : public edm::stream::EDProducer<> {
    public:
        explicit PatElectronFilter(const edm::ParameterSet&);
        ~PatElectronFilter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;

        bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
        bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
        bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 

        //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
        edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
        edm::EDGetTokenT<reco::BeamSpot> bsToken_;
        edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
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
PatElectronFilter::PatElectronFilter(const edm::ParameterSet& iConfig):
    elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions")))  
{
    produces<std::vector<pat::Electron>>("LooseElectrons");
    produces<std::vector<double>>("LooseElectronRelIso");
    produces<std::vector<pat::Electron>>("MediumElectrons");
    produces<std::vector<double>>("MediumElectronRelIso");
    produces<std::vector<pat::Electron>>("TightElectrons");
    produces<std::vector<double>>("TightElectronRelIso");

}


PatElectronFilter::~PatElectronFilter()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
PatElectronFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    Handle<std::vector<pat::Electron>> elecs;
    iEvent.getByToken(elecsToken_, elecs);
    Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();  
    std::unique_ptr<std::vector<pat::Electron>> filteredLooseElectrons;
    std::unique_ptr<std::vector<double>> filteredLooseElectronRelIso;
    std::unique_ptr<std::vector<pat::Electron>> filteredMediumElectrons;
    std::unique_ptr<std::vector<double>> filteredMediumElectronRelIso;
    std::unique_ptr<std::vector<pat::Electron>> filteredTightElectrons;
    std::unique_ptr<std::vector<double>> filteredTightElectronRelIso;
    std::vector<pat::Electron> looseVec;
    std::vector<double> looseIsoVec;
    std::vector<pat::Electron> mediumVec;
    std::vector<double> mediumIsoVec;
    std::vector<pat::Electron> tightVec;
    std::vector<double> tightIsoVec;
    for (size_t i = 0; i < elecs->size(); i++) {
        if (elecs->at(i).pt() < 10.) continue;
        if (fabs(elecs->at(i).eta()) > 3.) continue;

        bool isLoose = isLooseElec(elecs->at(i),conversions,beamspot);    
        bool isMedium = isMediumElec(elecs->at(i),conversions,beamspot);    
        bool isTight = isTightElec(elecs->at(i),conversions,beamspot);    

        double relIso = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();

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

    filteredLooseElectrons.reset(new std::vector<pat::Electron>(looseVec));
    filteredLooseElectronRelIso.reset(new std::vector<double>(looseIsoVec));
    filteredMediumElectrons.reset(new std::vector<pat::Electron>(mediumVec));
    filteredMediumElectronRelIso.reset(new std::vector<double>(mediumIsoVec));
    filteredTightElectrons.reset(new std::vector<pat::Electron>(tightVec));
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
PatElectronFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PatElectronFilter::endStream() {
}

// ------------ method check that an e passes loose ID ----------------------------------
    bool
PatElectronFilter::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
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
PatElectronFilter::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
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
PatElectronFilter::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
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

// ------------ method called when starting to processes a run  ------------
/*
   void
   PatElectronFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   PatElectronFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   PatElectronFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   PatElectronFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatElectronFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronFilter);
