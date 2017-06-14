#ifndef EgammaIsolationProducers_EgammaElectronTkIsolationProducer_h
#define EgammaIsolationProducers_EgammaElectronTkIsolationProducer_h

//*****************************************************************************
// File:      EgammaElectronTkIsolationProducer.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

class EgammaElectronTkIsolationProducer : public edm::stream::EDProducer<> {
 public:
  explicit EgammaElectronTkIsolationProducer(const edm::ParameterSet&);
  ~EgammaElectronTkIsolationProducer();
  
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::EDGetTokenT< reco::GsfElectronCollection> electronProducer_;
  edm::EDGetTokenT<reco::TrackCollection> trackProducer_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotProducer_;

  double ptMin_;
  double intRadiusBarrel_;
  double intRadiusEndcap_;
  double stripBarrel_;
  double stripEndcap_;
  double extRadius_;
  double maxVtxDist_;
  double drb_;
  
  edm::ParameterSet conf_;

};


#endif
