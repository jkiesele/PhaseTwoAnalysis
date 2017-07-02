// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/NTupler
// Class:      WeightCounter
// 
/**\class WeightCounter WeightCounter.cc PhaseTwoAnalysis/NTupler/plugins/WeightCounter.cc

Description: save the total number of events in then gen dataset

Implementation:
  the sum of the input weights is stored in a TH1F
*/
//
// Original Author:  Mirena Ivova Paneva
//         Created:  Mon, 16 Nov 2015 15:15:29 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"


//
// class declaration
//

class WeightCounter : public edm::EDAnalyzer {
  public:
    explicit WeightCounter(const edm::ParameterSet&);
    ~WeightCounter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    TH1F* weight;
    float T_Event_weight;

    //tokens
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;

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
WeightCounter::WeightCounter(const edm::ParameterSet& iConfig) : 
  genEventInfoToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
  //now do what ever initialization is needed

}


WeightCounter::~WeightCounter()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
WeightCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;



  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(genEventInfoToken_, genEvtInfo);

  T_Event_weight= genEvtInfo->weight();


  weight->Fill(0.,T_Event_weight);

}


// ------------ method called once each job just before starting event loop  ------------
  void 
WeightCounter::beginJob()
{
  edm::Service<TFileService> fs;  
  weight = fs->make<TH1F>("Event_weight",    ";Variation;Events", 1000, 0., 1000.); 


}

// ------------ method called once each job just after ending the event loop  ------------
  void 
WeightCounter::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WeightCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WeightCounter);
