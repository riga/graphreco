// -*- C++ -*-
//
// Package:    RecoHGCal/WindowNTupler
// Class:      WindowNTupler
//
/**\class WindowNTupler WindowNTupler.cc RecoHGCal/GraphReco/plugins/WindowNTupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jan Kieseler
//         Created:  Tue, 27 Aug 2019 17:25:47 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "../interface/Window.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class WindowNTupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit WindowNTupler(const edm::ParameterSet&);
      ~WindowNTupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;
      edm::EDGetTokenT<reco::CaloClusterCollection> layerClusters_;
      edm::EDGetTokenT<std::vector<SimCluster>> simClusters_;
      std::vector<edm::EDGetTokenT<HGCRecHitCollection> > rechitsTokens_;

      std::vector<Window> windows_;

      edm::Service<TFileService> fs_;
      TTree * outTree_;

    std::vector<std::vector<float> > * rechitFeatures_;
    std::vector<std::vector<float> > * layerClusterFeatures_;
    std::vector<std::vector<float> > * truthFractions_;
    std::vector<std::vector<int> > * truthIDs_;
    std::vector<std::vector<float> > * truthEnergies_;

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
WindowNTupler::WindowNTupler(const edm::ParameterSet& config)
 :
  tracksToken_(consumes<TrackCollection>(config.getUntrackedParameter<edm::InputTag>("tracks"))),
  layerClusters_(consumes<reco::CaloClusterCollection>(config.getUntrackedParameter<edm::InputTag>("layerClusters"))),
  simClusters_(consumes<std::vector<SimCluster>>(config.getUntrackedParameter<edm::InputTag>("simClusters"))),
  outTree_(nullptr),
  rechitFeatures_(new std::vector<std::vector<float> >()),
  layerClusterFeatures_(new std::vector<std::vector<float> >()),
  truthFractions_(new std::vector<std::vector<float> >()),
  truthIDs_(new std::vector<std::vector<int> >()),
  truthEnergies_(new std::vector<std::vector<float> >())

/* ... */

{
    for (edm::InputTag& recHitCollection : config.getParameter<
            std::vector<edm::InputTag> >("recHitCollections")) {
        rechitsTokens_.push_back(
                consumes<HGCRecHitCollection>(recHitCollection));
    }

    //DEBUG INFO: has checks built in
    windows_ = Window::createWindows(
            (size_t) config.getParameter<uint32_t>("nPhiSegments"),
            (size_t) config.getParameter<uint32_t>("nEtaSegments"),
            config.getParameter<double>("minEta"),
            config.getParameter<double>("maxEta"),
            config.getParameter<double>("etaFrameWidth"),
            config.getParameter<double>("phiFrameWidth"));

}


WindowNTupler::~WindowNTupler()
{
    delete rechitFeatures_;
    delete layerClusterFeatures_;
    delete truthFractions_;
    delete truthIDs_;
    delete truthEnergies_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
WindowNTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::vector<TrackWithHGCalPos> trackswithpos;
   HGCalTrackPropagator trackprop(iSetup);
   for(const auto& track : iEvent.get(tracksToken_) ) {
       auto propTrack = trackprop.propagateTrack(track);
       trackswithpos.push_back(propTrack);
   }

/// and fill the windows
}


// ------------ method called once each job just before starting event loop  ------------
void WindowNTupler::beginJob() {

    if (!fs_) {
        throw edm::Exception(edm::errors::Configuration,
                "TFile Service is not registered in cfg file");
    }

    outTree_ = fs_->make<TTree>("tree", "tree");
    outTree_->Branch("rechitFeatures",&rechitFeatures_);
    outTree_->Branch("layerClusterFeatures",&layerClusterFeatures_);
    outTree_->Branch("truthFractions",&truthFractions_);
    outTree_->Branch("truthIDs",&truthIDs_);
    outTree_->Branch("truthEnergies",&truthEnergies_);


}

// ------------ method called once each job just after ending the event loop  ------------
void
WindowNTupler::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WindowNTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WindowNTupler);
