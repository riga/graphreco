/*
 * Simple window structure with some convenience methods to determine if a RecHit matches the
 * windows eta and phi range, and to bookkeep input and output tensors for inference.
 *
 * Note 1: startPhi and startEta are inclusive, endPhi and endEta are exclusive (as in histograms)
 * Note 2: currently only positive eta ranges are supported for simplicity
 *
 * Author: Marcel Rieger <marcel.rieger@cern.ch>
 */

#ifndef RECOHGCAL_GRAPHRECO_INTERFACE_WINDOW_H_
#define RECOHGCAL_GRAPHRECO_INTERFACE_WINDOW_H_
#include <string>
#include "TTree.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include <vector>

#include "HGCalTrackPropagator.h"

struct HGCRecHitWithPos{
    const HGCRecHit * hit;
    const GlobalPoint  pos;
};


class Window {
public:

    //first prop all tracks then associate windows
    //geometry and propagator should be initialised before


    enum mode {
        useRechits, useLayerClusters
    };

    /*
     * startEta > endEta
     * innerRegionDEta, innerRegionDPhi > 0
     */
    Window(float centerEta, float centerPhi, float outerRegionDEta, float outerRegionDPhi,
            float innerRegionDEta, float innerRegionDPhi);

    ~Window() ;

    void setupTFInterface(size_t padSize, size_t nFeatures, bool batchedModel,
            const std::string& inputTensorName,
            const std::string& outputTensorName);

    void setMode(mode m){
        mode_=m;
    }
    mode getMode()const{
        return mode_;
    }

    inline bool accept(float phi, float eta) const {
        return fabs(reco::deltaPhi(phi, centerPhi_)) < outerRegionDPhi_
        && fabs(eta - centerEta_) < outerRegionDEta_;
    }

    bool maybeAddTrack(const TrackWithHGCalPos& t) {
        //potential cuts here!
        if (accept((float)t.pos.phi(), (float)t.pos.eta())
                && t.track->pt()>1) {
            tracks_.push_back(t);
            return true;
        }
        return false;
    }

    inline bool maybeAddRecHit(const HGCRecHit& recHit,
            const GlobalPoint& pos) {
        //potential cuts here!
        if (accept((float) pos.phi(), (float) pos.eta())
                && recHit.energy() > 0.01) {
            HGCRecHitWithPos trh = {
                    &recHit,
                    pos
            };
            recHits.push_back(trh);
            return true;
        }
        return false;
    }

    inline bool maybeAddLayerCluster(
            const reco::CaloCluster * layerCluster) {
        //potential cuts here!
        if (accept(layerCluster->phi(), layerCluster->eta())) {
            layerClusters_.push_back(layerCluster);
            return true;
        }
        return false;
    }

    inline bool maybeAddSimCluster(const SimCluster& sc){
        //potential cuts here!
        if (accept(sc.phi(),sc.eta())){
            simClusters_.push_back(&sc);
            return true;
        }
        return false;
    }

    inline bool isInner(const float& eta, const float& phi) {
        return fabs(reco::deltaPhi(phi, centerPhi_)) < innerRegionDPhi_
                && fabs(eta - centerEta_) < innerRegionDEta_;
    }

    inline size_t getNRecHits() const {
        return recHits.size();
    }


    void clear();

    void evaluate(tensorflow::Session* sess);

    //for output

    //gets size checks and adjustments before calling private fill
    void fillTTreeRechitFeatures(std::vector<std::vector<float > > * array) const;
    void fillTTreeLayerClusterFeatures(std::vector<std::vector<float > > * array) const;
    void fillTTreeTruthFractions(std::vector<std::vector<float > > * array) const;
    void fillTTreeTruthIDs(std::vector<std::vector<int > > * array) const; //one-hot
    void fillTTreeTruthEnergies(std::vector<std::vector<float > > * array) const; //"one-hot"

    static std::vector<Window> createWindows(size_t nSegmentsPhi,
            size_t nSegmentsEta, double minEta, double maxEta,
            double frameWidthEta, double frameWidthPhi);

    //debug functions

    void printDebug()const;

    const float& getCenterEta() const {
        return centerEta_;
    }

    const float& getCenterPhi() const {
        return centerPhi_;
    }

    const float& getOuterRegionDEta() const {
        return outerRegionDEta_;
    }

    const float& getOuterRegionDPhi() const {
        return outerRegionDPhi_;
    }

    const float& getInnerRegionDEta() const {
        return innerRegionDEta_;
    }

    const float& getInnerRegionDPhi() const {
        return innerRegionDPhi_;
    }


private:
    //for one rechit
    void fillRecHitFeatures(float*& data, const HGCRecHitWithPos * ) const;
    //for one layer cluster
    void fillLayerClusterFeatures(float*& data, const reco::CaloCluster * ) const;


    mode mode_;

    float centerEta_;
    float centerPhi_;

    float outerRegionDEta_;
    float outerRegionDPhi_;

    float innerRegionDEta_;
    float innerRegionDPhi_;

    std::vector<TrackWithHGCalPos > tracks_;
    std::vector<HGCRecHitWithPos> recHits;
    std::vector<const reco::CaloCluster * > layerClusters_;
    std::vector<const SimCluster*> simClusters_;

    tensorflow::Tensor inputTensor;
    tensorflow::NamedTensorList inputTensorList;
    tensorflow::Tensor outputTensor;
    std::string outputTensorName_;

    static const size_t nRechitFeatures_;
    static const size_t nLayerClusterFeatures_;
};

#endif  // RECOHGCAL_GRAPHRECO_INTERFACE_WINDOW_H_
