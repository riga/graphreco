/*
 * Window.cpp
 *
 *  Created on: 26 Aug 2019
 *      Author: jkiesele
 */



#include "FWCore/Utilities/interface/Exception.h"
#include "../interface/Window.h"

#define DEBUGPRINT(x) {std::cout << #x << ": " << x << std::endl;}


Window::Window(float centerEta, float centerPhi, float outerRegionDEta, float outerRegionDPhi,
        float innerRegionDEta, float innerRegionDPhi) :

        mode_(useLayerClusters),
        centerEta_(centerEta),centerPhi_(centerPhi),
        outerRegionDEta_(outerRegionDEta),outerRegionDPhi_(outerRegionDPhi),
        innerRegionDEta_(innerRegionDEta), innerRegionDPhi_(
                innerRegionDPhi) {

    //sanity checks FIXME: add more
    if (innerRegionDEta_ <= 0 || innerRegionDPhi_ <= 0) {
        throw cms::Exception("IncorrectWindowParameters")
                << "innerRegionDEta,innerRegionDPhi  must be > 0";
    }
    if (innerRegionDEta_ > outerRegionDEta_ || innerRegionDPhi_ > outerRegionDPhi_) {
        throw cms::Exception("IncorrectWindowParameters")
                << "innerRegionDEta,innerRegionDPhi  must be <= outerRegionDEta, outerRegionDPhi";
    }

}


Window::~Window() {
    clear();
}

void Window::setupTFInterface(size_t padSize, size_t nFeatures, bool batchedModel,
        const std::string& inputTensorName,
        const std::string& outputTensorName) {
    /*   tensorflow::TensorShape shape = { (int) padSize, (int) nFeatures };
    if (batchedModel) {
        shape.InsertDim(0, 1);
    }

    inputTensor = tensorflow::Tensor(tensorflow::DT_FLOAT, shape);
    inputTensorList = { {inputTensorName, inputTensor}};
    */
}



void Window::clear() {
    // this class does not own anything
    tracks_.clear();
    recHits.clear();
    layerClusters_.clear();
    simClusters_.clear();
}

void Window::evaluate(tensorflow::Session* sess) {
 /*   std::vector<tensorflow::Tensor> outputs;
            tensorflow::run(session_, inputTensorList, { outputTensorName_ },
                    &outputs);
            //do something with it
*/
}//FIXME

void Window::fillTTreeTrackFeatures(std::vector<std::vector<float > > * array) const{
    array->clear();
    for(const auto& t:tracks_){
        std::vector<float > features(nTrackFeatures_);
        auto data = &features.at(0);
        fillTrackFeatures(data, &t);
        array->push_back(features);
    }
}

void Window::fillTTreeRechitFeatures(std::vector<std::vector<float > > * array) const{
    array->clear();
    for(const auto& rh : recHits){
        std::vector<float > features(nRechitFeatures_);
        auto data = &features.at(0);
        fillRecHitFeatures(data, &rh);
        array->push_back(features);
    }
}
void Window::fillTTreeLayerClusterFeatures(std::vector<std::vector<float > > * array) const{
    array->clear();
    std::vector<float > features;

    //FILL

    array->push_back(features);
}
void Window::fillTTreeTruthFractions(std::vector<std::vector<float > > * array) const{
    array->clear();
    //for ...
    std::vector<float > fractions;

    //FILL

    array->push_back(fractions);
}
void Window::fillTTreeTruthIDs(std::vector<std::vector<int > > * array) const{
    //one-hot
    array->clear();
    //for ...
    std::vector<int > ids;

    //FILL

    array->push_back(ids);
}
void Window::fillTTreeTruthEnergies(std::vector<std::vector<float > > * array) const{
    //"one-hot"
    array->clear();
    //for ...
    std::vector<float > energies;

    //FILL

    array->push_back(energies);
}

//// private ////

const size_t Window::nTrackFeatures_=5;
void Window::fillTrackFeatures(float*& data, const TrackWithHGCalPos * ) const {
    //uses tracjs
    //creates inputTensorList
    // float* data
    // same for filling the output vector for root tuples
}//FIXME

const size_t Window::nRechitFeatures_=10;
void Window::fillRecHitFeatures(float*& data, const HGCRecHitWithPos * recHit) const {
    *(data++) = recHit->hit->energy();
    *(data++) = recHit->pos.eta();
    *(data++) = recHit->pos.phi();
    *(data++) = recHit->pos.theta();
    *(data++) = recHit->pos.mag();
    *(data++) = recHit->pos.x();
    *(data++) = recHit->pos.y();
    *(data++) = recHit->pos.z();
    *(data++) = (float)recHit->hit->detid();
    *(data++) = recHit->hit->time();
}

const size_t Window::nLayerClusterFeatures_=15;
void Window::fillLayerClusterFeatures(float*& data, const reco::CaloCluster * ) const {
    //uses recHits
    //creates inputTensorList
    // float* data
    // same for filling the output vector for root tuples
}//FIXME



// debug

void Window::printDebug()const{
     DEBUGPRINT(centerPhi_);
     DEBUGPRINT(centerEta_);
     DEBUGPRINT(outerRegionDEta_);
     DEBUGPRINT(outerRegionDPhi_);
     DEBUGPRINT(innerRegionDEta_);
     DEBUGPRINT(innerRegionDPhi_);
     std::cout << "coverage phi " << centerPhi_-outerRegionDPhi_ << "| " << centerPhi_-innerRegionDPhi_ << "[  :" <<
             centerPhi_ << ":  ]" << centerPhi_+innerRegionDPhi_ << " |" << centerPhi_ + outerRegionDPhi_ << std::endl;
     std::cout << "coverage eta " << centerEta_-outerRegionDEta_ << "| " << centerEta_-innerRegionDEta_ << "[  :" <<
             centerEta_ << ":  ]" << centerEta_+innerRegionDEta_ << " |" << centerEta_+outerRegionDEta_ << std::endl;
}



/// static

std::vector<Window> Window::createWindows(size_t nSegmentsPhi,
        size_t nSegmentsEta, double minEta, double maxEta, double frameWidthEta,
        double frameWidthPhi) {

    if (minEta <= 0 || maxEta <= 0 || minEta >= maxEta) {
        throw cms::Exception("IncorrectWindowCreationParameters")
        << "minEta, maxEta must be > 0 and maxEta > minEta (negative eta will be created automatically)";
    }
    if (frameWidthEta < 0 || frameWidthPhi < 0 ) {
        throw cms::Exception("IncorrectWindowCreationParameters")
        << "frameWidthEta, frameWidthPhi must be >= 0";
    }

    const float epsilon=1e-5;

    std::vector<Window> windows;
    float phiStep = 2. * M_PI / (float) nSegmentsPhi;
    float totalDPhi = (phiStep + 2. * frameWidthPhi)/2.;
    float etaStep = (maxEta - minEta) / (float) nSegmentsEta;
    float totalDEta = (etaStep + 2. * frameWidthEta)/2.;

    for (float phi_i = -M_PI; phi_i + epsilon < M_PI; phi_i += phiStep) {
        float phiCenter = phi_i + phiStep / 2.;
        for (float eta_j = minEta; eta_j + epsilon < maxEta; eta_j += etaStep) {
            float etaCenter = eta_j + etaStep / 2.;
            auto w = Window(etaCenter, phiCenter, totalDEta, totalDPhi,
                    etaStep / 2., phiStep / 2.);
            w.printDebug();
            windows.push_back(w);

        }



        float minEtaNeg = -maxEta;
        float maxEtaNeg = -minEta;
        for (float eta_j = minEtaNeg; eta_j + epsilon < maxEtaNeg; eta_j +=
                etaStep) {

            float etaCenter = eta_j + etaStep / 2.;

            auto w = Window(etaCenter, phiCenter, totalDEta, totalDPhi,
                    etaStep / 2., phiStep / 2.);
            w.printDebug();
            windows.push_back(w);
        }

    }
    return windows;
}















