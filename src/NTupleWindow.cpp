/*
 * NTupleWindow.cpp
 *
 *  Created on: 26 Sep 2019
 *      Author: jkiesele
 */


#include "../interface/NTupleWindow.h"
#include "DataFormats/Math/interface/deltaR.h"

std::vector<std::vector<float> > * NTupleWindow::sp_hitFeatures_=0;

std::vector<std::vector<float> > * NTupleWindow::sp_truthHitFractions_=0;
std::vector<int>                 * NTupleWindow::sp_truthHitAssignementIdx_=0;
std::vector<float>               * NTupleWindow::sp_truthHitAssignedEnergies_=0;
std::vector<float>               * NTupleWindow::sp_truthHitAssignedEtas_=0;
std::vector<float>               * NTupleWindow::sp_truthHitAssignedPhis_=0;
std::vector<std::vector<int> >   * NTupleWindow::sp_truthHitAssignedPIDs_=0;

std::vector<int>    * NTupleWindow::sp_truthSimclusterIdx_=0;
std::vector<std::vector<int> >   * NTupleWindow::sp_truthSimclusterPIDs_=0;
std::vector<float>  * NTupleWindow::sp_truthSimclusterEnergies_=0;
std::vector<float>  * NTupleWindow::sp_truthSimclusterEtas_=0;
std::vector<float>  * NTupleWindow::sp_truthSimclusterPhis_=0;

//static
std::vector<NTupleWindow> NTupleWindow::createWindows(size_t nSegmentsPhi,
            size_t nSegmentsEta, double minEta, double maxEta, double frameWidthEta,
            double frameWidthPhi){
    return WindowBase::createWindows<NTupleWindow>(nSegmentsPhi,nSegmentsEta,
            minEta,maxEta,frameWidthEta,frameWidthPhi);
}



void NTupleWindow::createTreeBranches(TTree* t){

   // NTupleWindow dummy;
   // dummy.assignTreePointers(); //so that the pointers are not null, maybe not needed? FIXME

    t->Branch("hitFeatures", &sp_hitFeatures_);

    t->Branch("truthHitFractions", &sp_truthHitFractions_);
    t->Branch("truthHitAssignementIdx", &sp_truthHitAssignementIdx_);
    t->Branch("truthHitAssignedEnergies", &sp_truthHitAssignedEnergies_);
    t->Branch("truthHitAssignedEtas", &sp_truthHitAssignedEtas_);
    t->Branch("truthHitAssignedPhis", &sp_truthHitAssignedPhis_);
    t->Branch("truthHitAssignedPIDs", &sp_truthHitAssignedPIDs_);

    t->Branch("truthSimclusterIdx",&sp_truthSimclusterIdx_);
    t->Branch("truthSimclusterPIDs",&sp_truthSimclusterPIDs_);
    t->Branch("truthSimclusterEnergies",&sp_truthSimclusterEnergies_);
    t->Branch("truthSimclusterEtas",&sp_truthSimclusterEtas_);
    t->Branch("truthSimclusterPhis",&sp_truthSimclusterPhis_);

}

NTupleWindow::NTupleWindow(float centerEta, float centerPhi,
        float outerRegionDEta, float outerRegionDPhi, float innerRegionDEta,
        float innerRegionDPhi) :
        WindowBase(centerEta, centerPhi, outerRegionDEta, outerRegionDPhi,
                innerRegionDEta, innerRegionDPhi),
                nSimclusters_(0),
                truthTotalEnergy_(0){
}


//static


void NTupleWindow::assignTreePointers()  {

    sp_hitFeatures_ = &hitFeatures_;

    sp_truthHitFractions_ = &truthHitFractions_;
    sp_truthHitAssignementIdx_ = &truthHitAssignementIdx_;
    sp_truthHitAssignedEnergies_ = &truthHitAssignedEnergies_;
    sp_truthHitAssignedEtas_ = &truthHitAssignedEtas_;
    sp_truthHitAssignedPhis_ = &truthHitAssignedPhis_;
    sp_truthHitAssignedPIDs_ = &truthHitAssignedPIDs_;

    sp_truthSimclusterIdx_ = &truthSimclusterIdx_;
    sp_truthSimclusterPIDs_ = &truthSimclusterPIDs_;
    sp_truthSimclusterEnergies_ = &truthSimclusterEnergies_;
    sp_truthSimclusterEtas_ = &truthSimclusterEtas_;
    sp_truthSimclusterPhis_ = &truthSimclusterPhis_;

}



void NTupleWindow::clear(){
    WindowBase::clear(); //clears rechits etc

    detIDHitAsso_.clear();

    hitFeatures_.clear();

    truthHitFractions_.clear();
    truthHitAssignementIdx_.clear();
    truthHitAssignedEnergies_.clear();
    truthHitAssignedEtas_.clear();
    truthHitAssignedPhis_.clear();
    truthHitAssignedPIDs_.clear();

    truthSimclusterIdx_.clear();
    truthSimclusterPIDs_.clear();
    truthSimclusterEnergies_.clear();
    truthSimclusterEtas_.clear();
    truthSimclusterPhis_.clear();
    //etc

    nSimclusters_=0;
    truthTotalEnergy_=0;
}

void NTupleWindow::fillFeatureArrays(){

    hitFeatures_.clear();
    if(getMode() == useRechits){
        for(const auto& rh:recHits){
            std::vector<float> feats(nRechitFeatures_);
            auto data = &feats.at(0);
            fillRecHitFeatures(data,rh);
            hitFeatures_.push_back(feats);
        }
    }
    else{
        for(const auto& lc: layerClusters_){
            std::vector<float> feats(nLayerClusterFeatures_);
            auto data = &feats.at(0);
            fillLayerClusterFeatures(data,lc);
            hitFeatures_.push_back(feats);
        }
    }
    //add tracks LAST!
    for(const auto& tr:tracks_){
        std::vector<float> feats(nTrackFeatures_);
        auto data = &feats.at(0);
        fillTrackFeatures(data,tr);
        hitFeatures_.push_back(feats);
    }

}

void NTupleWindow::createDetIDHitAssociation(){
    detIDHitAsso_.clear();

    for(size_t i=0;i<recHits.size();i++){
        detIDHitAsso_[recHits.at(i)->hit->detid()]={i,1.};
    }
    if(getMode() != useLayerClusters)
        return;

    //now map to layer clusters
    for(size_t i=0;i<layerClusters_.size();i++){
        auto hafs = layerClusters_.at(i)->hitsAndFractions();
        for(const auto& haf: hafs){
            auto pos = detIDHitAsso_.find(haf.first);
            if(pos == detIDHitAsso_.end()) //edges
                continue;
            pos->second = {i, haf.second};
        }
    }

}

void NTupleWindow::fillTruthArrays(){

    createDetIDHitAssociation();
    calculateSimclusterFeatures();//fills the simcluster properties
    calculateTruthFractions();//generates the truthHitFractions_ vector
    fillTruthAssignment();//fills the rest

}
void NTupleWindow::calculateSimclusterFeatures(){

    truthSimclusterIdx_.clear();
    truthSimclusterPIDs_.clear();
    truthSimclusterEnergies_.clear();
    truthSimclusterEtas_.clear();
    truthSimclusterPhis_.clear();

    nSimclusters_=simClusters_.size();
    truthTotalEnergy_=0;;

    for(size_t i=0;i<simClusters_.size();i++){
        truthSimclusterIdx_.push_back(i);
        truthSimclusterPIDs_.push_back(pdgToOneHot(simClusters_.at(i)->pdgId()));
        auto simCMomentum = simClusters_.at(i)->p4();
        float energy = simCMomentum.E();
        truthTotalEnergy_ += energy;
        truthSimclusterEnergies_.push_back(energy);
        truthSimclusterEtas_.push_back(simCMomentum.Eta());
        truthSimclusterPhis_.push_back(simCMomentum.Phi());
    }
}

//needs function to match truth tracks

void NTupleWindow::calculateTruthFractions(){

    truthHitFractions_.clear();
    truthHitFractions_.resize(hitFeatures_.size(),
            std::vector<float>(simClusters_.size(), 0)); //includes tracks

    for (size_t i_sc = 0; i_sc < simClusters_.size(); i_sc++) {
        const auto& hitsandfracs = simClusters_.at(i_sc)->hits_and_fractions();
        for(const auto& haf: hitsandfracs){
            auto pos = detIDHitAsso_.find(haf.first);
            if(pos == detIDHitAsso_.end()) //edges
                continue;
            size_t idx = pos->first;
            float totalfrac = pos->second.second * haf.second;
            truthHitFractions_.at(idx).at(i_sc) += totalfrac; //can be more than 1-1 for layer clusters
        }
    }

    //associate the tracks here, such that they look like hits, simple matching
    size_t trackStartIterator = recHits.size();
    if(getMode() == useLayerClusters)
        trackStartIterator = layerClusters_.size();

    //match, will be improved by direct truth matching in new simclusters on longer term
    //assumption: for every track there is charged simcluster
    std::vector<size_t> usedSimclusters;

    for(size_t i_t=0;i_t<tracks_.size();i_t++){

        const double momentumscaler = 0.1;
        double minDistance=0.1 + 0.1;

        size_t matchedSCIdx=simClusters_.size();
        for(size_t i_sc=0;i_sc<simClusters_.size();i_sc++){
            if(fabs(simClusters_.at(i_sc)->charge())<0.1)
                continue;
            if(std::find(usedSimclusters.begin(),usedSimclusters.end(),i_sc) != usedSimclusters.end())
                continue;
            double scEnergy = simClusters_.at(i_sc)->p4().E();
            double trackMomentum = tracks_.at(i_t)->track->p();
            double distance = reco::deltaR(simClusters_.at(i_sc)->eta(),
                    simClusters_.at(i_sc)->phi(),
                    (float)tracks_.at(i_t)->pos.eta(),
                    (float)tracks_.at(i_t)->pos.phi()) +
                            momentumscaler*std::abs(scEnergy - trackMomentum)/(scEnergy);

            if(distance<minDistance)
                matchedSCIdx=i_sc;
        }
        truthHitFractions_.at(i_t+trackStartIterator).at(matchedSCIdx) = 1.;
        usedSimclusters.push_back(matchedSCIdx);
    }
}


void NTupleWindow::fillTruthAssignment(){

    truthHitAssignementIdx_.resize(truthHitFractions_.size());
    truthHitAssignedEnergies_.resize(truthHitFractions_.size());
    truthHitAssignedEtas_.resize(truthHitFractions_.size());
    truthHitAssignedPhis_.resize(truthHitFractions_.size());
    truthHitAssignedPIDs_.resize(truthHitFractions_.size());

    for (size_t i_hit = 0; i_hit < truthHitFractions_.size(); i_hit++) {
        size_t maxfrac_idx = std::max_element(
                truthHitFractions_.at(i_hit).begin(),
                truthHitFractions_.at(i_hit).end())
                - truthHitFractions_.at(i_hit).begin();

        truthHitAssignementIdx_.at(i_hit) = maxfrac_idx;
        truthHitAssignedEnergies_.at(i_hit) = truthSimclusterEnergies_.at(
                maxfrac_idx);
        truthHitAssignedEtas_.at(i_hit) = truthSimclusterEtas_.at(maxfrac_idx);
        truthHitAssignedPhis_.at(i_hit) = truthSimclusterPhis_.at(maxfrac_idx);
        truthHitAssignedPIDs_.at(i_hit) = truthSimclusterPIDs_.at(maxfrac_idx);

    }
}

