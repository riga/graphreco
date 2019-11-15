/*
 * NTupleWindow.cpp
 *
 *  Created on: 26 Sep 2019
 *      Author: jkiesele
 */



#include "../interface/NTupleWindow.h"
#include "DataFormats/Math/interface/deltaR.h"

std::vector<std::vector<float>>* NTupleWindow::sp_hitFeatures_=0;
std::vector<float>* NTupleWindow::sp_recHitEnergy_;
std::vector<float>* NTupleWindow::sp_recHitEta_;
std::vector<float>* NTupleWindow::sp_recHitRelPhi_;
std::vector<float>* NTupleWindow::sp_recHitTheta_;
std::vector<float>* NTupleWindow::sp_recHitMag_;
std::vector<float>* NTupleWindow::sp_recHitX_;
std::vector<float>* NTupleWindow::sp_recHitY_;
std::vector<float>* NTupleWindow::sp_recHitZ_;
std::vector<float>* NTupleWindow::sp_recHitDetID_;
std::vector<float>* NTupleWindow::sp_recHitTime_;
std::vector<float>* NTupleWindow::sp_recHitID_;
std::vector<float>* NTupleWindow::sp_recHitPad_;
//std::vector<float>* NTupleWindow::sp_trackFeatures_=0;

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


float * NTupleWindow::sp_windowEta_=0;
float * NTupleWindow::sp_windowPhi_=0;

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

    if (flattenRechitFeatures_) {
        t->Branch("recHitEnergy", &sp_recHitEnergy_);
        t->Branch("recHitEta", &sp_recHitEta_);
        t->Branch("recHitRelPhi", &sp_recHitRelPhi_);
        t->Branch("recHitTheta", &sp_recHitTheta_);
        t->Branch("recHitMag", &sp_recHitMag_);
        t->Branch("recHitX", &sp_recHitX_);
        t->Branch("recHitY", &sp_recHitY_);
        t->Branch("recHitZ", &sp_recHitZ_);
        t->Branch("recHitDetID", &sp_recHitDetID_);
        t->Branch("recHitTime", &sp_recHitTime_);
        t->Branch("recHitID", &sp_recHitID_);
        t->Branch("recHitPad", &sp_recHitPad_);
    }
    else
        t->Branch("recHitFeatures", &sp_hitFeatures_);
    //t->Branch("trackFeatures", &sp_trackFeatures_);

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

    t->Branch("windowEta",sp_windowEta_);
    t->Branch("windowPhi",sp_windowPhi_);

}

NTupleWindow::NTupleWindow(float centerEta, float centerPhi,
        float outerRegionDEta, float outerRegionDPhi, float innerRegionDEta,
        float innerRegionDPhi) :
        WindowBase(centerEta, centerPhi, outerRegionDEta, outerRegionDPhi,
                innerRegionDEta, innerRegionDPhi) {

    windowEta_ = centerEta;
    windowPhi_ = centerPhi;
}


//static

// This is not a very clear/efficient way to do this, but it's the simplest way given
// the structure we already have
void NTupleWindow::flattenRechitFeatures() {
    for (size_t i = 0; i < hitFeatures_.size(); i++) {
        recHitEnergy_.push_back(hitFeatures_[i][kEnergy]);
        recHitEta_.push_back(hitFeatures_[i][kEta]);
        recHitRelPhi_.push_back(hitFeatures_[i][kRelPhi]);
        recHitTheta_.push_back(hitFeatures_[i][kTheta]);
        recHitMag_.push_back(hitFeatures_[i][kMag]);
        recHitX_.push_back(hitFeatures_[i][kx]);
        recHitY_.push_back(hitFeatures_[i][ky]);
        recHitZ_.push_back(hitFeatures_[i][kz]);
        recHitDetID_.push_back(hitFeatures_[i][kDetid]);
        recHitTime_.push_back(hitFeatures_[i][kTime]);
        recHitID_.push_back(hitFeatures_[i][kId]);
        recHitPad_.push_back(hitFeatures_[i][kPad]);
    }
}


void NTupleWindow::assignTreePointers()  {

    sp_hitFeatures_ = &hitFeatures_;
    //sp_trackFeatures_ = &hitFeatures_;
    if (flattenRechitFeatures_) {
        sp_recHitEnergy_  = &recHitEnergy_;
        sp_recHitEta_  = &recHitEta_;
        sp_recHitRelPhi_  = &recHitRelPhi_;
        sp_recHitTheta_  = &recHitTheta_;
        sp_recHitMag_  = &recHitMag_;
        sp_recHitX_  = &recHitX_;
        sp_recHitY_  = &recHitY_;
        sp_recHitZ_  = &recHitZ_;
        sp_recHitDetID_  = &recHitDetID_;
        sp_recHitTime_  = &recHitTime_;
        sp_recHitID_  = &recHitID_;
        sp_recHitPad_  = &recHitPad_;
    }

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

    sp_windowEta_ = &windowEta_;
    sp_windowPhi_ = &windowPhi_;

}



void NTupleWindow::clear(){
    WindowBase::clear(); //clears rechits etc

    detIDHitAsso_.clear();

    hitFeatures_.clear();
    recHitEnergy_.clear();
    recHitEta_.clear();
    recHitRelPhi_.clear();
    recHitTheta_.clear();
    recHitMag_.clear();
    recHitX_.clear();
    recHitY_.clear();
    recHitZ_.clear();
    recHitDetID_.clear();
    recHitTime_.clear();
    recHitID_.clear();
    recHitPad_.clear();

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



}

void NTupleWindow::fillFeatureArrays(){
    //NO CUTS HERE!

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
    //for(const auto& tr:tracks_){
    //    std::vector<float> feats(nTrackFeatures_);
    //    auto data = &feats.at(0);
    //    fillTrackFeatures(data,tr);
    //    hitFeatures_.push_back(feats);
    //}

}

void NTupleWindow::fillTruthArrays(){

    createDetIDHitAssociation();
    calculateSimclusterFeatures();
    calculateTruthFractions();
    fillTruthAssignment();

    DEBUGPRINT(getCenterEta());
    DEBUGPRINT(getCenterPhi());
    DEBUGPRINT(truthHitFractions_.size());
    DEBUGPRINT(hitFeatures_.size());
    DEBUGPRINT(truthSimclusterIdx_.size());
}

void NTupleWindow::createDetIDHitAssociation(){
    detIDHitAsso_.clear();

    if(getMode() == useRechits){
        for(size_t i=0;i<recHits.size();i++){
            detIDHitAsso_[recHits.at(i)->hit->detid()]={i,1.};
        }
    }
    else{
        for(size_t i=0;i<layerClusters_.size();i++){
            for(const auto& haf: layerClusters_.at(i)->hitsAndFractions()){
                detIDHitAsso_[haf.first] = {i, haf.second};
            }
        }
    }
}

void NTupleWindow::calculateSimclusterFeatures(){

    truthSimclusterIdx_.clear();
    truthSimclusterPIDs_.clear();
    truthSimclusterEnergies_.clear();
    truthSimclusterEtas_.clear();
    truthSimclusterPhis_.clear();


    for(size_t i=0;i<simClusters_.size();i++){
        truthSimclusterIdx_.push_back(i);
        truthSimclusterPIDs_.push_back(pdgToOneHot(simClusters_.at(i)->pdgId()));
        const auto& simCMomentum = simClusters_.at(i)->p4();
        truthSimclusterEnergies_.push_back(simCMomentum.E());
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
            if(pos == detIDHitAsso_.end()) //edges or not included in layer clusters
                continue;
            size_t idx = pos->second.first;
            float totalfrac = pos->second.second * haf.second;
            truthHitFractions_.at(idx).at(i_sc) += totalfrac; //can be more than 1-1 for layer clusters
        }
    }

    //associate the tracks here, such that they look like hits, simple matching
    //size_t trackStartIterator = recHits.size();
    //if(getMode() == useLayerClusters)
    //    trackStartIterator = layerClusters_.size();

    ////match, will be improved by direct truth matching in new simclusters on longer term
    ////assumption: for every track there is charged simcluster
    ///*
    // *
    // * this is just a temporary solution until a proper simcluster-simtrack integration exists
    // *
    // */
    //std::vector<size_t> usedSimclusters;

    //float debug_ntrackwithnoSC=0;

    //for(size_t i_t=0;i_t<tracks_.size();i_t++){

    //    const double momentumscaler = 0.0001;
    //    double minDistance=0.1 + 0.1;

    //    size_t matchedSCIdx=simClusters_.size();
    //    double distance = 0;
    //    for(size_t i_sc=0;i_sc<simClusters_.size();i_sc++){
    //        if(fabs(simClusters_.at(i_sc)->charge())<0.1)
    //            continue;
    //        if(std::find(usedSimclusters.begin(),usedSimclusters.end(),i_sc) != usedSimclusters.end())
    //            continue;
    //        double scEnergy = simClusters_.at(i_sc)->p4().E();
    //        double trackMomentum = tracks_.at(i_t)->track->p();
    //        distance = reco::deltaR(simClusters_.at(i_sc)->eta(),
    //                simClusters_.at(i_sc)->phi(),
    //                (float)tracks_.at(i_t)->pos.eta(),
    //                (float)tracks_.at(i_t)->pos.phi()) +
    //                        momentumscaler*std::abs(scEnergy - trackMomentum)/(scEnergy);

    //        if(distance<minDistance){
    //            matchedSCIdx=i_sc;
    //            minDistance=distance;
    //        }
    //    }
    //    if(matchedSCIdx<simClusters_.size()){
    //        truthHitFractions_.at(i_t+trackStartIterator).at(matchedSCIdx) = 1.;
    //        usedSimclusters.push_back(matchedSCIdx);
    //    }
    //    else{
    //        debug_ntrackwithnoSC++;
    //        DEBUGPRINT(distance);
    //        DEBUGPRINT(tracks_.at(i_t)->track->p());
    //        DEBUGPRINT(tracks_.at(i_t)->track->eta());
    //    }
    //}
    //DEBUGPRINT(debug_ntrackwithnoSC);
    //DEBUGPRINT(debug_ntrackwithnoSC/(float)tracks_.size());
}


void NTupleWindow::fillTruthAssignment(){

    truthHitAssignementIdx_.resize(truthHitFractions_.size());
    truthHitAssignedEnergies_.resize(truthHitFractions_.size());
    truthHitAssignedEtas_.resize(truthHitFractions_.size());
    truthHitAssignedPhis_.resize(truthHitFractions_.size());
    truthHitAssignedPIDs_.resize(truthHitFractions_.size());

    bool nosim = simClusters_.size() < 1;

    for (size_t i_hit = 0; i_hit < truthHitFractions_.size(); i_hit++) {

        bool allzero = std::all_of(truthHitFractions_.at(i_hit).begin(),
                truthHitFractions_.at(i_hit).end(), [](float i) {return i==0;});

        if(allzero || nosim){
            truthHitAssignementIdx_.at(i_hit) = -1;
            continue;
        }
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

