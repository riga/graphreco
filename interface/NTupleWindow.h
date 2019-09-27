/*
 * NTupleWindow.h
 *
 *  Created on: 26 Sep 2019
 *      Author: jkiesele
 */

#ifndef SRC_RECOHGCAL_GRAPHRECO_INTERFACE_NTUPLEWINDOW_H_
#define SRC_RECOHGCAL_GRAPHRECO_INTERFACE_NTUPLEWINDOW_H_

#include "../interface/WindowBase.h"
#include <algorithm>
#include <unordered_map>

class NTupleWindow: public WindowBase {
public:

    NTupleWindow(float centerEta, float centerPhi, float outerRegionDEta, float outerRegionDPhi,
            float innerRegionDEta, float innerRegionDPhi);

    static void createTreeBranches(TTree* t);

    //0 associate all rechits etc

    void fillFeatureArrays();
    //1
    void fillTruthArrays();
    //2
    void assignTreePointers() ;

    //3: tree->Fill();

    //4
    void clear();

    static std::vector<NTupleWindow> createWindows(size_t nSegmentsPhi,
            size_t nSegmentsEta, double minEta, double maxEta, double frameWidthEta,
            double frameWidthPhi);



private:
    NTupleWindow(){}

    void createDetIDHitAssociation();
    void calculateSimclusterFeatures();
    void calculateTruthFractions();
    void fillTruthAssignment();

    //temporary for (layer cluster) fraction calculation. detID to hit index and fraction
    std::unordered_map<DetId, std::pair<size_t, float>> detIDHitAsso_;


    //can be layer clusters or rechits according to mode
    std::vector<std::vector<float> >  hitFeatures_; //this includes tracks!

    std::vector<std::vector<float> >  truthHitFractions_;


    //hit (rechit or layer cluster) assigned to exactly one simcluster with index
    std::vector<int>                 truthHitAssignementIdx_;
    //truth energy of assigned simcluster
    std::vector<float>               truthHitAssignedEnergies_;
    //truth eta of assigned simcluster
    std::vector<float>               truthHitAssignedEtas_;
    //truth delta-phi to window center of assigned simcluster
    std::vector<float>               truthHitAssignedPhis_;
    //truth energy of assigned simcluster
    std::vector<std::vector<int> >   truthHitAssignedPIDs_;


    int nSimclusters_;
    float truthTotalEnergy_;
    std::vector<int>     truthSimclusterIdx_;
    std::vector<std::vector<int> >    truthSimclusterPIDs_;
    std::vector<float>   truthSimclusterEnergies_;
    std::vector<float>   truthSimclusterEtas_;
    std::vector<float>   truthSimclusterPhis_;


    //static pointers to create branches and fill tree
    static std::vector<std::vector<float> > * sp_hitFeatures_;

    static std::vector<std::vector<float> > * sp_truthHitFractions_;
    static std::vector<int>                 * sp_truthHitAssignementIdx_;
    static std::vector<float>               * sp_truthHitAssignedEnergies_;
    static std::vector<float>               * sp_truthHitAssignedEtas_;
    static std::vector<float>               * sp_truthHitAssignedPhis_;
    static std::vector<std::vector<int> >   * sp_truthHitAssignedPIDs_;

    static std::vector<int>    * sp_truthSimclusterIdx_;
    static std::vector<std::vector<int> >   * sp_truthSimclusterPIDs_;
    static std::vector<float>  * sp_truthSimclusterEnergies_;
    static std::vector<float>  * sp_truthSimclusterEtas_;
    static std::vector<float>  * sp_truthSimclusterPhis_;



};




#endif /* SRC_RECOHGCAL_GRAPHRECO_INTERFACE_NTUPLEWINDOW_H_ */
