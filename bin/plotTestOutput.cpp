/*
 * plotOutput.C
 *
 *  Created on: 27 Sep 2019
 *      Author: jkiesele
 */
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <iostream>


/*
 *
 * The meaning of the hit features can be found in WindowBase.cpp
 * It is different for rechits, layer clusters, and tracks
 *
 */
#define id_dx 10
#define id_rechit ((float)0)
#define id_layercluster ((float)1)
#define id_track ((float)-1)
enum rechitfeatures{rh_energy=0, rh_eta=1, rh_dphi=2}; // etc
enum trackfeatures{tr_momentum=0, tr_eta=1, tr_dphi=2};

void mergeOverflow(TH1D* h){
    auto bc = h->GetBinContent(h->GetNbinsX()+1);
    h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX())+ bc);
}


int main(){

    TFile f("testout.root","READ");
    TTree * tree = (TTree *)f.Get("WindowNTupler/tree");
    if(!tree || tree->IsZombie())
        return -1;
    int nentries = tree->GetEntries();

    std::cout << "nentries " << nentries << std::endl;

    std::vector<std::vector<float> > * hitFeatures=0, * truthHitFractions=0;
    tree->SetBranchAddress("hitFeatures",&hitFeatures);
    tree->SetBranchAddress("truthHitFractions",&truthHitFractions);
    std::vector<int> * truthHitAssignementIdx=0;
    tree->SetBranchAddress("truthHitAssignementIdx",&truthHitAssignementIdx);

    TH1D noisefractions("noisefractions","noisefractions",10,0,1);
    TH2D noisefractions_perenergy("noisefractions_perenergy","noisefractions_perenergy",5,0,1,5,0,100);

    TH1D nhits("nhits","nhits",10,0,10000);

    for(int entry=0;entry<nentries; entry++){
        tree->GetEntry(entry);
        double noiseenergy=0;
        double totalenergy=0;
        for(size_t i=0;i<hitFeatures->size();i++){
            totalenergy+=hitFeatures->at(i).at(rh_energy);
            if(truthHitAssignementIdx->at(i) < 0){
                noiseenergy+=hitFeatures->at(i).at(rh_energy);
            }


        }
        noisefractions.Fill(noiseenergy/totalenergy);
        noisefractions_perenergy.Fill(noiseenergy/totalenergy, totalenergy);
        nhits.Fill(hitFeatures->size());

    }

    mergeOverflow(&noisefractions);
    mergeOverflow(&nhits);

    TCanvas cv;
    noisefractions.Draw();
    cv.Print("noisefractions.pdf");
    noisefractions_perenergy.Draw("colz");
    cv.Print("noisefractions_perenergy.pdf");
    nhits.Draw();
    cv.Print("nhits.pdf");
}


void plotOutput(){
    main();
}
