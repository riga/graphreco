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


    for(int entry=0;entry<nentries; entry++){
        tree->GetEntry(entry);
        std::cout << hitFeatures->size() <<std::endl;

    }



}


void plotOutput(){
    main();
}
