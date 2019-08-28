/*
 * HGCalTrackPropagator.h
 *
 *  Created on: 28 Aug 2019
 *      Author: jkiesele
 */

#ifndef SRC_RECOHGCAL_GRAPHRECO_INTERFACE_HGCALTRACKPROPAGATOR_H_
#define SRC_RECOHGCAL_GRAPHRECO_INTERFACE_HGCALTRACKPROPAGATOR_H_

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/EventSetup.h"


struct TrackWithHGCalPos{
    const reco::Track * track;
    const GlobalPoint  pos;
};


class HGCalTrackPropagator{
public:
    enum zpos{ negZ=0, posZ=1};

    HGCalTrackPropagator(const edm::EventSetup &es); //sets up geometry etc.

    TrackWithHGCalPos propagateTrack(const reco::Track&)const;

private:
    edm::ESHandle<MagneticField> bField_;
    edm::ESHandle<Propagator> propagator_;
    std::unique_ptr<GeomDet> frontFaces_[2];
};



#endif /* SRC_RECOHGCAL_GRAPHRECO_INTERFACE_HGCALTRACKPROPAGATOR_H_ */
