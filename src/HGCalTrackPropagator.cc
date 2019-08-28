/*
 * HGCalTrackPropagator.cc
 *
 *  Created on: 28 Aug 2019
 *      Author: jkiesele
 */

#include "../interface/HGCalTrackPropagator.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

HGCalTrackPropagator::HGCalTrackPropagator(const edm::EventSetup &es){


    //get the propagator
    es.get<IdealMagneticFieldRecord>().get(bField_);
    es.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagator_);


    //create the hgcal inner surface for both z
    edm::ESHandle<HGCalDDDConstants> hdc;
    es.get<IdealGeometryRecord>().get("HGCalEESensitive", hdc);

    double frontZ = hdc.product()->waferZ(1, true);
    auto frontradii = hdc.product()->rangeR(frontZ, true);

    frontFaces_[posZ] = std::make_unique < GeomDet
            > (Disk::build(Disk::PositionType(0, 0, frontZ),
                    Disk::RotationType(),
                    SimpleDiskBounds(frontradii.first, frontradii.second,
                            frontZ - 0.5, frontZ + 0.5)).get());

    frontFaces_[negZ] = std::make_unique < GeomDet
            > (Disk::build(Disk::PositionType(0, 0, -frontZ),
                    Disk::RotationType(),
                    SimpleDiskBounds(frontradii.first, frontradii.second,
                            -frontZ - 0.5, -frontZ + 0.5)).get());


}

TrackWithHGCalPos HGCalTrackPropagator::propagateTrack(const reco::Track& t)const{

    zpos trackz = posZ;
    if(t.eta()<0) trackz = negZ;

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(t, bField_.product());

    TrajectoryStateOnSurface tsos = (*propagator_).propagate(fts, frontFaces_[trackz]->surface());
    if (tsos.isValid())
        return TrackWithHGCalPos{&t, tsos.globalPosition()};
    return TrackWithHGCalPos{&t, GlobalPoint(1,1,0)};

}

