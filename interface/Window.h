/*
 * Simple window structure with some convenience methods to determine if a RecHit matches the
 * windows eta and phi range, and to bookkeep input and output tensors for inference.
 *
 * Note 1: startPhi and startEta are inclusive, endPhi and endEta are exclusive (as in histograms)
 * Note 2: currently only positive eta ranges are supported for simplicity
 *
 * Author: Marcel Rieger <marcel.rieger@cern.ch>
 */

#ifndef RECOHGCAL_GRAPHRECO_WINDOW_H
#define RECOHGCAL_GRAPHRECO_WINDOW_H

struct Window
{
    Window(double startPhi, double endPhi, double startEta, double endEta)
        : startPhi(startPhi)
        , endPhi(endPhi)
        , startEta(startEta)
        , endEta(endEta)
    {
    }

    Window(double startPhi, double endPhi, double startEta, double endEta, size_t padSize,
        size_t nFeatures, bool batchedModel, std::string& inputTensorName)
        : Window(startPhi, endPhi, startEta, endEta)
    {
        setupTensors(padSize, nFeatures, batchedModel, inputTensorName);
    }

    void setupTensors(size_t padSize, size_t nFeatures, bool batchedModel,
        std::string& inputTensorName)
    {
        tensorflow::TensorShape shape = { (int)padSize, (int)nFeatures };
        if (batchedModel)
        {
            shape.InsertDim(0, 1);
        }

        inputTensor = tensorflow::Tensor(tensorflow::DT_FLOAT, shape);
        inputTensorList = { { inputTensorName, inputTensor } };
    }

    inline void addRecHit(const HGCRecHit& recHit)
    {
        recHits.push_back(&recHit);
    }

    inline bool acceptsRecHit(const HGCRecHit& recHit, float phi, float eta) const
    {
        return phi >= startPhi && phi < endPhi && eta >= startEta && eta < endEta;
    }

    inline bool maybeAddRecHit(const HGCRecHit& recHit, float phi, float eta)
    {
        if (acceptsRecHit(recHit, phi, eta))
        {
            addRecHit(recHit);
            return true;
        }
        return false;
    }

    inline size_t getNRecHits() const
    {
        return recHits.size();
    }

    inline void clear()
    {
        recHits.clear();
    }

    double startPhi;
    double endPhi;
    double startEta;
    double endEta;
    std::vector<const HGCRecHit*> recHits;
    tensorflow::Tensor inputTensor;
    tensorflow::NamedTensorList inputTensorList;
    tensorflow::Tensor outputTensor;
};

#endif // RECOHGCAL_GRAPHRECO_WINDOW_H
