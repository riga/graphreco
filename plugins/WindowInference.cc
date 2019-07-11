/*
 * CMSSW plugin that performs a Window-based inference of networks using RecHits.
 *
 * Author: Marcel Rieger <marcel.rieger@cern.ch>
 */

#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "RecoHGCal/GraphReco/interface/Window.h"

// macros for simplified logs
// message logger disabled for the moment
// #define INFO edm::LogInfo("WindowInference")
// #define WARNING edm::LogWarning("WindowInference")
// #define ERROR edm::LogError("WindowInference")
#define INFO std::cout << "WindowInference INFO   : "
#define WARNING std::cout << "WindowInference WARNING: "
#define ERROR std::cout << "WindowInference ERROR  : "

// datastructure hold by edm::GlobalCache
struct WindowInferenceCache
{
    WindowInferenceCache(const edm::ParameterSet& config)
        : graphDef(nullptr)
    {
    }

    std::atomic<tensorflow::GraphDef*> graphDef;
};

class WindowInference : public edm::stream::EDAnalyzer<edm::GlobalCache<WindowInferenceCache> >
{
public:
    explicit WindowInference(const edm::ParameterSet&, const WindowInferenceCache*);
    ~WindowInference();

    // methods for handling the global cache
    static std::unique_ptr<WindowInferenceCache> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const WindowInferenceCache*);

private:
    void beginStream(edm::StreamID);
    void endStream();
    void analyze(const edm::Event&, const edm::EventSetup&);

    void createWindows();
    void fillWindows(const edm::Event&);
    void evaluateWindow(Window*);
    void fillRecHitFeatures(const HGCRecHit*, float*);

    // dummy function for the moment
    void reconstructShowers();

    // options
    std::vector<edm::InputTag> recHitCollections_;
    double minPhi_;
    double maxPhi_;
    double minEta_;
    double maxEta_;
    double deltaPhi_;
    double deltaEta_;
    double overlapPhi_;
    double overlapEta_;
    std::string inputTensorName_;
    std::string outputTensorName_;
    bool batchedModel_;
    size_t padSize_;

    // tokens
    std::vector<edm::EDGetTokenT<HGCRecHitCollection> > recHitTokens_;

    // rechit tools
    hgcal::RecHitTools recHitTools_;

    // windows
    std::vector<Window*> windows_;

    // the tensorflow session
    tensorflow::Session* session_;

    // hardcoded values
    size_t nFeatures_;
    double epsilon_;
};

std::unique_ptr<WindowInferenceCache> WindowInference::initializeGlobalCache(
    const edm::ParameterSet& config)
{
    // this method is supposed to create, initialize and return a WindowInferenceCache instance
    WindowInferenceCache* windowInferenceCache = new WindowInferenceCache(config);

    // load the graph def and save it
    std::string graphPath = config.getParameter<std::string>("graphPath");
    INFO << "loading graph from " << graphPath << std::endl;
    windowInferenceCache->graphDef = tensorflow::loadGraphDef(graphPath);

    // set some global configs, such as the TF log level
    tensorflow::setLogging("0");

    return std::unique_ptr<WindowInferenceCache>(windowInferenceCache);
}

void WindowInference::globalEndJob(const WindowInferenceCache* windowInferenceCache)
{
    // reset the graphDef
    if (windowInferenceCache->graphDef != nullptr)
    {
        delete windowInferenceCache->graphDef;
    }
}

WindowInference::WindowInference(const edm::ParameterSet& config,
    const WindowInferenceCache* windowInferenceCache)
    : recHitCollections_(config.getParameter<std::vector<edm::InputTag> >("recHitCollections"))
    , minPhi_(config.getParameter<double>("minPhi"))
    , maxPhi_(config.getParameter<double>("maxPhi"))
    , minEta_(config.getParameter<double>("minEta"))
    , maxEta_(config.getParameter<double>("maxEta"))
    , deltaPhi_(config.getParameter<double>("deltaPhi"))
    , deltaEta_(config.getParameter<double>("deltaEta"))
    , overlapPhi_(config.getParameter<double>("overlapPhi"))
    , overlapEta_(config.getParameter<double>("overlapEta"))
    , inputTensorName_(config.getParameter<std::string>("inputTensorName"))
    , outputTensorName_(config.getParameter<std::string>("outputTensorName"))
    , batchedModel_(config.getParameter<bool>("batchedModel"))
    , padSize_((size_t)config.getParameter<uint32_t>("padSize"))
    , session_(nullptr)
    , nFeatures_(10)
    , epsilon_(1e-5)
{
    // sanity checks for sliding windows
    if (deltaPhi_ <= 0 || deltaEta_ <= 0 || overlapPhi_ <= 0 || overlapEta_ <= 0)
    {
        throw cms::Exception("IncorrectWindowParameters")
            << "deltaPhi, deltaE, overlapPhi and overlapEta must be > 0";
    }
    else if (deltaPhi_ <= overlapPhi_)
    {
        throw cms::Exception("IncorrectWindowParameters")
            << "deltaPhi must be larger than overlapPhi";
    }
    else if (deltaEta_ <= overlapEta_)
    {
        throw cms::Exception("IncorrectWindowParameters")
            << "deltaEta must be larger than overlapEta";
    }

    // get tokens
    for (edm::InputTag& recHitCollection : recHitCollections_)
    {
        recHitTokens_.push_back(consumes<HGCRecHitCollection>(recHitCollection));
    }

    // mount the graphDef stored in windowInferenceCache onto the session
    session_ = tensorflow::createSession(windowInferenceCache->graphDef);
}

WindowInference::~WindowInference()
{
}

void WindowInference::beginStream(edm::StreamID streamId)
{
    createWindows();
}

void WindowInference::endStream()
{
    // close the session
    tensorflow::closeSession(session_);
    session_ = nullptr;

    // delete windows
    for (Window*& window : windows_)
    {
        delete window;
        window = nullptr;
    }
    windows_.clear();
}

void WindowInference::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    recHitTools_.getEventSetup(setup);

    // fill rechits into windows
    fillWindows(event);

    // run the evaluation per window
    for (Window* window : windows_)
    {
        evaluateWindow(window);
        // std::cout << window->outputTensor.shape().DebugString() << std::endl;
    }

    // reconstruct showers using all windows and put them into the event
    reconstructShowers();

    // clear all windows
    for (Window* window : windows_)
    {
        window->clear();
    }
}

void WindowInference::createWindows()
{
    for (float phi = minPhi_; phi + epsilon_ < maxPhi_; phi += deltaPhi_ - overlapPhi_)
    {
        for (float eta = minEta_; eta + epsilon_ < maxEta_; eta += deltaEta_ - overlapEta_)
        {
            windows_.push_back(new Window(phi, phi + deltaPhi_, eta, eta + deltaEta_,
                padSize_, nFeatures_, batchedModel_, inputTensorName_));
        }
    }

    INFO << "built " << windows_.size() << " window(s)" << std::endl;
}

void WindowInference::fillWindows(const edm::Event& event)
{
    // read rechits from all collections and store them in appropriate windows
    for (edm::EDGetTokenT<HGCRecHitCollection>& token : recHitTokens_)
    {
        edm::Handle<HGCRecHitCollection> handle;
        event.getByToken(token, handle);
        for (const HGCRecHit& recHit : *handle)
        {
            // TODO: right now, all windows are checked per rechit which might be stopped earlier
            // e.g. in case in window rejects a rechit due to a too small eta value, and the
            // subsequent windows have even higher eta ranges, but since overlap rules might be
            // somewhat complex, go for the brute force approch for now
            const GlobalPoint position = recHitTools_.getPosition(recHit.detid());
            float phi = position.phi();
            float eta = position.eta();
            for (Window* window : windows_)
            {
                window->maybeAddRecHit(recHit, phi, eta);
            }
        }
    }
}

void WindowInference::evaluateWindow(Window* window)
{
    // fill rechit features
    float* data = window->inputTensor.flat<float>().data();
    size_t nFilled = std::min<size_t>(window->getNRecHits(), padSize_);
    for (size_t i = 0; i < nFilled; i++)
    {
        const HGCRecHit* recHit = window->recHits.at(i);
        fillRecHitFeatures(recHit, data);
    }

    // zero-padding of unfilled rechits
    if (nFilled < padSize_)
    {
        for (size_t i = 0; i < (padSize_ - nFilled) * nFeatures_; i++)
        {
            *(data++) = 0.;
        }
    }

    // define the output and run
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session_, window->inputTensorList, { outputTensorName_ }, &outputs);

    // store the output in the window
    window->outputTensor = outputs[0];
}

void WindowInference::fillRecHitFeatures(const HGCRecHit* recHit, float* data)
{
    // fill rechit features: energy, eta, phi, theta, r, x, y, z, detId, time
    // all features _must_ be float types, or otherwise the float pointer arithmetic will break
    // most features are extracted from the GlobalPoint of the sensor which already uses float types
    // (https://github.com/cms-sw/cmssw/blob/master/DataFormats/GeometryVector/interface/GlobalPoint.h#L7)

    const GlobalPoint position = recHitTools_.getPosition(recHit->detid());

    *(data++) = recHit->energy();
    *(data++) = position.eta();
    *(data++) = position.phi();
    *(data++) = position.theta();
    *(data++) = position.mag();
    *(data++) = position.x();
    *(data++) = position.y();
    *(data++) = position.z();
    *(data++) = (float)recHit->detid();
    *(data++) = recHit->time();
}

void WindowInference::reconstructShowers()
{
    // this is where the stitching magic happens
}

DEFINE_FWK_MODULE(WindowInference);
