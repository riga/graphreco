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

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

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
    void beginJob();
    void analyze(const edm::Event&, const edm::EventSetup&);
    void endJob();

    // per module copy, only the session is stored
    tensorflow::Session* session_;

    // options
    std::string inputTensorName_;
    std::string outputTensorName_;
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
    : session_(nullptr)
    , inputTensorName_(config.getParameter<std::string>("inputTensorName"))
    , outputTensorName_(config.getParameter<std::string>("outputTensorName"))
{
    // mount the graphDef stored in windowInferenceCache onto the session
    session_ = tensorflow::createSession(windowInferenceCache->graphDef);
}

WindowInference::~WindowInference()
{
}

void WindowInference::beginJob()
{
}

void WindowInference::endJob()
{
    // close the session
    tensorflow::closeSession(session_);
    session_ = nullptr;
}

void WindowInference::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, 100, 10 });

    // test: fill all values with 0.01
    float* d = input.flat<float>().data();
    for (size_t i = 0; i < 100; i++, d++)
    {
        for (size_t j = 0; j < 10; j++, d++)
        {
            *d = 0.01;
        }
    }

    // define the output and run
    std::vector<tensorflow::Tensor> outputs;
    INFO << "evaluate" << std::endl;
    tensorflow::run(session_, { { inputTensorName_, input } }, { outputTensorName_ }, &outputs);

    // check and print the output
    INFO << outputs[0].DebugString() << std::endl;
}

DEFINE_FWK_MODULE(WindowInference);
