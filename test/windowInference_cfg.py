# coding: utf-8

"""
Test config to run the WindowInference plugin.
"""


import os
import subprocess

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


# helper to determine the location if _this_ file
def get_this_dir():
    if "__file__" in globals():
        return os.path.dirname(os.path.abspath(__file__))
    else:
        return os.path.expandvars("$CMSSW_BASE/src/RecoHGCal/GraphReco/test")


# ensure that the graph exists
# if not, call the create_dummy_graph.py script in a subprocess since tensorflow complains
# when its loaded twice (once here in python, once in c++)
graph_path = os.path.abspath("graph.pb")
if not os.path.exists(graph_path):
    script_path = os.path.join(get_this_dir(), "create_dummy_graph.py")
    code = subprocess.call(["python", script_path, graph_path])
    if code != 0:
        raise Exception("create_dummy_graph.py failed")

# setup minimal options
options = VarParsing("python")
options.setDefault("inputFiles", "file:///eos/cms/store/cmst3/group/hgcal/CMG_studies/mrieger/hgcalsim/RecoTask/closeby_5.0To100.0_idsmix_dR0.4_n10_rnd1_s1/dev3_countRecHits/reco_0_n10.root")  # noqa: E501
options.parseArguments()

# define the process to run
process = cms.Process("HGR")

# minimal configuration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))

# process options
process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)

# load and configure the windowInference module
from RecoHGCal.GraphReco.windowInference_cfi import windowInference

process.windowInference = windowInference.clone(
    graphPath=cms.string(graph_path),
    inputTensorName=cms.string("input"),
    outputTensorName=cms.string("output"),
)

# define the path to run
process.p = cms.Path(process.windowInference)
