# coding: utf-8

"""
Initialization file for the WindowInference module.
"""


__all__ = ["windowInference"]


import FWCore.ParameterSet.Config as cms


windowInference = cms.EDAnalyzer("WindowInference",
    graphPath=cms.string("graph.pb"),
    inputTensorName=cms.string("input"),
    outputTensorName=cms.string("output"),
)
