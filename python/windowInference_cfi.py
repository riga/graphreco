# coding: utf-8

"""
Initialization file for the WindowInference module.
"""


__all__ = ["windowInference"]


import math

import FWCore.ParameterSet.Config as cms


windowInference = cms.EDAnalyzer("WindowInference",
    # the collections of rechits to use
    recHitCollections=cms.VInputTag(
        cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    ),
    # phi range of rechits to use
    minPhi=cms.double(0.0),
    maxPhi=cms.double(2 * math.pi),
    # eta range of rechits to use
    minEta=cms.double(1.6),
    maxEta=cms.double(3.0),
    # window size in phi and eta
    deltaPhi=cms.double(0.4),
    deltaEta=cms.double(0.4),
    # overlap in phi and eta
    overlapPhi=cms.double(0.2),
    overlapEta=cms.double(0.2),
    # graph to the trained model
    graphPath=cms.string("graph.pb"),
    # names of the input and output tensors
    inputTensorName=cms.string("input"),
    outputTensorName=cms.string("output"),
    # whether or not the model in the graph expects a batch dimension
    batchedModel=cms.bool(True),
    # dimension of the padding of the second dimension, i.e., the rec hits themselves
    padSize=cms.uint32(100),
)
