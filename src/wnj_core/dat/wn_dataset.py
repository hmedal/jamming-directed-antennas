'''
Created on Jun 14, 2014

@author: hmedal
'''

import lxml.etree as etree
import numpy as np
from gurobipy import *
import itertools
from multiprocessing import Pool
import time
import networkx as nx
#from pulp import *
import argparse
import ast

class WN_Dataset(object):
    '''
    classdocs
    '''
    ids = None
    numNodes = 0
    coors = None
    commodities = {}
    nodeOnlyGraph = None
    pos = None
    
    def __init__(self, path):
        '''
        Constructor
        '''
        print "path", path
        self.readInDataset_and_CreateNodesAndCommodities(path)
    
    def readInDataset_and_CreateNodesAndCommodities(self, path):
        global ids, numNodes, coors, commodities, nodeOnlyGraph
        #print "path: ", path
        #global numFacs, demPtWts, numDemPts, capacities, pairsDistMatrix
        d = etree.parse(open(path))
        #print "d", d
        #facility information
        self.ids = [int(i) for i in d.xpath('//node/id[1]/text()')]
        self.numNodes = len(self.ids)
        self.coors = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
        self.distToClosest = [ast.literal_eval(i) for i in d.xpath('//node/distToClosest[1]/text()')]
        coorsTuples = [(i[0], i[1]) for i in self.coors]
        #print "coors", self.coors
        #print "coorsTuples", coorsTuples
        self.commodities = {}
        self.nodeOnlyGraph = nx.Graph()
        for index in range(len(self.ids)):
            self.nodeOnlyGraph.add_node(self.ids[index], coor = coorsTuples[index], distToClosest = self.distToClosest[index])
        for commod in d.xpath('//odPair'):
            id = int(commod.xpath('./id/text()')[0])
            self.commodities[id] = {}
            origin = int(commod.xpath('./origin/text()')[0])
            destination = int(commod.xpath('./destination/text()')[0])
            demand = float(commod.xpath('./demand/text()')[0])
            self.commodities[id]['odPair'] = (origin, destination)
            self.commodities[id]['demand'] = demand 