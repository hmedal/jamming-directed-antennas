'''
Created on Jun 14, 2014

@author: hmedal
'''

import lxml.etree as etree
import numpy as np
import itertools
import ast
import networkx as nx
import math

FUZZ = 0.0001

class WN_Instance(object):
    '''
    classdocs
    '''
    #params
    commRange = 0
    multRadiosPerNode = False
    numChannels = 0
    maxNumHops = 0
    gridSize = 0
    # objects
    myDataset = None
    nodeOnlyGraphWithAttr = None
    graphFinal = None
    CGraph = None
    interferenceGraph = None
    #stuff from dataset
    ids = None
    numNodes = 0
    coors = None
    commodities = {}
    
    def __init__(self, path, dataset):
        '''
        Constructor
        '''
        #print "path", path
        #new from dataset
        self.ids = dataset.ids
        self.coors = dataset.coors
        self.numNodes = dataset.numNodes
        self.commodities = dataset.commodities        
        self.myDataset = dataset
        #instance parameters
        self.readInExperimentData(path)
        self.addAttributesAndObjectsTo_NodesOnlyGraph_FromInstance(dataset)
        self.createConnectivityGraph_FromDatasetAndInstance()
        self.createInterferenceGraph_FromDatasetAndInstance()
        
    def readInExperimentData(self, path):
        global numJamLocs, commRange, infRange, multRadiosPerNode, numChannels # instance params
        d = etree.parse(open(path))
        self.commRange = float(d.xpath('//instance/commRange[1]/text()')[0])
        self.infRange = float(d.xpath('//instance/infRange[1]/text()')[0])
        self.multRadiosPerNode = bool(int(d.xpath('//instance/multRadiosPerNode[1]/text()')[0]))
        self.numChannels = int(d.xpath('//instance/numChannels[1]/text()')[0])
        maxNumHopsMult = float(d.xpath('//instance/maxNumHopsMult[1]/text()')[0])
        self.linkCap = 1.0
        self.interfModelType = '802.11-MAC-protocol'
        self.modelType = d.xpath('//model/modelType[1]/text()')[0]
        #derived values
        self.gridSize = int(math.sqrt(self.numNodes))
        minPossibleHops = (self.gridSize - 1) * 2 + 2 # +2 is to include edge from origin to first real node and edge from last real node to destination
        self.maxNumHops = math.ceil(maxNumHopsMult * minPossibleHops)
    
    def euclidDist(self, coor1, coor2):
        return np.linalg.norm(np.array(coor1) - np.array(coor2))
    
    def addAttributesAndObjectsTo_NodesOnlyGraph_FromInstance(self, dataset):
        global nodeOnlyGraphWithAttr
        #gridSize = int(math.sqrt(dataset.numNodes))
        print "nodeOnlyGraph has ", len(dataset.nodeOnlyGraph.nodes()), "nodes and ", len(dataset.nodeOnlyGraph.edges()), "edges"
        self.nodeOnlyGraphWithAttr = dataset.nodeOnlyGraph.copy()
        for i in self.nodeOnlyGraphWithAttr.nodes():
            self.nodeOnlyGraphWithAttr.node[i]['power']= 1.0
            self.nodeOnlyGraphWithAttr.node[i]['commRange']= self.commRange
            self.nodeOnlyGraphWithAttr.node[i]['interferenceRange']= self.infRange
        print "nodes", self.nodeOnlyGraphWithAttr.nodes(data = True)
        counter = 0
        for headNode in self.nodeOnlyGraphWithAttr.nodes():
            for tailNode in self.nodeOnlyGraphWithAttr.nodes():
                coor1 = self.nodeOnlyGraphWithAttr.node[headNode]['coor']
                coor2 = self.nodeOnlyGraphWithAttr.node[tailNode]['coor']
                distCalc = self.euclidDist(coor1, coor2)
                if((self.interfModelType == 'simple-protocol') or (self.interfModelType == '802.11-MAC-protocol') or (self.interfModelType == 'none')):
                    if(headNode != tailNode and distCalc <= (self.nodeOnlyGraphWithAttr.node[headNode]['commRange'] + FUZZ)):
                        self.nodeOnlyGraphWithAttr.add_edge(headNode, tailNode, dist = distCalc, number = counter, capacity = self.linkCap)
                        counter += 1
        #print "nodeOnlyGraphWithAttr edges", self.nodeOnlyGraphWithAttr.edges()
        print "nodeOnlyGraphWithAttr has ", len(self.nodeOnlyGraphWithAttr.nodes()), "nodes and ", len(self.nodeOnlyGraphWithAttr.edges()), "edges"
        self.graphFinal = nx.MultiDiGraph()
        #print "nodes", G.nodes()
        #posNodeOnlyGraphEnhanced = {}
        counter = 0
        if(self.multRadiosPerNode):
            numRadiosPerNode = self.numChannels
        else:
            numRadiosPerNode = 1
        #print "numRadiosPerNode", numRadiosPerNode
        for radioTypeIndex in range(numRadiosPerNode):
            for n in sorted(self.nodeOnlyGraphWithAttr.nodes()):
                self.graphFinal.add_node((n, radioTypeIndex), power = self.nodeOnlyGraphWithAttr.node[n]['power'], 
                                         commRange = self.nodeOnlyGraphWithAttr.node[n]['commRange'], 
                                               coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], 
                             interferenceRange = self.nodeOnlyGraphWithAttr.node[n]['interferenceRange'], radioType = radioTypeIndex, index = counter)
                counter += 1
                self.graphFinal.add_node((n), coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], index = counter)#add super node
                self.graphFinal.add_edge((n), (n, radioTypeIndex), radioTypeIndex, edgeType = 'virtual')
                #posNodeOnlyGraphEnhanced[(n, radioTypeIndex)] = posG[n]
                #posNodeOnlyGraphEnhanced[(n)] = posG[n]
                counter += 1
        for n in sorted(self.nodeOnlyGraphWithAttr.nodes()):
            self.graphFinal.add_node((n), coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], interferenceRange = 0.0, commRange = float('inf'))#add super node
            #posNodeOnlyGraphEnhanced[(n)] = posG[n]
            for radioTypeIndex in range(numRadiosPerNode):
                self.graphFinal.add_edge((n), (n, radioTypeIndex), radioTypeIndex, dist = 0.0, radioType = radioTypeIndex, 
                                         capacity = 1000.0, channel = -1, origNum = -1, edgeType = 'virtual')
                self.graphFinal.add_edge((n, radioTypeIndex), (n), radioTypeIndex, dist = 0.0, radioType = radioTypeIndex, 
                                         capacity = 1000.0, channel = -1, origNum = -1, edgeType = 'virtual')
        #print "edges1", self.nodeOnlyGraphWithAttr.edges(data = True)
        for e in sorted(self.nodeOnlyGraphWithAttr.edges()):
            distBetween = self.nodeOnlyGraphWithAttr.edge[e[0]][e[1]]['dist']
            myNum = self.nodeOnlyGraphWithAttr.edge[e[0]][e[1]]['number']
            for radioTypeIndex in range(numRadiosPerNode):
                if(self.multRadiosPerNode):
                    self.graphFinal.add_edge((e[0], radioTypeIndex), (e[1], radioTypeIndex), radioType = radioTypeIndex, 
                                             channel = radioTypeIndex, dist = distBetween, capacity = self.linkCap, origNum = myNum, edgeType = 'real')
                    self.graphFinal.add_edge((e[1], radioTypeIndex), (e[0], radioTypeIndex), radioType = radioTypeIndex, 
                                             channel = radioTypeIndex, dist = distBetween, capacity = self.linkCap, origNum = myNum, edgeType = 'real')
                    #print "add", e, radioTypeIndex, radioTypeIndex, MDG.edge[e[0]][e[1]][radioTypeIndex]
                else:
                    for channelIndex in range(self.numChannels):
                        self.graphFinal.add_edge((e[0], radioTypeIndex), (e[1], radioTypeIndex), radioType = radioTypeIndex, 
                                                 channel = channelIndex, dist = distBetween, capacity = self.linkCap, 
                                                 origNum = myNum, edgeType = 'real')
                        self.graphFinal.add_edge((e[1], radioTypeIndex), (e[0], radioTypeIndex), radioType = radioTypeIndex, 
                                                 channel = channelIndex, dist = distBetween, capacity = self.linkCap, 
                                                 origNum = myNum, edgeType = 'real')
        print "graphFinal has ", len(self.graphFinal.nodes()), "nodes and ", len(self.graphFinal.edges()), "edges"
            
    def createConnectivityGraph_FromDatasetAndInstance(self):
        global CGraph
        self.CGraph = nx.MultiDiGraph()
        #CGraph = nx.convert.convert_to_directed(CGraph)
        #CGraph.add_nodes_from(originalGraph)
        #print "nodes", originalGraph.nodes()
        
        for i in self.graphFinal.nodes():
            #print "node", i, originalGraph.node[i]
            if isinstance(i, int):
                self.CGraph.add_node(i, type = 'super', interferenceRange = 0.0, commRange = float('inf'), coor = self.graphFinal.node[i]['coor'])
            else:
                self.CGraph.add_node(i, power = self.graphFinal.node[i]['power'], commRange = self.graphFinal.node[i]['commRange'], 
                        interferenceRange = self.graphFinal.node[i]['interferenceRange'], 
                        radioType = self.graphFinal.node[i]['radioType'], coor = self.graphFinal.node[i]['coor'])
        for edge in sorted(self.graphFinal.edges(data = True)):
            edgeAttr = edge[2]
            if(self.edgeExistsInConnectivityGraph(self.graphFinal, edge, self.interfModelType, edgeAttr['radioType'], edgeAttr['channel'])):
                #edgeIndex = edgeAttr.
                myDist = edgeAttr['dist']
                myNum = edgeAttr['origNum']
                myCap = edgeAttr['capacity']
                #print "added", edge, edgeIndex, radioTypeIndex, channelIndex, myDist, myNum, myCap
                self.CGraph.add_edge(edge[0], edge[1], radioType = edgeAttr['radioType'], channel = edgeAttr['channel'], 
                                     dist = myDist, origNum = myNum, capacity = myCap, edgeType = edgeAttr['edgeType'])
        print "CGraph has ", len(self.CGraph.nodes()), "nodes and ", len(self.CGraph.edges()), "edges"
            
    def createInterferenceGraph_FromDatasetAndInstance(self):
        None
        
    def edgeExistsInConnectivityGraph(self, G, edgeInfoAll, interfModelType, radioTypeIndex, channelIndex):
        edge = (edgeInfoAll[0], edgeInfoAll[1])
        edgeInfo = edgeInfoAll[2]
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol') or (interfModelType == 'none')):
            if(edge[0] != edge[1] and edgeInfo['dist'] <= (G.node[edge[0]]['commRange'] + FUZZ)):
                return True
        
    def createConflictGraph(self, G):
        None