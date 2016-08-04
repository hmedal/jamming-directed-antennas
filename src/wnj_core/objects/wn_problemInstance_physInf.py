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
from src.wnj_core.objects.wn_problemInstance import WN_Instance

FUZZ = 0.0001

class WN_Instance_PhysInf(WN_Instance):
    '''
    classdocs
    '''
    #params
    ambientNoise = 0
    snrThresh = 0
    proportionConstant = 0
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
        global numJamLocs,multRadiosPerNode, numChannels, snrThresh, proportionConstant # instance params
        d = etree.parse(open(path))
        self.ambientNoise = bool(int(d.xpath('//instance/ambientNoise[1]/text()')[0]))
        self.snrThresh = bool(int(d.xpath('//instance/snrThresh[1]/text()')[0]))
        self.proportionConstant = bool(int(d.xpath('//instance/proportionConstant[1]/text()')[0]))
        self.multRadiosPerNode = bool(int(d.xpath('//instance/multRadiosPerNode[1]/text()')[0]))
        self.numChannels = int(d.xpath('//instance/numChannels[1]/text()')[0])
        maxNumHopsMult = float(d.xpath('//instance/maxNumHopsMult[1]/text()')[0])
        self.linkCap = 1.0
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
        self.graphFinal = nx.MultiDiGraph()
        counter = 0
        if(self.multRadiosPerNode):
            numRadiosPerNode = self.numChannels
        else:
            numRadiosPerNode = 1
        #print "numRadiosPerNode", numRadiosPerNode
        for radioTypeIndex in range(numRadiosPerNode):
            for n in sorted(self.nodeOnlyGraphWithAttr.nodes()):
                self.graphFinal.add_node((n, radioTypeIndex), power = self.nodeOnlyGraphWithAttr.node[n]['power'], 
                                               coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], radioType = radioTypeIndex, index = counter)
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
        global interferenceGraph
        if(self.interfModelType == 'simple-protocol' or self.interfModelType == '802.11-MAC-protocol'):
            self.interferenceGraph = self.createConflictGraph_Protocol(self.CGraph, self.interfModelType)
        elif(self.interfModelType == 'simple-physical'):
            self.interferenceGraph = self.createConflictGraph_PhysicalModel(self.CGraph, self.interfModelType)
        #print "infGraph nodes 2", self.interferenceGraph.nodes()
        
    def edgeExistsInConnectivityGraph(self, G, edgeInfoAll, interfModelType, radioTypeIndex, channelIndex):
        edge = (edgeInfoAll[0], edgeInfoAll[1])
        edgeInfo = edgeInfoAll[2]
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol') or (interfModelType == 'none')):
            if(edge[0] != edge[1] and edgeInfo['dist'] <= (G.node[edge[0]]['commRange'] + FUZZ)):
                return True
        elif interfModelType == 'physical':
        	if (edge[0] != edge[1]):
        		if self.getSNR(self, G, edge, includeNoiseFromOtherNodes = False) <= self.snrThresh:
				 	return True
			else:
				return False
		return False
        
    def createConflictGraph_Protocol(self, G, interfModelType):
        #print "createConflictGraph_Protocol", G.edges(data=True)
        interferenceGraph = nx.DiGraph()
        #print "reading", len(G.edges()), "edges"
        for edgeInfo1 in G.edges(data = True):
            for edgeInfo2 in G.edges(data = True):
                #print edgeInfo1, edgeInfo2
                if(edgeInfo1 != edgeInfo2):
                    if(self.cannotBeActiveSimultaneously(G, edgeInfo1, edgeInfo2, interfModelType) is True):
                        edge1 = (edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['channel'])
                        edge2 = (edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['channel'])
                        #print "add:", edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['radioType'], edgeInfo1[2]['channel'], "***", edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['radioType'], edgeInfo2[2]['channel']
                        interferenceGraph.add_edge(edge1, edge2)
                    #else:
                        #print "don't add:", edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['radioType'], edgeInfo1[2]['channel'], "***", edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['radioType'], edgeInfo2[2]['channel']
        #nx.draw(ConflictGraph)
        #plt.show() # display
        for edgeInfo in G.edges(data = True):
            #print "edgeInfo", edgeInfo
            if edgeInfo[2]['edgeType'] == 'virtual':
                interferenceGraph.add_node((edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']))
        #for node in interferenceGraph.nodes(data = True):
        #    print "infGraph node", node
        #for edge in interferenceGraph.edges(data = True):
        #    print "infGraph edge", edge
        return interferenceGraph
    
    def createConflictGraph_PhysicalModel(self, G, interfModelType):
        None
        #print G.edges(data=True)
        #ConflictGraph = nx.DiGraph()
        #for edge1 in G.edges():
        #    for edge2 in G.edges():
        #        if(edge1 != edge2):
        #            if(interfModelType == 'simple-physical'):
        #                ConflictGraph.add_edge(edge1, edge2, weight = getFractionOfMaxPermissibleNoise(G, edge1, edge2))
        #nx.draw(ConflictGraph)
        #plt.show() # display
        #return ConflictGraph
        
    def cannotBeActiveSimultaneously(self, G, edge1, edge2, interfModelType):
        if(interfModelType == 'simple-protocol'):
            return self.cannotBeActiveSimultaneously_simpleProtocolModel(G, edge1, edge2)
        elif(interfModelType == '802.11-MAC-protocol'):
            #print "is 802-11"
            return self.cannotBeActiveSimultaneously_802_11_ProtocolModel(G, edge1, edge2)
        
    def cannotBeActiveSimultaneously_simpleProtocolModel(self, G, edge1, edge2):
        i=edge1[0]
        j=edge1[1]
        p=edge2[0]
        q=edge2[1]
        #print "cannotBeActiveSimultaneously_simpleProtocolModel", i, j, p, q
        #print i, j, p, q
        if(G.has_edge(i,q) and G.edge[i][q]['dist'] <= G.node[i]['interferenceRange']):
            return True
        if(G.has_edge(p,j) and G.edge[p][j]['dist'] <= G.node[p]['interferenceRange']):
            return True
    
    def cannotBeActiveSimultaneously_802_11_ProtocolModel(self, G, edgeInfo1, edgeInfo2):
        i=edgeInfo1[0]
        j=edgeInfo1[1]
        radioType1 = edgeInfo1[2]['radioType']
        channel1 = edgeInfo1[2]['channel']
        edgeType1 = edgeInfo1[2]['edgeType']
        p=edgeInfo2[0]
        q=edgeInfo2[1]
        radioType2 = edgeInfo2[2]['radioType']
        channel2 = edgeInfo2[2]['channel']
        edgeType2 = edgeInfo2[2]['edgeType']
        #print "cannotBeActiveSimultaneously_802_11_ProtocolModel", edge1, edge2
        #print i, q, G.edge[i][q]['dist']
        #print i, j, radioType1, channel1, "...", p, q, radioType2, channel2
        tuplesOfNodePairs = [(i,q), (q,i), (i,p), (p,i), (j,p), (p,j), (j,q), (q,j)]
        #print "tuplesOfNodePairs", i,j, p, q, tuplesOfNodePairs
        #print "types", edgeInfo1, edgeInfo2
        
        if(edgeType1 == 'virtual' or edgeType2 == 'virtual'):
            return False
        if(i == p or j == q):#same head and/or tail (interference due to same radio)
            if(radioType1 == radioType2):# cannot have same radio if same head or tail
                return True
            else:
                return False
        else:
            counter = 0
            for (a,b) in tuplesOfNodePairs:
                dist = self.euclidDist(G.node[a]['coor'], G.node[b]['coor'])
                #print i,j, a, b, G.node[a]['coor'], G.node[b]['coor'], "dist", dist
                #if(G.has_edge(a,b)):
                    #print G.edge[a][b][0], G.node[a]
                    #dist = self.euclidDist(G.node[a]['coor'], G.node[b]['coor'])
                    #print i,j, a, b, G.node[a]['coor'], G.node[b]['coor'], "dist", dist
                if(dist <= G.node[a]['interferenceRange'] and channel1 == channel2):
                    return True
                counter += 1
    
    def getSignalStrength(self, distance):
        return self.proportionConstant/distance**2
    
    def getSNR(self, G, edge, includeNoiseFromOtherNodes = False):
        noiseFromOtherNodes = 0
        if includeNoiseFromOtherNodes:
            for otherNode in G.nodes():
                if otherNode != edge[0]:
                    noiseFromOtherNodes += self.getSignalStrength(self, G.edge[otherNode][edge[1]]['dist'])
        return self.getSignalStrength(G.edge[edge[0]][edge[1]]['dist'])/(self.ambientNoise + noiseFromOtherNodes)

#     def getFractionOfMaxPermissibleNoise(self, G, edge1, edge2):
#         i=edge1[0]
#         j=edge1[1]
#         p=edge2[0]
#         if(p != j):
#             #print "getFractionOfMaxPermissibleNoise", edge1, edge2
#             #print G.edges()
#             #print p, j, G.edge[p][j]['dist']
#             if((p,j) in G.edges()):
#                 distance = G.edge[p][j]['dist']
#             else:
#                 distance = np.linalg.norm(np.array(originalGraph.node[p]['coor']) - np.array(originalGraph.node[j]['coor']))
#         else:
#             distance = 0.01
#         ss_pj = self.getSignalStrength(G.node[p]['power'],distance)
#         ss_ij = self.getSignalStrength(G.node[i]['power'], G.edge[i][j]['dist'])
#         #print edge1, edge2, distance, ss_pj, ss_pj / ((ss_ij/snr_threshold) - ambient_noise)
#         return  ss_pj / ((ss_ij/snr_threshold) - ambient_noise)