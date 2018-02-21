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
import csv

FUZZ = 0.0001


class Instance(object):
    '''
    classdocs
    '''
    #params
    numJamLocs = 0
    commRange = 0
    infRange = 0
    multRadiosPerNode = False
    numChannels = 0
    jamBudget = 0
    jamRange = 0
    maxNumHops = 0
    gridSize = 0
    # objects
    myDataset = None
    nodeOnlyGraphWithAttr = None
    graphFinal = None
    CGraph = None
    interferenceGraph = None
    interferenceGraph_nodesOnly = None
    jamGraph = None
    #stuff from dataset
    tcurr = 0
    trec = 0
    battCap = 0
    ids = None
    numNodes = 0
    coors = None
    commodities = {}
    transmissiondistance = 0
    
    
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
        self.tcurr = dataset.tcurr
        self.trec = dataset.trec
        self.battCap = dataset.battCap     
        self.myDataset = dataset
        #instance parameters
        self.readInExperimentData(path)
        print "experiment data read into instance"
        self.addAttributesAndObjectsTo_NodesOnlyGraph_FromInstance(dataset)
        print "nodes only graph created"
        self.createConnectivityGraph_FromDatasetAndInstance()
        print "connectivity graph created"
        self.createInterferenceGraph_FromDatasetAndInstance()
        print "inf graph created"
        self.createJammingGraph_FromDatasetAndInstance()
        print "jamming graph created"
        
    def readInExperimentData(self, path):
        global numJamLocs, commRange, infRange, multRadiosPerNode, numChannels, jamBudget, jamRange, infRangeMult # instance params
        d = etree.parse(open(path))
        self.numJamLocs = float(d.xpath('//instance/numJamLocs[1]/text()')[0])
        self.commRange = float(d.xpath('//instance/commRange[1]/text()')[0])
        self.infRangeMult = float(d.xpath('//instance/infRangeMult[1]/text()')[0])
        self.infRange = self.commRange * self.infRangeMult#print "boolTest", bool('F')
        #print "multRadiosStr", d.xpath('//instance/multRadiosPerNode[1]/text()')[0], bool(str(d.xpath('//instance/multRadiosPerNode[1]/text()')[0]))
        self.multRadiosPerNode = bool(int(d.xpath('//instance/multRadiosPerNode[1]/text()')[0]))
        self.numChannels = int(d.xpath('//instance/numChannels[1]/text()')[0])
        self.jamBudget= float(d.xpath('//instance/jamBudget[1]/text()')[0])
        self.jamRange = float(d.xpath('//instance/jamRange[1]/text()')[0])
        #self.jamRangeMult = float(d.xpath('//instance/jamRangeMult[1]/text()')[0])
        maxNumHopsMult = float(d.xpath('//instance/maxNumHopsMult[1]/text()')[0])
        self.linkCap = 1.0
        self.interfModelType = '802.11-MAC-protocol'
        self.modelType = d.xpath('//model/modelType[1]/text()')[0]
        #derived values
        self.gridSize = int(math.sqrt(self.numNodes))
        minPossibleHops = (self.gridSize - 1) * 2 + 2 # +2 is to include edge from origin to first real node and edge from last real node to destination
        self.maxNumHops = math.ceil(maxNumHopsMult * minPossibleHops)
        #print "maxNumHopsMult", maxNumHopsMult
        #print "self.maxNumHops", self.maxNumHops
    
    def euclidDist(self, coor1, coor2):
        return np.linalg.norm(np.array(coor1) - np.array(coor2))
    
    def addAttributesAndObjectsTo_NodesOnlyGraph_FromInstance(self, dataset):
        global nodeOnlyGraphWithAttr
        #gridSize = int(math.sqrt(dataset.numNodes))
        print "nodeOnlyGraph has ", len(dataset.nodeOnlyGraph.nodes()), "nodes and ", len(dataset.nodeOnlyGraph.edges()), "edges"
        self.nodeOnlyGraphWithAttr = dataset.nodeOnlyGraph.copy()
        for i in self.nodeOnlyGraphWithAttr.nodes():
            #self.nodeOnlyGraphWithAttr.node[i]['tcurr']= self.tcurr
            #self.nodeOnlyGraphWithAttr.node[i]['trec']= self.trec
            #self.nodeOnlyGraphWithAttr.node[i]['battCap']= self.battCap
            self.nodeOnlyGraphWithAttr.node[i]['commRange']= self.commRange
            self.nodeOnlyGraphWithAttr.node[i]['interferenceRange']= self.infRange
        print "nodes", self.nodeOnlyGraphWithAttr.nodes(data = True)
        counter = 0
#         for headNode in self.nodeOnlyGraphWithAttr.nodes():
#             print "headnode", headNode
#             for tailNode in self.nodeOnlyGraphWithAttr.nodes():
#                  print "tailnode", tailnode
#                  coor1 = self.nodeOnlyGraphWithAttr.node[headNode]['coor']
#                  coor2 = self.nodeOnlyGraphWithAttr.node[tailNode]['coor']
                #print "coors", np.array(coor1), np.array(coor2)
                #print "coorsDiff", np.array(coor1) - np.array(coor2)
                #distCalc = np.linalg.norm(np.array(self.nodeOnlyGraphWithAttr.node[headNode]['coor']) - 
                #                      np.array(self.nodeOnlyGraphWithAttr.node[tailNode]['coor']))
                #distCalc = self.euclidDist(coor1, coor2)
                #print "test", headNode, tailNode, coor1, coor2, distCalc, self.nodeOnlyGraphWithAttr.node[headNode]['commRange'], self.interfModelType
#                 if((self.interfModelType == 'simple-protocol') or (self.interfModelType == '802.11-MAC-protocol') or (self.interfModelType == 'none')):
#                     #print "possibly add edge", headNode != tailNode, distCalc <= (self.nodeOnlyGraphWithAttr.node[headNode]['commRange'] + FUZZ)
#                     if(headNode != tailNode and distCalc <= (self.nodeOnlyGraphWithAttr.node[headNode]['commRange'] + FUZZ)):
#                         #print "add edge"
#                         self.nodeOnlyGraphWithAttr.add_edge(headNode, tailNode, dist = distCalc, number = counter, capacity = self.linkCap)
#                         counter += 1
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
                #print sorted(self.nodeOnlyGraphWithAttr.nodes())
                #print "n", n
                
                self.graphFinal.add_node((n, radioTypeIndex), tcurr = self.nodeOnlyGraphWithAttr.node[n]['tcurr'], trec = self.nodeOnlyGraphWithAttr.node[n]['trec'], battCap = self.nodeOnlyGraphWithAttr.node[n]['battCap'],  
                                         commRange = self.nodeOnlyGraphWithAttr.node[n]['commRange'], 
                                               coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], 
                             interferenceRange = self.nodeOnlyGraphWithAttr.node[n]['interferenceRange'], radioType = radioTypeIndex, index = counter)
                counter += 1
                self.graphFinal.add_node((n), coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], tcurr = self.nodeOnlyGraphWithAttr.node[n]['tcurr'], trec = self.nodeOnlyGraphWithAttr.node[n]['trec'], battCap = self.nodeOnlyGraphWithAttr.node[n]['battCap'], index = counter)#add super node
                self.graphFinal.add_edge((n), (n, radioTypeIndex), radioTypeIndex, edgeType = 'virtual')
                #posNodeOnlyGraphEnhanced[(n, radioTypeIndex)] = posG[n]
                #posNodeOnlyGraphEnhanced[(n)] = posG[n]
                counter += 1
        for n in sorted(self.nodeOnlyGraphWithAttr.nodes()):
            self.graphFinal.add_node((n), coor = self.nodeOnlyGraphWithAttr.node[n]['coor'], tcurr = self.nodeOnlyGraphWithAttr.node[n]['tcurr'], trec = self.nodeOnlyGraphWithAttr.node[n]['trec'], battCap = self.nodeOnlyGraphWithAttr.node[n]['battCap'], interferenceRange = 0.0, commRange = float('inf'))#add super node
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
                self.CGraph.add_node(i, type = 'super', interferenceRange = 0.0, commRange = float('inf'), coor = self.graphFinal.node[i]['coor'], tcurr = self.nodeOnlyGraphWithAttr.node[i]['tcurr'], trec = self.nodeOnlyGraphWithAttr.node[i]['trec'], battCap = self.nodeOnlyGraphWithAttr.node[i]['battCap'])
            else:
                self.CGraph.add_node(i, tcurr = self.graphFinal.node[i]['tcurr'], trec = self.graphFinal.node[i]['trec'], battCap = self.graphFinal.node[i]['battCap'], commRange = self.graphFinal.node[i]['commRange'], 
                        interferenceRange = self.graphFinal.node[i]['interferenceRange'], 
                        radioType = self.graphFinal.node[i]['radioType'], coor = self.graphFinal.node[i]['coor'])
        for edge in sorted(self.graphFinal.edges(data = True)):
            edgeAttr = edge[2]
            #if(self.edgeExistsInConnectivityGraph(self.graphFinal, edge, self.interfModelType, edgeAttr['radioType'], edgeAttr['channel'])):
            edgeIndex = edgeAttr
        #if self.nodeOnlyGraphWithAttr.has_edge(edge[0],edge[1]):
            myDist = edgeAttr['dist']
            myNum = edgeAttr['origNum']
            myCap = edgeAttr['capacity']
            #print "added", edge, edgeIndex, radioTypeIndex, channelIndex, myDist, myNum, myCap
            self.CGraph.add_edge(edge[0], edge[1], radioType = edgeAttr['radioType'], channel = edgeAttr['channel'], 
                                 dist = myDist, origNum = myNum, capacity = myCap, edgeType = edgeAttr['edgeType'])
        print "CGraph has ", len(self.CGraph.nodes()), "nodes and ", len(self.CGraph.edges()), "edges"
        #print "CGraph edges", self.CGraph.edges()
        #for node in self.CGraph.nodes(data = True):
        #    print "CGraph node", node
        #for edge in self.CGraph.edges(data = True):
        #    print "CGraph edge", edge
            
    def createInterferenceGraph_FromDatasetAndInstance(self):
        global interferenceGraph
        if(self.interfModelType == 'simple-protocol' or self.interfModelType == '802.11-MAC-protocol'):
            self.interferenceGraph = self.createConflictGraph_Protocol(self.CGraph, self.interfModelType)
        elif(self.interfModelType == 'simple-physical'):
            self.interferenceGraph = self.createConflictGraph_PhysicalModel(self.CGraph, self.interfModelType)
        #print "infGraph nodes 2", self.interferenceGraph.nodes()
        
    def createJammingGraph_FromDatasetAndInstance(self):
        global jamGraph, numRadiosPerNode
        gridSizeForJamming = int(math.sqrt(self.numJamLocs))
        self.jamGraph = nx.Graph()
        pos = {}
        if(self.multRadiosPerNode):
            numRadiosPerNode = self.numChannels
        else:
            numRadiosPerNode = 1
        #print "numRadiosPerNode", numRadiosPerNode
        counter = 0
        newnumbercounter = 0
        fixList = []
        #comment
        for radioTypeIndex in range(numRadiosPerNode):
            for i in np.arange(0, gridSizeForJamming):
                for j in np.arange(0, gridSizeForJamming):
                    coordinates = (i/float(gridSizeForJamming - 1), j/float(gridSizeForJamming - 1))
                    #print "NOTE: coordinates are", coordinates[0], coordinates[1]
                    #newcoors={}
                    #newcoors = coordinates
                    #comment
                    #print "NOTE: newcoors are", coordinates[0], coordinates[1]
                    #while new_counter <= gridSizeForJamming:
                    fixList.append([coordinates[0], coordinates[1]])
        print "fixlist is", fixList
                    #Add sys.exit back in to see results at end of run sequence
                    
        #while newnumbercounter < 360:
            
        for radioTypeIndex in range(numRadiosPerNode):
            for i in np.arange(0, gridSizeForJamming):
                for j in np.arange(0, gridSizeForJamming):
                    coordinates = (i/float(gridSizeForJamming - 1), j/float(gridSizeForJamming - 1))
                        #new_counter = new_counter + 1
                    self.jamGraph.add_node((coordinates[0], coordinates[1], radioTypeIndex), name = counter, selected = 0.0, power = 1.0, 
                                      cost = 1.0, interferenceRange = self.jamRange, radioType = radioTypeIndex)
                    #print "add", i, j, counter, radioTypeIndex, g_JamGraph.nodes(data=True)
                    pos[coordinates] = coordinates
                    #newcounter += 1
                    counter += 1
                print "NOTE: jam coors are", coordinates[0], coordinates[1]        
            print "gridSizeForJamming is", gridSizeForJamming
            #print "The numRadiosPerNode is", numRadiosPerNode, "and radioTypeIndex is", radioTypeIndex
            #print "fixList is", fixList
            #import sys
            #sys.exit()
        #import sys
        #sys.exit()
        
    def edgeExistsInConnectivityGraph(self, G, edgeInfoAll, interfModelType, radioTypeIndex, channelIndex):
        edge = (edgeInfoAll[0], edgeInfoAll[1])
        #print "edge[0]", G.node[edge[0]]['coor']
        #print"edge[1]", G.node[edge[1]]['coor']
        edgeInfo = edgeInfoAll[2]
        #edgeIndex = radioTypeIndex * numChannels + channelIndex
        #nodeIndex = edge[0] * numRadiosPerNode + radioTypeIndex
        #print "indices", edgeIndex, nodeIndex
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol') or (interfModelType == 'none')):
            if (edge[0] != edge[1]): #and edgeInfo['dist'] <= (G.node[edge[0]]['commRange'] + FUZZ):
                #return True
            
                with open('/Users/wbl62/Desktop/directed-code/Transmitter_directed_new_high.csv', 'rU') as f:
                    reader1 = csv.reader(f)
                    mycsvlist1 = list(reader1)
                    degreenumber = [x[0] for x in mycsvlist1]
                    distancenumber = [x[1] for x in mycsvlist1]
#if G.has_edge(a,b):
                    #dist = self.euclidDist(edge[0]['coor'], edge[1]['coor'])
                    #print "single COORD", G.node[a]['coor'][0]
                    x1 = G.node[edge[0]]['coor'][0]
                    y1 = G.node[edge[0]]['coor'][1]
                    x2 = G.node[edge[1]]['coor'][0]
                    y2 = G.node[edge[1]]['coor'][1]
                    
                    dist_a_b_new = self.newdist(x1,x2,y1,y2)
                    #print "dist_a_b_new:", dist_a_b_new
                    #sys.exit()
                    #my code edits from 1-15-17 below
                    xdist = x2 - x1
                    ydist = y2 - y1
                    #print "xdist", xdist, "ydist", ydist
                    #check the above
                    #print "xdist", xdist
                    #print "ydist", ydist
                    #if xdist == 0 and ydist == 0:
                        #transmissiondistance = 0
                    if xdist == 0 and ydist != 0:
                        if y2 > y1:
                            anglefound = math.pi/2
                        elif y2 < y1:
                            anglefound = (3/2) * math.pi
                    elif ydist == 0 and xdist != 0:
                        if x2 > x1:
                            anglefound = 0
                        elif x2 < x1:
                            anglefound = math.pi
                    elif xdist == 0 and ydist == 0:
                            anglefound = 0.000001 #placeholder value to deal with node compared to itself (i.e., same x and y coordinates)
                    else:
                        if ydist > 0 and xdist > 0:
                            anglefound = np.arctan((ydist / xdist))
                        elif ydist > 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + math.pi
                        elif ydist < 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + math.pi
                        elif ydist < 0 and xdist > 0: 
                            anglefound = (np.arctan((ydist / xdist))) + (2*math.pi)
                        else:
                            print "Error"
                    #print "anglefound in radians:", anglefound, ",", "anglefound in degrees", int(math.degrees(anglefound))
                    #print "dist_a_b",dist_a_b, "dist_a_b_new", dist_a_b_new
                    #if dist_a_b <= Transmitter[Nodes_custom[count,3]-1,1]:#verifying the transmitter is transmitting far enough
                    #if Nodes_custom[count,2]>=Nodes_custom[dist_b_count,2]-1: #and Nodes_custom[count,2]<=Nodes_custom[dist_b_count,2]+1:
                    #print 'ok' #verifying that along the z axis there is line of sight with the antenna within some tolerance
                    if anglefound == 0.000001: #0.000001 is a placeholder value to deal with having to compare each node with itself
                        transmissiondistance = 0
                        #print "anglefound is 0000001"
                    else: #int(degreenumber[int(anglefound-1)]) == int(math.degrees(anglefound)): #math.degrees(math.int(anglefound)):
                                #print "counter", counter
                        transmissiondistance = distancenumber[int(math.degrees(anglefound))-1]
                        #print "transmissiondistance2222222", transmissiondistance
                                #print "transmissiondistance", transmissiondistance
                
                    
                    print "transmission distance is ", transmissiondistance, "and dist_a_b_new is ", dist_a_b_new
                    #print "transdist", transmissiondistance
                    #print "mutliplication", float(transmissiondistance)*self.infRangeMult
                    if float(transmissiondistance) == 0:
                           
                            #dist_b_count = dist_b_count + 1
                            #print "same node", count, "=", dist_b_count
                            print "no interference"
                    elif (float(transmissiondistance)>=dist_a_b_new):
                        print "they can connect"
                        return True
                

#else:
#                 print "no edge"
#         elif((interfModelType == 'simple-physical')):
#             print"simple-physical"
#         #    if(edge[0] != edge[1] and getSNR(G, edge) >= snr_threshold):
#         #        return True
        
    def createConflictGraph_Protocol(self, G, interfModelType):
        #NOTE: THIS IS WHERE EDGES BECOME NODES IN THE INTERFERENCE GRAPH
        #print "createConflictGraph_Protocol", G.edges(data=True)
        global interferenceGraph_nodesOnly
        interferenceGraph = nx.DiGraph()
        interferenceGraph_nodesOnly = nx.DiGraph()
        #print "reading", len(G.edges()), "edges"
        for edgeInfo1 in G.edges(data = True):
            interferenceGraph_nodesOnly.add_node(edgeInfo1[0])
            for edgeInfo2 in G.edges(data = True):
                #print edgeInfo1, edgeInfo2
                if(edgeInfo1 != edgeInfo2):
                    if(self.cannotBeActiveSimultaneously(G, edgeInfo1, edgeInfo2, interfModelType) is True):
                        edge1 = (edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['channel'])
                        #print "edge1[0], edge1[1], edge1[2]['channel'] is", edge1
                        edge2 = (edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['channel'])
                        #print "edge2[0], edge2[1], edge2[2]['channel'] is", edge2
                        #comment
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
        #    print "infGraph node", interferenceGraph.nodes()
        #for edge in interferenceGraph.edges(data = True):
        #    print "infGraph edge", edge
        print "interferenceGraph", interferenceGraph.nodes()
        #import sys
        #sys.exit()
        return interferenceGraph
    
    def createConflictGraph_PhysicalModel(self, G, interfModelType):
        global filler
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
    
    def newdist(self,x1,x2,y1,y2):
            return np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 ) 
    
    def cannotBeActiveSimultaneously_802_11_ProtocolModel(self, G, edgeInfo1, edgeInfo2):
        global transmissiondistance
        transmissiondistance = 0
        
        with open('/Users/wbl62/Desktop/directed-code/Transmitter_directed_new_high.csv', 'rU') as f:
            reader1 = csv.reader(f)
            mycsvlist1 = list(reader1)
            degreenumber = [x[0] for x in mycsvlist1]
            distancenumber = [x[1] for x in mycsvlist1]
            #print "degreenumber:", degreenumber, "distancenumber:", distancenumber
            #print "degreenumber[180]", degreenumber[180]
            
            #print "distance_number", distancenumber
            
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
        
        #tuplesOfNodePairs = [(i,j), (p,q)]

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
                if G.has_edge(a,b):
                    dist = self.euclidDist(G.node[a]['coor'], G.node[b]['coor'])
                    #print "single COORD", G.node[a]['coor'][0]
                    x1 = G.node[a]['coor'][0]
                    y1 = G.node[a]['coor'][1]
                    x2 = G.node[b]['coor'][0]
                    y2 = G.node[b]['coor'][1]
                    
                    dist_a_b_new = self.newdist(x1,x2,y1,y2)
                    #print "dist_a_b_new:", dist_a_b_new
                    #sys.exit()
                    #my code edits from 1-15-17 below
                    xdist = x2 - x1
                    ydist = y2 - y1
                    #print "xdist", xdist, "ydist", ydist
                    #check the above
                    #print "xdist", xdist
                    #print "ydist", ydist
                    #if xdist == 0 and ydist == 0:
                        #transmissiondistance = 0
                    if xdist == 0 and ydist != 0:
                        if y2 > y1:
                            anglefound = math.pi/2
                        elif y2 < y1:
                            anglefound = (3/2) * math.pi
                    elif ydist == 0 and xdist != 0:
                        if x2 > x1:
                            anglefound = 0
                        elif x2 < x1:
                            anglefound = math.pi
                    elif xdist == 0 and ydist == 0:
                            anglefound = 0.000001 #placeholder value to deal with node compared to itself (i.e., same x and y coordinates)
                    else:
                        if ydist > 0 and xdist > 0:
                            anglefound = np.arctan((ydist / xdist))
                        elif ydist > 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + math.pi
                        elif ydist < 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + math.pi
                        elif ydist < 0 and xdist > 0: 
                            anglefound = (np.arctan((ydist / xdist))) + (2*math.pi)
                        else:
                            print "Error"
                    #print "anglefound in radians:", anglefound, ",", "anglefound in degrees", int(math.degrees(anglefound))
                    #print "dist_a_b",dist_a_b, "dist_a_b_new", dist_a_b_new
                    #if dist_a_b <= Transmitter[Nodes_custom[count,3]-1,1]:#verifying the transmitter is transmitting far enough
                    #if Nodes_custom[count,2]>=Nodes_custom[dist_b_count,2]-1: #and Nodes_custom[count,2]<=Nodes_custom[dist_b_count,2]+1:
                    #print 'ok' #verifying that along the z axis there is line of sight with the antenna within some tolerance
                    if anglefound == 0.000001: #0.000001 is a placeholder value to deal with having to compare each node with itself
                        transmissiondistance = 0
                        #print "anglefound is 0000001"
                    else: #int(degreenumber[int(anglefound-1)]) == int(math.degrees(anglefound)): #math.degrees(math.int(anglefound)):
                                #print "counter", counter
                        transmissiondistance = distancenumber[int(math.degrees(anglefound))-1]
                                #print "transmissiondistance2222222", transmissiondistance
                                #print "transmissiondistance", transmissiondistance
                
                    
                    #print "transmission distance is ", transmissiondistance, "and dist_a_b_new is ", dist_a_b_new
                    #print "transdist", transmissiondistance
                    #print "mutliplication", float(transmissiondistance)*self.infRangeMult
                    if float(transmissiondistance) == 0:
                           
                            #dist_b_count = dist_b_count + 1
                            #print "same node", count, "=", dist_b_count
                            print "no interference"
                    elif (float(transmissiondistance)*self.infRangeMult >= dist_a_b_new  and channel1 == channel2):
                        #print "cannot be active simultaneously"
                        return True
                        
#                             self.nodeOnlyGraph.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
#                             self.nodeOnlyGraph.add_node(self.ids[dist_b_count], coor = coorsTuples[dist_b_count], distToClosest = self.distToClosest[dist_b_count], tcurr = self.tcurr[dist_b_count], trec = self.trec[dist_b_count], battCap = self.battCap[dist_b_count])
#                             self.nodeOnlyGraph.add_edge(count, dist_b_count, dist = dist_a_b_new, capacity = 1, number = count)
                        
#                             
#                             G_dirAnt.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
#                             
#                             G_dirAnt.add_edge(count, dist_b_count)
                            
#                     elif (float(transmissiondistance)*self.infRangeMult < dist_a_b_new):
                    else:
                        #print "CAN BE ACTIVE TOGETHER"
                        return False
                        
                            #dist_b_count = dist_b_count + 1
                            #print "failure with edge between", count, "and ", dist_b_count
                        
                    
                    
                    #print i,j, a, b, G.node[a]['coor'], G.node[b]['coor'], "dist", dist
                    #if(G.has_edge(a,b)):
                        #print G.edge[a][b][0], G.node[a]
                        #dist = self.euclidDist(G.node[a]['coor'], G.node[b]['coor'])
                        #print i,j, a, b, G.node[a]['coor'], G.node[b]['coor'], "dist", dist
                    #if(dist <= G.node[a]['interferenceRange'] and channel1 == channel2):
                        #print "infRange", G.node[a]['interferenceRange']
                        #print "   added b/c of dist:", a, b, channel1, channel2, counter, 'dist', G.edge[a][b][0]['dist']
                    #  return True
                counter += 1
    
    def getSignalStrength(self, power, distance):
        return power/distance**2

#     def getSNR(self, G, edge):
#         return self.getSignalStrength(G.node[edge[0]]['power'], G.edge[edge[0]][edge[1]]['dist'])/ambient_noise
#             
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