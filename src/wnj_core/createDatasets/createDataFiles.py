import math
from math import sqrt
import os
try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx
import numpy as np
import lxml.etree as etree
import xml.etree.cElementTree as ET
from ast import literal_eval as make_tuple

G = None
commodities = {}
nodeNames = []
nodeXCoor = []
nodeYCoor = []
normalizedXCoors = []
normalizedYCoors = []
numNodes = 0

dataFilePath = '/home/hmedal/Documents/2_msu/research_manager/data/wirelessNet'

def read_and_convert_MIT_net():
    global nodeNames, nodeXCoor, nodeYCoor
    nodeNames = []
    nodeXCoor = []
    nodeYCoor = []
    f = open('/home/hmedal/Documents/2_msu/research_manager/data/wirelessNet/mit-roofnet/mit_roofnet_nodes.txt', 'r')
    
    for line in f:
        lineArray = line.split()
        #print "lineArray", lineArray, len(lineArray)
        nodeNames.append(lineArray[0])
        nodeXCoor.append(float(lineArray[1]))
        nodeYCoor.append(float(lineArray[2]))
    createGraphForCustomDataset(drawGraph = True, graphName = 'mit')
    f.close()
        
def read_and_convert_Berkeley_net():
    global nodeNames, nodeXCoor, nodeYCoor
    nodeNames = []
    nodeXCoor = []
    nodeYCoor = []
    f = open('/home/hmedal/Documents/2_msu/research_manager/data/wirelessNet/berkeley-net/mote_locs.txt', 'r')
    for line in f:
        lineArray = line.split()
        #print "lineArray", lineArray, len(lineArray)
        nodeNames.append(lineArray[0])
        nodeXCoor.append(float(lineArray[1]))
        nodeYCoor.append(float(lineArray[2]))
    createGraphForCustomDataset(drawGraph = True, graphName = 'berkeley')
    f.close()

def createGraphForCustomDataset(drawGraph = True, graphName = 'noGraphName', includeCommRange = True, commRange = 80000/6.0):
    global G, normalizedXCoors, normalizedYCoors
    G = nx.Graph()
    xCoorRange = max(nodeXCoor) - min(nodeXCoor)
    yCoorRange = max(nodeYCoor) - min(nodeYCoor)
    normalizedXCoors = []
    normalizedYCoors = []
    #print 'ranges', xCoorRange, yCoorRange
    posG = {}
    for nodeIndex in range(len(nodeNames)):
        normalizedXCoor = (nodeXCoor[nodeIndex] - min(nodeXCoor)) / xCoorRange
        normalizedYCoor = (nodeYCoor[nodeIndex] - min(nodeYCoor)) / yCoorRange
        #print 'normalized', normalizedXCoor, normalizedYCoor
        G.add_node(nodeNames[nodeIndex], coor = (normalizedXCoor, normalizedYCoor))
        posG[nodeNames[nodeIndex]] = (normalizedXCoor, normalizedYCoor)
        normalizedXCoors.append(normalizedXCoor)
        normalizedYCoors.append(normalizedYCoor)
    for nodeIndex in range(len(nodeNames)):
        G.node[nodeNames[nodeIndex]]['distToClosest'] = getMaxClosestDistanceFromNode(nodeIndex)
    #print "normalizedXCoor", normalizedXCoors
    print G.nodes(data = True)
    if drawGraph:
        nx.draw_networkx_nodes(G, posG, node_size = commRange, with_labels = False, alpha = 0.2)
        nx.draw(G, posG, node_size = 500, font_size = 14, node_color = 'white', with_labels = True) 
        plt.savefig("/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/figures/" + graphName + "_nodesOnly.pdf") # save as jpg
        plt.show() # display
        jamGraph, posJam = createJammingGraph(25, 1/3.0)
        print(jamGraph.nodes(data = True))
        nx.draw_networkx_nodes(jamGraph, posJam, node_size = commRange*2, with_labels = False, alpha = 0.2, node_color = 'blue')
        nx.draw_networkx_nodes(jamGraph, posJam, node_size = 500, with_labels = False, node_color = 'blue')
        nx.draw(G, posG, node_size = 500, font_size = 14, node_color = 'white', with_labels = True)
        plt.savefig("/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/figures/" + graphName + "_jamLocsAndArcs.pdf") # save as jpg
        plt.show() # display
    getMaxClosestDistance()
    
def createJammingGraph(numJamLocs, jamRange):
    gridSizeForJamming = int(math.sqrt(numJamLocs))
    jamGraph = nx.Graph()
    pos = {}
    #print "numRadiosPerNode", numRadiosPerNode
    counter = 0
    for i in np.arange(0, gridSizeForJamming):
        for j in np.arange(0, gridSizeForJamming):
            coordinates = (i/float(gridSizeForJamming - 1), j/float(gridSizeForJamming - 1))
            jamGraph.add_node(counter)
            #print "add", i, j, counter, jamGraph.nodes(data=True)
            pos[counter] = coordinates
            counter += 1
    print "jam graph", jamGraph.nodes(data = True)
    return [jamGraph, pos]
                    
def getMaxClosestDistanceFromNode(node):
    closestDist = float('inf')
    for nodeIndex2 in range(len(nodeNames)):
        if node != nodeIndex2:
            dist = sqrt((normalizedXCoors[nodeIndex2] - normalizedXCoors[node])**2 + (normalizedYCoors[nodeIndex2] - normalizedYCoors[node])**2)
            #print "dist", dist
            if dist <= closestDist:
                closestDist = dist
    return closestDist

def getDistToClosestNode_NonRealGraph(graph, node):
    closestDist = float('inf')
    for nodeIndex2 in graph:
        if node != nodeIndex2:
            #print make_tuple(str(graph.node[node]['coor']))
            #print node, nodeIndex2, graph.node[node]['coor'], graph.node[nodeIndex2]['coor']
            x1 = make_tuple(str(graph.node[node]['coor']))[0]
            x2 = make_tuple(str(graph.node[nodeIndex2]['coor']))[0]
            y1 = make_tuple(str(graph.node[node]['coor']))[1]
            y2 = make_tuple(str(graph.node[nodeIndex2]['coor']))[1]
            #print x1, x2, y1, y2
            dist = sqrt((x2 - x1)**2 + (y2 - y1)**2)
            #print node, nodeIndex2, "dist", dist
            if dist <= closestDist:
                closestDist = dist
    return closestDist

def getMaxClosestDistance():
    maxClosest = 0
    for nodeIndex1 in range(len(nodeNames)):
        closestDist = getMaxClosestDistanceFromNode(nodeIndex1)
        if closestDist >= maxClosest:
            maxClosest = closestDist
            #print "dist", dist
    #print "numNodes", len(nodeNames)
    print "maxClosest", maxClosest
        
def generate_2D_Grid_Graph(gridSize):
    global G
    G = None
    G = nx.grid_2d_graph(gridSize, gridSize)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted", label_attribute='coor')
    # normalize coordinates
    for n in G.nodes():
        #print "coorOrig", G.node[n]['coor'], gridSize
        prevCoorValue = G.node[n]['coor']
        G.node[n]['coor'] = (prevCoorValue[0]/float(gridSize - 1), prevCoorValue[1]/float(gridSize - 1))
    for n in G.nodes():
        #print n, "coor", G.node[n]['coor']
        G.node[n]['distToClosest'] = getDistToClosestNode_NonRealGraph(G, n)
        #print "coor", G.node[n]['coor']

def generateRandomGraph(n, cubeSize = 1, dim = 2, pos = None):
    global G
    G = None
    G=nx.Graph()
    G.name="Random Graph"
    G.add_nodes_from(range(n))
    #print "nodes", G.nodes()
    if pos is None:
        # random positions 
        for n in G:
            G.node[n]['coor'] = [cubeSize * np.random.random() for i in range(0, dim)]
            #print "coor", G.node[n]['coor']
    else:
        nx.set_node_attributes(G, 'coor', pos)

def generateGraph(numNodes, topology = 'grid', commodLocType = 'corners', commodDemandType = 'rand01'):
    #print "generateGraph", topology
    if topology == 'grid':
        generate_2D_Grid_Graph(int(math.sqrt(numNodes)))
    elif topology == 'rand':
        generateRandomGraph(numNodes, 1)
    else:
        print "ERROR"
        
def readRealGraph(topology = 'mit', commodLocType = 'rand', commodDemandType = 'rand01'):
    if topology == 'mit':
        read_and_convert_MIT_net()
    elif topology == 'berkeley':
        read_and_convert_Berkeley_net()
    else:
        print "ERROR"

def generateCommodityDemand(commodDemandType):
    if commodDemandType == 'inf':
        return float('inf')
    elif commodDemandType == 'unif2':
        return 2.0
    elif commodDemandType == 'unif3':
        return 3.0
    elif commodDemandType == 'unif5':
        return 5.0
    elif commodDemandType == 'unif1':
        return 1.0
    elif commodDemandType == 'rand01':
        return np.random.rand()
    elif commodDemandType == 'rand02':
        return 2.0 * np.random.rand()
    else:
        raise Exception("invalid commodity demand type", commodDemandType)

def addBidirectionalCommodityPair(myList, source, dest, addReverse = True):
    myList.append((source, dest))
    if addReverse:
        myList.append((dest, source))
    
def createCommodities_Corners(gridSize, numCommodities, commodDemandType):
    global commodities
    commodities = {}
    bottomLeft = 0
    bottomLeftShifted = gridSize + 1
    topLeft = gridSize - 1
    topLeftShifted = gridSize * 2 - 2
    topRight = gridSize * gridSize - 1
    topRightShifted = gridSize * (gridSize - 1) - 2
    bottomRight = gridSize * (gridSize - 1)
    bottomRightShifted = gridSize * (gridSize - 2) + 1
    
    halfOfGridSize = gridSize / 2
    leftMid = halfOfGridSize
    topMid = (halfOfGridSize + 1) * gridSize - 1
    rightMid = (gridSize - 1) * gridSize + halfOfGridSize
    bottomMid = halfOfGridSize * gridSize
    #print "corners", bottomLeft, topLeft, topRight, bottomRight
    odPairs = []
    shifted = False
    if(shifted):
        raise Exception("under construction")
#         odPairs.append((bottomLeftShifted, topRightShifted))
#         odPairs.append((topLeftShifted, bottomRightShifted))
#         odPairs.append((bottomLeftShifted, topLeftShifted))
#         odPairs.append((bottomLeftShifted, bottomRightShifted))
#         odPairs.append((topLeftShifted, topRightShifted))
#         odPairs.append((topRightShifted, bottomRightShifted))
    else:
        addBidirectionalCommodityPair(odPairs, bottomLeft, topRight)#8 total
        addBidirectionalCommodityPair(odPairs, topLeft, bottomRight)
        addBidirectionalCommodityPair(odPairs, bottomLeft, topLeft)
        addBidirectionalCommodityPair(odPairs, bottomLeft, bottomRight)
        addBidirectionalCommodityPair(odPairs, topLeft, topRight)
        addBidirectionalCommodityPair(odPairs, topRight, bottomRight)
        addBidirectionalCommodityPair(odPairs, leftMid, rightMid)
        addBidirectionalCommodityPair(odPairs, topMid, bottomMid)
    counter = 0
    print "numODPairs", len(odPairs), "numCommods", numCommodities, odPairs
    for index in range(numCommodities):
        commodities[counter] = {}
        commodities[counter]['odPair'] = (odPairs[index][0], odPairs[index][1])
        commodities[counter]['demand'] = generateCommodityDemand(commodDemandType)
        counter += 1
    print "commodities", commodities
    
def createCommodities_RandomLoc(gridSize, numCommodities, commodDemandType):
    global commodities
    commodities = {}
    odPairs = []
    print "numCommodities", numCommodities
    #print "values", gridSize**2, numCommodities*2
    numCommodPairs = int(numCommodities / 2.0)  
    #print "numCommodPairs", numCommodPairs
    if numCommodPairs == 0:
        sample = np.random.choice(gridSize**2, 2, replace = False)
        print "sample", sample
        odPairs.append((sample[0], sample[1]))
    elif numCommodPairs > 0:
        sample = np.random.choice(gridSize**2, numCommodPairs*2, replace = False)
        origins = sample[0:numCommodPairs]
        destinations = sample[numCommodPairs:numCommodPairs*2]
        for commod in range(numCommodPairs):
            odPairs.append((origins[commod], destinations[commod]))
            odPairs.append((destinations[commod], origins[commod]))
    #print "odPairs", odPairs
    counter = 0
    for index in range(numCommodities):
        #for radioIndex in range(numRadiosPerNode):
        commodities[counter] = {}
        commodities[counter]['odPair'] = (odPairs[index][0], odPairs[index][1])
        commodities[counter]['demand'] = generateCommodityDemand(commodDemandType)
        counter += 1
    #print "commodities", commodities
    
def createDataset(isReal, numNodesToGen, datasetIndex, topology = 'grid', numCommodArg = 1, commodLocType = 'corners', commodDemandTypeArg = 'rand01', name = 'grid'):
    #print "createDataset", topology
    global numNodes
    if not isReal:
        generateGraph(numNodesToGen, topology)
    numNodes = len(G.nodes())
    #print "numNodes", numNodes
    gridSize = int(math.sqrt(numNodes))
    if commodLocType == 'corners':
        createCommodities_Corners(gridSize, numCommodArg, commodDemandTypeArg)
    elif commodLocType == 'rand':
        createCommodities_RandomLoc(gridSize, numCommodArg, commodDemandTypeArg)
    #create XML dataset
    root = ET.Element("data")
    
    nodes = ET.SubElement(root, "nodes")
    nodeCounter = 0
    for i in G.nodes():
        nodeElem = ET.SubElement(nodes, "node")
        ET.SubElement(nodeElem, "name").text = str(i)
        ET.SubElement(nodeElem, "id").text = str(nodeCounter)
        ET.SubElement(nodeElem, "coor").text = str(G.node[i]['coor'])
        ET.SubElement(nodeElem, "distToClosest").text = str(G.node[i]['distToClosest'])
        nodeCounter += 1
        
    edgesElem = ET.SubElement(root, "odPairs")
    for commodIndex in commodities.keys():
        edgeElem = ET.SubElement(edgesElem, "odPair")
        ET.SubElement(edgeElem, "id").text = str(commodIndex)
        ET.SubElement(edgeElem, "origin").text = str(commodities[commodIndex]['odPair'][0])
        ET.SubElement(edgeElem, "destination").text = str(commodities[commodIndex]['odPair'][1])
        ET.SubElement(edgeElem, "demand").text = str(commodities[commodIndex]['demand'])
    
    tree = ET.ElementTree(root)
    if not isReal:
        if topology is 'grid':
            datasetName = "grid-" + str(gridSize) + "x" + str(gridSize)
        elif topology is 'rand':
            "rand-" + str(gridSize**2)
    else:
        datasetName = name
    headString = datasetName + '_index-' + str(datasetIndex) + "_commods-" + str(numCommodArg) \
        + "_commodLoc-" + str(commodLocType) + "_commodDemand-" + str(commodDemandTypeArg)
    exprFile = dataFilePath + '/' + headString + '.xml'
    tree.write(exprFile)

def clearOld(clear = False):
    if(clear):
        bashCommand1 = "for i in /home/hmedal/Documents/2_msu/research_manager/data/wirelessNet/*.xml; do rm -f $i; done"
        os.system(bashCommand1)

if __name__ == "__main__":
    print "START"
    clearOld(True)
    #create deterministic datasets
    gridSizes = [3, 4, 5, 6, 7,8,9]
    commodNums = [1,2,4,8,16]
    demTypes = ['unif1', 'rand02', 'unif5']
    numReps = 2
    createGridDatasets = True
    if createGridDatasets:
        for gridSize in gridSizes:
            for numCommod in commodNums:
                for commodDemandType in demTypes:
                    createDataset(False, int(gridSize**2), 1, topology = 'grid', numCommodArg = numCommod, commodLocType = 'corners', 
                                  commodDemandTypeArg = commodDemandType)
        #create datasets with random commodity demand
        for gridSize in gridSizes:
            for numCommod in commodNums:
                for datasetIndex in range(1, numReps + 1):
                        for commodDemandType in ['rand01']:
                            createDataset(False, int(gridSize**2), datasetIndex, topology = 'grid', numCommodArg = numCommod, commodLocType = 'corners', 
                                          commodDemandTypeArg = commodDemandType)
    #create randomly-generated datasets
    # for gridSize in [3, 5, 6, 9, 12]:
    #     for numCommod in range(1, 17):
    #         for datasetIndex in range(1, 10):
    #                 for commodDemandType in ['rand01']:
    #                     createDataset(False, int(gridSize**2), datasetIndex, topology = 'rand', numCommodArg = numCommod, commodLocType = 'rand', 
    #                                   commodDemandTypeArg = commodDemandType)
    #create MIT and Berkeley datasets
    
    for top in ['berkeley']:
        print "top", top
        if top == 'mit':
            read_and_convert_MIT_net()
        elif top == 'berkeley':
            read_and_convert_Berkeley_net()
        for numCommod in commodNums:
            for commodDemandType in demTypes:
                for datasetIndex in range(1, numReps + 1):
                    #print "createDataset", top, numCommod, datasetIndex
                    createDataset(True, None, datasetIndex, topology = top, numCommodArg = numCommod, commodLocType = 'rand', 
                                  commodDemandTypeArg = commodDemandType, name = top)
    print "FINISHED"