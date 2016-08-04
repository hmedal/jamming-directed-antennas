
import networkx as nx
import numpy as np
from graph_tool.all import *
from random import sample
import random

colorsList = ["red", "blue", "green"]
dashStylesList = [[0.1,0.1,0], [0.05,0.05,0], [0.1,0.1,0]]

def getArcsInPath(path):
    arcsSet = []
    #print "pathList", path
    for index in range(len(path)-1):
        arc = (path[index], path[index+1])
        #print "arc", arc
        arcsSet.append(arc)
    return arcsSet

def getResidualCapOfPath(G, path):
    arcsOnPath = getArcsInPath(path)
    minResid = float('inf')
    for arc in arcsOnPath:
        arcCap = G[arc[0]][arc[1]]['capacity']
        arcFlow = G[arc[0]][arc[1]]['flow']
        resid = arcCap - arcFlow
        if minResid > resid:
            minResid = resid
    return minResid

def augment(G, path, amount):
    arcsOnPath = getArcsInPath(path)
    for arc in arcsOnPath:
        G[arc[0]][arc[1]]['flow'] = G[arc[0]][arc[1]]['flow'] + amount
    
def generateNetworkXGraph(gridDim1 = 2, gridDim2 = 2, type = 'none', addCrossArcsArg = True, addSrcAndSink = True, arcsToReverse = [], 
                          generateRandomCapacities = True, flows = '0'):
    G=nx.grid_2d_graph(gridDim1, gridDim2)
    if addCrossArcsArg is True:
        addCrossArcs(G, gridDim1, gridDim2)
    G=nx.convert.convert_to_directed(G)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted", label_attribute='coor')
    print "edges", G.edges()
    if type is 'acyclic':
        makeAcyclic(G)
    elif type is 'randdir':
        makeRandDir(G)
    elif type is 'none':
        None
    else:
        raise Exception("invalid type")
    if addSrcAndSink:
        addSourceAndSink(G, gridDim1, gridDim2)
    print "nodes", G.nodes()
    print "edges", G.edges()
    for edgePair in arcsToReverse:
        G.remove_edge(edgePair[0], edgePair[1])
        G.add_edge(edgePair[1], edgePair[0])
    for e in G.edges():
        randCap = np.random.randint(1, 10)
        randWeight = np.random.randint(1, 10)
        G[e[0]][e[1]]['weight'] = randWeight
        G[e[0]][e[1]]['capacity'] = randCap
        G[e[0]][e[1]]['flow'] = 0
    if flows is 'rand':
        allPaths = list(nx.all_simple_paths(G, 's', 't'))
        for pathIndex in [0,1,2]:
            path = allPaths[pathIndex]
            print "path", getArcsInPath(path)
            resid = getResidualCapOfPath(G, path)
            augmentAmt = np.random.randint(0, resid + 1)
            augment(G, path, augmentAmt)
    for e in G.edges():
        print "edge from nx", e, G[e[0]][e[1]]['flow'], G[e[0]][e[1]]['capacity']
    #filename = "/home/hmedal/Documents/2_msu/teacher/current_courses/net_flows/course_materials/course_pack/img/networkXDot.dot"
    #nx.write_dot(G, filename)
    return G

def addSourceAndSink(G, gridDim1, gridDim2):
    G.add_node('s', coor = (-1, 0.5 * (gridDim2 - 1)))
    G.add_node('t', coor = (gridDim1, 0.5 * (gridDim2 - 1)))
    for head in range(gridDim2):
        print "head", head
        G.add_edge('s', head)
    for tail in range(gridDim2 * (gridDim1 - 1), gridDim1*gridDim2):
        G.add_edge(tail, 't')
        
def addCrossArcs(G, gridDim1, gridDim2):
    for i in range(gridDim1 - 1):
            for j in range(gridDim2 - 1):
                G.add_edge((i, j), (i + 1, j + 1))
    for i in range(1, gridDim1):
        for j in range(gridDim2 - 1):
            G.add_edge((i, j), (i - 1, j + 1))
                
def makeAcyclic(G):
    for edgePair in G.edges():
        tail = edgePair[0]
        head = edgePair[1]
        if head < tail:
            G.remove_edge(tail, head)
            
def makeRandDir(G):
    toRemove = []
    for edgePair in G.edges():
        tail = edgePair[0]
        head = edgePair[1]
        print tail, head
        keepDir = np.random.choice([True, False])
        if keepDir:
            toRemove.append((head, tail))
        else:
            toRemove.append((tail, head))
    for edgePair in toRemove:
        G.remove_edge(edgePair[0], edgePair[1])

def drawGraphToolNetwork(networkXGraph, gridDim1 = 3, gridDim2 = 2, randLabels=True, includeFlows = False, includeEdges = True, includeEdgeLabels = True, 
                         suffix = '', arcsToColorMap = {}):
    #G= generateNetworkXGraph(gridDim1, gridDim2, type = 'acyclic', arcsToReverse = []) #sweep
    #G= generateNetworkXGraph(gridDim1, gridDim2, type = 'acyclic', arcsToReverse = [(2,5), (1, 4)]) # dijkstra
    G= networkXGraph
    numNodes = gridDim1*gridDim2+2
    g = Graph()
    vertexLabels = g.new_vertex_property("string")
    edgeLabels = g.new_edge_property("string")
    edgeColors = g.new_edge_property("string")
    edgeDashStyles = g.new_edge_property("vector<float>")
    vertexPos = g.new_vertex_property("vector<double>")
    g.vertex_properties["label"] = vertexLabels
    g.edge_properties["label"] = edgeLabels
    verticesMap = {}
    if randLabels is True:
        labelsArray = sample(xrange(numNodes), numNodes)
    index = 0
    for n in G.nodes(data = True):
        print "n", n
        nodeIndex = n[0]
        #print "nodeIndex", nodeIndex
        #print "G[n]", G.node[nodeIndex]
        myPos = G.node[nodeIndex]['coor']
        if randLabels is False:
            if isinstance(nodeIndex, int):
                myLabel = nodeIndex + 1
            else:
                myLabel = nodeIndex
        else:
            myLabel = labelsArray[index]
        #print "myLabel", nodeIndex, myLabel
        v = g.add_vertex()
        vertexLabels[v] = myLabel
        verticesMap[nodeIndex] = v
        vertexPos[v] = [myPos[0], myPos[1]]
        #print "vertexPos[v]", v, vertexPos[v]
        index += 1
    print "verticesMap", verticesMap
    if includeEdges:
        for edgePair in G.edges():
            print "edge", edgePair, G[edgePair[0]][edgePair[1]]['capacity']
            flow = G[edgePair[0]][edgePair[1]]['flow']
            cap = G[edgePair[0]][edgePair[1]]['capacity']
            e = g.add_edge(verticesMap[edgePair[0]], verticesMap[edgePair[1]])
            print "e", e
            flowStr = ""
            if includeFlows:
                flowStr = str(flow)
            if includeEdgeLabels:
                edgeLabels[e] = flowStr + "/" + str(cap)
            else:
                edgeLabels[e] = ""
            if edgePair in arcsToColorMap:
                edgeColors[e] = arcsToColorMap[edgePair]['color']
                edgeDashStyles[e] = arcsToColorMap[edgePair]['dash_style']
            else:
                edgeColors[e] = "gray"
                edgeDashStyles[e] = []
    filename = "~/Documents/2_msu/1_MSU_Projects/PAPER_JammingSpatialInterference/figures/wiredVsWireless" + suffix + ".pdf"
    saveFig = True
    if saveFig:
        outputMode = filename
    else:
        outputMode = None
    graph_draw(g, pos=vertexPos,  
               vertex_text=vertexLabels, vertex_color = 'black', vertex_fill_color = 'white',vertex_font_size=24, 
               edge_font_size=18, edge_pen_width=3, edge_color = edgeColors, edge_text = edgeLabels, edge_marker_size = 20, edge_end_marker = 'arrow', edge_dash_style = edgeDashStyles,
               output = outputMode)
    #g.save("~/Documents/2_msu/teacher/current_courses/net_flows/course_materials/course_pack/img/maxFlow_residual.dot")

if __name__=='__main__':
    #g = price_network(1500)
    #deg = g.degree_property_map("in")
    #print "deg", deg
    #print "start"
    gridDim1 = 4
    gridDim2 = 4
    myType = 'none'
    addCrossArcsArg = False
    #pathToDraw = {}
    #pathToDraw[1] = {'path' : [(0,1), (1,2), (2,3), (3, 7), (7,11), (11,15)], 'color' : 'red', 'dash_style' : [0.05,0.05,0]}
    #pathToDraw[2] = {'path' : [(0,4), (4,8), (8,12), (12,13), (13,14), (14,15)], 'color' : 'red', 'dash_style' : [0.05,0.05,0]}
    arcsToColorMap = {}
    for edge in [(2,3), (11,15), (0,4), (12,13)]:
        arcsToColorMap[edge] = {'color' : 'green', 'dash_style' : [0.02,0.02,0]}
    for edge in [(1,2), (4,8), (7,11), (13,14)]:
        arcsToColorMap[edge] = {'color' : 'red', 'dash_style' : [0.05,0.05,0]}
    for edge in [(0,1), (3,7), (8,12), (14,15)]:
        arcsToColorMap[edge] = {'color' : 'blue', 'dash_style' : [0.15,0.02,0]}
    networkXGraph = generateNetworkXGraph(gridDim1, gridDim2, type = myType, addCrossArcsArg = addCrossArcsArg, addSrcAndSink = False, arcsToReverse = []) # F-W
    drawGraphToolNetwork(networkXGraph, gridDim1 = gridDim1, gridDim2 = gridDim2, randLabels=False, includeFlows = True, includeEdges = True, includeEdgeLabels = False, suffix = '_wireless', arcsToColorMap = arcsToColorMap)
    #drawGraphToolNetwork(networkXGraph, gridDim1 = gridDim1, gridDim2 = gridDim2, randLabels=False, includeFlows = True, includeEdges = False, includeEdgeLabels = False, suffix = '_wireless')
    #drawGraphToolNetwork(networkXGraph, gridDim1 = gridDim1, gridDim2 = gridDim2, randLabels=False, includeFlows = False, includeEdges = True, suffix = '_noFlow')
    #drawGraphToolNetwork(networkXGraph, gridDim1 = gridDim1, gridDim2 = gridDim2, randLabels=False, includeFlows = False, includeEdges = False, suffix = '_noEdges')