try:
    import matplotlib.pyplot as plt
except:
    raise

import networkx as nx
import numpy as np
import gurobipy
import pylab
import time
import lxml.etree as etree
from scipy import stats
import itertools
import logging
from wnj_core.objects import executableModel, wnj_problemInstance
from wnj_core.data import wnj_dataset
import wnj_core.output.dbUtil as dbUtil
import wnj_core.myUtil as myutil
from wnj_projects.wnjJamInf.models import wnj_model
from wnj_projects.wnjJamInf.models import wnj_interference

def test_showGraph(interference, gridSize):
    global interfModelType, g_ConflictGraph, commodities, Paths, ISets, maxNumHops
    interfModelType = 'none'
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    #commodities[1] = (bottomLeft, topLeft)
    #commodities[2] = (bottomLeft, bottomRight)
    #commodities[3] = (topLeft, topRight)
    #commodities[4] = (topLeft, bottomRight)
    #commodities[5] = (topRight, bottomRight)
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid', False, 
                                                 wnj_interference.numRadiosPerNode, wnj_interference.numChannels)
    if(interference):
        nx.draw(G, wnj_interference.pos, edgelist = [])
    else:
        nx.draw(G, wnj_interference.pos)
    #print G.nodes()
    nx.draw(G, wnj_interference.pos, edgelist = [], nodelist = [0, 15], node_color = '0.9')
    if(interference):
        nx.draw_networkx_nodes(G, wnj_interference.pos, alpha = 0.1, node_size = 50000, node_color='purple', node_shape='o', nodelist = G.nodes(), with_labels=False)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln.jpg") # save as jpg
    plt.show() # display
    
def test_showWiredSolution_withJamming():
    print "runProtocolWithJamming_FullMIP_RowGenCallback_PathBased_Test"
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 5
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    #commodities[1] = (bottomLeft, topLeft)
    #commodities[2] = (bottomLeft, bottomRight)
    #commodities[3] = (topLeft, topRight)
    #commodities[4] = (topLeft, bottomRight)
    #commodities[5] = (topRight, bottomRight)
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid', False)
    print G.edges(data=True)
    
    
    g_ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    jammingSolution[39] = 1
    print "jammingSolution", jammingSolution
    jammedNodes = []
    jammedPos = {}
    counter = 0
    for node in JammingGraph.nodes():
        print node, counter
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        if(jammingSolution[counter] == 1.0):
            jammedNodes.append(node)
            jammedPos[node] = node
        counter += 1
    print "jammedNodes", jammedNodes
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    g_JamGraph = JammingGraph
    
    throughput, flowSoln = wnj_interference.getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", throughput
    print "wiredFlows", flowSoln
    flowLinks = []
    for key in flowSoln.keys():
        for edge in flowSoln[key].keys():
            flowLinks.append(edge)
    print "flowLinks", flowLinks
    nx.draw(G, wnj_interference.pos)
    nx.draw_networkx_nodes(G, jammedPos, alpha = 0.1, node_size = 25000, node_color='blue', node_shape='o', nodelist = jammedNodes, with_labels=False)
    nx.draw_networkx_edges(G, wnj_interference.pos, edge_color='red', width=3.0, edgelist = flowLinks)
    nx.draw_networkx_nodes(G, jammedPos, node_size = 100, node_color='green', node_shape='s', nodelist = jammedNodes, with_labels=False)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln_withJamming.jpg") # save as jpg
    plt.show() # display

def test_showWiredSolution_revised(interdictedNodes):
    global interfModelType, g_ConflictGraph, commodities, Paths, ISets, maxNumHops
    interfModelType = 'none'
    gridSize = 4
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    #commodities[1] = (bottomLeft, topLeft)
    #commodities[2] = (bottomLeft, bottomRight)
    #commodities[3] = (topLeft, topRight)
    #commodities[4] = (topLeft, bottomRight)
    #commodities[5] = (topRight, bottomRight)
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid', False)
    print "edges", G.edges(data=True)
    
    
    #g_ConflictGraph = createConflictGraph(G, interfModelType)
    #print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    for node in interdictedNodes:
        jammingSolution[node] = 1
    print "jammingSolution", jammingSolution
    
    jammedNodes = []
    jammedPos = {}
    counter = 0
    for node in JammingGraph.nodes():
        print node, counter
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        if(jammingSolution[counter] == 1.0):
            jammedNodes.append(node)
            jammedPos[node] = node
        counter += 1
    print "jammedNodes", jammedNodes
    
    ISets = []
    maxNumHops = gridSize * 2
    
    PathsList = nx.all_simple_paths(G, 0, gridSize * gridSize - 1)
    Paths = {}
    for commod in commodities.keys():
        Paths[commod] = []
    edgesJammed = wnj_interference.getEdgesJammedByJammers(G, jammedNodes, interfModelType)
    print "edgesJammed", edgesJammed
    flowVarSoln, isetUsageSoln = wnj_interference.solveThroughputProblem_Pricing_CormicanPathBased(G, edgesJammed, maxNumHops, 'none', JammingGraph)
    print "flowVarSoln", flowVarSoln

    pathCtr = 0
    pathColors = ['red', 'green', 'blue', 'orange', 'yellow', 'purple', 'black']
    pathStyle = ['solid', 'dotted', 'dashdot', 'dashed']
    #nx.draw_networkx_nodes(G, pos, nodelist = G.nodes(), with_labels=True)
    nx.draw(G, wnj_interference.pos)
    nx.draw(G, wnj_interference.pos, edgelist = [], nodelist = [0, 15], node_color = '0.9')
    for commod in flowVarSoln.keys():
        print "numPaths", len(Paths[commod])
        for path in flowVarSoln[commod].keys():
            flowAmt = flowVarSoln[commod][path]
            if(flowVarSoln[commod][path] > wnj_interference.FUZZ):
                print "flowAmt", flowAmt, Paths[commod][path][1]
                for edge in Paths[commod][path][1]:
                    flowLinks = []
                    flowLinks.append(edge)
                    nx.draw_networkx_edges(G, wnj_interference.pos,edge_color = pathColors[0], width= flowAmt * 10.0, style = pathStyle[pathCtr], edgelist=flowLinks)
                print "pathCtr", pathCtr
                pathCtr += 1
    print "flowLinks", flowLinks
    nx.draw_networkx_nodes(G, jammedPos, alpha = 0.1, node_size = 90000, node_color='blue', node_shape='o', nodelist = jammedNodes, with_labels=False)
    nx.draw_networkx_nodes(G, jammedPos, node_size = 100, node_color='black', node_shape='s', nodelist = jammedNodes, with_labels=False)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln.jpg") # save as jpg
    plt.show() # display

def test_showWiredSolution():
    print "runProtocolWithJamming_FullMIP_RowGenCallback_PathBased_Test"
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 5
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    #commodities[1] = (bottomLeft, topLeft)
    #commodities[2] = (bottomLeft, bottomRight)
    #commodities[3] = (topLeft, topRight)
    #commodities[4] = (topLeft, bottomRight)
    #commodities[5] = (topRight, bottomRight)
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid', False)
    print G.edges(data=True)
    
    
    g_ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    jammingSolution[15] = 1
    print "jammingSolution", jammingSolution
    
    throughput, flowSoln = wnj_interference.getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", throughput
    print "wiredFlows", flowSoln
    flowLinks = []
    for key in flowSoln.keys():
        for edge in flowSoln[key].keys():
            flowLinks.append(edge)
    nx.draw(G, wnj_interference.pos)
    nx.draw_networkx_edges(G, wnj_interference.pos,edge_color='red', width=3.0, edgelist=flowLinks)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln.jpg") # save as jpg
    plt.show() # display

def test_showInterferenceSolution(interdictedNodes):
    global interfModelType, g_ConflictGraph, commodities, Paths, ISets, maxNumHops
    interfModelType = '802.11-MAC-protocol'
    gridSize = 4
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    #commodities[1] = (bottomLeft, topLeft)
    #commodities[2] = (bottomLeft, bottomRight)
    #commodities[3] = (topLeft, topRight)
    #commodities[4] = (topLeft, bottomRight)
    #commodities[5] = (topRight, bottomRight)
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid', False)
    print G.edges(data=True)
    
    
    g_ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    for node in interdictedNodes:
        jammingSolution[node] = 1
    print "jammingSolution", jammingSolution
    
    jammedNodes = []
    jammedPos = {}
    counter = 0
    for node in JammingGraph.nodes():
        print node, counter
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        if(jammingSolution[counter] == 1.0):
            jammedNodes.append(node)
            jammedPos[node] = node
        counter += 1
    print "jammedNodes", jammedNodes
    
    ISets = []
    maxNumHops = gridSize * 2
    
    PathsList = nx.all_simple_paths(G, 0, gridSize * gridSize - 1)
    Paths = {}
    for commod in commodities.keys():
        Paths[commod] = []
    edgesJammed = wnj_interference.getEdgesJammedByJammers(G, jammedNodes, interfModelType)
    print "edgesJammed", edgesJammed
    
    flowVarSoln, isetUsageSoln = wnj_interference.solveThroughputProblem_Pricing_CormicanPathBased(G, edgesJammed, maxNumHops, interfModelType, JammingGraph)
    print "flowVarSoln", flowVarSoln
    print "isetUsageSoln", isetUsageSoln
    ISetsUsed = [ISets[k] for k in range(len(ISets)) if isetUsageSoln[k] > wnj_interference.FUZZ]
    
    #throughput, flowSoln = getNetworkMaxTotalThroughput(G, commodities, interfModelType)
    #print "wiredNetworkFlow:", throughput
    print "ISetsUsed", ISetsUsed
    pathCtr = 0
    pathColors = ['red', 'green', 'blue', 'orange', 'yellow', 'purple']
    pathStyle = ['solid', 'dotted', 'dashdot', 'dashed']
    #nx.draw_networkx(G, pos, style = 'dashed')
    nx.draw(G, wnj_interference.pos, edgelist = [])
    nx.draw(G, wnj_interference.pos, edgelist = [], nodelist = [0, 15], node_color = '0.9')
    for commod in flowVarSoln.keys():
        print "numPaths", len(Paths[commod])
        for path in flowVarSoln[commod].keys():
            flowAmt = flowVarSoln[commod][path]
            if(flowVarSoln[commod][path] > wnj_interference.FUZZ):
                print "flowAmt", flowAmt
                for edge in Paths[commod][path][1]:
                    flowLinks = []
                    flowLinks.append(edge)
                    isetIndex = wnj_interference.getISetIndexForEdge(edge, ISetsUsed)
                    print edge, "isetIndex", isetIndex
                    nx.draw_networkx_edges(G, wnj_interference.pos,edge_color = pathColors[isetIndex], width= flowAmt * 25.0, style = pathStyle[pathCtr], edgelist=flowLinks)
                #print "pathCtr", pathCtr
                pathCtr += 1
    #print "flowLinks", flowLinks
    nx.draw_networkx_nodes(G, jammedPos, alpha = 0.1, node_size = 90000, node_color='blue', node_shape='o', nodelist = jammedNodes, with_labels=False)
    nx.draw_networkx_nodes(G, jammedPos, node_size = 100, node_color='black', node_shape='s', nodelist = jammedNodes, with_labels=False)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln.jpg") # save as jpg
    plt.show() # display
    
def showStuffForPaper():
    test_showGraph(True)
    #test_showWiredSolution_revised([])
    #test_showWiredSolution_revised([11])
    #test_showInterferenceSolution([])
    test_showInterferenceSolution([11])