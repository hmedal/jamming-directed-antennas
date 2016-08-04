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
from wnj_projects.wnjJamInf.models import wnj_interference

def runPhysicalWithJamming_Test():
    global interfModelType, g_ConflictGraph
    interfModelType = 'simple-physical'
    gridSize = 3
    source = 0
    dest = gridSize * gridSize - 1
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = nx.max_flow(G, source, dest)
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    #jammingSolution[18] = 1
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    jammingDuals = dict([(node, 1.0) for node in JammingGraph.nodes()])
    count = 0
    maxWeight, selected = wnj_interference.maxWtIndepSet(g_ConflictGraph, capacityDuals, interfModelType)
    ISets.append(selected)
    wnj_interference.createThroughputModel(G, source, dest, [selected], JammingGraph, interfModelType)
    while True:
        print "ITERATION", count
        nodeWeights = wnj_interference.getWeightsForMaxIndSet(G, JammingGraph, capacityDuals, jammingDuals, interfModelType)
        print "node weights", nodeWeights
        maxWeight, selected = wnj_interference.maxWtIndepSet(g_ConflictGraph, nodeWeights, interfModelType)
        ISets.append(selected)
        throughput, capacityDuals, jammingDuals, usageDual = wnj_interference.addISetAsCol_AndSolveTHProb(G, JammingGraph, ISets, selected, count, interfModelType)
        print "selected ", selected
        if((maxWeight - usageDual) <= 0.0001):
            break
        count += 1
        
def testIsEdgeJammedByJammer():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 5
    commodities[0] = (0, gridSize * gridSize - 1)
    commodities[1] = ((gridSize - 1) * gridSize, gridSize - 1)
    commodities[2] = (gridSize, 3*gridSize - 1)
    print commodities
    #source = 0
    #dest = gridSize * gridSize - 1
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = wnj_interference.getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = wnj_interference.createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    for e in G.edges():
        for loc in JammingGraph.nodes():
            print e, loc, wnj_interference.isEdgeJammedByJammer_Protocol(G, e, loc, interfModelType)
    
def runPhysicalTest():
    interfModelType = 'simple-physical'
    gridSize = 3
    source = 0
    dest = gridSize * gridSize - 1
    G = wnj_interference.createConnectivityGraph(gridSize, interfModelType)
    wiredNetworkFlow = nx.max_flow(G, source, dest)
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    ConflictGraph = wnj_interference.createConflictGraph(G, interfModelType)
    #print ConflictGraph.edges(data=True)
    
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    count = 0
    maxWeight, selected = wnj_interference.maxWtIndepSet(ConflictGraph, capacityDuals, interfModelType)
    ISets.append(selected)
    #print "selected ", selected
    interdictionFeasibleRegion = [(i+0.5, j+0.5) for i in range(gridSize - 1) for j in range(gridSize - 1)]
    interdictionVars = dict([(tuple, 1.0) for tuple in interdictionFeasibleRegion])
    wnj_interference.createThroughputModel(G, source, dest, [selected], interdictionVars)
    while True:
        print "ITERATION", count
        maxWeight, selected = wnj_interference.maxWtIndepSet(ConflictGraph, capacityDuals, interfModelType)
        ISets.append(selected)
        print "selected ", selected
        throughput, capacityDuals, jammingDuals, usageDual = wnj_interference.addISetAsCol_AndSolveTHProb(G, interdictionVars, ISets, selected, count)
        if((maxWeight - (usageDual + sum([jammingDuals[key] for key in jammingDuals.keys()]))) <= 0.0001):
            break
        count += 1
        
def runProtocolTest():
    interfModelType = '802.11-protocol'
    gridSize = 3
    source = 0
    dest = gridSize * gridSize - 1
    JainExample(gridSize, source, dest, interfModelType)
    
def runProtocolWithJamming_Test():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 3
    commodities[0] = (0, gridSize * gridSize - 1)
    #source = 0
    #dest = gridSize * gridSize - 1
    G = createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    jammingSolution[15] = 1
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    solveThroughputProblem_Pricing(G, JammingGraph)
    print "FINISHED"
    
def runProtocolWithJamming_FullMIP_Test():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 3
    commodities[0] = (0, gridSize * gridSize - 1)
    #source = 0
    #dest = gridSize * gridSize - 1
    G = createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    #jammingDuals = dict([(node, dict([(edge, 0.0) for edge in G.edges()])) for node in JammingGraph.nodes()])
    jammingDuals = dict([(edge, 0.0) for edge in G.edges()])
    maxWeight, selected = maxWtIndepSet(g_ConflictGraph, capacityDuals, interfModelType)
    ISets.append(selected)
    ISets.append([(3, 6)])
    ISets.append([(6, 7)])
    ISets.append([(0, 3), (7, 8)])
    ISets.append([(0, 1), (5, 8)])
    ISets.append([(1, 2), (7, 8)])
    ISets.append([(2, 5), (3, 6)])
    ISets.append([(1, 2), (6, 7)])
    ISets.append([(2, 5), (0, 3)])
    
    #createThroughputModel(G, commodities, ISets, JammingGraph, interfModelType)
    create_SingleLevelMIP_Cormican_ArcBased(G, commodities, ISets, JammingGraph, interfModelType)
    print "FINISHED"

def runProtocolWithJamming_FullMIP_RowGen_Test():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 4
    commodities[0] = (0, gridSize * gridSize - 1)
    commodities[1] = ((gridSize - 1) * gridSize, gridSize - 1)
    commodities[2] = (gridSize, 3*gridSize - 1)
    #commodities[3] = (1, 5*gridSize+1)
    #commodities[4] = (23, 33)
    print "commodities", commodities
    #source = 0
    #dest = gridSize * gridSize - 1
    G = createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    nodeWeights = capacityDuals
    #jammingDuals = dict([(node, dict([(edge, 0.0) for edge in G.edges()])) for node in JammingGraph.nodes()])
    #jammingDuals = dict([(edge, 0.0) for edge in G.edges()])
    maxWeight, selected = maxWtIndepSet(g_ConflictGraph, nodeWeights, interfModelType)
    ISets.append(selected)
    create_SingleLevelMIP_Cormican_ArcBased(G, commodities, ISets, JammingGraph, interfModelType)
    count = 0
    while True:
        print "ITERATION", count
        print "nodeWeights", nodeWeights
        maxWeight, selected = maxWtIndepSet(g_ConflictGraph, nodeWeights, interfModelType)
        ISets.append(selected)
        print "selected ", selected
        
        throughput, capacityDuals, usageDual, interdicted = addISetAsRow_AndSolveFullMIP(G, selected, count, interfModelType)
        nodeWeights = dict([(e, capacityDuals[e] * G.edge[e[0]][e[1]]['capacity']) for e in G.edges()])
        if((maxWeight - usageDual) <= 0.0001):
            break
        count += 1
    print "FINISHED"
    
    
def test_InterferenceGraph(interference, gridSize):
    global interfModelType, g_ConflictGraph, commodities, Paths, ISets, maxNumHops
    interfModelType = '802.11-MAC-protocol'
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
    #G = createConnectivityGraph(gridSize, interfModelType, 'grid', False, multipleRadiosPerNode, numChannels)
    #print "G edges", G.edges(data = True)
    #g_ConflictGraph = createConflictGraph(G, interfModelType)
    
def test_JamGraph(interference, gridSize):
    global interfModelType, g_ConflictGraph, commodities, Paths, ISets, maxNumHops
    interfModelType = '802.11-MAC-protocol'
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
    #G = createConnectivityGraph(gridSize, interfModelType, 'grid', False, multipleRadiosPerNode, numChannels)
    #print "G edges", G.edges(data = True)
    #g_ConflictGraph = createConflictGraph(G, interfModelType)
    print "g_ConflictGraph edges", g_ConflictGraph.edges()
    JammingGraph = createJamGraph(gridSize)
    
    
def test_MultiGraph(interference, gridSize, numRadiosPerNode, numChannels):
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
    G = createConnectivityGraph(gridSize, interfModelType, 'grid', False, numRadiosPerNode, numChannels)
    if(interference):
        nx.draw(G, pos, edgelist = [])
    else:
        nx.draw(G, pos)
    #print G.nodes()
    nx.draw(G, pos, edgelist = [], nodelist = [0, 15], node_color = '0.9')
    if(interference):
        nx.draw_networkx_nodes(G, pos, alpha = 0.1, node_size = 50000, node_color='purple', node_shape='o', nodelist = G.nodes(), with_labels=False)
    plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/wiredSoln.jpg") # save as jpg
    plt.show() # display
    
def runProtocolWithJamming_Cormican_PathBased_SingleMIP_Test():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 6
    bottomLeft = 0
    topLeft = gridSize - 1
    topRight = gridSize * gridSize - 1
    bottomRight = gridSize * (gridSize - 1)
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    commodities[0] = (bottomLeft, topRight)
    commodities[1] = (bottomLeft, topLeft)
    commodities[2] = (bottomLeft, bottomRight)
    commodities[3] = (topLeft, topRight)
    commodities[4] = (topLeft, bottomRight)
    commodities[5] = (topRight, bottomRight)
    G = createConnectivityGraph(gridSize, interfModelType, 'grid')
    wiredNetworkFlow = getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = createConflictGraph(G, interfModelType)

    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    ISets = []
    max_hop_length = gridSize * 3
    Paths = {}
    solve_FullMIP_DelayedRowGen_Cormican_PathBased(G, JammingGraph, Paths, ISets, max_hop_length)
    print "FINISHED"
    
def runProtocolWithJamming_FullMIP_RowGenCallback_Test():
    global interfModelType, g_ConflictGraph, commodities
    interfModelType = '802.11-MAC-protocol'
    gridSize = 3
    commodities[0] = (0, gridSize * gridSize - 1)
    #commodities[1] = ((gridSize - 1) * gridSize, gridSize - 1)
    #commodities[2] = (gridSize, 3*gridSize - 1)
    #commodities[3] = (1, 5*gridSize+1)
    #commodities[4] = (23, 33)
    print "commodities", commodities
    #source = 0
    #dest = gridSize * gridSize - 1
    G = createConnectivityGraph(gridSize, interfModelType, 'grid')
    print G.edges(data=True)
    wiredNetworkFlow = getNetworkMaxTotalThroughput(G, commodities, 'none')
    print "wiredNetworkFlow:", wiredNetworkFlow
    
    g_ConflictGraph = createConflictGraph(G, interfModelType)
    print "CG edges", g_ConflictGraph.edges(data=True)
    
    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2)**2
    print "jammingSolution", jammingSolution
    
    counter = 0
    for node in JammingGraph.nodes():
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    nodeWeights = capacityDuals
    #jammingDuals = dict([(node, dict([(edge, 0.0) for edge in G.edges()])) for node in JammingGraph.nodes()])
    #jammingDuals = dict([(edge, 0.0) for edge in G.edges()])
    maxWeight, selected = maxWtIndepSet(g_ConflictGraph, nodeWeights, interfModelType)
    ISets.append(selected)
    throughput = solveFullMIP_DelayedRowGen_withCallbacks(G, JammingGraph, ISets, interfModelType)
    print "throughput", throughput
    print "FINISHED"