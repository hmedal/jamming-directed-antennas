from collections import OrderedDict
import itertools
import logging
import math
import sys
import time
import operator

import gurobipy
import pylab
from scipy import stats

import lxml.etree as etree
import networkx as n
import numpy as np
from wnj_core.dat import wnj_dataset
import wnj_core.myIO.databaseUtil as dbUtil
from wnj_core.objects import executableModel
from wnj_proj.wnj_jamInf.mod import wnj_model


try:
    import matplotlib.pyplot as pltg
except:
    raise

try:
    from wnj_core.objects import wnj_problemInstance
except:
    raise Exception("cannot import from src.wnj_core.objects import wnj_problemInstance")

FUZZ = 0.000001

#models
maxWtIndSetModel = None
gurobiModel = None
masterProbLP = None
gurobiThroughputModel = None
shortestPathModels = {}

pos = None
#commodities = {}
pathNumForCommod = {}
interfModelType = None

#constraints
flowBalanceConstrs = None
capConstraints = None
usageConstr = None
commodDemConstr = None
jamConstrs = None
throughputEqualsValueConstr = None
trustRegionConstr = None
slackUBConstr = None

#variables
edgeInterdict = None
jamPlaced = None
jamPlacedLP = None
flowVars = None
inSet = None
iSetUsage = None
findPathflowVars = None
returnFlow = None
slackVars = None

#dual variables
capDuals = None
usageDual = None
flowVarsAsDual = None
isetUsageVarsAsDual = None
demForCommodDuals = None
g_capacityDuals = None # for throughput problem

# parameters
slackVariableMultipliers = None
slackVarUBs = None

#sets
Paths = {}
ISets = []
edgeTriples = []
edgeTriplesInterdictable = []
jammersThatCanJamEdge = {}
edgesToIndexInterdictDict = {}
nodesToIndexInterdictDict = {}
edgesToIndexDict = {}
commodToIndexDict = {}

#constants
maxNumHops = 0
maxNumHopsMult = 0
snr_threshold = 1.5
ambient_noise= 0.01
interdictBudget = 0

#flags
writeToFile = False
debug = False

previous_mip_gap = 0.0
bestBoundGlobal = 0.0
bestObjGlobal = 0.0

#counts
multipleRadiosPerNode = False
numRadiosPerNode = 1
numChannels = 1

#counters
numMaxWtIndSetProbsSolved = 0
callbackCtr = 0
numTypeICuts = 0
numNodesExplored = 0

#limits
numTypeICutsLimit = float('inf')
trust_region_max_iter = 10

#paths
dataFilePath = None

#algorithm params
time_limit = 0
numThreads = 1

#strings
modelType = None
jamVarsType = None
flowVarsType = None

#other
inst = None
clusterName = None
bestMIP_upperBound = 0
smallestSubprobThroughput = float('inf')
sviTypeI_Cuts = []
sviTypeI_NodesInLHS = []
branchAndBoundVariations = ['branch-and-cut', 'b-and-c-reg']
bendersVariations = ['benders', 'benders-tr', 'benders-ki', 'benders-pareto', 'benders-svi', 'benders-pareto-tr', 'benders-ki-pareto', 
                     'benders-ki-pareto-tr', 'benders-ki-pareto-svi', 'benders-pareto-svi-tr', 'benders-ki-pareto-svi-tr',
                     'benders-callback', 'benders-callback-ki', 'benders-callback-pareto', 'benders-callback-ki-pareto']

includeKnapsack = False
includeSuperValid = False
includeParetoOptimalBendersCuts = False
includeTrustRegion = False

def getPathFromEdgeList(edgeList, origin, destination):
    path = []
    currentNode = origin
    path.append(origin)
    remainingEdges = edgeList
    edgeToRemove = None
    while(True):
        #print "remainingEdges", remainingEdges
        for edge in remainingEdges:
            if(edge[0] == currentNode):
                currentNode = edge[1]
                path.append(currentNode)
                edgeToRemove = edge
                break
        remainingEdges.remove(edgeToRemove)
        if(currentNode == destination):
            break
    return path
    
def getEdgesJammedByJammers(G, jammedNodes, interfModelType):
    edgesJammed = []
    for edge in edgeTriplesInterdictable:
        for jammedNode in jammedNodes:
            if(isEdgeJammedByJammer_Protocol(G, edge, jammedNode, interfModelType)):
                edgesJammed.append(edge)
                break
    return edgesJammed        
            
def isEdgeJammedByJammer_Protocol(G, edgeTriple, jamLocWithRadioType, interfModelType):
    #print "jamLoc", jamLocWithRadioType, "edgeTriple", edgeTriple, "infRange", instance.jamGraph.node[jamLocWithRadioType]['interferenceRange']
    #print g_JamGraph.nodes(data=True)
    #print g_JamGraph.node[jamLoc]
    jamLoc = (jamLocWithRadioType[0], jamLocWithRadioType[1])
    jamRadioType = jamLocWithRadioType[2]
    edge = (edgeTriple[0], edgeTriple[1])
    edgeChannel = edgeTriple[2]
    #jamLoc = (jamLocWithRadioType[0], jamLocWithRadioType[1])
    if(jamRadioType == edgeChannel):
        distToReceiver = np.linalg.norm(np.array(jamLoc) - np.array(instance.CGraph.node[edge[1]]['coor']))
        #print "distToReceiver", distToReceiver
        if(distToReceiver <= instance.jamGraph.node[jamLocWithRadioType]['interferenceRange']):
            return True
        distToSender = np.linalg.norm(np.array(jamLoc) - np.array(instance.CGraph.node[edge[0]]['coor']))
        #print "distToSender", distToSender
        if(distToSender <= instance.jamGraph.node[jamLocWithRadioType]['interferenceRange']):
            return True
    else:
        return False
    return False

def isEdgeJammedByJammers(G, edge, jamLocs, interfModelType):
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        for jamLoc in jamLocs:
            #print "edge", edge, "loc", jamLoc, "canBeJammed", isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)
            if isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType):
                #print "edge ", edge, " is jammed by  ", jamLoc
                return True
        return False
    
def isEdgeJammedByJammers_JamGraphSelectedValue(G, JammingGraph, edge, jamLocs, interfModelType):
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        for jamLoc in jamLocs:
            #print "edge", edge, "loc", jamLoc, "canBeJammed", isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)
            if(isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType) & (JammingGraph.node[jamLoc]['selected'] >= (1-FUZZ))):
                #print "edge ", edge, " is jammed by  ", jamLoc
                return True
        return False
    
def canEdgeBeJammedByJammers(G, edge, jamLocs, interfModelType):
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        for jamLoc in jamLocs:
            #print "edge", edge, "loc", jamLoc, "canBeJammed", isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)
            if(isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)):
                return True
        return False
    
def isEdgeSetJammedByJammer(G, edges, jamLoc, interfModelType):
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        for edge in edges:
            if(isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)):
                return True
        return False
    
def getJammersThatCanJamEdge(G, JammingGraph, edge, interfModelType):
    jammersList = []
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        for jamLoc in JammingGraph.nodes():
            #print "isJammed", edge, jamLoc, isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)
            if(isEdgeJammedByJammer_Protocol(G, edge, jamLoc, interfModelType)):
                jammersList.append(jamLoc)
    print "jammerlist", jammersList
    
    return jammersList

def createJammingMasterProblem(G, JammingGraph, interfModelType):
    global jamPlaced, edgeInterdict, approxVar, g_JamGraph
    g_JamGraph= JammingGraph
    jamMasterProb = gurobipy.Model("Jamming Master Problem")
    try:
        # Create variables
        jamPlaced = dict([(n, jamMasterProb.addVar(0, vtype=gurobipy.GRB.BINARY, name="jamVar_"+str(n))) for n in JammingGraph.nodes()])
        approxVar = jamMasterProb.addVar(0, vtype=gurobipy.GRB.CONTINUOUS, name="approxVar")
        jamMasterProb.update() # Integrate new variables
        # Set objective
        jamMasterProb.setObjective(approxVar, gurobipy.GRB.MINIMIZE)
        # constraints
        jamMasterProb.addConstr(sum([JammingGraph.node[l]['cost'] * jamPlaced[l] for l in JammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        if modelType == 'cormican-path-based-jam-vars':
            None
        elif modelType == 'cormican-path-based-jam-and-interdict-vars':
            edgeInterdict = dict([(e, jamMasterProb.addVar(vtype=gurobipy.GRB.BINARY, name="edgeInterdict_e"+str(e))) for e in G.edges()])
            [jamMasterProb.addConstr(edgeInterdict[e] <= sum([jamPlaced[l] for l in JammingGraph.nodes() 
                                        if (int(isEdgeJammedByJammer_Protocol(G, e, l, interfModelType)) == 1)]), "arcInterdict_e"+str(e)) for e in G.edges()]
        jamMasterProb.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            jamMasterProb.write("/tmp/JamMasterProb.lp")
    except gurobipy.GurobiError as e:
        print str(e)
    return jamMasterProb

def solveJammingMasterProblem_Benders(rel_mip_gap_tol = 0.0001):
    createSubProbMods()
    jamMasterProb = createJammingMasterProblem(instance.CGraph, instance.jamGraph, interfModelType)
    jamMasterProb._thetaVars = approxVar
    jamMasterProb._jamPlacedVars = [jamPlaced[n] for n in instance.jamGraph.nodes()]
    jamMasterProb.params.LazyConstraints = 1
    jamMasterProb.Params.MIPGap = rel_mip_gap_tol
    if modelType == 'cormican-path-based-jam-vars':
        jamMasterProb.optimize(bendersCutCallback_pathBased_JamVars)
    elif modelType == 'cormican-path-based-jam-and-interdict-vars':
        jamMasterProb._edgeInterdictVars = [edgeInterdict[e] for e in instance.CGraph.edges()]
        jamMasterProb.optimize(bendersCutCallback_pathBased_JamAndInterdictVars)
    print "OPTIMIZATION FINISHED"
    print "jamLocs", [n for n in instance.jamGraph.nodes() if jamPlaced[n].X > FUZZ]
    return jamMasterProb.objVal

def bendersCutCallback_pathBased_JamVars(model, where):
    if where == gurobipy.GRB.callback.MIPSOL:
        resetSets()
        print "***bendersCutCallback***"
        print "approxVar value", model.cbGetSolution(model._thetaVars)
        try:
            edgesToIndexDict = {}
            count = 0
            for e in instance.CGraph.edges():
                edgesToIndexDict[e] = count
                count += 1
            jamPlacedValues = model.cbGetSolution(model._jamPlacedVars)
            print "jamPlacedValues", jamPlacedValues
            jamPlacedValuesDict = {}
            jamLocsToIndexDict = {}
            jamLocVarsDict = {}
            count = 0
            for n in g_JamGraph.nodes():
                jamPlacedValuesDict[n] = jamPlacedValues[count]
                jamLocVarsDict[n] = model._jamPlacedVars[count]
                jamLocsToIndexDict[n] = count
                count += 1
            print "jamPlacedValuesDict", [(n, jamPlacedValuesDict[n]) for n in jamPlacedValuesDict.keys() if jamPlacedValuesDict[n] > 0.0]
            jamLocs = [n for n in jamPlacedValuesDict.keys() if jamPlacedValuesDict[n] > 0.0]
            print "jamLocs", jamLocs
            for n in g_JamGraph.nodes():
                g_JamGraph.node[n]['selected'] = jamPlacedValuesDict[n]
            flowVarSoln, isetUsageSoln, throughput = solveThroughputProblem_Pricing_CormicanPathBased(instance.CGraph, None, jamLocs, interfModelType)
            print "after flowVarSoln"
            expr = getConstraintExpressionForBendersCut_JamVarsOnly(g_JamGraph, flowVarSoln, jamLocVarsDict, model._thetaVars)
            print "expr", expr
            model.cbLazy(expr)
            dummyExpr = 0.5 - 0.25 * model._jamPlacedVars[0] <= model._thetaVars
            #print "dummyExpr", dummyExpr
            #model.cbLazy(dummyExpr)
        except gurobipy.GurobiError as e:
            print "bendersCutCallback error: ", e.message

def getJamLocsFinal():
    #print "jamPlaced", jamPlaced
    jamLocs = [n for n in instance.jamGraph.nodes() if jamPlaced[n].X > FUZZ]
    return jamLocs

def bendersCutCallback_arcBased_JamVars(model, where):
    global bestMIP_upperBound, smallestSubprobThroughput
    if where == gurobipy.GRB.callback.MIP:
        bestMIP_upperBound = model.cbGet(gurobipy.GRB.callback.MIP_OBJBST)
    if where == gurobipy.GRB.callback.MIPSOL:
        #resetSets()
        print "***bendersCutCallback***"
        obj     = model.cbGet(gurobipy.GRB.callback.MIPSOL_OBJ)
        myTime = time.time()
        print "obj", obj, "nodeCnt", myTime
        #objbnd = model.cbGet(gurobipy.GRB.callback.MIP_OBJBND)
        #print "ub", objbnd
        print "approxVar value", model.cbGetSolution(model._thetaVars)
        try:
            edgesToIndexDict = {}
            count = 0
            for e in edgeTriples:
                edgesToIndexDict[e] = count
                count += 1
            #print "edgesToIndexDict", edgesToIndexDict
            jamPlacedValues = model.cbGetSolution(model._jamPlaced)
            #print "jamPlacedValues", jamPlacedValues
            jamPlacedValuesDict = {}
            jamLocsToIndexDict = {}
            jamLocVarsDict = {}
            count = 0
            for n in instance.jamGraph.nodes():
                jamPlacedValuesDict[n] = jamPlacedValues[count]
                jamLocVarsDict[n] = model._jamPlaced[count]
                jamLocsToIndexDict[n] = count
                count += 1
            jamLocs = [n for n in jamPlacedValuesDict.keys() if jamPlacedValuesDict[n] > 0.0]
            print "jamLocs", jamLocs
            for n in instance.jamGraph.nodes():
                instance.jamGraph.node[n]['selected'] = jamPlacedValuesDict[n]
            resetObjectiveForThroughputProblem(instance.CGraph, jamLocs, instance.commodities)
            flowVarSolnForBendersCut, isetUsageSoln, throughput, throughputWithoutPenaltyForBendersCut = solveThroughputProblem_Pricing_CormicanArcBased(instance.CGraph, None, interfModelType, g_capacityDuals)
            #bendersExpr = None
            #bendersExpr = getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarSoln, throughputWithoutPenalty, jamLocVarsDict, model._thetaVars)
            #print "bendersExpr", bendersExpr
            if throughput < smallestSubprobThroughput:
                smallestSubprobThroughput = throughput
            if includeParetoOptimalBendersCuts is True:
                #print "pareto cut"
                corePoint = getCorePoint_Fractional()
                #print "corePoint", corePoint
                resetObjectiveForThroughputProblem_JamVector(instance.CGraph, corePoint, instance.commodities)
                addConstraintForThroughputValue(instance.CGraph, jamLocs, instance.commodities, throughput)
                flowVarSolnForBendersCut, isetUsageSolnPareto, throughputPareto, throughputWithoutPenaltyForBendersCut = solveThroughputProblem_Pricing_CormicanArcBased(instance.CGraph, None, interfModelType, g_capacityDuals, suffix = 'pareto' + str(myTime))
                bendersExpr = getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarSolnForBendersCut, throughputWithoutPenaltyForBendersCut, jamLocVarsDict, model._thetaVars)
                print "paretoCutExpr", bendersExpr
                gurobiThroughputModel.remove(throughputEqualsValueConstr)
            else:
                bendersExpr = getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarSolnForBendersCut, throughputWithoutPenaltyForBendersCut, jamLocVarsDict, model._thetaVars)
                print "bendersExpr", bendersExpr
            model.cbLazy(bendersExpr)
            if includeParetoOptimalBendersCuts:
                masterProbLP.addConstr(bendersExpr)
            if includeKnapsack is True:
                knapsackExpr = getConstraintExpressionFor_KI_Cut_JamVarsOnly_ArcBased(bestMIP_upperBound, flowVarSolnForBendersCut, throughputWithoutPenaltyForBendersCut, jamLocVarsDict, model._thetaVars)
                #print "rhs", bestMIP_upperBound, throughputWithoutPenaltyForBendersCut, math.floor(bestMIP_upperBound - throughputWithoutPenaltyForBendersCut)
                #print "knapsackExpr", knapsackExpr
                model.cbLazy(knapsackExpr)
                if includeParetoOptimalBendersCuts:
                    masterProbLP.addConstr(knapsackExpr)
        except gurobipy.GurobiError as e:
            print "bendersCutCallback error: ", e.message
            
def bendersCutCallback_pathBased_JamAndInterdictVars(model, where):
    if where == gurobipy.GRB.callback.MIPSOL:
        print "***bendersCutCallback***"
        try:
            edgeInterdictValues = model.cbGetSolution(model._edgeInterdictVars)
            edgeInterdictValuesDict = {}
            edgeInterdictVarsDict = {}
            edgesToIndexDict = {}
            count = 0
            for e in instance.CGraph.edges():
                edgeInterdictValuesDict[e] = edgeInterdictValues[count]
                edgeInterdictVarsDict[e] = model._edgeInterdictVars[count]
                edgesToIndexDict[e] = count
                count += 1
            print "edgeInterdictVarsDict", edgeInterdictVarsDict
            jamPlacedValues = model.cbGetSolution(model._jamPlacedVars)
            print "jamPlacedValues", jamPlacedValues
            jamPlacedValuesDict = {}
            jamLocsToIndexDict = {}
            count = 0
            for n in g_JamGraph.nodes():
                jamPlacedValuesDict[n] = jamPlacedValues[count]
                jamLocsToIndexDict[n] = count
                count += 1
            print "jamPlacedValuesDict", [(n, jamPlacedValuesDict[n]) for n in jamPlacedValuesDict.keys() if jamPlacedValuesDict[n] > 0.0]
            print "edgeInterdictValuesDict", [(e, edgeInterdictValuesDict[e]) for e in edgeInterdictValuesDict.keys() if edgeInterdictValuesDict[e] > 0.0]
            edgesJammed = [e for e in edgeInterdictValuesDict.keys() if edgeInterdictValuesDict[e] > 0.0]
            print "edgesJammed", edgesJammed
            for n in g_JamGraph.nodes():
                g_JamGraph.node[n]['selected'] = jamPlacedValuesDict[n]
            flowVarSoln, isetUsageSoln = solveThroughputProblem_Pricing_CormicanPathBased(instance.CGraph, edgesJammed, maxNumHops, interfModelType)
            #print "edgeInterdictValuesDict", [(e, edgeInterdictValuesDict[e]) for e in edgeInterdictValuesDict.keys() if edgeInterdictValuesDict[e] > 0.0]
            #for commod in commodities.keys():
            #    print "flowVarSoln", [(commod, path, flowVarSoln[commod][path], Paths[commod][path]) for path in flowVars[commod].keys() if flowVarSoln[commod][path] > FUZZ]
            print "after flowVarSoln"
            expr = getConstraintExpressionForBendersCut_JamAndEdgeInterdictVars(g_JamGraph, flowVarSoln, edgeInterdictVarsDict, model._thetaVars)
            print "expr", expr
            model.cbLazy(expr)
        except gurobipy.GurobiError as e:
            print "bendersCutCallback error: ", e.message
            
# def solveFullMIP_DelayedRowGen_withCallbacks(G, JammingGraph, ISets, interfModelType):
#     create_SingleLevelMIP_Cormican_ArcBased(G, instance.commodities, ISets, JammingGraph, interfModelType)
#     gurobiModel._usageDual = usageDual
#     gurobiModel._capDuals = [capDuals[e] for e in G.edges()]
#     gurobiModel.params.LazyConstraints = 1
#     gurobiModel.optimize(newISCutCallback)
#     return gurobiModel.objVal

def getInitialEdgeWeights():
    weights = {}
    for e in edgeTriples:
        weights[e] = 1.0
    return weights

def getInitialNodeWeights():
    initialNodeWeights = {}
    for edgeTriple in edgeTriples:
        attrIndex = getAttrIndex(edgeTriple)
        initialNodeWeights[edgeTriple] = instance.CGraph.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity']
    return initialNodeWeights

def setVars(modelType, useVirtualEdges = True):
    gurobiModel._usageDual = usageDual
    if useVirtualEdges:
        edgeSetToIterateOver = edgeTriples
    else:
        edgeSetToIterateOver = edgeTriplesInterdictable
    gurobiModel._capDuals = [capDuals[triple] for triple in edgeSetToIterateOver]
    gurobiModel._jamPlaced = [jamPlaced[node] for node in instance.jamGraph.nodes()]
    gurobiModel._demForCommodDuals = [demForCommodDuals[commod] for commod in instance.commodities.keys()]
    if modelType == 'regular':
        gurobiModel._jamUBDual = []
        for edgeInfo in instance.CGraph.edges(data = True):
            edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
            for node in jammersThatCanJamEdge[edgeTriple]:
                gurobiModel._jamUBDual.append(jamUBDual[edgeTriple][node])
    if jamVarsType == 'jam-and-interdict-vars':
        gurobiModel._edgeInterdict = [edgeInterdict[triple] for triple in edgeTriplesInterdictable]
        
def setVars_Benders(modelType, useVirtualEdges = True):
    gurobiModel._usageDual = usageDual
    if useVirtualEdges:
        edgeSetToIterateOver = edgeTriples
    else:
        edgeSetToIterateOver = edgeTriplesInterdictable
    gurobiModel._thetaVars = approxVar
    gurobiModel._jamPlaced = [jamPlaced[node] for node in instance.jamGraph.nodes()]
    if modelType == 'regular':
        gurobiModel._jamUBDual = []
        for edgeInfo in instance.CGraph.edges(data = True):
            edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
            for node in jammersThatCanJamEdge[edgeTriple]:
                gurobiModel._jamUBDual.append(jamUBDual[edgeTriple][node])
    if jamVarsType == 'jam-and-interdict-vars':
        gurobiModel._edgeInterdict = [edgeInterdict[triple] for triple in edgeTriplesInterdictable]
        
def setVars_Benders_Classic(modelType, useVirtualEdges = True):
    None
    
def setParams(rel_mip_gap_tol, lazy = True):
    if lazy:
        gurobiModel.params.LazyConstraints = 1
    if algType == 'b-and-c-reg':
        gurobiModel.Params.Cuts = 0
    gurobiModel.Params.MIPGap = rel_mip_gap_tol
    gurobiModel.Params.Threads = numThreads
    gurobiModel.Params.MIPFocus = 2
    gurobiModel.Params.TimeLimit = time_limit
    
def createSubProbMods():
    global shortestPathModels, findPathflowVars
    if flowVarsType == 'path-based':
        findPathflowVars = {}
        shortestPathModels = {}
        initialEdgeWeights = getInitialNodeWeights()
        for commod in instance.commodities.keys():
            create_ShortestPathProb(instance.CGraph, initialEdgeWeights, maxNumHops, commod)   
    create_MaxWtIndepSet_Model(instance.interferenceGraph, getInitialNodeWeights(), interfModelType)
    createThroughputModel_ContinuousJamming_Cormican(instance.CGraph, instance.commodities, ISets, [], interfModelType)
    
def solveFullMIP_DelayedRowGen_withCallbacks(G, JammingGraph, Paths, ISets, max_hop_length, interfModelType, rel_mip_gap_tol = 0.0001):
    #G is the connectivity graph
    print "solveFullMIP_DelayedRowGen_withCallbacks maxNumHops", max_hop_length
    createSubProbMods()
    if flowVarsType == 'path-based' and modelType == 'cormican':
        create_SingleLevelMIP_Cormican_PathBased(G, instance.jamGraph, instance.commodities, Paths, ISets, interfModelType)
    elif flowVarsType == 'path-based' and modelType == 'regular':
        create_SingleLevelMIP_Regular_PathBased(G, instance.jamGraph, instance.commodities, Paths, ISets, interfModelType)
    elif flowVarsType == 'arc-based' and modelType == 'cormican': 
        create_SingleLevelMIP_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType) # use this one
    elif flowVarsType == 'arc-based' and modelType == 'regular':
        create_SingleLevelMIP_Regular_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType) 
    setVars(instance.modelType)
    setParams(rel_mip_gap_tol)
    if flowVarsType == 'path-based' and modelType == 'cormican':
        gurobiModel.optimize(new_PathAndIS_CutCallback)
    elif flowVarsType == 'path-based' and modelType == 'regular':
        gurobiModel.optimize(new_PathAndIS_CutCallback)
    elif flowVarsType == 'arc-based' and modelType == 'cormican':
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
            gurobiModel.optimize(newISCutCallback) # use this one; this solve model (9) using callbacks with model (10) as the separation problem (newISCutCallback is the callback function)
        else:
            gurobiModel.optimize()
    elif flowVarsType == 'arc-based' and modelType == 'regular':
        gurobiModel.optimize(newISCutCallback)
    return gurobiModel.objVal

def runBenders_Callback(G, JammingGraph, Paths, ISets, max_hop_length, interfModelType, rel_gap_tol = 0.0001):
    print "solveFullMIP_DelayedRowGen_withCallbacks maxNumHops", max_hop_length
    createSubProbMods()
    if flowVarsType == 'arc-based' and modelType == 'cormican':
        create_BendersMasterProb_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
        if includeParetoOptimalBendersCuts:
            create_BendersMasterProb_LP_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
    elif flowVarsType == 'arc-based' and modelType == 'regular':
        raise Exception("configuration must be flowVarsType == 'arc-based' and modelType == 'cormican'")
    setVars_Benders(instance.modelType)
    setParams(rel_gap_tol)
    if flowVarsType == 'arc-based' and modelType == 'cormican':
        gurobiModel.optimize(bendersCutCallback_arcBased_JamVars)
    if includeSuperValid is True:
        #print "print pareto"
        for cut in sviTypeI_Cuts:
            print "cut", cut
        for nodesSet in sviTypeI_NodesInLHS:
            print "nodesSet", nodesSet
        sviTypeI_NodesInLHS_Set = []
        for nodesSet in sviTypeI_NodesInLHS:
            sviTypeI_NodesInLHS_Set.append(set(nodesSet))
        print "intersect", set.intersection(*sviTypeI_NodesInLHS_Set)
        print "status", gurobiModel.status
    
    lb = gurobiModel.objVal
    ub = gurobiModel.getAttr('ObjBound')
    jamLocs = [n for n in jamPlaced.keys() if jamPlaced[n].X > 0.0]
    return lb, ub, jamLocs

def runBenders_Classic(G, JammingGraph, Paths, ISets, max_hop_length, interfModelType, rel_gap_tol = 0.0001):
    global numNodesExplored
    print "run benders classic"
    createSubProbMods()
    if flowVarsType == 'arc-based' and modelType == 'cormican':
        create_BendersMasterProb_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
        if includeParetoOptimalBendersCuts:
            create_BendersMasterProb_LP_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
    elif flowVarsType == 'arc-based' and modelType == 'regular':
        raise Exception("configuration must be flowVarsType == 'arc-based' and modelType == 'cormican'")
    setVars_Benders_Classic(instance.modelType)
    setParams(rel_gap_tol, lazy = False)
    lb = 0
    ub = float('inf')
    iteration = 1
    startTime = time.time()
    while not hasTerminated(lb, ub, rel_gap_tol, time.time() - startTime):
        timeRemaining = time_limit - (time.time() - startTime)
        gurobiModel.Params.TimeLimit = timeRemaining
        gurobiModel.optimize()
        jamLocs = [n for n in jamPlaced.keys() if jamPlaced[n].X > 0.0]
        numNodesExplored += gurobiModel.NodeCount
        lb = gurobiModel.objVal
        resetObjectiveForThroughputProblem(instance.CGraph, jamLocs, instance.commodities)
        flow, isSetUsage, th, thNoPen = solveThroughputProblem_Pricing_CormicanArcBased(instance.CGraph, None, interfModelType, g_capacityDuals)
        ub = min(th, ub)
        addCuts_BendersClassic(gurobiModel, jamLocs, flow, th, thNoPen, jamPlaced, approxVar, iteration, ub)
        printBendersClassicInfo(iteration, lb, ub, jamLocs)
        iteration += 1
    return lb, ub, jamLocs

def computeThroughputForJammingSolution(G, JammingGraph, jamLocs, Paths, ISets, max_hop_length, interfModelType, rel_gap_tol = 0.0001):
    global numNodesExplored
    createSubProbMods()
    if flowVarsType == 'arc-based' and modelType == 'cormican':
        create_BendersMasterProb_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
        if includeParetoOptimalBendersCuts:
            create_BendersMasterProb_LP_Cormican_ArcBased(G, instance.jamGraph, instance.commodities, ISets, interfModelType)
    elif flowVarsType == 'arc-based' and modelType == 'regular':
        raise Exception("configuration must be flowVarsType == 'arc-based' and modelType == 'cormican'")
    setParams(rel_gap_tol, lazy = False)
    resetObjectiveForThroughputProblem(instance.CGraph, jamLocs, instance.commodities)
    flow, isSetUsage, th, thNoPen = solveThroughputProblem_Pricing_CormicanArcBased(instance.CGraph, None, interfModelType, g_capacityDuals)
    return th

def printBendersClassicInfo(iteration, lb, ub, jamLocs):
    print "BENDERS ITERATION", iteration, lb, ub, "jamLocs", jamLocs

def hasTerminated(lb, ub, rel_gap_tol, time):
    if (ub - lb) <= lb * rel_gap_tol:
        return True
    elif time >= time_limit:
        return True
    else:
        return False
    
def addCuts_BendersClassic(myModel, jamLocs, flowSoln, throughput, thNoPen, jamLocVarsDict, approxVar, iteration, bestUB):
    global smallestSubprobThroughput
    flowVarSolnForBendersCut = flowSoln
    throughputNoPenForBendersCut = thNoPen
    jamLocs = [n for n in jamLocVarsDict.keys() if jamLocVarsDict[n].X > 0.0]
    if includeTrustRegion is True:
        if iteration <= trust_region_max_iter:
            threshold = max(2, instance.jamBudget)
            setTrustRegionConstraint(myModel, jamLocs, threshold, iteration)
        else:
            myModel.remove(trustRegionConstr)
    if throughput < smallestSubprobThroughput:
        smallestSubprobThroughput = throughput
    if includeParetoOptimalBendersCuts is True:
        corePoint = getCorePoint_Fractional()
        resetObjectiveForThroughputProblem_JamVector(instance.CGraph, corePoint, instance.commodities)
        addConstraintForThroughputValue(instance.CGraph, jamLocs, instance.commodities, throughput)
        flowVarSolnForBendersCut, isetUsageSolnPareto, throughputPareto, throughputNoPenForBendersCut = solveThroughputProblem_Pricing_CormicanArcBased(instance.CGraph, None, interfModelType, g_capacityDuals, suffix = 'classic-pareto')
        bendersExpr = getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarSolnForBendersCut, throughputNoPenForBendersCut, jamLocVarsDict, approxVar)
        print "paretoCutExpr", bendersExpr
        gurobiThroughputModel.remove(throughputEqualsValueConstr)
    else:
        bendersExpr = getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarSolnForBendersCut, throughputNoPenForBendersCut, jamLocVarsDict, approxVar)
        print "bendersExpr", bendersExpr
    myModel.addConstr(bendersExpr)
    if includeParetoOptimalBendersCuts:
        masterProbLP.addConstr(bendersExpr)
    if includeKnapsack is True:
        knapsackExpr = getConstraintExpressionFor_KI_Cut_JamVarsOnly_ArcBased(bestUB, flowVarSolnForBendersCut, throughputNoPenForBendersCut, jamLocVarsDict, approxVar)
        #print "rhs", bestMIP_upperBound, throughputWithoutPenaltyForBendersCut, math.floor(bestMIP_upperBound - throughputWithoutPenaltyForBendersCut)
        print "knapsackExpr", knapsackExpr
        myModel.addConstr(knapsackExpr)
        if includeParetoOptimalBendersCuts:
            masterProbLP.addConstr(knapsackExpr)
    if includeSuperValid is True:
        print "smallestSubprobThroughput", smallestSubprobThroughput
        addAndUpdateSupervalidCuts_Classic(myModel, bestUB, smallestSubprobThroughput, flowVarSolnForBendersCut, throughputNoPenForBendersCut, jamLocVarsDict)
        
def setTrustRegionConstraint(myModel, jamLocs, trustRegionThreshold, iteration):
    global trustRegionConstr
    # remove old
    if iteration > 1:
        myModel.remove(trustRegionConstr)
        myModel.update()
    trustRegionConstr = None
    lhsSum = 0
    print "jammers placed", jamLocs 
    for n in instance.jamGraph.nodes():
        if n in jamLocs:
            lhsSum += (1 - jamPlaced[n])
        else:
            lhsSum += jamPlaced[n]
    trustRegionConstr = myModel.addConstr(lhsSum <= trustRegionThreshold)
    myModel.update()
    print "trustRegionConstr", lhsSum <= trustRegionThreshold
       
def getNonZeroValues(myDict):
    return [(key, myDict[key]) for key in myDict.keys() if myDict[key] > FUZZ]

def createJamPlacedDicts(jamPlacedValues):
    jamPlacedDict = {}
    count = 0
    for node in instance.jamGraph.nodes():
        jamPlacedDict[node] = jamPlacedValues[count]
        count += 1
    return jamPlacedDict

def createNodesToIndexDicts():
    global nodesToIndexInterdictDict
    nodesToIndexInterdictDict = {}
    count = 0
    for node in instance.jamGraph.nodes():
        nodesToIndexInterdictDict[node] = count
        count += 1
    return nodesToIndexInterdictDict

def createDemForCommodDualsDict(demForCommodDualsValues):
    demForCommodDualsValuesDict = {}
    count = 0
    for commod in instance.commodities.keys():
        demForCommodDualsValuesDict[commod] = demForCommodDualsValues[count]
        count += 1
    return demForCommodDualsValuesDict

def createCommodToIndexDict():
    global commodToIndexDict
    commodToIndexDict = {}
    count = 0
    for commod in instance.commodities.keys():
        commodToIndexDict[commod] = count
        count += 1
    return commodToIndexDict

def createEdgeInterdictDict(edgeInterdictValues):
    edgeInterdictsDict = {}
    count = 0
    for e in edgeTriplesInterdictable:
        edgeInterdictsDict[e] = edgeInterdictValues[count]
        count += 1
    return edgeInterdictsDict

def OLD_createEdgeInterdictAndEdgesToIndexDicts(edgeInterdictValues):
    edgeInterdictsDict = {}
    edgesToIndexInterdictDict = {}
    count = 0
    for e in edgeTriplesInterdictable:
        edgeInterdictsDict[e] = edgeInterdictValues[count]
        edgesToIndexInterdictDict[e] = count
        count += 1
    return edgeInterdictsDict, edgesToIndexInterdictDict

def getShortestPathEdgeWeights(capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, jamLocs, 
                               nodesToIndexInterdictDict, modelType, jamUBDualValues):
    if jamVarsType == 'jam-vars' and modelType == 'regular':
        jamUBDualsValuesDict = getJamUBDUalValuesDict(jamUBDualValues)
    #print "edgesToIndexInterdictDict", edgesToIndexInterdictDict
    weights = {}
    for e in edgeTriples:
        value = capDualValues[edgesToIndexDict[e]]
        if jamVarsType == 'jam-and-interdict-vars' and modelType == 'cormican':
            if e in edgeTriplesInterdictable:
                value += edgeInterdictValues[edgesToIndexInterdictDict[e]]
        elif jamVarsType == 'jam-vars' and modelType == 'cormican':
            if e in edgeTriplesInterdictable:
                value += int(isEdgeJammedByJammers(instance.CGraph, e, jamLocs, interfModelType))
        elif jamVarsType == 'jam-vars' and modelType == 'regular':
            if e in edgeTriplesInterdictable:
                for node in getJammersThatCanJamEdge(instance.CGraph, instance.jamGraph, e, interfModelType):
                    value += jamUBDualsValuesDict[e][node]
        weights[e] = value
    return weights

# def getShortestPathEdgeWeights_JamVars(capDualValues, edgesToIndexDict, jamPlacedValues, nodesToIndexInterdictDict, modelType):
#     weights = {}
#     for e in edgeTriples:
#         value = capDualValues[edgesToIndexDict[e]]
#         if modelType == 'cormican-path-based-jam-vars':
#             for node in jammersThatCanJamEdge[e]:
#                 value += jamPlacedValues[nodesToIndexInterdictDict[node]]
#         weights[e] = value
#     return weights

def addISetToModel(model, capDualValuesDict, edgesToIndexDict, usageDualValue, gap_for_ISetProb = 0.05): # this solves model (10) and adds a constraint to the master problem (model (9))
    global ISets
    nodeWeights = {}
    for edgeTriple in edgeTriples:
        attrIndex = getAttrIndex(edgeTriple)
        nodeWeights[edgeTriple] = capDualValuesDict[edgeTriple] * instance.CGraph.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity']
    # FIRST, TRY SOLVING IT HEURISTICALLY (I started this but couldn't get it to work without some effort)
#     maxWeights, solnsSet = modifyAndSolve_maxWtIndepSet_ViaHeuristic(instance.interferenceGraph, nodeWeights, interfModelType, callbackCtr, gap_for_ISetProb)
#     threshold = 1.0 * maxWeights[0] #only the best are added
#     numberAllowedToAdd = 1
#     for index in range(len(maxWeights)):
#         if maxWeights[index] >= max(usageDualValue, threshold) and index < numberAllowedToAdd :
#             print "TO ADD NEW HEURISTIC CUT"
#             ISets.append(solnsSet[index])
#             attrIndex = getAttrIndex(edgeTriple)
#             expr = sum([-instance.CGraph.edge[e[0]][e[1]][attrIndex]['capacity'] * model._capDuals[edgesToIndexDict[e]] for e in solnsSet[index]]) + model._usageDual >= 0
#             model.cbLazy(expr)
#             print "   add new heuristic cut for ISEt:", solnsSet[index], expr
#             return
    # IF THAT FAILS, SOLVE TO OPTIMALITY
    print "Solve to optimality"
    maxWeights, solnsSet = modifyAndSolve_maxWtIndepSet(instance.interferenceGraph, nodeWeights, interfModelType, callbackCtr, gap_for_ISetProb) # changes beta values and resolves (10)
    #print "node weights: ", [(edgeTriple, capDualValuesDict[edgeTriple]) for edgeTriple in edgeTriples if capDualValuesDict[edgeTriple] > FUZZ]
    print "CALLBACK (add ISet): maxWt", maxWeights, "usageDual", usageDualValue
    threshold = 1.0 * maxWeights[0] #only the best are added
    numberAllowedToAdd = 1
    for index in range(len(maxWeights)): # adds new independent sets back to (9)
        if maxWeights[index] >= max(usageDualValue, threshold) and index < numberAllowedToAdd :
            ISets.append(solnsSet[index])
            attrIndex = getAttrIndex(edgeTriple)
            expr = sum([-instance.CGraph.edge[e[0]][e[1]][attrIndex]['capacity'] * model._capDuals[edgesToIndexDict[e]] for e in solnsSet[index]]) + model._usageDual >= 0
            model.cbLazy(expr)
            print "   add new cut for ISEt:", solnsSet[index], expr
            
def newISCutCallback(model, where): # this solves model (10) and adds a new constraint to model (9)
    global previous_mip_gap, bestBoundGlobal, bestObjGlobal, callbackCtr
    if where == gurobipy.GRB.callback.MIP:
        bestBoundGlobal = model.cbGet(gurobipy.GRB.callback.MIP_OBJBND)
        bestObjGlobal = model.cbGet(gurobipy.GRB.callback.MIP_OBJBST)
    if where == gurobipy.GRB.callback.MIPSOL:
        print "newISCutCallback"
        try:
            #get variables values and create data structures
            runTime = model.cbGet(gurobipy.GRB.callback.RUNTIME)
            bestObjLocal = model.cbGet(gurobipy.GRB.callback.MIPSOL_OBJBST)
            bestBoundLocal = model.cbGet(gurobipy.GRB.callback.MIPSOL_OBJBND)
            if(bestObjGlobal != 0.0 and bestObjGlobal < 1e+100):
                mip_gap = (bestObjGlobal - bestBoundGlobal)/ bestObjGlobal
            else:
                mip_gap = float("inf")
            node_count = model.cbGet(gurobipy.GRB.callback.MIPSOL_NODCNT)
            gap_change = previous_mip_gap - mip_gap
            if(gap_change < 0.01 and mip_gap < float("inf")):
                gap_for_ISetProb = 0.0001
            else:
                gap_for_ISetProb = 0.10
            usageDualValue = model.cbGetSolution(model._usageDual)
            capDualValues = model.cbGetSolution(model._capDuals)
            if jamVarsType == 'jam-and-interdict-vars':
                edgeInterdictValues = model.cbGetSolution(model._edgeInterdict)
                edgeInterdictsDict = createEdgeInterdictDict(edgeInterdictValues)
            jamPlacedValues = model.cbGetSolution(model._jamPlaced)
            capDualValuesDict = {}
            count = 0
            for e in edgeTriples:
                capDualValuesDict[e] = capDualValues[count]
                count += 1
            jamPlacedDict = createJamPlacedDicts(jamPlacedValues)
            addISetToModel(model, capDualValuesDict, edgesToIndexDict, usageDualValue, gap_for_ISetProb) # look at this; this solves model (10) 
            previous_mip_gap = mip_gap
            callbackCtr += 1
        except gurobipy.GurobiError as e:
            print "newISCutCallback error", e.message
            
def doMIPStuff(model):
    bestBoundGlobal = model.cbGet(gurobipy.GRB.callback.MIP_OBJBND)
    bestObjGlobal = model.cbGet(gurobipy.GRB.callback.MIP_OBJBST)
    
def doMIPNodeStuff(model):
    #print "where == gurobipy.GRB.callback.MIPNODE", where == gurobipy.GRB.callback.MIPNODE
    if jamVarsType == 'jam-and-interdict-vars':
        edgeInterdictsDict = createEdgeInterdictDict(model.cbGetNodeRel(model._edgeInterdict))
        fracEdgeInterdict = getNonZeroValues(edgeInterdictsDict)
        #if len(fracEdgeInterdict) > 0:
        #    print "   edgeInterdictDict_frac", fracEdgeInterdict
    jamPlacedDict = createJamPlacedDicts(model.cbGetNodeRel(model._jamPlaced))
    fracJamPlaced = getNonZeroValues(jamPlacedDict)
    
def getGapForISetProb(model):
    runTime = model.cbGet(gurobipy.GRB.callback.RUNTIME)
    bestObjLocal = model.cbGet(gurobipy.GRB.callback.MIPSOL_OBJBST)
    bestBoundLocal = model.cbGet(gurobipy.GRB.callback.MIPSOL_OBJBND)
    if(bestObjGlobal != 0.0 and bestObjGlobal < 1e+100):
        mip_gap = (bestObjGlobal - bestBoundGlobal)/ bestObjGlobal
    else:
        mip_gap = float("inf")
    node_count = model.cbGet(gurobipy.GRB.callback.MIPSOL_NODCNT)
    gap_change = previous_mip_gap - mip_gap
    if(gap_change < 0.01 and mip_gap < float("inf")):
        gap_for_ISetProb = 0.0001
    else:
        gap_for_ISetProb = 0.10
    return gap_for_ISetProb

def getPathsToAdd(model, modelType, jamLocs, capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, 
                  nodesToIndexInterdictDict, demForCommodDualsValuesDict, edgesJammed):
    NewPaths = []
    #startTime = time.time()
    jamUBDualValues = None
    if modelType == 'regular':
        jamUBDualValues = model.cbGetSolution(model._jamUBDual)
    for commod in instance.commodities.keys():
        weights = getShortestPathEdgeWeights(capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, jamLocs, nodesToIndexInterdictDict, instance.modelType, jamUBDualValues)
        #print "weights", weights
        pathLength, edgesOnNewPath = modifyAndSolve_ShortestPathProb(instance.CGraph, weights, maxNumHops, commod, callbackCtr)
        if newPathPricesOut(model, modelType, jamLocs, demForCommodDualsValuesDict, commod, edgesOnNewPath, pathLength, edgesJammed, jamUBDualValues):
            #print "new path: ", pathLength
            NewPaths.append((commod, pathNumForCommod[commod], edgesOnNewPath, pathLength))
            Paths[commod].append(edgesOnNewPath)
            pathNumForCommod[commod] = pathNumForCommod[commod] + 1
    return NewPaths

# def getPathsToAdd_OLD(model, modelType, jamLocs, capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, jamPlacedValues, 
#                   nodesToIndexInterdictDict, demForCommodDualsValuesDict, edgesJammed):
#     NewPaths = []
#     #startTime = time.time()
#     jamUBDualValues = None
#     if modelType == 'regular':
#         jamUBDualValues = model.cbGetSolution(model._jamUBDual)
#     for commod in instance.commodities.keys():
#         weights = getShortestPathEdgeWeights(capDualValues, edgesToIndexDict, jamPlacedValues, nodesToIndexInterdictDict, instance.modelType, jamUBDualValues)
#         #print "weights", weights
#         pathLength = 0
#         pathsToExclude = []
#         index = 0
#         edgesOnNewPath = None
#         while index < 1 and newPathPricesOut(model, modelType, jamLocs, demForCommodDualsValuesDict, commod, edgesOnNewPath, pathLength, edgesJammed, jamUBDualValues):# 3 based on limited pilot experiments
#             pathLength, edgesOnNewPath = modifyAndSolve_ShortestPathProb(instance.CGraph, weights, maxNumHops, commod, callbackCtr, pathsToExclude)
#             #print "pathLength, edgesOnNewPath", pathLength, edgesOnNewPath
#             pathsToExclude.append(edgesOnNewPath)
#             if newPathPricesOut(model, modelType, jamLocs, demForCommodDualsValuesDict, commod, edgesOnNewPath, pathLength, edgesJammed, jamUBDualValues):
#                 print "new path: ", pathLength
#                 #if(pathLength < shortestPathLength):
#                 #    shortestPathLength = pathLength
#                 #    shortestPathIndex = shortestPathLengthCtr
#                 NewPaths.append((commod, pathNumForCommod[commod], edgesOnNewPath, pathLength))
#                 Paths[commod].append(edgesOnNewPath)
#                 pathNumForCommod[commod] = pathNumForCommod[commod] + 1
#                 #shortestPathLengthCtr += 1
#             index += 1
#     return NewPaths

def getCapDualsDict(capDualValues):
    capDualValuesDict = {}
    count = 0
    for e in edgeTriples:
        capDualValuesDict[e] = capDualValues[count]
        count += 1
    return capDualValuesDict

def getJamUBDUalToIndexDict():
    jamUBDualToIndexDict = {}
    count = 0
    for e in edgeTriples:
        jamUBDualToIndexDict[e] = {}
        for node in jammersThatCanJamEdge[e]:
            jamUBDualToIndexDict[e][node] = count
            count += 1
    #print "jamUBDualToIndexDict", jamUBDualToIndexDict
    return jamUBDualToIndexDict

def getJamUBDUalValuesDict(jamUBDualValues):
    #print "jamUBDualValues", jamUBDualValues, edgeTriples
    jamUBDualToIndexDict = {}
    count = 0
    for e in edgeTriples:
        jamUBDualToIndexDict[e] = {}
        for node in jammersThatCanJamEdge[e]:
            jamUBDualToIndexDict[e][node] = jamUBDualValues[count]
            count += 1
    return jamUBDualToIndexDict

def getSumTwo(model, modelType, edgesOnPath, edgesToIndexInterdictDict, nodesToIndexInterdictDict):
    if jamVarsType == 'jam-and-interdict-vars' and modelType == 'cormican':
        interdictableEdgesOnPath = list(set(edgesOnPath) & set(edgeTriplesInterdictable))
        sum2 = sum([model._edgeInterdict[edgesToIndexInterdictDict[e]] for e in interdictableEdgesOnPath])
    elif jamVarsType == 'jam-vars' and modelType == 'cormican':
        sum2 = sum([model._jamPlaced[nodesToIndexInterdictDict[node]] for node in getNodesThatCanInterdictEdgeSet(edgesOnPath)])
    elif jamVarsType == 'jam-vars' and modelType == 'regular':
        #print "sumTwo", edgesOnPath
        sum2 = 0
        myIndex = 0
        jamUBDualToIndexDict = getJamUBDUalToIndexDict()
        #print "jamUBDualToIndexDict", jamUBDualToIndexDict
        for edgeTriple in edgesOnPath:
            #print "edgeTriple", edgeTriple, jammersThatCanJamEdge[edgeTriple]
            for node in jammersThatCanJamEdge[edgeTriple]:
                #print "node", node, model._jamUBDual[myIndex]
                sum2 += model._jamUBDual[jamUBDualToIndexDict[edgeTriple][node]]
                myIndex += 1
    return sum2

def addNewPathCutsToModel(model, modelType, NewPaths, edgesToIndexDict, edgesToIndexInterdictDict, nodesToIndexInterdictDict, commodsToIndexDict):
    for item in NewPaths:
        commod = item[0]
        edgesOnPath = item[2]
        sum1 = sum([model._capDuals[edgesToIndexDict[e]] for e in edgesOnPath])
        sum2 = getSumTwo(model, modelType, edgesOnPath, edgesToIndexInterdictDict, nodesToIndexInterdictDict)
        model.cbLazy(sum1 + sum2 + model._demForCommodDuals[commodsToIndexDict[commod]] >= 1)
        print "ADD PATH", item[0], item[1], item[3], sum1 + sum2 >= 1


def new_PathAndIS_CutCallback(model, where):
    global previous_mip_gap, bestBoundGlobal, bestObjGlobal, callbackCtr, Paths
    if where == gurobipy.GRB.callback.MIP:
        doMIPStuff(model)
    elif where == gurobipy.GRB.callback.MIPNODE:
        doMIPNodeStuff(model)
    elif where == gurobipy.GRB.callback.MIPSOL:
        print "***callback***"
        try:
            gap_for_ISetProb = getGapForISetProb(model)
            usageDualValue = model.cbGetSolution(model._usageDual)
            print "usageDualValue", usageDualValue
            capDualValues = model.cbGetSolution(model._capDuals)
            demForCommodDualsValues = model.cbGetSolution(model._demForCommodDuals)
            edgesJammed = None
            edgeInterdictValues = None
            edgeInterdictsDict = None
            if jamVarsType == 'jam-and-interdict-vars':
                edgeInterdictValues = model.cbGetSolution(model._edgeInterdict)
                edgeInterdictsDict = createEdgeInterdictDict(edgeInterdictValues)
                edgesJammed = [key for key in edgeInterdictsDict.keys() if edgeInterdictsDict[key] > FUZZ]
            jamPlacedValues = model.cbGetSolution(model._jamPlaced)
            capDualValuesDict = getCapDualsDict(capDualValues)
            jamPlacedDict = createJamPlacedDicts(jamPlacedValues)
            demForCommodDualsValuesDict = createDemForCommodDualsDict(demForCommodDualsValues)
            jamLocs = [key for key in jamPlacedDict.keys() if jamPlacedDict[key] > FUZZ]
            print "jamPlacedDict", [(key, jamPlacedDict[key]) for key in jamPlacedDict.keys() if jamPlacedDict[key] > FUZZ]
            print "capDualValuesDict", [(key, capDualValuesDict[key]) for key in capDualValuesDict.keys() if capDualValuesDict[key] > FUZZ]
            NewPaths = getPathsToAdd(model, modelType, jamLocs, capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, 
                  nodesToIndexInterdictDict, demForCommodDualsValuesDict, edgesJammed)
            newPathAvailable = (len(NewPaths) > 0)
            if(newPathAvailable):# new paths added
                #print "new path avail"
                addNewPathCutsToModel(model, modelType, NewPaths, edgesToIndexDict, edgesToIndexInterdictDict, nodesToIndexInterdictDict, commodToIndexDict)
            else:#no new paths added
                print "no new paths"
                addISetToModel(model, capDualValuesDict, edgesToIndexDict, usageDualValue, gap_for_ISetProb)
            callbackCtr += 1
        except gurobipy.GurobiError as e:
            print "new_PathAndIS_CutCallback error", e.message
            
# def new_PathAndIS_CutCallback_Regular(model, where):
#     global previous_mip_gap, bestBoundGlobal, bestObjGlobal, callbackCtr
#     if where == gurobipy.GRB.callback.MIP:
#         doMIPStuff(model)
#     elif where == gurobipy.GRB.callback.MIPNODE:
#         doMIPNodeStuff(model)
#     elif where == gurobipy.GRB.callback.MIPSOL:
#         print "***callback***"
#         try:
#             #get variables values and create data structures
#             gap_for_ISetProb = getGapForISetProb(model)
#             usageDualValue = model.cbGetSolution(model._usageDual)
#             capDualValues = model.cbGetSolution(model._capDuals)
#             #print "capDualValues", capDualValues
#             demForCommodDualsValues = model.cbGetSolution(model._demForCommodDuals)
#             if jamVarsType == 'jam-and-interdict-vars':
#                 edgeInterdictValues = model.cbGetSolution(model._edgeInterdict)
#                 edgeInterdictsDict, edgesToIndexInterdictDict = createEdgeInterdictAndEdgesToIndexDicts(edgeInterdictValues)
#             jamPlacedValues = model.cbGetSolution(model._jamPlaced)
#             capDualValuesDict = {}
#             edgesToIndexDict = {}
#             count = 0
#             for e in edgeTriples:
#                 capDualValuesDict[e] = capDualValues[count]
#                 edgesToIndexDict[e] = count
#                 count += 1
#             jamPlacedDict, nodesToIndexInterdictDict = createJamPlacedAndNodesToIndexDicts(jamPlacedValues)
#             demForCommodDualsValuesDict, commodsToIndexDict = createCommodsToIndexDict(demForCommodDualsValues)
#             NewPaths = []
#             shortestPathLength = 1000000
#             shortestPathIndex = -1
#             shortestPathLengthCtr = 0
#             startTime = time.time()
#             for commod in instance.commodities.keys():
#                 #print "before shortest path"
#                 if jamVarsType == 'jam-and-interdict-vars':
#                     weights = getShortestPathEdgeWeights_JamAndInterdictVars(capDualValues, edgesToIndexDict, edgeInterdictValues, edgesToIndexInterdictDict, jamPlacedValues, nodesToIndexInterdictDict, instance.modelType)
#                 elif jamVarsType == 'jam-vars':
#                     weights = getShortestPathEdgeWeights_JamVars(capDualValues, edgesToIndexDict, jamPlacedValues, nodesToIndexInterdictDict, instance.modelType)
#                 pathLength = 0
#                 pathsToExclude = []
#                 index = 0
#                 while index < 1 and pathLength < (1.0 - FUZZ):# 3 based on limited pilot experiments
#                     pathLength, newPath = modifyAndSolve_ShortestPathProb(instance.CGraph, weights, maxNumHops, commod, callbackCtr, pathsToExclude)
#                     pathsToExclude.append(newPath)
#                     #print index, "pathLength", pathLength, newPath
#                     if(pathLength < (1.0 - FUZZ)):
#                         if(pathLength < shortestPathLength):
#                             shortestPathLength = pathLength
#                             shortestPathIndex = shortestPathLengthCtr
#                         NewPaths.append((commod, pathNumForCommod[commod], newPath, pathLength))
#                         pathNumForCommod[commod] = pathNumForCommod[commod] + 1
#                         shortestPathLengthCtr += 1
#                     index += 1
#             newPathAvailable = (len(NewPaths) > 0)
#             if(newPathAvailable):# new paths added
#                 #print "new path available"
#                 item = NewPaths[shortestPathIndex]
#                 edgesOnPath = item[2]
#                 commod = item[0]
#                 sum1 = sum([model._capDuals[edgesToIndexDict[e]] for e in edgesOnPath])
#                 if jamVarsType == 'jam-and-interdict-vars':
#                     interdictableEdgesOnPath = list(set(edgesOnPath) & set(edgeTriplesInterdictable))
#                     sum2 = sum([model._edgeInterdict[edgesToIndexInterdictDict[e]] for e in interdictableEdgesOnPath])
#                 elif jamVarsType == 'jam-vars' and modelType == 'cormican':
#                     sum2 = sum([model._jamPlaced[nodesToIndexInterdictDict[node]] for node in getNodesThatCanInterdictEdgeSet(edgesOnPath)])
#                 elif jamVarsType == 'jam-vars' and modelType == 'regular':
#                     sum2 = 0
#                     myIndex = 0
#                     for edgeInfo in instance.CGraph.edges(data = True):
#                         edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
#                         #print "edgeTriple", edgeTriple, jammersThatCanJamEdge[edgeTriple]
#                         for node in jammersThatCanJamEdge[edgeTriple]:
#                             #print "node", node, model._jamUBDual[myIndex]
#                             sum2 += model._jamUBDual[myIndex]
#                             myIndex += 1
#                 model.cbLazy(sum1 + sum2 + model._demForCommodDuals[commodsToIndexDict[commod]] >= 1)
#                 print "ADD PATH", item[0], item[1], item[3], sum1 + sum2 >= 1
#             else:#no new paths added
#                 #solve maxWtIndSet separation problem and add cut if needed
#                 addISetToModel(model, capDualValuesDict, edgesToIndexDict, usageDualValue, gap_for_ISetProb)
#                 #jamPlacedDict, nodesToIndexInterdictDict = createJamPlacedAndNodesToIndexDicts(jamPlacedValues)
#             callbackCtr += 1
#         except gurobipy.GurobiError as e:
#             print "new_PathAndIS_CutCallback error", e.message

def getConstraintExpressionForBendersCut_JamVarsOnly(JammingGraph, flowVarValues, jamLocVarsDict, approxVarArg):
    firstSum = sum([flowVarValues[commod][p] for commod in instance.commodities.keys() for p in flowVars[commod].keys()])
    #print "firstSum", firstSum
    secondSum = 0
    for commod in instance.commodities.keys():
        for p in flowVars[commod].keys():
            if flowVarValues[commod][p] > FUZZ:
                secondSum += sum(flowVarValues[commod][p] * jamLocVarsDict[node] for node in getNodesThatCanInterdictEdgeSet(Paths[commod][p][1]))
    return firstSum - secondSum <= approxVarArg

def getConstraintExpressionFor_KI_Cut_JamVarsOnly_ArcBased(bestMIP_upperBound, flowVarValues, throughput, jamLocVarsDict, approxVarArg):
    lhsSum = 0
    for commod in instance.commodities.keys():
        for e in flowVars[commod].keys():
            if flowVarValues[commod][e] > FUZZ:
                lhsSum += sum(math.floor(-flowVarValues[commod][e]) * jamLocVarsDict[node] for node in jammersThatCanJamEdge[e])
    return lhsSum <= math.floor(bestMIP_upperBound - throughput)

def getExprForNew_Supervalid_TypeI_Cut_JamVarsOnly_ArcBased(bestMIP_upperBound, smallestSubprobThroughput, coefs, throughput, jamLocVarsDict):
    lhsSum = 0
    nodes = []
    print "getExprForNew_Supervalid_TypeI_Cut_JamVarsOnly_ArcBased", bestMIP_upperBound, smallestSubprobThroughput, throughput
    for node in instance.jamGraph.nodes():
        if (throughput + coefs[node]) <= smallestSubprobThroughput + FUZZ:
            lhsSum += jamLocVarsDict[node]
            nodes.append(node)
    return lhsSum >= 1, nodes

def OLD_getExprForNew_Supervalid_TypeI_Cut_JamVarsOnly_ArcBased(bestMIP_upperBound, smallestSubprobThroughput, coefs, throughput, jamLocVarsDict):
    lhsSum = 0
    for node in instance.jamGraph.nodes():
        if (smallestSubprobThroughput - throughput + coefs[node]) >= FUZZ:
            lhsSum += jamLocVarsDict[node]
    return lhsSum >= min(instance.jamBudget, get_Tightened_RHS_for_TypeICut(smallestSubprobThroughput, throughput, coefs))

def updateOldSVIs(myModel, bestMIP_upperBound, flowVarValues, throughput, jamLocVarsDict):
    None

def get_Tightened_RHS_for_TypeICut(smallestSubprobThroughput, throughput, coefs):
    coefsOrderedDict = OrderedDict(sorted(coefs.items(), key=lambda t: t[1]))
    mySum = throughput
    rhs = 0
    for node in instance.jamGraph.nodes():
        if mySum <= smallestSubprobThroughput:
            break
        else:
            rhs += 1
        mySum += coefsOrderedDict[node]
    return rhs
    
    
def addAndUpdateSupervalidCuts_Classic(myModel, bestMIP_upperBound, smallestSubprobThroughput, flowVarValues, throughput, jamLocVarsDict):
    global sviTypeI_Cuts, numTypeICuts
    # add new cut
    coefs = getCoefsForJamVars(flowVarValues)
    newSVI_TypeI, nodes = getExprForNew_Supervalid_TypeI_Cut_JamVarsOnly_ArcBased(bestMIP_upperBound, smallestSubprobThroughput, coefs, throughput, jamLocVarsDict)
    #print "throughput", throughput
    print "coefs", coefs
    print "newSVI_TypeI expr", numTypeICuts, newSVI_TypeI
    if numTypeICuts <= numTypeICutsLimit:
        print "add SVI cut"
        myModel.addConstr(newSVI_TypeI)
        if includeParetoOptimalBendersCuts:
            masterProbLP.addConstr(newSVI_TypeI)
        sviTypeI_NodesInLHS.append(nodes)
        sviTypeI_Cuts.append(newSVI_TypeI)
        numTypeICuts += 1
    # update old cuts
    #updateOldSVIs(myModel, bestMIP_upperBound, flowVarValues, throughput, jamLocVarsDict)

def getCoefsForJamVars(flowVarValues):
    coefs = {}
    for node in instance.jamGraph.nodes():
        coefs[node] = 0
    for commod in instance.commodities.keys():
        for e in flowVars[commod].keys():
            for node in jammersThatCanJamEdge[e]:
                coefs[node] -= flowVarValues[commod][e]
    return coefs

def getConstraintExpressionForBendersCut_JamVarsOnly_ArcBased(flowVarValues, throughput, jamLocVarsDict, approxVarArg):
    secondSum = 0
    for commod in instance.commodities.keys():
        for e in flowVars[commod].keys():
            if flowVarValues[commod][e] > FUZZ:
                secondSum += sum(flowVarValues[commod][e] * jamLocVarsDict[node] for node in jammersThatCanJamEdge[e])
    return throughput - secondSum <= approxVarArg
        
def getConstraintExpressionForBendersCut_JamAndEdgeInterdictVars(JammingGraph, flowVarValues, interdictVarsDict, approxVarArg):
    firstSum = sum([flowVarValues[commod][p] for commod in instance.commodities.keys() for p in flowVars[commod].keys()])
    secondSum = sum([flowVarValues[commod][p] * interdictVarsDict[e] for commod in instance.commodities.keys() for p in flowVars[commod].keys() for e in Paths[commod][p][1] if flowVarValues[commod][p] > FUZZ])
    #print "secondSum", secondSum
    return firstSum - secondSum <= approxVarArg
    #return 0 <= approxVarArg

def create_MaxWtIndepSet_Model(ConflictGraph, nodeWeights, interfModelType):
    global maxWtIndSetModel, inSet
    print "create_MaxWtIndepSet_Model", interfModelType
    maxWtIndSetModel = gurobipy.Model("max wt indep set")
    try:
        # Create variables
        inSet = dict([(node, maxWtIndSetModel.addVar(0, 1, vtype = gurobipy.GRB.BINARY, name="x_"+str(node))) 
                      for node in ConflictGraph.nodes()])
        maxWtIndSetModel.update() # Integrate new variables
        # Set objective
        maxWtIndSetModel.setObjective(sum([nodeWeights[node] * inSet[node] for node in ConflictGraph.nodes()]), gurobipy.GRB.MAXIMIZE)
        # adjacency constraints
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            [maxWtIndSetModel.addConstr( inSet[node] + inSet[adjNode] <= 1, "adjConstr_"+str(node)+","+str(adjNode))
             for node in ConflictGraph.nodes() for adjNode in ConflictGraph.neighbors(node) if adjNode != node]
        elif(interfModelType == 'simple-physical'):
            M = 100000
            [maxWtIndSetModel.addConstr(sum([ConflictGraph.edge[node][otherNode]['weight'] * inSet[otherNode] for otherNode in ConflictGraph.nodes() 
                                        if otherNode != node]) <= 1.0 * inSet[node] + M*(1.0 - inSet[node]),
                                    "snr-constr_"+str(node))
                                        for node in ConflictGraph.nodes()]
        maxWtIndSetModel.update() # integrate objective and constraints
        maxWtIndSetModel.setParam('OutputFlag', False ) #turn output off
        #maxWtIndSetModel.setParam('CliqueCuts', 2) #turn output off
        if(writeToFile):
            maxWtIndSetModel.write("/tmp/maxWtIndepSet_initial.lp")
    except gurobipy.GurobiError as e:
        print "maxWt error: ", str(e)

def getNodeToRoundUp(node0, node1, fracSoln, nodeWeights, roundingHeuristic):
    #print "getNodeToRoundUp", node0, node1
    if roundingHeuristic == 'random':
        denom =  nodeWeights[node0] *  fracSoln[node0] + nodeWeights[node1] * fracSoln[node1]
        #print "denom", denom, nodeWeights[node0], fracSoln[node0], nodeWeights[node1], fracSoln[node1]
        if denom != 0:
            prob = ( nodeWeights[node1] *  fracSoln[node1] ) / ( nodeWeights[node0] *  fracSoln[node0] + nodeWeights[node1] * fracSoln[node1])
        else:
            prob = 0.5
        #print "prob", prob
        if np.random.binomial(1, prob) == 0.0:
            return node0
        else:
            return node1
    elif roundingHeuristic == 'greedy':
        if nodeWeights[node0] *  fracSoln[node0] > nodeWeights[node1] * fracSoln[node1]:
            return node0
        else:
            return node1

def is01Fractional(value):
    if value > FUZZ and value < (1.0 - FUZZ):
        return True
    else:
        return False
    
def roundFractionalSoln(fracSoln, nodeWeights, roundingHeuristic = 'random'):
    #print "roundFractionalSoln"
    solnMap = {}
    for key in fracSoln.keys():
        solnMap[key] = fracSoln[key]
    for node in instance.interferenceGraph.nodes():
        for adjNode in instance.interferenceGraph.neighbors(node):
            if adjNode != node:
                #print "info", node, solnMap[node], adjNode, solnMap[adjNode]
                if (is01Fractional(solnMap[node]) or is01Fractional(solnMap[adjNode])) or solnMap[node] + solnMap[adjNode] > 1:
                    nodeToRound = getNodeToRoundUp(node, adjNode, fracSoln, nodeWeights, roundingHeuristic)
                    #print "nodeToRound", nodeToRound
                    if nodeToRound == node:
                        solnMap[node] = 1
                        solnMap[adjNode] = 0
                    elif nodeToRound == adjNode:
                        solnMap[adjNode] = 1
                        solnMap[node] = 0
    return solnMap

def modifyAndSolve_maxWtIndepSet(ConflictGraph, nodeWeights, interfModelType, iteration, rel_mip_gap_tol = 0.0001, setsToExclude = []):
    global numMaxWtIndSetProbsSolved
    numMaxWtIndSetProbsSolved += 1
    modifiedNodeWeights = {node : nodeWeights[node] for node in ConflictGraph.nodes()}
    try:
        for node in ConflictGraph.nodes():
            if modifiedNodeWeights[node] == 0.0:
                modifiedNodeWeights[node] = FUZZ
        for node in ConflictGraph.nodes():
            inSet[node].setAttr(gurobipy.GRB.attr.Obj, modifiedNodeWeights[node])
        maxWtIndSetModel.update() # integrate objective and constraints
        constrAdded = []
        for setToExclude in setsToExclude:
            constrAdded.append(maxWtIndSetModel.addConstr(sum([inSet[node] for node in setToExclude]) <= len(setToExclude) - 1, "cutoff"))
        maxWtIndSetModel.update()
        maxWtIndSetModel.Params.MIPGap = rel_mip_gap_tol
        maxWtIndSetModel.optimize() # this is telling Gurobi to solve it
        if(writeToFile):
            maxWtIndSetModel.write("/tmp/maxWtIndepSet_"+str(iteration) +".lp")
        solnCount = maxWtIndSetModel.SolCount
        numBinVars = maxWtIndSetModel.NumBinVars
        solnSet = []
        altObjVals = []
        if numBinVars > 0:
            for index in range(min(solnCount, 5)):
                maxWtIndSetModel.Params.SolutionNumber = index
                solnSet.append([key for key in inSet.keys() if inSet[key].Xn > 0.0001])
                altObjVals.append(sum([nodeWeights[key] * inSet[key].Xn for key in inSet.keys()]))
        else:
            solnMap = roundFractionalSoln({key : inSet[key].X for key in inSet.keys()}, nodeWeights, 'random')
            solnSet.append([key for key in solnMap.keys()])
            altObjVals.append(sum([nodeWeights[key] * solnMap[key] for key in solnMap.keys()]))
        for constr in constrAdded:
            maxWtIndSetModel.remove(constr)
    except gurobipy.GurobiError as e:
        print "maxWt error: ", str(e)
    return altObjVals, solnSet

def modifyAndSolve_maxWtIndepSet_ViaHeuristic(ConflictGraph, nodeWeights, interfModelType, iteration, rel_mip_gap_tol = 0.0001, setsToExclude = []):
    #print "modifyAndSolve_maxWtIndepSet_ViaHeuristic", iteration, "nodeWts", nodeWeights
    tempNodeWts = {node : nodeWeights[node] for node in ConflictGraph.nodes()}
    #print "tempNodeWts", tempNodeWts
    try:
        solnSet = []
        altObjVals = []
        for index in range(1): # number of times we run the heuristic
            soln = []
            sorted_tempNodeWts = sorted(tempNodeWts.items(), key=operator.itemgetter(1))
            #print "sorted_tempNodeWts", sorted_tempNodeWts
            while len(sorted_tempNodeWts) > 0:
                bestNode = sorted_tempNodeWts.pop(0)[0]
                print "bestNode", bestNode
                # does an edge exist between the existing nodes in the solution and the best node? if so, add it to the solution. if not, terminate the while loop
                print "n.all_neighbors(ConflictGraph, bestNode)", set(n.all_neighbors(ConflictGraph, bestNode))
                intersectionList = list(set(soln) & set(n.all_neighbors(ConflictGraph, bestNode)))
                print "intersectionList", intersectionList
                if len(intersectionList) == 0:
                    soln.append(bestNode)
                    print "added"
                print "soln", soln
            print "FINAL SOLN", soln
            solnSet.append(soln)
            altObjVals.append(sum([nodeWeights[node] for node in soln]))
        print "altObjVals", altObjVals
    except gurobipy.GurobiError as e:
        print "maxWt error: ", str(e)
    return altObjVals, solnSet

def maxWtIndepSet(ConflictGraph, nodeWeights, interfModelType, iteration, rel_mip_gap_tol = 0.0001):
    gurobiModel = gurobipy.Model("myLP")
    #print "maxWtIndSet", ConflictGraph.nodes()
    #print "nodeWeights", nodeWeights
    try:
        # Create variables
        #print "ConflictGraph nodes", ConflictGraph.nodes()
        inSet = dict([(node, gurobiModel.addVar(0, vtype = gurobipy.GRB.BINARY, name="x_"+str(node))) 
                      for node in ConflictGraph.nodes()])
        #print "created vars"
        gurobiModel.update() # Integrate new variables
        # Set objective
        
        gurobiModel.setObjective(sum([nodeWeights[node] * inSet[node] for node in ConflictGraph.nodes()]), gurobipy.GRB.MAXIMIZE)
        #print "created obj"
        # adjacency constraints
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            [gurobiModel.addConstr( inSet[node] + inSet[adjNode] <= 1, "adjConstr_"+str(node)+","+str(adjNode))
             for node in ConflictGraph.nodes() for adjNode in ConflictGraph.neighbors(node) if adjNode != node]
        elif(interfModelType == 'simple-physical'):
            M = 100000
            #print "nodes", G.nodes()
            #print "edges", G.edges(data=True)
            [gurobiModel.addConstr(sum([ConflictGraph.edge[node][otherNode]['weight'] * inSet[otherNode] for otherNode in ConflictGraph.nodes() 
                                        if otherNode != node]) <= 1.0 * inSet[node] + M*(1.0 - inSet[node]),
                                    "snr-constr_"+str(node))
                                        for node in ConflictGraph.nodes()]
        gurobiModel.update() # integrate objective and constraints
        #print "model created"
        gurobiModel.setParam('OutputFlag', False ) #turn output off
        gurobiModel.Params.MIPGap = rel_mip_gap_tol
        gurobiModel.optimize()
        if(writeToFile):
            gurobiModel.write("/tmp/maxWtIndepSet_"+str(iteration) +".lp")
        maxWeight = gurobiModel.objVal
        #print "maxWeight:", maxWeight
        selected = dict([(key, inSet[key].X) for key in inSet.keys() if inSet[key].X > 0.0001])
        #print "selected:", selected
    except gurobipy.GurobiError as e:
        print "maxWt error: ", str(e)
    return maxWeight, [key for key in selected.keys()]

def create_ShortestPathProb(G, edgeWeights, maxNumHops, commod, useVirtualEdges = True):
    global shortestPathModels, findPathflowVars
    print "maxNumHops", maxNumHops
    shortestPathModels[commod] = gurobipy.Model("leastCostPathLP_" + str(commod))
    try:
        # Create variables
        findPathflowVars[commod] = dict([(edge, shortestPathModels[commod].addVar(0, 1, vtype=gurobipy.GRB.BINARY, name="y_"+str(edge))) for edge in edgeTriples])
        dummyArc = shortestPathModels[commod].addVar(0, 1, vtype=gurobipy.GRB.BINARY, name="dummyArc")
        shortestPathModels[commod].update() # Integrate new variables
        # Set objective
        shortestPathModels[commod].setObjective(sum([edgeWeights[edge] * findPathflowVars[commod][edge] for edge in edgeTriples]) + 2.0 * dummyArc, gurobipy.GRB.MINIMIZE)
        # flow balance constraints
        for node in G.nodes():
            outEdgeTriples = []
            for edgeInfo in G.edges([node], data = True):
                if useVirtualEdges or edgeInfo[2]['edgeType'] == 'real':
                    outEdgeTriples.append((edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']))
            inEdgeTriples = []
            for edgeInfo in G.in_edges([node], data = True):
                if useVirtualEdges or edgeInfo[2]['edgeType'] == 'real':
                    inEdgeTriples.append((edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']))
            shortestPathModels[commod].addConstr(sum([findPathflowVars[commod][edge] for edge in outEdgeTriples])
                               - sum([findPathflowVars[commod][edge] for edge in inEdgeTriples])
                               + int(node == instance.commodities[commod]['odPair'][0]) * dummyArc
                               - int(node == instance.commodities[commod]['odPair'][1]) * dummyArc
                               == int(node == instance.commodities[commod]['odPair'][0]) - int(node == instance.commodities[commod]['odPair'][1]), 
                               "flowBal_"+str(node))
            shortestPathModels[commod].addConstr(sum([findPathflowVars[commod][edge] for edge in outEdgeTriples]) <= 1, "outFlowLeq1_"+str(node))
        # max hop constraint
        shortestPathModels[commod].addConstr(sum([findPathflowVars[commod][edge] for edge in edgeTriples]) <= instance.maxNumHops, "max_num_hops")
        shortestPathModels[commod].update() # integrate objective and constraints
        shortestPathModels[commod].setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            shortestPathModels[commod].write("/tmp/leastCostPath_initial_" + str(commod)+ ".lp")
    except gurobipy.GurobiError as e:
        print "findLeastCostPath error", str(e)
        
def modifyAndSolve_ShortestPathProb(G, edgeWeights, maxNumHops, commod, iteration, pathsToExclude = []):
    try:
        for edge in edgeTriples:
            findPathflowVars[commod][edge].setAttr(gurobipy.GRB.attr.Obj, edgeWeights[edge])
        shortestPathModels[commod].update() # integrate objective and constraints
        constrAdded = []
        for pathToExclude in pathsToExclude:
            constrAdded.append(shortestPathModels[commod].addConstr(sum([findPathflowVars[commod][edge] for edge in pathToExclude]) <= len(pathToExclude) - 1, "cutoff"))
        shortestPathModels[commod].update()
        shortestPathModels[commod].optimize()
        pathLength = shortestPathModels[commod].objVal
        #print "runtime:", shortestPathModels[commod].Runtime
        selected = dict([(key, findPathflowVars[commod][key].X) for key in findPathflowVars[commod].keys() if findPathflowVars[commod][key].X > FUZZ])
        for constr in constrAdded:
            shortestPathModels[commod].remove(constr)
    except gurobipy.GurobiError as e:
        print "findLeastCostPath error", str(e)
    return pathLength, [key for key in selected.keys()]

def findLeastCostPath(G, edgeWeights, maxNumHops, interdictedEdges, commod, iteration):
    #print "findLeastCostPath", iteration, interdictedEdges, edgeWeights
    #print "edges", G.edges()
    gurobiModel = gurobipy.Model("leastCostPathLP")
    try:
        # Create variables
        findPathflowVarsLocal = dict([(edge, gurobiModel.addVar(0, 1, vtype=gurobipy.GRB.BINARY, name="y_"+str(edge))) for edge in edgeTriples])
        dummyArc = gurobiModel.addVar(0, 1, vtype=gurobipy.GRB.BINARY, name="dummyArc")
        gurobiModel.update() # Integrate new variables
        #print "edgeWts", iteration, edgeWeights
        # Set objective
        #print "sum", sum([edgeWeights[edge] * findPathflowVars[edge] for edge in edgeTriples])
        gurobiModel.setObjective(sum([edgeWeights[edge] * findPathflowVarsLocal[edge] for edge in edgeTriples]) + 2.0 * dummyArc, gurobipy.GRB.MINIMIZE)
        # flow balance constraints
        for node in G.nodes():
            outEdgeTriples = [(edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']) for edgeInfo in G.edges([node], data = True)]
            inEdgeTriples = [(edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']) for edgeInfo in G.in_edges([node], data = True)]
            #print "node", node, commod, commodities[commod]['odPair'][0], commodities[commod]['odPair'][1], outEdgeTriples, inEdgeTriples
            gurobiModel.addConstr(sum([findPathflowVarsLocal[edge] for edge in outEdgeTriples])
                               - sum([findPathflowVarsLocal[edge] for edge in inEdgeTriples])
                               + int(node == instance.commodities[commod]['odPair'][0]) * dummyArc
                               - int(node == instance.commodities[commod]['odPair'][1]) * dummyArc
                               == int(node == instance.commodities[commod]['odPair'][0]) - int(node == instance.commodities[commod]['odPair'][1]), 
                               "flowBal_"+str(node))
        # max hop constraint
        gurobiModel.addConstr(sum([findPathflowVarsLocal[edge] for edge in edgeTriples]) <= instance.maxNumHops, "max_num_hops")
        # constraint: cannot use interdicted edges
        #print "interdictedEdges before constraints", interdictedEdges
        for edgeTriple in interdictedEdges:
            gurobiModel.addConstr(findPathflowVarsLocal[edgeTriple] <= 0, "cannot use " + str(edgeTriple))
        #print "myListInterdictedNotUse", myListInterdictedNotUse
        gurobiModel.update() # integrate objective and constraints
        gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/leastCostPath2_"+str(iteration)+"_"+ str(commod)+ ".lp")
        gurobiModel.optimize()
        pathLength = gurobiModel.objVal
        #print "dummyArc flow:", dummyArc.X
        selected = dict([(key, findPathflowVarsLocal[key].X) for key in findPathflowVarsLocal.keys() if findPathflowVarsLocal[key].X > 0.0001])
        #print "arcs selected:", selected
    except gurobipy.GurobiError as e:
        print "findLeastCostPath error", str(e)
    return pathLength, [key for key in selected.keys()]

def createThroughputModel(G, commodities, ISets, jammingGraph, interferenceModelType, jamType='continuous', interdictModelType = 'cormican'):
    #if(jamType == 'continuous' and interdictModelType == 'regular'):
    #    createThroughputModel_ContinuousJamming_Regular(G, commodities, ISets, jammingGraph, interferenceModelType)
    if(jamType == 'continuous' and interdictModelType == 'cormican'):
        createThroughputModel_ContinuousJamming_Cormican(G, commodities, ISets, jammingGraph, interferenceModelType)

# def createThroughputModel_ContinuousJamming_Regular(G, commodities, ISets, jammingGraph, interfModelType):
#     global gurobiModel, iSetUsage, capConstraints, usageConstr, jamConstrs
#     numISets = len(ISets)
#     gurobiModel = gurobipy.Model("myLP")
#     isJammedList = {}
#     for edge in G.edges():
#         if(sum([int(isEdgeJammedByJammer_Protocol(G, edge, node, interfModelType))*float(jammingGraph.node[node]['selected']) 
#                 for node in jammingGraph.nodes()]) >= 1):
#             isJammedList[edge] = 1
#         else:
#             isJammedList[edge] = 0
#     print "createThroughputModel_ContinuousJamming_Regular", "isJammed", [isJammedList[i] for i in isJammedList.keys() if isJammedList[i] > 0.0001]
#     try:
#         # Create variables
#         flowVars = dict([(commodity, dict([(edge, gurobiModel.addVar(0, name="flow_"+str(edge[0])+","+str(edge[1]))) for edge in G.edges()]))
#             for commodity in commodities.keys()])
#         iSetUsage = [gurobiModel.addVar(0, 1, name="lamda_"+str(k)) for k in range(numISets)]
#         gurobiModel.update() # Integrate new variables
#         # Set objective
#         gurobiModel.setObjective(sum([flowVars[commodity][edge] for edge in G.edges([commodities[commodity]['odPair'][0]]) for commodity in commodities.keys()]), 
#                                  gurobipy.GRB.MAXIMIZE)
#         # flow balance constraints
#         #print "src", source, "dest", destination
#         #print "nodes: ", [node for node in G.nodes() if node != source and node != destination]
#         [gurobiModel.addConstr(sum([flowVars[commodity][edge] for edge in G.edges([node])])
#                                == sum([flowVars[commodity][edge] for edge in G.in_edges([node])]), "flowBal_"+str(node))
#                                     for node in G.nodes() if (node != commodities[commodity]['odPair'][0] and node != commodities[commodity]['odPair'][1]) 
#                                     for commodity in commodities.keys()]
#         # flow into source
#         [gurobiModel.addConstr(sum([flowVars[commodity][edge] for edge in G.in_edges([commodities[commodity]['odPair'][0]])]) == 0.0, "flowIntoSource") 
#             for commodity in commodities.keys()]
#         # flow out of destination
#         [gurobiModel.addConstr(sum([flowVars[commodity][edge] for edge in G.out_edges([commodities[commodity]['odPair'][1]])]) == 0.0, "flowOutOfDest")
#             for commodity in commodities.keys()]
#         # capacity
#         capConstraints = dict([(edge, gurobiModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) <= 
#                                sum([iSetUsage[k] * int(edge in ISets[k]) * G.edge[edge[0]][edge[1]]['capacity'] for k in range(numISets)]), 
#                                    "capacity_"+str(edge[0])+","+str(edge[1]))) for edge in G.edges()])
#         # usage
#         usageConstr = gurobiModel.addConstr(sum([iSetUsage[k] for k in range(numISets)]) <= 1, "iSetUsage")
#         # jamming
#         if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
#             #jamConstrs = {}
#             #for node in jammingGraph.nodes():
#                 #jamConstrs[node] = dict([(edge, gurobiModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) <= 
#                     #G.edge[edge[0]][edge[1]]['capacity'] * 
#                     #(1 - float(jammingGraph.node[node]['selected'])), "jammed_"+str(edge)+","+str(node)))
#                     #           for edge in G.edges() if int(isEdgeJammedByJammer_Protocol(G, edge, node, modelType)) == 1])
#             jamConstrs = dict([(edge, gurobiModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) <= 
#                 G.edge[edge[0]][edge[1]]['capacity'] * 
#                 (1 - isJammedList[edge]), "jammed_"+str(edge)))
#                            for edge in G.edges()])
#             
#         elif(interfModelType == 'simple-physical'):
#             raise Exception("not implemented yet")
#         
#         gurobiModel.update() # integrate objective and constraints
#         gurobiModel.setParam('OutputFlag', False ) #turn output off
#         if(writeToFile):
#             gurobiModel.write("/home/hmedal/Documents/Temp/JainInitial_"+interfModelType + ".lp")
#     except gurobipy.GurobiError as e:
#         print "createThroughputModel_ContinuousJamming", str(e)

def resetObjectiveForThroughputProblem(G, jamLocs, commodities):
    gurobiThroughputModel.setObjective(getObjective_ThroughputModel_Cormican(G, jamLocs, commodities))
    
def resetObjectiveForThroughputProblem_JamVector(G, jammingVarValuesMap, commodities):
    gurobiThroughputModel.setObjective(getObjective_ThroughputModel_Cormican_JamVector(G, jammingVarValuesMap, commodities))

def getObjective_ThroughputModel_Cormican_JamVector(G, jammingVarValuesMap, commodities):
    sum1 = sum([returnFlow[commod] for commod in commodities.keys()])
    sum2 = 0
    for commodity in commodities.keys():
        for edge in edgeTriples:
            for node in getJammersThatCanJamEdge(G, instance.jamGraph, edge, interfModelType):
                sum2 += flowVars[commodity][edge] * jammingVarValuesMap[node]
    #print "getObjective_ThroughputModel_Cormican_JamVector", sum2
    return sum1 - sum2

def getObjective_ThroughputModel_Cormican(G, jamLocs, commodities):
    isJammedList = {}
    for edge in edgeTriples:
        if(sum([int(isEdgeJammedByJammer_Protocol(G, edge, node, interfModelType)) for node in jamLocs]) >= 1):
            isJammedList[edge] = 1
        else:
            isJammedList[edge] = 0
    #print "isJammed", [(i,isJammedList[i]) for i in isJammedList.keys() if isJammedList[i] > 0.0001]
    sum1 = sum([returnFlow[commod] for commod in commodities.keys()])
    #print "isJammedList", isJammedList
    sum2 = sum([flowVars[commodity][edge] * isJammedList[edge] for edge in edgeTriples for commodity in commodities.keys()])
    return sum1 - sum2

def addConstraintForThroughputValue(G, jamLocs, commodities, throughput):
    global throughputEqualsValueConstr
    throughputEqualsValueConstr = gurobiThroughputModel.addConstr(getObjective_ThroughputModel_Cormican(G, jamLocs, commodities) == throughput, "throughputEqualsValue")

def getObjective_Penalty_Term(G, commodities, slackVars, slackVariableMultipliers):
    penaltyTerm = 0
    constrName = 'flowBal'
    for commod in commodities.keys():
        for node in G.nodes():
            penaltyTerm -= slackVariableMultipliers['-'][constrName][commod][node]*slackVars['-'][constrName][commod][node]
            penaltyTerm += slackVariableMultipliers['+'][constrName][commod][node]*slackVars['+'][constrName][commod][node]
    constrName = 'cap'
    for edge in edgeTriples:
        penaltyTerm -= slackVariableMultipliers['-'][constrName][edge]*slackVars['-'][constrName][edge]
        penaltyTerm += slackVariableMultipliers['+'][constrName][edge]*slackVars['+'][constrName][edge]
    constrName = 'commodDemand'
    for commod in commodities.keys():
        penaltyTerm -= slackVariableMultipliers['-'][constrName][commod]*slackVars['-'][constrName][commod]
        penaltyTerm += slackVariableMultipliers['+'][constrName][commod]*slackVars['+'][constrName][commod]
    constrName = 'usage'
    penaltyTerm -= slackVariableMultipliers['-'][constrName]*slackVars['-'][constrName]
    penaltyTerm += slackVariableMultipliers['+'][constrName]*slackVars['+'][constrName]
    return penaltyTerm

def createSlackVariables(G, commodities):
    slackVars = {}
    for sign in ['-', '+']:
        slackVars[sign] = {}
        constrName = 'flowBal'
        slackVars[sign][constrName] = {}
        for commod in commodities.keys():
            slackVars[sign][constrName][commod] = {}
            for node in G.nodes():
                slackVars[sign][constrName][commod][node] = gurobiThroughputModel.addVar(0, float('inf'), name= "slackVars_"+sign+"_"+constrName+"_"+str(commod)+"_"+str(node))
        constrName = 'cap'
        slackVars[sign][constrName] = {}
        for edge in edgeTriples:
            slackVars[sign][constrName][edge] = gurobiThroughputModel.addVar(0, float('inf'), name= "slackVars_"+sign+"_"+constrName+"_"+str(edge))
        constrName = 'commodDemand'
        slackVars[sign][constrName] = {}
        for commod in commodities.keys():
            slackVars[sign][constrName][commod] = gurobiThroughputModel.addVar(0, float('inf'), name= "slackVars_"+sign+"_"+constrName+"_"+str(commod))
        constrName = 'usage'
        slackVars[sign][constrName] = gurobiThroughputModel.addVar(0, float('inf'), name= "slackVars_"+sign+"_"+constrName)
    return slackVars

def createSlackVariableMultipliers(G, commodities, defaultValue, halfWidth = 0.001):
    slackVarMultipliers = {}
    addend = {'-' : -halfWidth, '+' : halfWidth}
    for sign in ['-', '+']:
        slackVarMultipliers[sign] = {}
        constrName = 'flowBal'
        slackVarMultipliers[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarMultipliers[sign][constrName][commod] = {}
            for node in G.nodes():
                slackVarMultipliers[sign][constrName][commod][node] = defaultValue + addend[sign]
        constrName = 'cap'
        slackVarMultipliers[sign][constrName] = {}
        for edge in edgeTriples:
            slackVarMultipliers[sign][constrName][edge] = defaultValue + addend[sign]
        constrName = 'commodDemand'
        slackVarMultipliers[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarMultipliers[sign][constrName][commod] = defaultValue + addend[sign]
        constrName = 'usage'
        slackVarMultipliers[sign][constrName] = defaultValue + addend[sign]
    return slackVarMultipliers

def createSlackVariableUBs(G, commodities, defaultValue):
    slackVarUBs = {}
    for sign in ['-', '+']:
        slackVarUBs[sign] = {}
        constrName = 'flowBal'
        slackVarUBs[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarUBs[sign][constrName][commod] = {}
            for node in G.nodes():
                slackVarUBs[sign][constrName][commod][node] = defaultValue
        constrName = 'cap'
        slackVarUBs[sign][constrName] = {}
        for edge in edgeTriples:
            slackVarUBs[sign][constrName][edge] = defaultValue
        constrName = 'commodDemand'
        slackVarUBs[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarUBs[sign][constrName][commod] = defaultValue
        constrName = 'usage'
        slackVarUBs[sign][constrName] = defaultValue
    return slackVarUBs
    
def createThroughputModel_ContinuousJamming_Cormican(G, commodities, ISets, jamLocs, interfModelType, stabilize = False):
    global gurobiThroughputModel, iSetUsage, flowBalanceConstrs, capConstraints, usageConstr, jamConstrs, flowVars, commodDemConstr, returnFlow, slackVars, slackUBConstr, slackVariableMultipliers, slackVarUBs
    numISets = len(ISets)
    gurobiThroughputModel = gurobipy.Model("createThroughputModel_ContinuousJamming_Cormican")
    try:
        ## CREATE VARIABLES
        # Flow variables
        flowVars = dict([
            (commodity,
            dict([(edge, gurobiThroughputModel.addVar(0, name="flow_"+str(edge[0])+","+str(edge[1]))) for edge in edgeTriples])
            )
            for commodity in commodities.keys()])
        dummyVar = gurobiThroughputModel.addVar(0, obj = 0.0, name="dummyVar")
        
        # Slack variables
        if stabilize:
            slackVars = createSlackVariables(G, commodities)
        
        # Other
        returnFlow = dict([(commod, gurobiThroughputModel.addVar(0, name="returnFlow_"+str(commod))) for commod in commodities.keys()])
        iSetUsage = [gurobiThroughputModel.addVar(0, 1, name="lamda_"+str(k)) for k in range(numISets)]
        gurobiThroughputModel.update() # Integrate new variables
        
        ## SET OBJECTIVE
        regularTerm = getObjective_ThroughputModel_Cormican(G, jamLocs, commodities)
        #print "regularTerm", regularTerm
        if stabilize:
            if slackVariableMultipliers is None:
                slackVariableMultipliers = createSlackVariableMultipliers(G, commodities, 0.0)
            penaltyTerm = getObjective_Penalty_Term(G, commodities, slackVars, slackVariableMultipliers)
            gurobiThroughputModel.setObjective(regularTerm + penaltyTerm, gurobipy.GRB.MAXIMIZE) #####  NEED TO CHANGE THIS WHEN FINISHED (add penalty term)
        else:
            gurobiThroughputModel.setObjective(regularTerm, gurobipy.GRB.MAXIMIZE)
        
        ## CREATE CONSTRAINTS
        # Flow balance constraints
        flowBalanceConstrs = {}
        for commod in commodities.keys():
            flowBalanceConstrs[commod] = {}
            for node in G.nodes():
                indicator = -int(node is commodities[commod]['odPair'][0]) * returnFlow[commod] + int((node is commodities[commod]['odPair'][1])) * returnFlow[commod]
                sumOut = 0
                for edgeInfo in G.edges([node], data = True):
                    edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                    sumOut += flowVars[commod][edgeTriple]
                sumIn = 0
                for edgeInfo in G.in_edges([node], data = True):
                    edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                    sumIn += flowVars[commod][edgeTriple]
                if stabilize is False:
                    flowBalanceConstrs[commod][node] = gurobiThroughputModel.addConstr(indicator + sumOut == sumIn, "flowBal_"+str(node))
                else:
                    flowBalanceConstrs[commod][node] = gurobiThroughputModel.addConstr(indicator + sumOut -slackVars['-']['flowBal'][commod][node] + \
                                                                                   slackVars['+']['flowBal'][commod][node] == sumIn, "flowBal_"+str(node))
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
            # Capacity
            capConstraints = {}
            for edge in edgeTriples:
                if stabilize is False:
                    capConstraints[edge] = gurobiThroughputModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) <= 
                                       sum([iSetUsage[k] * int(edge in ISets[k]) * G.edge[edge[0]][edge[1]]['capacity'] for k in range(numISets)]), 
                                           "capacity_"+str(edge[0])+","+str(edge[1]))
                else:
                    capConstraints[edge] = gurobiThroughputModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) -slackVars['-']['cap'][edge] + \
                                                                                       slackVars['+']['cap'][edge] <= 
                                       sum([iSetUsage[k] * int(edge in ISets[k]) * G.edge[edge[0]][edge[1]]['capacity'] for k in range(numISets)]), 
                                           "capacity_"+str(edge[0])+","+str(edge[1]))
            # Usage
            if numISets < 1:
                if stabilize is False:
                    usageConstr = gurobiThroughputModel.addConstr(dummyVar <= 1, "iSetUsage")
                else:
                    usageConstr = gurobiThroughputModel.addConstr(dummyVar - slackVars['-']['usage'] + slackVars['+']['usage'] <= 1, "iSetUsage")
            else:
                if stabilize is False:
                    usageConstr = gurobiThroughputModel.addConstr(sum([iSetUsage[k] for k in range(numISets)]) <= 1, "iSetUsage")
                else:
                    usageConstr = gurobiThroughputModel.addConstr(sum([iSetUsage[k] for k in range(numISets)]) -slackVars['-']['usage'] + slackVars['+']['usage'] <= 1, "iSetUsage")
        elif(interfModelType == 'none'):
            capConstraints = {}
            for edge in edgeTriples:
                name = "capacity_"+str(edge[0])+","+str(edge[1])
                if stabilize is False:
                    capConstraints[edge] = gurobiThroughputModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) <= \
                                                                           G.edge[edge[0]][edge[1]][0]['capacity'], name)
                else:
                    capConstraints[edge] = gurobiThroughputModel.addConstr(sum([flowVars[commodity][edge] for commodity in commodities.keys()]) -slackVars['-']['cap'][edge] + \
                                                                                   slackVars['+']['cap'][edge] <= G.edge[edge[0]][edge[1]][0]['capacity'], name)
        
        # Demand
        commodDemConstr = {}
        for commod in commodities.keys():
            if stabilize is False:
                commodDemConstr[commod] = gurobiThroughputModel.addConstr(returnFlow[commod] <= commodities[commod]['demand'], "demandForCommod_"+str(commod))
            else:
                commodDemConstr[commod] = gurobiThroughputModel.addConstr(returnFlow[commod] -slackVars['-']['commodDemand'][commod] + \
                                                                                   slackVars['+']['commodDemand'][commod] <= commodities[commod]['demand'], "demandForCommod_"+str(commod))
            
        # Slack var upper bound
        if stabilize is True:
            if slackVarUBs is None:
                slackVarUBs = createSlackVariableUBs(G, commodities, 0.0)
            slackUBConstr = {}
            for sign in ['-', '+']:
                slackUBConstr[sign] = {}
                constrName = 'flowBal'
                slackUBConstr[sign][constrName] = {}
                for commod in commodities.keys():
                    slackUBConstr[sign][constrName][commod] = {}
                    for node in G.nodes():
                        slackUBConstr[sign][constrName][commod][node] = gurobiThroughputModel.addConstr(slackVars[sign][constrName][commod][node] <= slackVarUBs[sign][constrName][commod][node], "slackUB_"+sign+"_"+constrName+"_"+str(commod)+"_"+str(node))
                constrName = 'cap'
                slackUBConstr[sign][constrName] = {}
                for edge in edgeTriples:
                    slackUBConstr[sign][constrName][edge] = gurobiThroughputModel.addConstr(slackVars[sign][constrName][edge] <= slackVarUBs[sign][constrName][edge], "slackUB_"+sign+"_"+constrName+"_"+str(edge))
                constrName = 'commodDemand'
                slackUBConstr[sign][constrName] = {}
                for commod in commodities.keys():
                    slackUBConstr[sign][constrName][commod] = gurobiThroughputModel.addConstr(slackVars[sign][constrName][commod] <= slackVarUBs[sign][constrName][commod], "slackUB_"+sign+"_"+constrName+"_"+str(commod))
                constrName = 'usage'
                slackUBConstr[sign][constrName] = gurobiThroughputModel.addConstr(slackVars[sign][constrName] <= slackVarUBs[sign][constrName], "slackUB_"+sign+"_"+constrName)
        gurobiThroughputModel.update() # integrate objective and constraints
        gurobiThroughputModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiThroughputModel.write("/tmp/ThroughputInitial_Cormican_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print "createThroughputModel_ContinuousJamming", str(e)

def updateSlackStuffInModel(G, commodities):
    global slackUBConstr, slackVars
    # update upper bounds
    for sign in ['-', '+']:
        constrName = 'flowBal'
        for commod in commodities.keys():
            for node in G.nodes():
                slackUBConstr[sign][constrName][commod][node].setAttr("rhs", slackVarUBs[sign][constrName][commod][node])
        constrName = 'cap'
        for edge in edgeTriples:
            slackUBConstr[sign][constrName][edge].setAttr("rhs", slackVarUBs[sign][constrName][edge])
        constrName = 'commodDemand'
        for commod in commodities.keys():
            slackUBConstr[sign][constrName][commod].setAttr("rhs", slackVarUBs[sign][constrName][commod])
        constrName = 'usage'
        slackUBConstr[sign][constrName].setAttr("rhs", slackVarUBs[sign][constrName])
    gurobiThroughputModel.update()
    # update objective coefficients
    constrName = 'flowBal'
    for commod in commodities.keys():
        for node in G.nodes():
            slackVars['-'][constrName][commod][node].obj = slackVariableMultipliers['-'][constrName][commod][node]
            slackVars['+'][constrName][commod][node].obj = slackVariableMultipliers['+'][constrName][commod][node]
    constrName = 'cap'
    for edge in edgeTriples:
        slackVars['-'][constrName][edge].obj = slackVariableMultipliers['-'][constrName][edge]
        slackVars['+'][constrName][edge].obj = slackVariableMultipliers['+'][constrName][edge]
    constrName = 'commodDemand'
    for commod in commodities.keys():
        slackVars['-'][constrName][commod].obj = slackVariableMultipliers['-'][constrName][commod]
        slackVars['+'][constrName][commod].obj = slackVariableMultipliers['+'][constrName][commod]
    constrName = 'usage'
    slackVars['-'][constrName].obj = slackVariableMultipliers['-'][constrName]
    slackVars['+'][constrName].obj = slackVariableMultipliers['+'][constrName]
    gurobiThroughputModel.update()
                
def getEdgesJammedList(edgesJammed, jammingGraph = None):
    isJammedList = {}
    if(edgesJammed is not None):
        for edge in instance.CGraph.edges():
            if(edge in edgesJammed):
                isJammedList[edge] = 1
            else:
                isJammedList[edge] = 0
    else:
        if(jammingGraph is not None):
            for edge in instance.CGraph.edges():
                if(isEdgeJammedByJammers(instance.CGraph, jammingGraph, edge, jammingGraph.nodes(), interfModelType)):
                    isJammedList[edge] = 1
                else:
                    isJammedList[edge] = 0
    #print "isJammed", isJammedList
    return isJammedList

def createCapConstraints(G, commodities, numISets, dummyVar):
    global gurobiThroughputModel
    capConstraints = {}
    for edge in edgeTriples:
        flowVarsSet = [flowVars[commod][p] for commod in commodities.keys() for p in flowVars[commod].keys() if edge in Paths[commod][p][1]]
        #print "edge", G.edge[edge[0]][edge[1]]
        ISetVarsSet = [int(edge in ISets[k]) for k in range(numISets)]
        attrIndex = getAttrIndex(edge)
        ISetVarsExprs = [iSetUsage[k] * int(edge in ISets[k]) * G.edge[edge[0]][edge[1]][attrIndex]['capacity'] for k in range(numISets)]
        if(len(flowVarsSet) > 0):
            #print "if"
            capConstraints[edge] = gurobiThroughputModel.addConstr(sum(flowVarsSet) <= sum(ISetVarsExprs), "capacityReg_"+str(edge[0])+","+str(edge[1]))
        else:
            #print "else"
            if(sum(ISetVarsSet) > 0):
                #print "if2"
                capConstraints[edge] = gurobiThroughputModel.addConstr(dummyVar <= sum(ISetVarsExprs), "cap2_"+str(edge[0])+","+str(edge[1]))
            else:
                #print "else2", interfModelType
                if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
                    #print "if3"
                    capConstraints[edge] = gurobiThroughputModel.addConstr(dummyVar <= sum(ISetVarsExprs), "capacityEmpty_"+str(edge[0])+","+str(edge[1]))
                elif(interfModelType == 'none'):
                    print "none"
                    capConstraints[edge] = gurobiThroughputModel.addConstr(dummyVar <= G.edge[edge[0]][edge[1]][attrIndex]['capacity'], "capacityEmpty_"+str(edge[0])+","+str(edge[1]))
    return capConstraints

def createFlowVars(commodities):
    global gurobiThroughputModel
    flowVars = {}
    for commodity in commodities.keys():
        flowVars[commodity] = {}
        for pathInfo in Paths[commodity]:
            flowVars[commodity][pathInfo[0]] = gurobiThroughputModel.addVar(0, name="flow_"+str(commodity)+","+str(pathInfo[0]))
    return flowVars

def createThroughputModel_ContinuousJamming_PathBased_Cormican(G, edgesJammed, jamLocs, commodities, interfModelType, jammingGraph = None):
    global gurobiThroughputModel, iSetUsage, flowBalanceConstrs, capConstraints, commodDemConstr, usageConstr, jamConstrs, flowVars
    numISets = len(ISets)
    gurobiThroughputModel = gurobipy.Model("myLP")
    if modelType == 'cormican-path-based-jam-and-interdict-vars':
        isJammedList = getEdgesJammedList(edgesJammed, jammingGraph)
    try:
        # Create variables
        flowVars = createFlowVars(commodities)
        dummyVar = gurobiThroughputModel.addVar(0, obj = 0.0, name="dummyVar")
        dummyISetVar = gurobiThroughputModel.addVar(0, obj = 0.0, name="dummyISetVar")
        iSetUsage = [gurobiThroughputModel.addVar(0, 1, name="lamda_"+str(k)) for k in range(numISets)]
        gurobiThroughputModel.update() # Integrate new variables
        #print "initial flowVars", flowVars
        # Set objective
        if jamVarsType == 'jam-vars':
            flowSum = sum([flowVars[commod][p] for commod in commodities.keys() for p in flowVars[commod].keys()])
            print "flowSum", flowSum
            jammedFlowSum = 0
            for commod in commodities.keys():
                for p in flowVars[commod].keys():
                    jamLocsThatAreOnPath = list(set(jamLocs) & set(getNodesThatCanInterdictEdgeSet(Paths[commod][p][1])))
                    print "jamLocsThatAreOnPath", jamLocsThatAreOnPath
                    jammedFlowSum -= len(jamLocsThatAreOnPath) * flowVars[commod][p] 
            print "jammedFlowSum", jammedFlowSum
            print "totalSum", flowSum - jammedFlowSum
            gurobiThroughputModel.setObjective(flowSum - jammedFlowSum, gurobipy.GRB.MAXIMIZE)
        elif jamVarsType == 'jam-and-interdict-vars':
            gurobiThroughputModel.setObjective(sum([flowVars[commod][p] for commod in commodities.keys() for p in flowVars[commod].keys()])
                - sum([flowVars[commod][p] * isJammedList[edge] for commod in commodities.keys() for p in flowVars[commod].keys() for edge in Paths[commod][p][1]]), 
                                     gurobipy.GRB.MAXIMIZE)
        #constraints
        capConstraints = createCapConstraints(G, commodities, numISets, dummyVar)
        gurobiThroughputModel.update()
        # demand for each commodity
        commodDemConstr = {}
        for commod in commodities.keys():
            flowVarsSet = [flowVars[commod][p] for p in flowVars[commod].keys()]
            if(len(flowVarsSet) > 0):
                commodDemConstr[commod] = gurobiThroughputModel.addConstr(sum(flowVarsSet) <= commodities[commod]['demand'], "demandForCommod_"+str(commod))
            else:
                commodDemConstr[commod] = gurobiThroughputModel.addConstr(dummyVar <= commodities[commod]['demand'], "demandForCommod_"+str(commod))
        gurobiThroughputModel.update()
        # usage
        if(numISets > 0):
            usageConstr = gurobiThroughputModel.addConstr(sum([iSetUsage[k] for k in range(numISets)]) <= 1, "iSetUsage")
        else:
            if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
                usageConstr = gurobiThroughputModel.addConstr(dummyISetVar <= 1, "iSetUsage")
        jamConstrs = dict([(edge, gurobiThroughputModel.addConstr(0, "jamConstr")) for edge in edgeTriples])
        gurobiThroughputModel.update() # integrate objective and constraints
        gurobiThroughputModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiThroughputModel.write("/tmp/JainInitialCormicanPathBased2_"+str(interfModelType) + ".lp")
    except gurobipy.GurobiError as e:
        print "createThroughputModel_ContinuousJamming error", str(e)

def create_flowVarRedCostConstraints_arcBased(G, commodities, flowBalanceDuals):
    global gurobiModel
    flowVarsAsDual = {}
    if jamVarsType == 'jam-and-interdict-vars':
        for commod in commodities.keys():
            flowVarsAsDual[commod] = {}
            for edgeInfo in G.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                mySum = flowBalanceDuals[commod][edgeTriple[0]] - flowBalanceDuals[commod][edgeTriple[1]] + capDuals[edgeTriple] + edgeInterdict[edgeTriple]
                flowVarsAsDual[commod][edgeTriple] = gurobiModel.addConstr(mySum >= 0, "flowVarRedCost_"+str(commod)+","+str(edgeTriple))
    elif jamVarsType == 'jam-vars':
        for commod in commodities.keys():
            #print "flowBalanceDuals", flowBalanceDuals[commod]
            flowVarsAsDual[commod] = {}
            for edgeInfo in G.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                #print 'edgeTriple', edgeTriple, flowBalanceDuals[commod][edgeTriple[0]], flowBalanceDuals[commod][edgeTriple[1]]
                #print "sum", sum(jamPlaced[node] for node in jammersThatCanJamEdge[edgeTriple])
                mySum = flowBalanceDuals[commod][edgeTriple[0]] - flowBalanceDuals[commod][edgeTriple[1]] + capDuals[edgeTriple]
                mySum += sum(jamPlaced[node] for node in jammersThatCanJamEdge[edgeTriple])
                flowVarsAsDual[commod][edgeTriple] = gurobiModel.addConstr(mySum >= 0, "flowVarRedCost_"+str(commod)+","+str(edgeTriple))
    return flowVarsAsDual

def create_flowVarRedCostConstraints_Regular_ArcBased(G, commodities, flowBalanceDuals):
    global gurobiModel
    flowVarsAsDual = {}
    if jamVarsType == 'jam-vars':
        for commod in commodities.keys():
            #print "flowBalanceDuals", flowBalanceDuals[commod]
            flowVarsAsDual[commod] = {}
            for edgeInfo in G.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                #print 'edgeTriple', edgeTriple, flowBalanceDuals[commod][edgeTriple[0]], flowBalanceDuals[commod][edgeTriple[1]]
                #print "sum", sum(jamPlaced[node] for node in jammersThatCanJamEdge[edgeTriple])
                mySum = flowBalanceDuals[commod][edgeTriple[0]] - flowBalanceDuals[commod][edgeTriple[1]] + capDuals[edgeTriple]
                mySum += sum(jamUBDual[edgeTriple][node] for node in jammersThatCanJamEdge[edgeTriple])
                flowVarsAsDual[commod][edgeTriple] = gurobiModel.addConstr(mySum >= 0, "flowVarRedCost_"+str(commod)+","+str(edgeTriple))
    return flowVarsAsDual

def getCapDualsVars(G, useVirtualEdges = True):
    capDuals = {}
    if useVirtualEdges:
        edgesToIterateOver = G.edges(data = True)
    else:
        edgesToIterateOver = edgeTriplesInterdictable
    for edgeInfo in edgesToIterateOver:
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        capDuals[edgeTriple] = gurobiModel.addVar(0.0, name="capD_e"+str(edgeTriple))
    print "capDuals", capDuals
    
    return capDuals

def getEdgeInterdictVars(G, jammingGraph, useVirtualEdges = False):
    edgeInterdict = {}
    if useVirtualEdges:
        edgesToIterateOver = G.edges(data = True)
    else:
        edgesToIterateOver = edgeTriplesInterdictable
    for edgeInfo in edgesToIterateOver:
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        edgeInterdict[edgeTriple] = gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="edgeInterdict_e"+str(edgeTriple))
    jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])

def createFlowBalanceDuals(G, commodities):
    flowBalanceDuals = {}
    for commod in commodities.keys():
        print str(commodities.keys())
        flowBalanceDuals[commod] = {}
        for node in G.nodes():
            #flowBalanceDuals[commod][node] = gurobiModel.addVar(-1, 1, name="flowBalD_i"+str(node)+",commod"+str(commod))
            flowBalanceDuals[commod][node] = gurobiModel.addVar(name="flowBalD_i"+str(node)+",commod"+str(commod))
            print "node", str(node)
            print "commod", str(commod)
    print "flowBalanceDuals", flowBalanceDuals
    
    return flowBalanceDuals

def printFlowBalanceDualsValues():
    for commod in instance.commodities.keys():
        for node in instance.CGraph.nodes():
            print commod, node, flowBalanceDuals[commod][node].X
    return flowBalanceDuals

def getRegularObjFn_FullMIP(G, useVirtualEdges = True):
    capDuals = {}
    if useVirtualEdges:
        edgesToIterateOver = G.edges(data = True)
    else:
        edgesToIterateOver = edgeTriplesInterdictable
    objSum = 0
    for edgeInfo in edgesToIterateOver:
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        attrIndex = getAttrIndex(edgeTriple)
        for node in jammersThatCanJamEdge[edgeTriple]:
            objSum += G.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity'] * (jamUBDual[edgeTriple][node] - jamPlacedUBDual[edgeTriple][node])
    objSum += usageDual
    sumForDemMetDuals = sum([instance.commodities[commod]['demand'] * demForCommodDuals[commod] for commod in instance.commodities.keys()])
    objSum += sumForDemMetDuals
    return objSum

def getCormicanObjFn_FullMIP(G, interfModelType):
    #print "demand", instance.commodities[0]['demand']
    sumForDemMetDuals = sum([instance.commodities[commod]['demand'] * demForCommodDuals[commod] for commod in instance.commodities.keys()])
    print "sumForDemMetDuals", sumForDemMetDuals
    sys.exit()
    if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
        return usageDual + sumForDemMetDuals
    else:
        sumForCapDuals = 0
        for edgeInfo in G.edges(data = True):
            edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
            attrIndex = getAttrIndex(edgeTriple)
            sumForCapDuals += G.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity'] * capDuals[edgeTriple]
        sumForDemMetDuals = sum([instance.commodities[commod]['demand'] * demForCommodDuals[commod] for commod in instance.commodities.keys()])
        return sumForCapDuals + sumForDemMetDuals

def getJamPlacedUBDualVar(G):
    jamPlacedUBDual = {}
    for edgeInfo in G.edges(data = True):
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        jamPlacedUBDual[edgeTriple] = {}
        for node in jammersThatCanJamEdge[edgeTriple]:
            jamPlacedUBDual[edgeTriple][node] = gurobiModel.addVar(0.0, 1.0, name = "jamPlacedTimesUBDual_e_" + str(edgeTriple) + "_l_"+str(node))
    return jamPlacedUBDual

def getJamUBDual(G):
    jamUBDual = {}
    for edgeInfo in G.edges(data = True):
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        jamUBDual[edgeTriple] = {}
        for node in jammersThatCanJamEdge[edgeTriple]:
            jamUBDual[edgeTriple][node] = gurobiModel.addVar(0.0, 1.0, name = "jamUBDual_e_" + str(edgeTriple) + "_l_"+str(node))
    return jamUBDual

def create_BottleneckConstraints(G, numISets):
    isetUsageVarsAsDual = []
    for m in range(numISets):
        firstSum = 0
        for edgeInfo in G.edges(data = True):
            edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
            attrIndex = getAttrIndex(edgeTriple)
            if edgeTriple in ISets[m]:
                firstSum += -G.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity'] * capDuals[edgeTriple]
        isetUsageVarsAsDual.append(gurobiModel.addConstr(firstSum + usageDual >= 0.0, "bottleneck_"+str(m)))
    return isetUsageVarsAsDual

def createCapDualsUBConstraints(G):
    for edgeInfo in G.edges(data = True):
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        gurobiModel.addConstr(capDuals[edgeTriple] <= 1, "capDualsUB_" + str(edgeTriple)) 

def createDemandForCommodDuals():
    return dict([(commod, gurobiModel.addVar(0.0, name="demForCommodDual_commod"+str(commod))) for commod in instance.commodities.keys()])
    
def create_SingleLevelMIP_Regular_ArcBased(G, jammingGraph, commodities, ISets, interfModelType):
    global gurobiModel, demForCommodDuals, capDuals, usageDual, edgeInterdict, jamPlaced, flowBalanceDuals, jamUBDual, jamPlacedUBDual
    numISets = len(ISets)
    gurobiModel = gurobipy.Model("SingleLevelMIP_Cormican_ArcBased")
    try:
        # Create variables
        flowBalanceDuals = createFlowBalanceDuals(G, commodities)
        capDuals = getCapDualsVars(G)
        usageDual = gurobiModel.addVar(0, name="usageD")
        jamUBDual = getJamUBDual(G)
        demForCommodDuals = createDemandForCommodDuals()
        jamPlacedUBDual = getJamPlacedUBDualVar(G)
        jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])
        gurobiModel.update() # Integrate new variables
        # Set objective
        gurobiModel.setObjective(getRegularObjFn_FullMIP(G), gurobipy.GRB.MINIMIZE)
        gurobiModel.update()
        #constraints
        isetUsageVarsAsDual = create_BottleneckConstraints(G, numISets)
        flowVarsAsDual = create_flowVarRedCostConstraints_Regular_ArcBased(G, commodities, flowBalanceDuals)
        for commod in commodities.keys():
            src = commodities[commod]['odPair'][0]
            sink = commodities[commod]['odPair'][1]
            gurobiModel.addConstr(demForCommodDuals[commod] + flowBalanceDuals[commod][sink] - flowBalanceDuals[commod][src] >= 1, "flowBalDualSTDiff"+str(commod))
        gurobiModel.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlaced[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        createCapDualsUBConstraints(G)
        createLinearizationConstraints(G)
        gurobiModel.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/FullMIP_Regular_ArcBased_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print "create_SingleLevelMIP_Cormican_ArcBased", str(e)
        
def create_SingleLevelMIP_Cormican_ArcBased(G, jammingGraph, commodities, ISets, interfModelType): # creates model from equation (9)
    global gurobiModel, demForCommodDuals, capDuals, usageDual, edgeInterdict, jamPlaced, flowBalanceDuals
    numISets = len(ISets)
    gurobiModel = gurobipy.Model("SingleLevelMIP_Cormican_ArcBased")
    try:
        # Create variables
        flowBalanceDuals = createFlowBalanceDuals(G, commodities) # alpha

        capDuals = getCapDualsVars(G) #beta
        usageDual = gurobiModel.addVar(0, name="usageD") #gamma
        demForCommodDuals = createDemandForCommodDuals() # psi_m
        if jamVarsType == 'jam-and-interdict-vars': # don't use this
            edgeInterdict = getEdgeInterdictVars(G, jammingGraph)
        jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])
        print "jammingGraph.nodes()", jammingGraph.nodes()
        
        gurobiModel.update() # Integrate new variables
        # Set objective
        gurobiModel.setObjective(getCormicanObjFn_FullMIP(G, interfModelType), gurobipy.GRB.MINIMIZE)
        gurobiModel.update()
        #constraints       
        flowVarsAsDual = create_flowVarRedCostConstraints_arcBased(G, commodities, flowBalanceDuals) #4b
        for commod in commodities.keys():
            gurobiModel.addConstr(demForCommodDuals[commod] + flowBalanceDuals[commod][commodities[commod]['odPair'][1]] - flowBalanceDuals[commod][commodities[commod]['odPair'][0]] >= 1, "flowBalDualSTDiff"+str(commod))
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')): # use 802.11 MAC
            isetUsageVarsAsDual = create_BottleneckConstraints(G, numISets) #4d
        gurobiModel.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlaced[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget") #2b
        if jamVarsType == 'jam-and-interdict-vars':
            for edgeInfo in G.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                gurobiModel.addConstr(edgeInterdict[edgeTriple] <= sum([jamPlaced[l] for l in jammingGraph.nodes() 
                        if (int(isEdgeJammedByJammer_Protocol(G, edgeTriple, l, interfModelType)) == 1)]), "arcInterdict_e"+str(edgeTriple))
        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
            createCapDualsUBConstraints(G) #may not need to do this anymore
        gurobiModel.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/FullMIPCormican_ArcBased_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print "create_SingleLevelMIP_Cormican_ArcBased", str(e)
        
def create_BendersMasterProb_Cormican_ArcBased(G, jammingGraph, commodities, ISets, interfModelType):
    global gurobiModel, demForCommodDuals, capDuals, usageDual, edgeInterdict, jamPlaced, approxVar
    gurobiModel = gurobipy.Model("BendersMasterProb_Cormican_ArcBased")
    try:
        # Create variables
        jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jam_l"+str(n[0])+","+str(n[1]))) for n in jammingGraph.nodes()])
        approxVar = gurobiModel.addVar(0, vtype=gurobipy.GRB.CONTINUOUS, name="approxVar")
        gurobiModel.update() # Integrate new variables
        # Set objective
        gurobiModel.setObjective(approxVar, gurobipy.GRB.MINIMIZE)
        gurobiModel.update()
        #constraints       
        gurobiModel.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlaced[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        gurobiModel.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/jamMasterProb_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print "create_BendersMasterProb_Cormican_ArcBased", str(e)
        
def create_BendersMasterProb_LP_Cormican_ArcBased(G, jammingGraph, commodities, ISets, interfModelType):
    global masterProbLP, jamPlacedLP
    masterProbLP = gurobipy.Model("BendersMasterProb_Cormican_ArcBased")
    try:
        # Create variables
        jamPlacedLP = dict([(n, masterProbLP.addVar(0, 1, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])
        approxVarLP = masterProbLP.addVar(0, vtype=gurobipy.GRB.CONTINUOUS, name="approxVar")
        masterProbLP.update() # Integrate new variables
        # Set objective
        masterProbLP.setObjective(approxVarLP, gurobipy.GRB.MINIMIZE)
        masterProbLP.update()
        #constraints       
        masterProbLP.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlacedLP[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        masterProbLP.update() # integrate objective and constraints
        masterProbLP.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            masterProbLP.write("/tmp/jamMasterProb_LP_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print "create_BendersMasterProb_Cormican_ArcBased", str(e)

def create_SingleLevelMIP_Cormican_PathBased(G, jammingGraph, commodities, Paths, ISets, interfModelType):
    global gurobiModel, usageConstr, jamConstrs, capDuals, demForCommodDuals, capDualsUB_duals, usageDual, edgeInterdict, jamPlaced, flowVarsAsDual, isetUsageVarsAsDual
    numISets = len(ISets)
    gurobiModel = gurobipy.Model("SingleLevelMIP_Cormican_PathBased")
    try:
        print "create variables"
        # Create variables
        capDuals = getCapDualsVars(G)
        usageDual = gurobiModel.addVar(0, name="usageD")
        demForCommodDuals = createDemandForCommodDuals()
        #binary variables
        if jamVarsType == 'jam-and-interdict-vars':
            edgeInterdict = getEdgeInterdictVars(G, jammingGraph)
        jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])
        gurobiModel.update() # Integrate new variables
        # Set objective
        print "create objective"
        gurobiModel.setObjective(getCormicanObjFn_FullMIP(G), gurobipy.GRB.MINIMIZE)
        gurobiModel.update()
        #constraints
        print "create constraints"
        if jamVarsType == 'jam-and-interdict-vars':
            flowVarsAsDual = {}
            for commod in commodities.keys():
                flowVarsAsDual[commod] = {}
                for pathInfo in Paths[commod]:
                    capDualsSum = sum([capDuals[e] for e in pathInfo[1]])
                    edgeInterdictsSum = sum([edgeInterdict[e] for e in pathInfo[1]])
                    flowVarsAsDual[commod][pathInfo[0]] = gurobiModel.addConstr(capDualsSum + edgeInterdictsSum + demForCommodDuals[commod] >= 1, 
                                                                                "flowVarRedCost_"+str(commod)+","+str(pathInfo[0]))
        elif jamVarsType == 'jam-vars':
            flowVarsAsDual = {}
            for commod in commodities.keys():
                flowVarsAsDual[commod] = {}
                for pathInfo in Paths[commod]:
                    capDualsSum = sum([capDuals[e] for e in pathInfo[1]])
                    nodesThatCanInterdictPath = getNodesThatCanInterdictEdgeSet(pathInfo[1])
                    nodeInterdictSum = sum([jamPlaced[node] for node in nodesThatCanInterdictPath])
                    flowVarsAsDual[commod][pathInfo[0]] = gurobiModel.addConstr(capDualsSum + nodeInterdictSum + demForCommodDuals[commod] >= 1, 
                                   "flowVarRedCost_"+str(commod)+","+str(pathInfo[0]))
        isetUsageVarsAsDual = create_BottleneckConstraints(G, numISets)
        gurobiModel.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlaced[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        if jamVarsType == 'jam-and-interdict-vars':
            for edgeInfo in G.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                gurobiModel.addConstr(edgeInterdict[edgeTriple] <= sum([jamPlaced[l] for l in jammingGraph.nodes() 
                        if (int(isEdgeJammedByJammer_Protocol(G, edgeTriple, l, interfModelType)) == 1)]), "arcInterdict_e"+str(edgeTriple))
        createCapDualsUBConstraints(G)
        gurobiModel.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/FullMIPCormican_PathBased_"+interfModelType + "_" + jamVarsType +".lp")
    except gurobipy.GurobiError as e:
        print "ERROR: create_SingleLevelMIP_Cormican_PathBased", str(e)

def createLinearizationConstraints(G):
    for edgeInfo in G.edges(data = True):
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        for node in jammersThatCanJamEdge[edgeTriple]:
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] <= jamUBDual[edgeTriple][node], "lin1"+str(edgeTriple)+","+str(node))
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] <= jamPlaced[node], "lin2"+str(edgeTriple)+","+str(node))
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] >= jamUBDual[edgeTriple][node] - 1 + jamPlaced[node], "lin3"+str(edgeTriple)+","+str(node))
            
def createLinearizedVariables(G):
    for edgeInfo in G.edges(data = True):
        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
        for node in jammersThatCanJamEdge[edgeTriple]:
            print "values", edgeTriple, node, jamPlaced[node], jamUBDual[edgeTriple][node], jamPlacedUBDual[edgeTriple][node]
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] <= jamUBDual[edgeTriple][node], "lin1"+str(edgeTriple)+","+str(node))
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] <= jamPlaced[node], "lin2"+str(edgeTriple)+","+str(node))
            gurobiModel.addConstr(jamPlacedUBDual[edgeTriple][node] <= jamUBDual[edgeTriple][node] - 1 + jamPlaced[node], "lin3"+str(edgeTriple)+","+str(node))
                
def create_SingleLevelMIP_Regular_PathBased(G, jammingGraph, commodities, Paths, ISets, interfModelType):
    global gurobiModel, demForCommodDuals, usageConstr, jamConstrs, capDuals, capDualsUB_duals, usageDual, edgeInterdict, jamPlaced, flowVarsAsDual, isetUsageVarsAsDual, jamUBDual, jamPlacedUBDual
    numISets = len(ISets)
    gurobiModel = gurobipy.Model("SingleLevelMIP_Regular_PathBased")
    try:
        # Create variables
        capDuals = getCapDualsVars(G)
        usageDual = gurobiModel.addVar(0, name="usageD")
        demForCommodDuals = createDemandForCommodDuals()
        jamUBDual = getJamUBDual(G)
        jamPlacedUBDual = getJamPlacedUBDualVar(G)
        jamPlaced = dict([(n, gurobiModel.addVar(vtype=gurobipy.GRB.BINARY, name="jamPlaced_l"+str(n))) for n in jammingGraph.nodes()])
        gurobiModel.update() # Integrate new variables
        # Set objective
        gurobiModel.setObjective(getRegularObjFn_FullMIP(G), gurobipy.GRB.MINIMIZE)
        gurobiModel.update()
        #constraints
        if jamVarsType == 'jam-vars':
            flowVarsAsDual = {}
            for commod in commodities.keys():
                flowVarsAsDual[commod] = {}
                for pathInfo in Paths[commod]:
                    nodesThatCanInterdictPath = getNodesThatCanInterdictEdgeSet(pathInfo[1])
                    firstSum = sum([capDuals[e] for e in pathInfo[1]])
                    secondSum = sum([jamPlaced[node] for node in nodesThatCanInterdictPath])
                    for edgeInfo in G.edges(data = True):
                        edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                        for node in jammersThatCanJamEdge[edgeTriple]:
                            secondSum += jamUBDual[edgeTriple][node]
                    flowVarsAsDual[commod][pathInfo[0]] = gurobiModel.addConstr(firstSum + secondSum + demForCommodDuals[commod] >= 1, 
                                                                                "flowVarRedCost_"+str(commod)+","+str(pathInfo[0]))
        else:
            raise Exception("jamVarsType " + jamVarsType + " is not appropriate for the regular path-based model")
        isetUsageVarsAsDual = create_BottleneckConstraints(G, numISets)
        gurobiModel.addConstr(sum([jammingGraph.node[l]['cost'] * jamPlaced[l] for l in jammingGraph.nodes()]) <= instance.jamBudget, "interdictBudget")
        createCapDualsUBConstraints(G)
        createLinearizationConstraints(G)
        gurobiModel.update() # integrate objective and constraints
        #gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/tmp/FullMIPCormican_PathBased_"+interfModelType + "_" + jamVarsType +".lp")
    except gurobipy.GurobiError as e:
        print "create_SingleLevelMIP_Cormican_PathBased", str(e)
        
def getNodesThatCanInterdictEdgeSet(edgeSet, includeRepeats = False):
    nodesList = []
    for edge in edgeSet:
        nodesList.extend(jammersThatCanJamEdge[edge])
    if includeRepeats:
        return nodesList
    else:
        return set(nodesList)
    
def getJammedEdgesInEdgeSet(jamLocs, edgeSet):
    edgeList = []
    for edge in edgeSet:
        if canEdgeBeJammedByJammers(instance.CGraph, edge, jamLocs, interfModelType):
            edgeList.append(edge)
    return edgeList

def solveThroughputProblem_Pricing(G, JammingGraph):
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    #jammingDuals = dict([(node, dict([(edge, 0.0) for edge in G.edges()])) for node in JammingGraph.nodes()])
    jammingDuals = dict([(edge, 0.0) for edge in G.edges()])
    count = 0
    maxWeight, selected = maxWtIndepSet(instance.CGraph, capacityDuals, interfModelType)
    ISets.append(selected)
    createThroughputModel(G, instance.commodities, [selected], JammingGraph, interfModelType)
    while True:
        print "ITERATION", count
        nodeWeights = getWeightsForMaxIndSet(G, JammingGraph, capacityDuals, jammingDuals, interfModelType)
        #print "node weights", [(e, nodeWeights[e]) for e in nodeWeights.keys() if nodeWeights[e] > 0]
        maxWeight, selected = maxWtIndepSet(instance.CGraph, nodeWeights, interfModelType)
        ISets.append(selected)
        #print "selected ", selected
        throughput, capacityDuals, jammingDuals, usageDual, throughputWithoutPenalty = addISetAsCol_AndSolveTHProb(G, JammingGraph, ISets, selected, count, interfModelType)
        if((maxWeight - usageDual) <= 0.0001):
            break
        count += 1
    # compute and return these: coefDualTerms, constantDualTerms
    
def solveThroughputProblem_WiredNetwork(G, JammingGraph):
    ISets = []
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    #jammingDuals = dict([(node, dict([(edge, 0.0) for edge in G.edges()])) for node in JammingGraph.nodes()])
    jammingDuals = dict([(edge, 0.0) for edge in G.edges()])
    count = 0
    maxWeight, selected = maxWtIndepSet(instance.CGraph, capacityDuals, interfModelType)
    ISets.append(selected)
    createThroughputModel(G, instance.commodities, [selected], JammingGraph, interfModelType)
    while True:
        print "ITERATION", count
        nodeWeights = getWeightsForMaxIndSet(G, JammingGraph, capacityDuals, jammingDuals, interfModelType)
        #print "node weights", [(e, nodeWeights[e]) for e in nodeWeights.keys() if nodeWeights[e] > 0]
        maxWeight, selected = maxWtIndepSet(instance.CGraph, nodeWeights, interfModelType)
        ISets.append(selected)
        #print "selected ", selected
        throughput, capacityDuals, jammingDuals, usageDual, throughputWithoutPenalty = addISetAsCol_AndSolveTHProb(G, JammingGraph, ISets, selected, count, interfModelType)
        if((maxWeight - usageDual) <= 0.0001):
            break
        count += 1

def printFlowSoln():
    flowVarSoln = dict([(commod, dict([(path, flowVars[commod][path].X) for path in flowVars[commod].keys()])) for commod in instance.commodities.keys()])
    print "FLOW"
    for commod in instance.commodities.keys():
        for path in flowVars[commod].keys():
            if flowVarSoln[commod][path] > 0.0:
                print commod, path, instance.commodities[commod]['odPair'],  flowVarSoln[commod][path]
    
def newPathPricesOut(model, modelType, jamLocs, demForCommodDualsValuesDict, commod, edgesOnNewPath, pathLength, edgesJammed, jamUBDualValues):
    #print "newPathPricesOut"
    if jamVarsType == 'jam-vars' and modelType == 'cormican':
        amountToSubtract = len(getJammedEdgesInEdgeSet(jamLocs, edgesOnNewPath))
        jammedLocsThatCanJamPath = list(set(jamLocs) & set(getNodesThatCanInterdictEdgeSet(edgesOnNewPath, includeRepeats=False)))
        amountToAdd = len(jammedLocsThatCanJamPath) - amountToSubtract
    elif jamVarsType == 'jam-and-interdict-vars' and modelType == 'cormican':
        amountToAdd = 0
    elif jamVarsType == 'jam-vars' and modelType == 'regular':
        amountToAdd = 0
    elif jamVarsType == 'jam-and-interdict-vars' and modelType == 'regular':
        raise Exception("this combination is invalid")
    print "pathLength=", pathLength, "demForCommodDualsValuesDict[" +str(commod)+ "]=", demForCommodDualsValuesDict[commod], "amountToAdd=", amountToAdd, "total=",(pathLength + demForCommodDualsValuesDict[commod] + amountToAdd), (1.0 - FUZZ)
    if (pathLength + demForCommodDualsValuesDict[commod] + amountToAdd) < (1.0 - FUZZ):
        return True
    else:
        return False
    
# def newPathPricesOut_OLD(model, modelType, jamLocs, demForCommodDualsValues, commod, edgesOnNewPath, pathLength, edgesJammed, jamUBDualValues):
#     #print "newPathPricesOut"
#     if jamVarsType == 'jam-vars' and modelType == 'cormican':
#         jammedLocsThatCanJamPath = list(set(jamLocs) & set(getNodesThatCanInterdictEdgeSet(edgesOnNewPath, includeRepeats=False)))
#         amountToSubtract = len(jammedLocsThatCanJamPath)
#     elif jamVarsType == 'jam-and-interdict-vars' and modelType == 'cormican':
#         amountToSubtract = len(set(edgesJammed) & set(edgesOnNewPath))
#     elif jamVarsType == 'jam-vars' and modelType == 'regular':
#         jamUBDualsValuesDict = getJamUBDUalValuesDict(jamUBDualValues)
#         #print "jamUBDualsValuesDict", jamUBDualsValuesDict
#         amountToSubtract = 0
#         for edge in edgesOnNewPath:
#             for node in jammersThatCanJamEdge[edge]:
#                 amountToSubtract += jamUBDualsValuesDict[edge][node]
#     elif jamVarsType == 'jam-and-interdict-vars' and modelType == 'regular':
#         raise Exception("this combination is invalid")
#     print "pathLength", pathLength, demForCommodDualsValues[commod], (1.0 - amountToSubtract - FUZZ)
#     if((pathLength + demForCommodDualsValues[commod]) < (1.0 - amountToSubtract - FUZZ)):
#         return True
#     else:
#         return False
  
def solveThroughputProblem_Pricing_CormicanPathBased(G, edgesJammed, jamLocs, interfModelType, JammingGraph = None):
    #print "solveThroughputProblem_Pricing_CormicanPathBased", "interfModelType", interfModelType
    capacityDuals = dict([(edge, 0.0) for edge in edgeTriples])
    jammingDuals = dict([(edge, 0.0) for edge in edgeTriples])
    pathNumForCommod = dict([(commod, 0) for commod in instance.commodities.keys()])
    demForCommodDuals = dict([(commod, 0.0) for commod in instance.commodities.keys()])
    count = 0
    createThroughputModel_ContinuousJamming_PathBased_Cormican(G, edgesJammed, jamLocs, instance.commodities, interfModelType)
    newPathAvailable = False
    newISetAvailable = False
    usageDual = 0
    throughput = 0
    while True:
        print "ITERATION", count
        NewPaths = []
        for commod in instance.commodities.keys():
            pathLength, edgesOnNewPath = modifyAndSolve_ShortestPathProb(G, capacityDuals, instance.maxNumHops, commod, count)
            #if modelType == 'cormican-path-based-jam-vars':
            #    jammedLocsThatCanJamPath = list(set(jamLocs) & set(getNodesThatCanInterdictEdgeSet(edgesOnNewPath)))
            #    amountToSubtract = len(jammedLocsThatCanJamPath)
            #elif modelType == 'cormican-path-based-jam-and-interdict-vars':
            #    amountToSubtract = len(set(edgesJammed) & set(edgesOnNewPath))
            #if((pathLength + demForCommodDuals[commod]) < (1.0 - amountToSubtract - FUZZ)):
            if new_PathAndIS_CutCallback(modelType, jamLocs, demForCommodDuals, commod, edgesOnNewPath, pathLength, edgesJammed):
                #print "path accepted: length < 1.0", pathLength, edgesOnNewPath
                NewPaths.append((commod, pathNumForCommod[commod], edgesOnNewPath))
                pathNumForCommod[commod] = pathNumForCommod[commod] + 1
        newPathAvailable = (len(NewPaths) > 0)
        if(newPathAvailable):# new paths added]
            throughput, capacityDuals, jammingDuals, usageDual, demForCommodDuals = addPathsAsColumns_AndSolveThroughputProb(G, edgesJammed, jamLocs, NewPaths, interfModelType, count)
            #print "throughput", throughput
            #printFlowSoln()
        elif(interfModelType != 'none'):#no new paths added
            nodeWeights = getWeightsForMaxIndSet(G, capacityDuals, jammingDuals, interfModelType)
            #print "node weights > 0", [(e, nodeWeights[e]) for e in nodeWeights.keys() if nodeWeights[e] > 0], "capDuals", capacityDuals
            maxWeights, solnsSet = modifyAndSolve_maxWtIndepSet(instance.interferenceGraph, nodeWeights, interfModelType, count)
            maxWeight = maxWeights[0]
            selected = solnsSet[0]
            print "maxWt", maxWeight, "usageDual", usageDual, selected                      
            newISetAvailable = ((maxWeight - usageDual) > FUZZ)
            if(newISetAvailable):
                ISets.append(selected)
                throughput, capacityDuals, jammingDuals, usageDual, demForCommodDuals, throughputWithoutPenalty = addISetAsCol_AndSolveTHProb(G, ISets, selected, count, interfModelType)
                #print "throughput", throughput
                #printFlowSoln()
        if((not newPathAvailable) and (not newISetAvailable)):
            print "throughput", throughput
            printFlowSoln()
            break
        count += 1
    flowVarSoln = dict([(commod, dict([(path, flowVars[commod][path].X) for path in flowVars[commod].keys()])) for commod in instance.commodities.keys()])
    isetUsageSoln = [iSetUsage[k].X for k in range(len(ISets))]
    return flowVarSoln, isetUsageSoln, throughput

def updateSlackVarUBs(G, commodities, updates):
    global slackVarUBs
    sumOfUBs = 0
    for sign in ['-', '+']:
        constrName = 'flowBal'
        for commod in commodities.keys():
            for node in G.nodes():
                slackVarUBs[sign][constrName][commod][node] *= 1/2.0**updates
                sumOfUBs += slackVarUBs[sign][constrName][commod][node]
        constrName = 'cap'
        for edge in edgeTriples:
            slackVarUBs[sign][constrName][edge] *= 1/2.0**updates
            sumOfUBs += slackVarUBs[sign][constrName][edge]
        constrName = 'commodDemand'
        for commod in commodities.keys():
            slackVarUBs[sign][constrName][commod] *= 1/2.0**updates
            sumOfUBs += slackVarUBs[sign][constrName][commod]
        constrName = 'usage'
        slackVarUBs[sign][constrName] *= 1/2.0**updates
        sumOfUBs += slackVarUBs[sign][constrName]
    print "sum of UBs", sumOfUBs
        

def updateSlackVarObjCoefsAroundDualSolutionAndKeepRange(G, commodities, flowBalDuals, capacityDuals, jammingDuals, usageDual, demForCommodDuals):
    global slackVarMultipliers
    constrName = 'flowBal'
    for commod in commodities.keys():
        for node in G.nodes():
            currentRange = slackVarMultipliers['+'][constrName][commod][node] - slackVarMultipliers['-'][constrName][commod][node]
            slackVarMultipliers['-'][constrName][commod][node] = flowBalDuals[commod][node] - currentRange/2.0
            slackVarMultipliers['+'][constrName][commod][node] = flowBalDuals[commod][node] + currentRange/2.0
    constrName = 'cap'
    for edge in edgeTriples:
        currentRange = slackVarMultipliers['+'][constrName][edge] - slackVarMultipliers['-'][constrName][edge]
        slackVarMultipliers['-'][constrName][edge] = capacityDuals[edge] - currentRange/2.0
        slackVarMultipliers['+'][constrName][edge] = capacityDuals[edge] + currentRange/2.0
    constrName = 'commodDemand'
    for commod in commodities.keys():
        currentRange = slackVarMultipliers['+'][constrName][commod] - slackVarMultipliers['-'][constrName][commod]
        slackVarMultipliers['-'][constrName][commod] = demForCommodDuals[commod] - currentRange/2.0
        slackVarMultipliers['+'][constrName][commod] = demForCommodDuals[commod] + currentRange/2.0
    constrName = 'usage'
    currentRange = slackVarMultipliers['+'][constrName] - slackVarMultipliers['-'][constrName]
    slackVarMultipliers['-'][constrName] = usageDual - currentRange/2.0
    slackVarMultipliers['+'][constrName] = usageDual + currentRange/2.0
    
def updateSlackVarObjCoefsAroundPrevDualSolutionAndIncreaseRange(G, commodities, prevFlowBalDuals, prevCapacityDuals, prevJammingDuals, 
                                                                 prevUsageDual, prevDemForCommodDuals, slackVarValues, rangeIncreaseMultiplier = 2.0):
    global slackVarMultiplier
    multipliers = {'-' : -1, '+' : 1}
    constrName = 'flowBal'
    for commod in commodities.keys():
        for node in G.nodes():
            theRange = slackVarMultipliers['+'][constrName][commod][node] - slackVarMultipliers['-'][constrName][commod][node]
            for sign in ['-', '+']:
                constraintSlack = slackVarUBs[sign][constrName][commod][node] - slackVarValues[sign][constrName][commod][node]
                if(constraintSlack < FUZZ):
                    theRange = theRange * rangeIncreaseMultiplier
                slackVarMultipliers[sign][constrName][commod][node] = prevFlowBalDuals[commod][node] + multipliers[sign] * (theRange/2.0)
    constrName = 'cap'
    for edge in edgeTriples:
        theRange = slackVarMultipliers['+'][constrName][edge] - slackVarMultipliers['-'][constrName][edge]
        for sign in ['-', '+']:
            constraintSlack = slackVarUBs[sign][constrName][edge] - slackVarValues[sign][constrName][edge]
            if(constraintSlack < FUZZ):
                theRange = theRange * rangeIncreaseMultiplier
            slackVarMultipliers[sign][constrName][edge] = prevCapacityDuals[edge] + multipliers[sign] * (theRange/2.0)
    constrName = 'commodDemand'
    for commod in commodities.keys():
        theRange = slackVarMultipliers['+'][constrName][commod] - slackVarMultipliers['-'][constrName][commod]
        for sign in ['-', '+']:
            constraintSlack = slackVarUBs[sign][constrName][commod] - slackVarValues[sign][constrName][commod]
            if(constraintSlack < FUZZ):
                theRange = theRange * rangeIncreaseMultiplier
            slackVarMultipliers[sign][constrName][commod] = prevDemForCommodDuals[commod] + multipliers[sign] * (theRange/2.0)
    constrName = 'usage'
    theRange = slackVarMultipliers['+'][constrName] - slackVarMultipliers['-'][constrName]
    for sign in ['-', '+']:
        constraintSlack = slackVarUBs[sign][constrName] - slackVarValues[sign][constrName]
        if(constraintSlack < FUZZ):
            theRange = theRange * rangeIncreaseMultiplier
        slackVarMultipliers[sign][constrName] = prevUsageDual + multipliers[sign] * (theRange/2.0)

def slackVarsAreZero(slackVarValues, G, commodities, stabilize):
    if stabilize is False or slackVarValues is None:
        return True
    else:
        sumOfValues = 0
        for sign in ['-', '+']:
            constrName = 'flowBal'
            for commod in commodities.keys():
                for node in G.nodes():
                    sumOfValues += slackVarValues[sign][constrName][commod][node]
            constrName = 'cap'
            for edge in edgeTriples:
                sumOfValues += slackVarValues[sign][constrName][edge]
            constrName = 'commodDemand'
            for commod in commodities.keys():
                sumOfValues += slackVarValues[sign][constrName][commod]
            constrName = 'usage'
            sumOfValues += slackVarValues[sign][constrName]
        #print "sumOfValues", sumOfValues
        return (sumOfValues < FUZZ)

def setInitialSlackVarUBValues(G, commodities, initialValues):
    global slackVarUBs
    slackVarUBs = createSlackVariableUBs(G, commodities, initialValues)

def setInitialSlackVarObjCoefValues(G, commodities, initialValues):
    global slackVarMultipliers
    slackVarMultipliers = createSlackVariableMultipliers(G, commodities, initialValues)

# obtain lower bound for problem during column generation algorithm
def getDualBound(flowBalDuals, capacityDuals, usageDual, demForCommodDuals, throughput, maxWeight):
    dualBound = 0
    dualBound += maxWeight # based on dual bound from Tebboth paper (don't need to include master problem objective)
    return dualBound
    
def dualSolutionIsBestSoFar(dualBound, bestDualBound):
    return dualBound < bestDualBound

def solveThroughputProblem_Pricing_CormicanArcBased(G, edgesJammed, interfModelType, capacityDualsInit, JammingGraph = None, suffix = '', stabilize = False):
    global g_capacityDuals, ISets
    flowBalDuals = {}
    slackVarValues = None
    capacityDuals = capacityDualsInit
    jammingDuals = dict([(edge, 0.0) for edge in edgeTriples])
    pathNumForCommod = dict([(commod, 0) for commod in instance.commodities.keys()])
    demForCommodDuals = dict([(commod, 0.0) for commod in instance.commodities.keys()])
    count = 0
    newISetAvailable = False
    usageDual = 0
    throughput = 0
    throughputWithoutPenalty = 0
    updates = 0
    bestDualBound = float('inf') #we are maximizing
    if stabilize:
        setInitialSlackVarUBValues(G, instance.commodities, 0.01)
        setInitialSlackVarObjCoefValues(G, instance.commodities, 0.1)
        updateSlackStuffInModel(G, instance.commodities)
    while True:
        print "ITERATION", count
        prevFlowBalDuals = flowBalDuals.copy()
        prevCapacityDuals = capacityDuals.copy()
        prevJammingDuals = jammingDuals.copy()
        prevUsageDual = usageDual
        prevDemForCommodDuals = demForCommodDuals.copy()
        throughput, capacityDuals, jammingDuals, usageDual, demForCommodDuals, throughputWithoutPenalty, slackVarValues, slackVarUBDuals, flowBalDuals = \
            solveThroughputProblem(G, count, interfModelType, suffix, stabilize)
        nodeWeights = getWeightsForMaxIndSet(G, capacityDuals, jammingDuals, interfModelType)
        maxWeights, solnsSet = modifyAndSolve_maxWtIndepSet(instance.interferenceGraph, nodeWeights, interfModelType, count)
        maxWeight = maxWeights[0]
        selected = solnsSet[0]
        print "throughput", throughput
        print "maxWt", maxWeight, "usageDual", usageDual #, selected
        newISetAvailable = ((maxWeight - usageDual) > FUZZ)
        if((newISetAvailable is False) and slackVarsAreZero(slackVarValues, G, instance.commodities, stabilize)):
            break
        else:
            if newISetAvailable:
                ISets.append(selected)
                addISetAsCol(G, ISets, selected, count, interfModelType, suffix, stabilize)
            if stabilize:
                if not newISetAvailable:
                    updateSlackVarUBs(G, instance.commodities, updates)
                    updates += 1
                dualBound = getDualBound(flowBalDuals, capacityDuals, usageDual, demForCommodDuals, throughput, maxWeight)
                isDualSolnBestSoFar = dualSolutionIsBestSoFar(dualBound, bestDualBound) # is dual solution the best known estimate of the optimal dual solution?
                bestDualBound =  min(dualBound, bestDualBound)
                print "best dual bound: ", bestDualBound, isDualSolnBestSoFar
                if isDualSolnBestSoFar:
                    updateSlackVarObjCoefsAroundDualSolutionAndKeepRange(G, instance.commodities, flowBalDuals, capacityDuals, jammingDuals, usageDual, demForCommodDuals)
                else:
                    updateSlackVarObjCoefsAroundPrevDualSolutionAndIncreaseRange(G, instance.commodities, prevFlowBalDuals, prevCapacityDuals, prevJammingDuals, prevUsageDual, \
                                                                                 prevDemForCommodDuals, slackVarValues, 2.0)
                updateSlackStuffInModel(G, instance.commodities)
        count += 1
    g_capacityDuals = capacityDuals
    flowVarSoln = dict([(commod, dict([(edge, flowVars[commod][edge].X) for edge in flowVars[commod].keys()])) for commod in instance.commodities.keys()])
    isetUsageSoln = [iSetUsage[k].X for k in range(len(ISets))]
    return flowVarSoln, isetUsageSoln, throughput, throughputWithoutPenalty

def OLD_getCorePoint_Rounded():
    masterProbLP.optimize()
    jamPlacedLPDict = {n : jamPlacedLP[n].X for n in instance.jamGraph.nodes()}
    jamPlacedLPDict_Sorted = OrderedDict(sorted(jamPlacedLPDict.items(), key=lambda t: t[1]))
    jamLocsCorePoint = []
    counter = 0
    for n in jamPlacedLPDict_Sorted.keys():
        if counter <= instance.jamBudget:
            jamLocsCorePoint.append(jamPlacedLPDict[n])
        counter += 1
    return jamLocsCorePoint

def getCorePoint_Fractional():
    masterProbLP.optimize()
    jamPlacedLPDict = {n : jamPlacedLP[n].X for n in instance.jamGraph.nodes()}
    return jamPlacedLPDict
        
def solve_FullMIP_DelayedRowGen_Cormican_PathBased(G, JammingGraph, InitialPaths, InitialISets, maxNumHops):
    #Paths = InitialPaths
    ISets = InitialISets
    interdictedEdges = dict([(edge, 0.0) for edge in G.edges()])
    capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
    pathNumForCommod = dict([(commod, 0) for commod in instance.commodities.keys()])
    count = 0
    create_SingleLevelMIP_Cormican_PathBased(G, instance.commodities, InitialPaths, ISets, interfModelType)
    newPathAvailable= False
    newISetAvailable= False
    usageDual = 0
    while True:
        print "ITERATION", count
        #numPaths = len(Paths)
        NewPaths=[]
        #print "capDualsGr0", [(edge, capacityDuals[edge]) for edge in capacityDuals.keys() if capacityDuals[edge] > 0.001]
        for commod in instance.commodities.keys():
            #print "interdictedEdges", interdictedEdges
            weights = dict([(e, capacityDuals[e] + interdictedEdges[e]) for e in G.edges()])
            pathLength, newPath = findLeastCostPath(G, weights, maxNumHops, [], commod, count)
            if(pathLength < (1.0 - FUZZ)):
                NewPaths.append((commod, pathNumForCommod[commod], newPath))
                pathNumForCommod[commod] = pathNumForCommod[commod] + 1
        #print "NewPaths", count, NewPaths
        #print "ISets ", ISets
        newPathAvailable = (len(NewPaths) > 0)
        if(newPathAvailable):# new paths added]
            throughput, capacityDuals, usageDual, interdictedEdges = addPathsAsRows_AndSolveFullMIP(G, JammingGraph, NewPaths, interfModelType, count)
        else: #no new paths added
            nodeWeights = dict([(e, G.edge[e[0]][e[1]]['capacity'] * capacityDuals[e]) for e in G.edges()])
            maxWeight, selected = maxWtIndepSet(instance.CGraph, nodeWeights, interfModelType)
            print "selected ", selected
            ISets.append(selected)
            #print "ISets ", ISets
            newISetAvailable = (maxWeight > usageDual)
            print "maxWeight, usageDual", maxWeight, usageDual, newISetAvailable
            if(newISetAvailable):
                throughput, capacityDuals, usageDual, interdictedEdges = addISetAsRow_AndSolveFullMIP(G, selected, count, interfModelType)
                #iSetUsageVals = [(ISets[k], iSetUsage[k].X) for k in range(len(ISets)) if iSetUsage[k].X > 0.001]
        print "jamPlaced", [n for n in g_JamGraph.nodes() if jamPlaced[n].X > FUZZ]
        #print "throughput", throughput
        #print "iSetUsageVals", iSetUsageVals
        #print "test", newPathAvailable, newISetAvailable, (not newPathAvailable) and (not newISetAvailable)
        if((not newPathAvailable) and (not newISetAvailable)):
            break
        count += 1
        
def getSlackVarUBDuals(G):
    slackVarUBDuals = {}
    for sign in ['-', '+']:
        slackVarUBDuals[sign] = {}
        constrName = 'flowBal'
        slackVarUBDuals[sign][constrName] = {}
        for commod in instance.commodities.keys():
            slackVarUBDuals[sign][constrName][commod] = {}
            for node in G.nodes():
                slackVarUBDuals[sign][constrName][commod][node] = slackUBConstr[sign][constrName][commod][node].Pi
        constrName = 'cap'
        slackVarUBDuals[sign][constrName] = {}
        for edge in edgeTriples:
            slackVarUBDuals[sign][constrName][edge] = slackUBConstr[sign][constrName][edge].Pi
        constrName = 'commodDemand'
        slackVarUBDuals[sign][constrName] = {}
        for commod in instance.commodities.keys():
            slackVarUBDuals[sign][constrName][commod] = slackUBConstr[sign][constrName][commod].Pi
        constrName = 'usage'
        slackVarUBDuals[sign][constrName] = slackUBConstr[sign][constrName].Pi
    return slackVarUBDuals
    
def getSlackVarValues(G, commodities):
    slackVarsValues = {}
    for sign in ['-', '+']:
        slackVarsValues[sign] = {}
        constrName = 'flowBal'
        slackVarsValues[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarsValues[sign][constrName][commod] = {}
            for node in G.nodes():
                slackVarsValues[sign][constrName][commod][node] = slackVars[sign][constrName][commod][node].X
        constrName = 'cap'
        slackVarsValues[sign][constrName] = {}
        for edge in edgeTriples:
            slackVarsValues[sign][constrName][edge] = slackVars[sign][constrName][edge].X
        constrName = 'commodDemand'
        slackVarsValues[sign][constrName] = {}
        for commod in commodities.keys():
            slackVarsValues[sign][constrName][commod] = slackVars[sign][constrName][commod].X
        constrName = 'usage'
        slackVarsValues[sign][constrName] = slackVars[sign][constrName].X
    return slackVarsValues

def solveThroughputProblem(G, iteration, interfModelType, suffix = '', stabilize = False):
    try:
        if(writeToFile):
            gurobiThroughputModel.write("/tmp/ThroughputProblem_"+ str(suffix) + "_" + str(iteration)+".lp")
        gurobiThroughputModel.optimize()
        throughputWithPenalty = gurobiThroughputModel.objVal
        throughputWithoutPenalty = sum([returnFlow[commod].X for commod in instance.commodities.keys()])
        
        flowBalDuals = {}
        for commod in instance.commodities.keys():
            flowBalDuals[commod] = {}
            for node in G.nodes():
                flowBalDuals[commod][node] = flowBalanceConstrs[commod][node].Pi
                                                                                   
        capacityDuals = dict([(edge, capConstraints[edge].Pi) for edge in edgeTriples])

        if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol')):
            usageDualVal = usageConstr.Pi

        elif(interfModelType == 'none'):
            usageDualVal = 0

        demForCommodDuals = dict([(commod, commodDemConstr[commod].Pi) for commod in instance.commodities.keys()])
        slackVarValues = None
        slackVarUBDuals = None
        if stabilize:
            slackVarValues = getSlackVarValues(G, instance.commodities)
            slackVarUBDuals = getSlackVarUBDuals(G)
    except gurobipy.GurobiError as e:
        print "solveThroughputProblem error: ", str(e)
    return throughputWithPenalty, capacityDuals, {}, usageDualVal, demForCommodDuals, throughputWithoutPenalty, slackVarValues, slackVarUBDuals, flowBalDuals

def solveFullMIP(G, iteration):
    #print "solveThroughputProblemForISets", len(ISets)
    try:
        if(writeToFile):
            gurobiModel.write("/home/hmedal/Documents/Temp/JainFullMIP_"+str(iteration)+".lp")
        gurobiModel.optimize()
        throughput = gurobiModel.objVal
        print "throughput:", throughput
        capacityDualVals = dict([(edge, capDuals[edge].X) for edge in G.edges()])
        print "capacityDuals:", [(key, capacityDualVals[key]) for key in capacityDualVals.keys() if capacityDualVals[key] > 0]
        print "usageDual:", usageDual.X
        interdicted = dict([(e, edgeInterdict[e].X) for e in G.edges()])
    except gurobipy.GurobiError as e:
        print str(e)
    return throughput, capacityDualVals, usageDual.X, interdicted

def getInterdictedEdges(G, JammingGraph, interfModelType):
    jamLocs = [loc for loc in JammingGraph.nodes()]
    #print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    #print [JammingGraph.node[loc]['selected'] for loc in jamLocs]
    interdictedEdges = []
    for e in G.edges():
        #print e, isEdgeJammedByJammers(G, JammingGraph, e, jamLocs, interfModelType)
        if(isEdgeJammedByJammers(G, JammingGraph, e, jamLocs, interfModelType)):
        #if(canEdgeBeJammedByJammers(G, e, jamLocs, interfModelType) & (sum([JammingGraph.node[loc]['selected'] for loc in jamLocs]) > 0)):
            interdictedEdges.append(e)
        #print interdictedEdges
    return interdictedEdges

def addPathsAsColumns_AndSolveThroughputProb(G, edgesJammed, jamLocs, NewPaths, interfModelType, iteration, JammingGraph = None):
    global Paths, flowVars
    #pathCounter = len(Paths)
    #interdictedEdges = getInterdictedEdges(G, JammingGraph, interfModelType)
    #print "addPathsAsColumns_AndSolveThroughputProb", "interdictedEdges", edgesJammed
    try:
        for item in NewPaths:
            #print "item", item
            commod = item[0]
            pathNum = item[1]
            #print "commod", commod, "pathNum", pathNum
            edgesOnPath = item[2]
            Paths[commod].append((pathNum, edgesOnPath))
            if modelType == 'cormican-path-based-jam-vars':
                jammedLocsThatCanJamPath = list(set(jamLocs) & set(getNodesThatCanInterdictEdgeSet(edgesOnPath)))
                amountToSubtract = len(jammedLocsThatCanJamPath)
            elif modelType == 'cormican-path-based-jam-and-interdict-vars':
                amountToSubtract = len(set(edgesJammed) & set(edgesOnPath))
            objCoef = 1 - amountToSubtract
            #print "objCoef", objCoef
            #print "capConstraints", capConstraints
            constrs = []
            coefs = []
            for edge in edgeTriples:
                constrs.append(capConstraints[edge])
                coefs.append(int(edge in edgesOnPath))
            for commodIndex in instance.commodities.keys():
                constrs.append(commodDemConstr[commodIndex])
                coefs.append(1.0)
            #coefs = [int(edge in edgesOnPath) for edge in edgeTriples]
            #print "add path", getPathFromEdgeList(edgesOnPath, commodities[commod]['odPair'][0], commodities[commod][1]), objCoef, coefs
            flowVars[commod][pathNum] = gurobiThroughputModel.addVar(0.0, obj = objCoef, name="flowAdd_"+str(commod)+","+str(pathNum), column = gurobipy.Column(coefs, constrs))
            gurobiThroughputModel.update()
            #print "var", commod, pathNum, flowVars[commod][pathNum]
            if(writeToFile):
                gurobiThroughputModel.write("/tmp/JainPathsUpdated_"+str(iteration)+".lp")
    except gurobipy.GurobiError as e:
        print "addISetAndSolveThroughputProb ERROR:", str(e), "numEdges", len(G.edges())
    return solveThroughputProblem(G, iteration, interfModelType)

def addPathsAsRows_AndSolveFullMIP(G, JammingGraph, NewPaths, interfModelType, iteration):
    #print "addPathsAsColumns_AndSolveThroughputProb", "interdictedEdges", interdictedEdges
    try:
        for item in NewPaths:
            #print "item", item
            commod = item[0]
            pathNum = item[1]
            edgesOnPath = item[2]
            gurobiModel.addConstr(sum([capDuals[e] for e in edgesOnPath]) + sum([edgeInterdict[e] for e in edgesOnPath]) >= 1, "new path commod" + str(commod)+"_path" + str(pathNum))
            gurobiModel.update()
            if(writeToFile):
                gurobiModel.write("/home/hmedal/Documents/Temp/JainPathsUpdated_RowGen_"+str(iteration)+".lp")
    except gurobipy.GurobiError as e:
        print "addISetAndSolveThroughputProb ERROR:", str(e), "numEdges", len(G.edges())
    return solveFullMIP(G, iteration)

def getAttrIndex(edgeTriple):
    if(numRadiosPerNode):
        attrIndex = 0
    else:
        attrIndex = edgeTriple[2]
    return attrIndex

def addISetAsCol(G, ISets, ISet, iteration, interfModelType, suffix = '', stabilize = False):
    print "interfModelType", interfModelType
    try:
        numISets = len(ISets)
        constrs = [capConstraints[edge] for edge in edgeTriples]
        coefs = []
        for edgeTriple in edgeTriples:
            attrIndex = getAttrIndex(edgeTriple)
            coefs.append(-int(edgeTriple in ISet) * G.edge[edgeTriple[0]][edgeTriple[1]][attrIndex]['capacity'])
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            None
        constrs.append(usageConstr)  
        coefs.append(1)
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            iSetUsage.append(gurobiThroughputModel.addVar(0.0, 1.0, vtype=gurobipy.GRB.CONTINUOUS, name = "lambdaAdd_"+str(numISets-1), column = gurobipy.Column(coefs, constrs)))
        gurobiThroughputModel.update()
        if(writeToFile):
            gurobiThroughputModel.write("/tmp/JainISetUpdated_"+str(iteration)+".lp")
    except gurobipy.GurobiError as e:
        print "addISet ERROR:", str(e), "numEdges", len(G.edges())

def addISetAsCol_AndSolveTHProb(G, ISets, ISet, iteration, interfModelType, suffix = '', stabilize = False):
    try:
        addISetAsCol(G, ISets, ISet, iteration, interfModelType, suffix, stabilize)
    except gurobipy.GurobiError as e:
        print "addISetAndSolveThroughputProb ERROR:", str(e), "numEdges", len(G.edges())
    return solveThroughputProblem(G, iteration, interfModelType, suffix, stabilize)

def addISetAsRow_AndSolveFullMIP(G, ISet, iteration, interfModelType):
    #print "addISetAndSolveThroughputProb", ISet, [-int(edge in ISet) for edge in G.edges()]
    try:
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            None
        gurobiModel.addConstr(sum([-G.edge[e[0]][e[1]]['capacity'] * capDuals[e] for e in ISet]) + usageDual >= 0, "new ISet row "+ str(iteration))
        gurobiModel.update()
    except gurobipy.GurobiError as e:
        print "addISetAndSolveThroughputProb ERROR:", str(e), "numEdges", len(G.edges())
    return solveFullMIP(G, iteration)
    
def getWeightsForMaxIndSet(G, capacityDuals, jammingDuals, interfModelType):
    #print "getWeightsForMaxIndSet"
    return dict([(edge, capacityDuals[edge]) for edge in edgeTriples])
    
def getNetworkMaxTotalThroughput(G, commodities, interfModelType):
    createThroughputModel_ContinuousJamming_Cormican(G, commodities, ISets, g_JamGraph, interfModelType)
    throughput, capacityDuals, jammingDuals, usageConstr = solveThroughputProblem(G, 0, interfModelType)
    flowOnArcs = dict([(commodity, dict([(edge, flowVars[commodity][edge].X) for edge in G.edges() if flowVars[commodity][edge].X > FUZZ])) for commodity in commodities.keys()])
    return throughput, flowOnArcs
    
def runProtocolWithJamming_CormicanPathBasedThroughput(gridSize, interdict_budget_arg, numCommodities):
    global interfModelType
    interdictBudget = interdict_budget_arg
    interfModelType = '802.11-MAC-protocol'
    
    capacityDuals = dict([(edge, 1.0) for edge in edgeTriples])
    nodeWeights = capacityDuals
    maxWeight, selected = maxWtIndepSet(instance.CGraph, nodeWeights, interfModelType)
    edgesJammed = []
    flowVarSoln, isetUsageSoln, throughput = solveThroughputProblem_Pricing_CormicanPathBased(instance.CGraph, edgesJammed, 
                                                                                              interfModelType, instance.jamGraph)
    print "THROUGHPUT", throughput
    print "FINISHED"
    
def runProtocolWithJamming_FullMIP_RowGenCallback():
    print "runProtocolWithJamming_FullMIP_RowGenCallback"
    #print "jamLocs", instance.jamGraph.nodes()
    global interfModelType
    interfModelType = '802.11-MAC-protocol'
    startTime = time.time()
    if flowVarsType == 'path-based':
        print "maxNumHops", instance.maxNumHops
    throughput = solveFullMIP_DelayedRowGen_withCallbacks(instance.CGraph, instance.jamGraph, Paths, ISets, 
                                                                instance.maxNumHops, interfModelType, 0.0001) # solves (9) via delayed row generation
    runTime = time.time() - startTime
    print "runTime", instance.gridSize, instance.maxNumHops, instance.jamBudget, runTime
    
    printOutVarValues()
    #printOutConstrValues()
    #print "ISets", ISets
    #for myISet in ISets:
    #    sumOfWts = sum([capDuals[e].X for e in myISet])
    #    print "set", sumOfWts, myISet
#     for commod in instance.commodities.keys():
#         orig, dest = instance.commodities[commod]['odPair'][0], instance.commodities[commod]['odPair'][1]
#         print "Commod", commod, orig, dest
#         for path in Paths[commod]:
#             sumOfWts = sum([capDuals[e].X for e in path])
#             varValueOfJamsThatCanInterdict = [jamPlaced[n].X for n in getNodesThatCanInterdictEdgeSet(path)]
#             totalInterdicts = sum(varValueOfJamsThatCanInterdict)
#             print "   pathInfo", sumOfWts + totalInterdicts, totalInterdicts, "      ", path
#             print "      ", [(n, jamPlaced[n].X) for n in getNodesThatCanInterdictEdgeSet(path)]
#             print "      ", [(e, capDuals[e].X) for e in path if capDuals[e].X > FUZZ]
    print "throughput", throughput
    numSolutionsFound = gurobiModel.SolCount
    numNodesExplored = gurobiModel.NodeCount
    print "numSolutionsFound", numSolutionsFound
    print "numNodesExplored", numNodesExplored
    print "numMaxWtIndSetProbsSolved", numMaxWtIndSetProbsSolved
    tableName = "CompResults"
    lb = gurobiModel.objVal
    ub = gurobiModel.getAttr('ObjBound')
    infrastructureInfo = [inst, clusterName]
    modelInfo = [flowVarsType]
    dataSetInfo = [datasetName, datasetIndex, numCommod, commodLocType, commodDemandType]
    instanceInfo = [instance.numJamLocs, instance.commRange, instance.infRange, instance.multRadiosPerNode, instance.numChannels, instance.jamBudget, instance.jamRange, instance.maxNumHops]
    algParams = [algType]
    algOutput = [runTime, numNodesExplored, numMaxWtIndSetProbsSolved, lb, ub]
    solnOutput = [throughput]
    #dbUtil.printResultsToDB(databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput)
    print databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput
    print "FINISHED"
    
def runProtocolWithJamming_Benders():
    global numNodesExplored
    print "runProtocolWithJamming_FullMIP_RowGenCallback"
    #print "jamLocs", instance.jamGraph.nodes()
    global interfModelType, ISets
    
    interfModelType = '802.11-MAC-protocol'
    startTime = time.time()
    if flowVarsType == 'path-based':
        raise Exception("Benders not implemented for path-based")
    if 'callback' in algType:
        lb, ub, jamLocs = runBenders_Callback(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    else:
        lb, ub, jamLocs = runBenders_Classic(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    runTime = time.time() - startTime
    if 'callback' in algType:
        numNodesExplored = gurobiModel.NodeCount
    print "runTime", runTime
    print "numNodesExplored", numNodesExplored
    print "numMaxWtIndSetProbsSolved", numMaxWtIndSetProbsSolved
    #printOutVarValues_Benders()
    #printOutConstrValues()
    print "throughput", ub
    tableName = "CompResults"
    infrastructureInfo = [inst, clusterName]
    modelInfo = [flowVarsType]
    dataSetInfo = [datasetName, datasetIndex, numCommod, commodLocType, commodDemandType]
    instanceInfo = [instance.numJamLocs, instance.commRange, instance.infRange, instance.multRadiosPerNode, instance.numChannels, instance.jamBudget, instance.jamRange, instance.maxNumHops]
    algParams = [algType]
    algOutput = [runTime, numNodesExplored, numMaxWtIndSetProbsSolved, lb, ub]
    solnOutput = [ub]
    dbUtil.printResultsToDB(databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput)
    print "FINISHED"
    
# not working
def runNoInterferenceHeuristic_Benders():
    global interfModelType, ISets
    interfModelType = 'none'
    startTime = time.time()
    print "\n", "RUNNING HEURISTIC"
    lb, heuristicThroughput, heuristicSolution = runBenders_Classic(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    runTimeOfHeuristic = time.time() - startTime

    interfModelType = '802.11-MAC-protocol'
    ISets = []
    print "\n", "COMPUTING OBJ VALUE FOR HEURISTIC SOLN"
    lb, actualThroughput, jamLocs = runBenders_Classic(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    
    print "\n", heuristicSolution, heuristicThroughput, actualThroughput, "\n"
    
    tableName = "NoInfHeuristic"
    infrastructureInfo = [inst, clusterName]
    modelInfo = [flowVarsType]
    dataSetInfo = [datasetName, datasetIndex, numCommod, commodLocType, commodDemandType]
    instanceInfo = [instance.numJamLocs, instance.commRange, instance.infRange, instance.multRadiosPerNode, instance.numChannels, instance.jamBudget, instance.jamRange, instance.maxNumHops]
    algParams = [algType]
    algOutput = [runTimeOfHeuristic, numMaxWtIndSetProbsSolved]
    solnOutput = [heuristicThroughput, actualThroughput]
    dbUtil.printResultsToDB(databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput)
    print "FINISHED"

def runNoInterferenceHeuristic_FullMIP_RowGenCallback():
    global interfModelType, ISets
    
    interfModelType = 'none'
    startTime = time.time()
    if flowVarsType == 'path-based':
        print "maxNumHops", instance.maxNumHops
    heuristicThroughput = solveFullMIP_DelayedRowGen_withCallbacks(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    runTime = time.time() - startTime
    heuristicSoln = [n for n in jamPlaced.keys() if jamPlaced[n].X > 0.0]
    printOutVarValues()
    print "runTime", instance.gridSize, instance.maxNumHops, instance.jamBudget, runTime, heuristicSoln, "\n"
    
    interfModelType = '802.11-MAC-protocol'
    ISets = []
    actualThroughput = computeThroughputForJammingSolution(instance.CGraph, instance.jamGraph, heuristicSoln, Paths, ISets, instance.maxNumHops, interfModelType, rel_gap_tol = 0.0001)
    
    interfModelType = '802.11-MAC-protocol'
    ISets = []
    bestThroughput = solveFullMIP_DelayedRowGen_withCallbacks(instance.CGraph, instance.jamGraph, Paths, ISets, instance.maxNumHops, interfModelType, 0.0001)
    bestSoln = [n for n in jamPlaced.keys() if jamPlaced[n].X > 0.0]
    
    print "\n", heuristicSoln, bestSoln, heuristicThroughput, actualThroughput, bestThroughput, "\n"
    
    tableName = "NoInterferenceHeuristic"
    infrastructureInfo = [inst, clusterName]
    modelInfo = [flowVarsType]
    dataSetInfo = [datasetName, datasetIndex, numCommod, commodLocType, commodDemandType]
    instanceInfo = [instance.numJamLocs, instance.commRange, instance.infRange, instance.multRadiosPerNode, instance.numChannels, instance.jamBudget, instance.jamRange, instance.maxNumHops]
    algParams = [algType]
    algOutput = [runTime]
    solnOutput = [str(heuristicSoln), str(bestSoln), heuristicThroughput, actualThroughput, bestThroughput]
    dbUtil.printResultsToDB(databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput)
    print "FINISHED"
    
def printOutVarValues():
    print "jamPlaced", [(n, jamPlaced[n].X) for n in instance.jamGraph.nodes() if jamPlaced[n].X > FUZZ]
    if jamVarsType == 'jam-and-interdict-vars':
        print "edgeInterdict", [(triple, edgeInterdict[triple].X) for triple in edgeTriplesInterdictable if edgeInterdict[triple].X > FUZZ]
    print "capDuals", [(e, capDuals[e].X) for e in edgeTriplesInterdictable if abs(capDuals[e].X) > FUZZ]
    printFlowBalanceDualsValues()
    #if flowVarsType == 'arc-based':
    #    for commod in instance.commodities.keys():
    #        for node in instance.CGraph.nodes():
     #           print "flowBalanceDuals", commod, node, flowBalanceDuals[commod][node].X
        
def printOutVarValues_Benders():
    print "jamPlaced", [(n, jamPlaced[n].X) for n in instance.jamGraph.nodes() if jamPlaced[n].X > FUZZ]
                
def printOutConstrValues():
    if flowVarsType == 'arc-based':
        for commod in instance.commodities.keys():
            for edgeInfo in instance.CGraph.edges(data = True):
                edgeTriple = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
                #print 'edgeTriple', edgeTriple, flowBalanceDuals[commod][edgeTriple[0]], flowBalanceDuals[commod][edgeTriple[1]]
                #print "sum", sum(jamPlaced[node] for node in jammersThatCanJamEdge[edgeTriple])
                #print "flowVarRedCost", commod, edgeTriple, flowBalanceDuals[commod][edgeTriple[0]].X, flowBalanceDuals[commod][edgeTriple[1]].X, capDuals[edgeTriple].X
    
def runProtocolWithJamming_FullMIP_RowGenCallback_PathBased_JamVarsOnly():
    print "runProtocolWithJamming_FullMIP_RowGenCallback_PathBased_Test"
    global interfModelType
    interdictBudget = instance.jamBudget
    interfModelType = '802.11-MAC-protocol'
    startTime = time.time()
    print "maxNumHops", instance.maxNumHops
    print "startTime", startTime
    throughput = solveFullMIP_DelayedRowGen_withCallbacks(instance.CGraph, instance.jamGraph, Paths, ISets, 
                                                                    instance.maxNumHops, interfModelType, 0.0001)
    runTime = time.time() - startTime
    print "runTime", instance.gridSize, instance.maxNumHops, instance.jamBudget, runTime
    print "jamPlaced", [(n, jamPlaced[n].X) for n in instance.jamGraph.nodes() if jamPlaced[n].X > FUZZ]
    print "capDuals", [(e, capDuals[e].X) for e in edgeTriples if abs(capDuals[e].X) > FUZZ]
    print "throughput", throughput
    print "FINISHED"
    
def getISetIndexForEdge(edge, ISets):
    counter = 0
    for ISet in ISets:
        if(edge in ISet):
            return counter
        else:
            counter += 1
    
def initializeSets():
    global Paths, ISets, edgeTriples, edgeTriplesInterdictable, g_capacityDuals
    print "initializeSets"
    
    Paths = {}
    for commod in instance.commodities.keys():
        Paths[commod] = []
    ISets = []
    edgeTriples = [(edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']) for edgeInfo in instance.CGraph.edges(data = True)]
    print "edgeTriples"
    print edgeTriples
    
    edgeTriplesInterdictable = [(edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel']) 
                                for edgeInfo in instance.CGraph.edges(data = True) if edgeInfo[2]['edgeType'] == 'real']
    
    print "edgeTriplesInt"
    print edgeTriplesInterdictable
    
    g_capacityDuals = dict([(edge, 0.1) for edge in edgeTriples])
    print "g_capacityDuals", g_capacityDuals
    #print "edgeTriplesInterdictable", edgeTriplesInterdictable
    
def resetSets():
    global Paths, ISets, edgeTriples, edgeTriplesInterdictable
    Paths = {}
    for commod in instance.commodities.keys():
        Paths[commod]=[]
    ISets = []
    
def run_Cormican_Benders_PathBased():
    global interfModelType
    interfModelType = '802.11-MAC-protocol'
    #print "jamPlaced", [n for n in g_JamGraph.nodes() if jamPlaced[n].X > FUZZ]
    throughput = solveJammingMasterProblem_Benders()
    #print "edgeInterdict", [e for e in G.edges() if edgeInterdict[e].X > FUZZ]
    print "throughput", throughput
    print "FINISHED"

def defineCutsPresent():
    global includeKnapsack, includeSuperValid, includeParetoOptimalBendersCuts, includeTrustRegion
    print "define cuts present"
    if algType == 'benders' or algType == 'benders-callback' or algType == 'branch-and-cut' or algType == 'b-and-c-reg' or algType == 'noInterferenceHeuristic':
        None
    elif algType == 'benders-ki' or algType == 'benders-callback-ki':
        includeKnapsack = True
    elif algType == 'benders-pareto' or algType == 'benders-callback-pareto':
        includeParetoOptimalBendersCuts = True
    elif algType == 'benders-svi':
        includeSuperValid = True
    elif algType == 'benders-tr':
        includeTrustRegion = True
    elif algType == 'benders-ki-pareto' or algType == 'benders-callback-ki-pareto':
        includeKnapsack = True
        includeParetoOptimalBendersCuts = True
    elif algType == 'benders-pareto-tr':
        includeParetoOptimalBendersCuts = True
        includeTrustRegion = True
    elif algType == 'benders-pareto-svi-tr':
        includeParetoOptimalBendersCuts = True
        includeSuperValid = True
        includeTrustRegion = True
    elif algType == 'benders-ki-pareto-svi':
        includeKnapsack = True
        includeSuperValid = True
        includeParetoOptimalBendersCuts = True
    elif algType == 'benders-ki-pareto-tr':
        includeKnapsack = True
        includeParetoOptimalBendersCuts = True
        includeTrustRegion = True
    elif algType == 'benders-ki-pareto-svi-tr':
        includeKnapsack = True
        includeSuperValid = True
        includeParetoOptimalBendersCuts = True
        includeTrustRegion = True
    else:
        raise Exception("algType:" + str(algType) + " invalid")
    
def doPrelimStuff():
    global exprFilePath, algType, instance, dataset
    exprFilePath, algType = executableModel.parseExprParamsFilePath()
    readInExperimentData(exprFilePath)
    defineCutsPresent()
    print "experiment data read"
    dataset = wnj_dataset.Wnj_Dataset(dataFilePath)
    print "dataset read"
    instance = wnj_problemInstance.Instance(exprFilePath, dataset)
    print "instance created"
    afterReadData()
    print "dominated", getLocationsDominatedBySingleLocation()

def afterReadData():
    global pathNumForCommod, jammersThatCanJamEdge
    initializeSets()
    #derived parameter values
    if(multipleRadiosPerNode):
        numRadiosPerNode = numChannels
    else:
        numRadiosPerNode = 1
    pathNumForCommod = dict([(commod, 0) for commod in instance.commodities.keys()])
    print "pathnumforCommod", pathNumForCommod
    
    jammersThatCanJamEdge = {}    
    for edgeInfo in instance.CGraph.edges(data = True):
        if jamVarsType == 'jam-vars':
            eTrip = (edgeInfo[0], edgeInfo[1], edgeInfo[2]['channel'])
            if edgeInfo[2]['edgeType'] == 'virtual':
                jammersThatCanJamEdge[eTrip] = []
            elif edgeInfo[2]['edgeType'] == 'real':
                jammersThatCanJamEdge[eTrip] = getJammersThatCanJamEdge(instance.CGraph, instance.jamGraph, eTrip, interfModelType)
    createShortestPathGraph()
    createEdgesInterdictableToIndexDict()
    createEdgesToIndexDict()
    createNodesToIndexDicts()
    createCommodToIndexDict()

def createShortestPathGraph():
    global shortestPathGraph
    shortestPathGraph = instance.CGraph.copy()
    for edgeInfo in shortestPathGraph.edges(data = True):
        edgeInfo[2]['weight'] = 1.0
    
def createEdgesInterdictableToIndexDict():
    global edgesToIndexInterdictDict
    edgesToIndexInterdictDict = {}
    count = 0
    for e in edgeTriplesInterdictable:
        edgesToIndexInterdictDict[e] = count
        count += 1
    return edgesToIndexInterdictDict

def createEdgesToIndexDict():
    global edgesToIndexDict
    edgesToIndexDict = {}
    count = 0
    for e in edgeTriples:
        edgesToIndexDict[e] = count
        count += 1
        
def getLocationsDominatedBySingleLocation():
    dominatedLocs = []
    for loc in instance.jamGraph.nodes():
        edgesJammedByLoc = getEdgesJammedByJammers(instance.CGraph, [loc], interfModelType)
        for loc2 in instance.jamGraph.nodes():
            edgesJammedByLoc2 = getEdgesJammedByJammers(instance.CGraph, [loc2], interfModelType)
            if loc2 != loc and set.issubset(set(edgesJammedByLoc2), set(edgesJammedByLoc)):
                #print "issubset", loc2, edgesJammedByLoc2, loc, edgesJammedByLoc
                if not loc2 in dominatedLocs:
                    dominatedLocs.append(loc2)
    return dominatedLocs
    
def readInExperimentData(path):
    global inst, clusterName #cluster stuff
    global modelType, jamVarsType, flowVarsType #model stuff
    global dataFilePath, datasetName, datasetIndex, numCommod, commodLocType, commodDemandType #dataset stuff
    global numJamLocs, commRange, infRange, multipleRadiosPerNode, numChannels, jamBudget, jamRange, linkCap, interfModelType, maxNumHopsMult, maxNumHops # instance params
    global time_limit, numThreads #algorithm params
    global numRadiosPerNode # derived params
    global databaseName
    d = etree.parse(open(path))
    modelType = str(d.xpath('//model/modelType[1]/text()')[0])
    jamVarsType = str(d.xpath('//model/jamVarsType[1]/text()')[0])
    flowVarsType = str(d.xpath('//model/flowVarsType[1]/text()')[0])
    
    # infrastructure stuff
    inst = str(d.xpath('//other/inst[1]/text()')[0])
    clusterName = str(d.xpath('//other/clusterName[1]/text()')[0])
    # dataset stuff
    dataFilePath = str(d.xpath('//dataset/path[1]/text()')[0])
    datasetName = str(d.xpath('//dataset/name[1]/text()')[0])
    datasetIndex = str(d.xpath('//dataset/datasetIndex[1]/text()')[0])
    numCommod = str(d.xpath('//dataset/numCommod[1]/text()')[0])
    commodLocType = str(d.xpath('//dataset/commodLocType[1]/text()')[0])
    commodDemandType = str(d.xpath('//dataset/commodDemandType[1]/text()')[0])

    #instance stuff
    numJamLocs = float(d.xpath('//instance/numJamLocs[1]/text()')[0])
    #commRange = float(d.xpath('//instance/commRange[1]/text()')[0])
    #infRange = float(d.xpath('//instance/infRange[1]/text()')[0])
    commRange = float(d.xpath('//instance/commRange[1]/text()')[0])
    infRange = float(d.xpath('//instance/infRange[1]/text()')[0])
    multipleRadiosPerNode = bool(int(d.xpath('//instance/multRadiosPerNode[1]/text()')[0]))
    numChannels = int(d.xpath('//instance/numChannels[1]/text()')[0])
    jamBudget= float(d.xpath('//instance/jamBudget[1]/text()')[0])
    jamRange = float(d.xpath('//instance/jamRange[1]/text()')[0])
    maxNumHopsMult = float(d.xpath('//instance/maxNumHopsMult[1]/text()')[0])
    
    databaseName = str(d.xpath('//other/databaseName[1]/text()')[0])
    linkCap = 1.0
    interfModelType = '802.11-MAC-protocol'
    #algorithm stuff
    time_limit = float(d.xpath('//algorithm/timeLimit[1]/text()')[0])
    numThreads = int(d.xpath('//algorithm/numThreads[1]/text()')[0])
    
if __name__ == "__main__":
    if(debug):
        print "!!!WARNING: DEBUG MODE!!!"
    doPrelimStuff()
    if algType in branchAndBoundVariations:
        runProtocolWithJamming_FullMIP_RowGenCallback() #START
    elif algType in bendersVariations:
        runProtocolWithJamming_Benders()
    elif algType == 'noInterferenceHeuristic':
        runNoInterferenceHeuristic_FullMIP_RowGenCallback()
    else:
        raise Exception("invalid algType", algType)