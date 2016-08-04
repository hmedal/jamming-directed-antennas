def createJamGraph(gridSize):
    JammingGraph = createPotentialJammingLocations(gridSize)
    jammingSolution = [0] * (gridSize*2*numRadiosPerNode)**2
    
    counter = 0
    for node in JammingGraph.nodes():
        #print "jamming node", node
        JammingGraph.node[node]['selected'] = jammingSolution[counter]
        counter += 1
    #print "jammed", [node for node in JammingGraph.nodes() if JammingGraph.node[node]['selected'] > 0.0]
    return JammingGraph

def createCommodities_Base(gridSize, numCommodities):
    global commodities
    bottomLeft = 0
    bottomLeftShifted = gridSize + 1
    topLeft = gridSize - 1
    topLeftShifted = gridSize * 2 - 2
    topRight = gridSize * gridSize - 1
    topRightShifted = gridSize * (gridSize - 1) - 2
    bottomRight = gridSize * (gridSize - 1)
    bottomRightShifted = gridSize * (gridSize - 2) + 1
    print "corners", bottomLeft, topLeft, topRight, bottomRight
    odPairs = []
    shifted = False
    if(shifted):
        odPairs.append((bottomLeftShifted, topRightShifted))
        odPairs.append((bottomLeftShifted, topLeftShifted))
        odPairs.append((bottomLeftShifted, bottomRightShifted))
        odPairs.append((topLeftShifted, topRightShifted))
        odPairs.append((topLeftShifted, bottomRightShifted))
        odPairs.append((topRightShifted, bottomRightShifted))
    else:
        odPairs.append((bottomLeft, topRight))
        odPairs.append((bottomLeft, topLeft))
        odPairs.append((bottomLeft, bottomRight))
        odPairs.append((topLeft, topRight))
        odPairs.append((topLeft, bottomRight))
        odPairs.append((topRight, bottomRight))
    counter = 0
    for index in range(numCommodities):
        #for radioIndex in range(numRadiosPerNode):
        commodities[counter] = {}
        commodities[counter]['odPair'] = (odPairs[index][0], odPairs[index][1])
        commodities[counter]['demand'] = 2.0
        counter += 1
    print "commodities", commodities
    
# def getSignalStrength(power, distance):
#     return power/distance**2
# 
# def getSNR(G, edge):
#     return getSignalStrength(G.node[edge[0]]['power'], G.edge[edge[0]][edge[1]]['dist'])/ambient_noise
# 
# def edgeExistsInConnectivityGraph(G, edgeInfoAll, interfModelType, radioTypeIndex, channelIndex):
#     edge = (edgeInfoAll[0], edgeInfoAll[1])
#     edgeInfo = edgeInfoAll[2]
#     #edgeIndex = radioTypeIndex * numChannels + channelIndex
#     #nodeIndex = edge[0] * numRadiosPerNode + radioTypeIndex
#     #print "indices", edgeIndex, nodeIndex
#     if((interfModelType == 'simple-protocol') or (interfModelType == '802.11-MAC-protocol') or (interfModelType == 'none')):
#         if(edge[0] != edge[1] and edgeInfo['dist'] <= G.node[edge[0]]['commRange']):
#             return True
#     elif((interfModelType == 'simple-physical')):
#         if(edge[0] != edge[1] and getSNR(G, edge) >= snr_threshold):
#             return True
#         
# def getFractionOfMaxPermissibleNoise(G, edge1, edge2):
#     i=edge1[0]
#     j=edge1[1]
#     p=edge2[0]
#     if(p != j):
#         #print "getFractionOfMaxPermissibleNoise", edge1, edge2
#         #print G.edges()
#         #print p, j, G.edge[p][j]['dist']
#         if((p,j) in G.edges()):
#             distance = G.edge[p][j]['dist']
#         else:
#             distance = np.linalg.norm(np.array(originalGraph.node[p]['coor']) - np.array(originalGraph.node[j]['coor']))
#     else:
#         distance = 0.01
#     ss_pj = getSignalStrength(G.node[p]['power'],distance)
#     ss_ij = getSignalStrength(G.node[i]['power'], G.edge[i][j]['dist'])
#     #print edge1, edge2, distance, ss_pj, ss_pj / ((ss_ij/snr_threshold) - ambient_noise)
#     return  ss_pj / ((ss_ij/snr_threshold) - ambient_noise)

# def getFractionOfMaxPermissibleNoise_DueToJammer(G, JammingGraph, edge1, jamNode):
#     i=edge1[0]
#     j=edge1[1]
#     #print "i", i, "j", j, G.nodes(data=True)
#     distance = np.linalg.norm(np.array(jamNode) - np.array(originalGraph.node[j]['coor']))
#     ss_pj = getSignalStrength(JammingGraph.node[jamNode]['power'],distance)
#     ss_ij = getSignalStrength(originalGraph.node[i]['power'], originalGraph.edge[i][j]['dist'])
#     #print edge1, edge2, distance, ss_pj, ss_pj / ((ss_ij/snr_threshold) - ambient_noise)
#     return  ss_pj / ((ss_ij/snr_threshold) - ambient_noise)
    
# def cannotBeActiveSimultaneously_simpleProtocolModel(G, edge1, edge2):
#     i=edge1[0]
#     j=edge1[1]
#     p=edge2[0]
#     q=edge2[1]
#     print "cannotBeActiveSimultaneously_simpleProtocolModel", i, j, p, q
#     #print i, j, p, q
#     if(G.has_edge(i,q) and G.edge[i][q]['dist'] <= G.node[i]['interferenceRange']):
#         return True
#     if(G.has_edge(p,j) and G.edge[p][j]['dist'] <= G.node[p]['interferenceRange']):
#         return True
#     
# def cannotBeActiveSimultaneously_802_11_ProtocolModel(G, edgeInfo1, edgeInfo2):
#     i=edgeInfo1[0]
#     j=edgeInfo1[1]
#     radioType1 = edgeInfo1[2]['radioType']
#     channel1 = edgeInfo1[2]['channel']
#     p=edgeInfo2[0]
#     q=edgeInfo2[1]
#     radioType2 = edgeInfo2[2]['radioType']
#     channel2 = edgeInfo2[2]['channel']
#     #print "cannotBeActiveSimultaneously_802_11_ProtocolModel", edge1, edge2
#     #print i, q, G.edge[i][q]['dist']
#     #print i, j, radioType1, channel1, "...", p, q, radioType2, channel2
#     tuplesOfNodePairs = [(i,q), (q,i), (i,p), (p,i), (j,p), (p,j), (j,q), (q,j)]
#     #print "tuplesOfNodePairs", i,j, p, q, tuplesOfNodePairs
#     if(i == p or j == q):#same head and/or tail (interference due to same radio)
#         if(radioType1 == radioType2):# cannot have same radio if same head or tail
#             return True
#         else:
#             return False
#     else:
#         counter = 0
#         for (a,b) in tuplesOfNodePairs:
#             if(G.has_edge(a,b)):
#                 print G.edge[a][b][0], G.node[a]
#                 if(G.edge[a][b][0]['dist'] <= G.node[a]['interferenceRange'] and channel1 == channel2):
#                     #print "   added b/c of dist:", a, b, channel1, channel2, counter, 'dist', G.edge[a][b][0]['dist']
#                     return True
#             counter += 1
            
# def cannotBeActiveSimultaneously(G, edge1, edge2, interfModelType):
#     if(interfModelType == 'simple-protocol'):
#         return cannotBeActiveSimultaneously_simpleProtocolModel(G, edge1, edge2)
#     elif(interfModelType == '802.11-MAC-protocol'):
#         #print "is 802-11"
#         return cannotBeActiveSimultaneously_802_11_ProtocolModel(G, edge1, edge2)

# def createConflictGraph_Protocol(G, interfModelType):
#     #print "createConflictGraph_Protocol", G.edges(data=True)
#     ConflictGraph = nx.DiGraph()
#     for edgeInfo1 in G.edges(data = True):
#         for edgeInfo2 in G.edges(data = True):
#             if(edgeInfo1 != edgeInfo2):
#                 if(cannotBeActiveSimultaneously(G, edgeInfo1, edgeInfo2, interfModelType) is True):
#                     edge1 = (edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['channel'])
#                     edge2 = (edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['channel'])
#                     #print "add:", edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['radioType'], edgeInfo1[2]['channel'], "***", edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['radioType'], edgeInfo2[2]['channel']
#                     ConflictGraph.add_edge(edge1, edge2)
#                 #else:
#                     #print "don't add:", edgeInfo1[0], edgeInfo1[1], edgeInfo1[2]['radioType'], edgeInfo1[2]['channel'], "***", edgeInfo2[0], edgeInfo2[1], edgeInfo2[2]['radioType'], edgeInfo2[2]['channel']
#     #nx.draw(ConflictGraph)
#     #plt.show() # display
#     return ConflictGraph
# 
# def createConflictGraph_PhysicalModel(G, interfModelType):
#     #print G.edges(data=True)
#     ConflictGraph = nx.DiGraph()
#     for edge1 in G.edges():
#         for edge2 in G.edges():
#             if(edge1 != edge2):
#                 if(interfModelType == 'simple-physical'):
#                     ConflictGraph.add_edge(edge1, edge2, weight = getFractionOfMaxPermissibleNoise(G, edge1, edge2))
#     #nx.draw(ConflictGraph)
#     #plt.show() # display
#     return ConflictGraph

# def createConflictGraph(G, interfModelType):
#     if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
#         return createConflictGraph_Protocol(G, interfModelType)
#     elif(interfModelType == 'simple-physical'):
#         return createConflictGraph_PhysicalModel(G, interfModelType)

        
def createThroughputModel_OLD(G, source, destination, ISets, jammingGraph, interfModelType):
    global gurobiModel, iSetUsage, capConstraints, usageConstr, jamConstrs
    numISets = len(ISets)
    gurobiModel = gurobipy.Model("myLP")
    try:
        # Create variables
        flowVars = dict([(edge, gurobiModel.addVar(0, name="x_"+str(edge[0])+","+str(edge[1]))) for edge in G.edges()])
        iSetUsage = [gurobiModel.addVar(0, 1, name="lamda_"+str(k)) for k in range(numISets)]
        gurobiModel.update() # Integrate new variables
        # Set objective
        gurobiModel.setObjective(sum([flowVars[edge] for edge in G.edges([source])]), gurobipy.GRB.MAXIMIZE)
        # flow balance constraints
        #print "src", source, "dest", destination
        #print "nodes: ", [node for node in G.nodes() if node != source and node != destination]
        [gurobiModel.addConstr(sum([flowVars[edge] for edge in G.edges([node])])
                               == sum([flowVars[edge] for edge in G.in_edges([node])]), "flowBal_"+str(node)) 
                                    for node in G.nodes() if (node != source and node != destination)]
        # flow into source
        gurobiModel.addConstr(sum([flowVars[edge] for edge in G.in_edges([source])]) == 0.0, "flowIntoSource")
        # flow out of destination
        gurobiModel.addConstr(sum([flowVars[edge] for edge in G.out_edges([destination])]) == 0.0, "flowOutOfDest")
        # capacity
        capConstraints = dict([(edge, gurobiModel.addConstr(flowVars[edge] <= 
                               sum([iSetUsage[k] * int(edge in ISets[k]) * G.edge[edge[0]][edge[1]]['capacity'] for k in range(numISets)]), 
                                   "capacity_"+str(edge[0])+","+str(edge[1]))) for edge in G.edges()])
        # usage
        usageConstr = gurobiModel.addConstr(sum([iSetUsage[k] for k in range(numISets)]) <= 1, "iSetUsage")
        # jamming
        if(interfModelType == 'simple-protocol' or interfModelType == '802.11-MAC-protocol'):
            jamConstrs = dict([(node, gurobiModel.addConstr(sum([int(isEdgeSetJammedByJammer(G, ISets[k], node, interfModelType)) * iSetUsage[k] 
                                                                 for k in range(numISets)]) <= 
                                                 1.0 - float(jammingGraph.node[node]['selected']), 
                                                     "interdict_" + str(node)))  for node in jammingGraph.nodes()])
        elif(interfModelType == 'simple-physical'):
            raise Exception("not implemented yet")
        gurobiModel.update() # integrate objective and constraints
        #print jamConstrs
        gurobiModel.setParam('OutputFlag', False ) #turn output off
        if(writeToFile):
            gurobiModel.write("/home/hmedal/Documents/Temp/JainInitial_"+interfModelType + ".lp")
    except gurobipy.GurobiError as e:
        print str(e)
        
def create2D_Grid_Graph(gridSize, show = False, multipleRadiosPerNodeArg = False, numChannels = 1, linkCap = 1.0, commRangeArg = 1.0, infRange = 1.0):
    global pos, numRadiosPerNode
    G = nx.grid_2d_graph(gridSize, gridSize)
    #print "edges", G.edges()
    G = nx.convert.convert_to_directed(G)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted", label_attribute='coor')
    posG = {}
    counter = 0
    for i in range(gridSize):
        for j in range(gridSize):
            posG[counter] = (i, j)
            counter += 1
    if(show):
        nx.draw(G, posG)
        plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/within_node_graph_v2.jpg") # save as jpg
        plt.show() # display
    for i in G.nodes():
        G.node[i]['power']= 1.0
        G.node[i]['commRange']= commRangeArg
        G.node[i]['interferenceRange']= infRange
    counter = 0
    for e in sorted(G.edges()):
        #print "dist", e[0], e[1], np.linalg.norm(np.array(G.node[e[0]]['coor']) - np.array(G.node[e[1]]['coor']))
        G.edge[e[0]][e[1]]['dist']= np.linalg.norm(np.array(G.node[e[0]]['coor']) - np.array(G.node[e[1]]['coor']))
        G.edge[e[0]][e[1]]['number']= counter
        G.edge[e[0]][e[1]]['capacity']=1.0
        counter += 1
    MDG = nx.MultiDiGraph()
    #print "nodes", G.nodes()
    posMDG = {}
    counter = 0
    if(multipleRadiosPerNodeArg):
        numRadiosPerNode = numChannels
    else:
        numRadiosPerNode = 1
    for radioTypeIndex in range(numRadiosPerNode):
        for n in sorted(G.nodes()):
            MDG.add_node((n, radioTypeIndex), power = G.node[n]['power'], commRange = G.node[n]['commRange'], coor = G.node[n]['coor'], 
                         interferenceRange = G.node[n]['interferenceRange'], radioType = radioTypeIndex, index = counter)
            counter += 1
            MDG.add_node((n), coor = G.node[n]['coor'], index = counter)#add super node
            MDG.add_edge((n), (n, radioTypeIndex), radioTypeIndex)
            posMDG[(n, radioTypeIndex)] = posG[n]
            posMDG[(n)] = posG[n]
            counter += 1
    for n in sorted(G.nodes()):
        MDG.add_node((n), coor = G.node[n]['coor'], interferenceRange = 0.0, commRange = float('inf'))#add super node
        posMDG[(n)] = posG[n]
        for radioTypeIndex in range(numRadiosPerNode):
            MDG.add_edge((n), (n, radioTypeIndex), radioTypeIndex, dist = 0.0, radioType = radioTypeIndex, capacity = 1.0, channel = -1, origNum = -1)
            MDG.add_edge((n, radioTypeIndex), (n), radioTypeIndex, dist = 0.0, radioType = radioTypeIndex, capacity = 1.0, channel = -1, origNum = -1)
    for e in sorted(G.edges()):
        distBetween = G.edge[e[0]][e[1]]['dist']
        myNum = G.edge[e[0]][e[1]]['number']
        for radioTypeIndex in range(numRadiosPerNode):
            if(multipleRadiosPerNodeArg):
                MDG.add_edge((e[0], radioTypeIndex), (e[1], radioTypeIndex), radioType = radioTypeIndex, channel = radioTypeIndex, dist = distBetween, capacity = linkCap, origNum = myNum)
                #print "add", e, radioTypeIndex, radioTypeIndex, MDG.edge[e[0]][e[1]][radioTypeIndex]
            else:
                for channelIndex in range(numChannels):
                    MDG.add_edge((e[0], radioTypeIndex), (e[1], radioTypeIndex), radioType = radioTypeIndex, channel = channelIndex, dist = distBetween, capacity = linkCap, origNum = myNum)
                    #print "add", e, radioTypeIndex, channelIndex, MDG.edge[(e[0], radioTypeIndex)][(e[1], radioTypeIndex)], MDG.edges(data = True)
    #print "MDG nodes", MDG.nodes(data = True)
    #print "MDG edges", MDG.edges(data = True)
    #for e in MDG.edges(data = True):
    #    print "MDG edge", e
    #nx.draw(MDG, posMDG, MDG.nodes())
    #plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/within_node_graph_v2.jpg") # save as jpg
    #plt.show() # display
    #for node in MDG.nodes(data = True):
    #    print "MDG node", node
    #print "MDG edges", MDG.edges(data = True)
    pos = posMDG
    return MDG

def createRandomGraph(n, commRange, intRange, cubeSize, dim = 2, pos = None, show = False, multRadiosPerNodeArg = False, numChannels = 1):
    G=nx.Graph()
    G.name="Random Graph"
    G.add_nodes_from(range(n))
    if pos is None:
        # random positions 
        for n in G:
            G.node[n]['pos']=[cubeSize * np.random.random() for i in range(0,dim)]
    else:
        nx.set_node_attributes(G,'pos',pos)
    
    for i in G.nodes():
        G.node[i]['power']= 1.0
        G.node[i]['commRange']= commRange
        G.node[i]['interferenceRange']= intRange
        
    nodes = G.nodes(data=True)
    while nodes:
        u,du = nodes.pop()
        pu = du['pos']
        for v,dv in nodes:
            pv = dv['pos']
            d = sum(((a-b)**2 for a,b in zip(pu,pv)))
            if d <= commRange**2:
                G.add_edge(u,v)
    
    print G.nodes(data=True)
    pos = dict([(n, G.node[n]['pos']) for n in G.nodes()])
    if(show):
        nx.draw(G, pos)
        plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/within_node_graph_v2.jpg") # save as jpg
        plt.show() # display
    return G

# def createConnectivityGraph(gridSize, interfModelType, graphType, show, multipleRadiosPerNodeArg = False, numChannelsArg = 1):
#     global originalGraph
#     #print "equals-createConnectivityGraph", (interfModelType == '802.11-protocol')
#     if(graphType == 'grid'):
#         originalGraph = create2D_Grid_Graph(gridSize, show, multipleRadiosPerNode, numChannelsArg)
#     elif(graphType == 'random'):
#         originalGraph = createRandomGraph(gridSize**2, 2.0, 2.0, gridSize, show, multipleRadiosPerNode, numChannelsArg)
#     CGraph = nx.MultiDiGraph()
#     #CGraph = nx.convert.convert_to_directed(CGraph)
#     #CGraph.add_nodes_from(originalGraph)
#     #print "nodes", originalGraph.nodes()
#     for i in originalGraph.nodes():
#         #print "node", i, originalGraph.node[i]
#         if isinstance(i, int):
#             CGraph.add_node(i, type = 'super', interferenceRange = 0.0, commRange = float('inf'))
#         else:
#             CGraph.add_node(i, power = originalGraph.node[i]['power'], commRange = originalGraph.node[i]['commRange'], 
#                              interferenceRange = originalGraph.node[i]['interferenceRange'], radioType = originalGraph.node[i]['radioType'])
#         
#         #CGraph.add_node(i, power = originalGraph.node[i]['power']
#         #CGraph.node[i]['commRange']= originalGraph.node[i]['commRange']
#         #CGraph.node[i]['interferenceRange']= originalGraph.node[i]['interferenceRange']
#        # CGraph.node[i]['radioType']= originalGraph.node[i]['radioType']
#     #print CGraph.nodes(data=True)
#     for edge in sorted(originalGraph.edges(data = True)):
#         edgeAttr = edge[2]
#         #print "edge", edge, edgeAttr
#         #for radioTypeIndex in range(numRadiosPerNode):
#             #for channelIndex in range(numChannels):
#             #print edge, edgeExistsInConnectivityGraph(G, edge, interfModelType)
#         if(edgeExistsInConnectivityGraph(originalGraph, edge, interfModelType, edgeAttr['radioType'], edgeAttr['channel'])):
#             #edgeIndex = edgeAttr.
#             myDist = edgeAttr['dist']
#             myNum = edgeAttr['origNum']
#             myCap = edgeAttr['capacity']
#             #print "added", edge, edgeIndex, radioTypeIndex, channelIndex, myDist, myNum, myCap
#             CGraph.add_edge(edge[0], edge[1], radioType = edgeAttr['radioType'], channel = edgeAttr['channel'], dist = myDist, origNum = myNum, capacity = myCap)
#             #print CGraph.edges(data = True)
#             #CGraph.edge[edge[0]][edge[1]]['dist'] = 
#             #CGraph.edge[edge[0]][edge[1]]['number'] = originalGraph.edge[edge[0]][edge[1]]['number']
#             #CGraph.edge[edge[0]][edge[1]]['capacity'] = originalGraph.edge[edge[0]][edge[1]]['capacity']
#     #print "CGraph nodes", CGraph.nodes(data=True)
#     #print "CGraph edges", CGraph.edges(data=True)
#     return CGraph

# def JainExample(gridSize, source, dest, interfModelType):
#     #G = createConnectivityGraph(gridSize, interfModelType)
#     wiredNetworkFlow = nx.max_flow(G, source, dest)
#     print "wiredNetworkFlow:", wiredNetworkFlow
#     
#     #ConflictGraph = createConflictGraph(G, interfModelType)
# 
#     ISets = []
#     #capacityDuals = dict([(edge, 1.0) for edge in G.edges()])
#     count = 0
#     #maxWeight, selected = maxWtIndepSet(ConflictGraph, capacityDuals)
#     ISets.append(selected)
#     #print "selected ", selected
#     jammerGraph = createPotentialJammingLocations()
#     createThroughputModel(G, source, dest, [selected], jammerGraph)
#     while True:
#         print "ITERATION", count
#         #maxWeight, selected = maxWtIndepSet(ConflictGraph, capacityDuals)
#         ISets.append(selected)
#         throughput, capacityDuals, usageDual = addISetAsCol_AndSolveTHProb(G, jammerGraph, ISets, selected, count)
#         print "selected ", selected
#         if((maxWeight - usageDual) <= 0.0001):
#             break
#         count += 1

def createPotentialJammingLocations(gridSize):
    global g_JamGraph, numRadiosPerNode
    g_JamGraph = nx.Graph()
    pos = {}
    if(multipleRadiosPerNode):
        numRadiosPerNode = numChannels
    else:
        numRadiosPerNode = 1
    #print "numRadiosPerNode", numRadiosPerNode
    counter = 0
    for radioTypeIndex in range(numRadiosPerNode):
        for i in np.arange(0, gridSize, 0.5):
            for j in np.arange(0, gridSize, 0.5):
                g_JamGraph.add_node((i, j, radioTypeIndex), name = counter, selected = 0.0, power = 1.0, cost = 1.0, interferenceRange = 1.0, 
                                    radioType = radioTypeIndex)
                #print "add", i, j, counter, radioTypeIndex, g_JamGraph.nodes(data=True)
                pos[(i, j)] = (i, j)
                counter += 1
    #for nodeInfo in g_JamGraph.nodes(data=True):
    #    print "nodeInfo", nodeInfo 
    #nx.draw(G, pos)
    #plt.savefig("/home/hmedal/Dropbox/2_msu/PROJECT_ncitec_safe_multimodal/figures/within_node_graph_v2.jpg") # save as jpg
    #plt.show() # display
    return g_JamGraph