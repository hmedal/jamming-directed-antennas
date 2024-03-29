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

#needed for directional antenna file parsing
import csv
import math
import matplotlib.pyplot as plt
pi = math.pi

class Wnj_Dataset(object):
    '''
    classdocs
    '''
    ids = None
    numNodes = 0
    coors = None
    commodities = {}
    nodeOnlyGraph = nx.DiGraph()
    tempGraph = nx.DiGraph()
    pos = None
    G_dirAnt = None
    
    #directional antennas additions
    dist_a_b_new = 0
    transmissiondistance = 0
    
    def __init__(self, path):
        '''
        Constructor
        '''
        print "path", path
        self.readInDataset_and_CreateNodesAndCommodities(path)
    
    def readInDataset_and_CreateNodesAndCommodities(self, path):
        
        
        with open('/Users/wbl62/Desktop/directed-code/Transmitter_directed_med.csv', 'rU') as f:
            reader1 = csv.reader(f)
            mycsvlist1 = list(reader1)
            degreenumber = [x[0] for x in mycsvlist1]
            distancenumber = [x[1] for x in mycsvlist1]
            print "degreenumber:", degreenumber, 
            print "distancenumber:", distancenumber[180]
        
        
        G_dirAnt = nx.Graph()
        
        d = etree.parse(open(path))
        #print "d", d
        #facility information
        self.ids = [int(i) for i in d.xpath('//node/id[1]/text()')]
        self.numNodes = len(self.ids)
        self.coors = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
        self.distToClosest = [ast.literal_eval(i) for i in d.xpath('//node/distToClosest[1]/text()')]
        self.tcurr = [ast.literal_eval(i) for i in d.xpath('//node/tcurr[1]/text()')]
        #tcurr = [ast.literal_eval(i) for i in d.xpath('//node/tcurr[1]/text()')]
        #print "initial tcurr value is ", self.tcurr
        #print "tcurr1 value is ", tcurr1
        #print "distToClosest", self.distToClosest
        #sys.exit()
        self.trec = [ast.literal_eval(i) for i in d.xpath('//node/trec[1]/text()')]
        self.battCap = [ast.literal_eval(i) for i in d.xpath('//node/battCap[1]/text()')]
        coorsTuples = [(i[0], i[1]) for i in self.coors]
        #print "coors", self.coors
        #commented out: print "coorsTuples are as follows", coorsTuples
        
        #wrong: tcurrTuples = [i for i in self.tcurr]
        #wrong: print "tcurrTuples are", tcurrTuples
        #test comment
        self.commodities = {}
        
        #Start directed antenna processing to determine the nodes and edges in the graph based on the antenna pattern for each node
        coordinate1 = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
        #for above: related syntax from wnj_dataset is coors = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
        print "c1", coordinate1

        List_of_lists = map(list, coordinate1)
        print "list of lists", List_of_lists
        #Nodes_custom = ()
        #Nodes_custom[0] = coordinate1
        #Nodes_custom = None
        Nodes_custom1 = [x[0] for x in coordinate1]
        #x coordinates of nodes
        Nodes_custom2 = [x[1] for x in coordinate1]
        #y coordinates of nodes
        print "NodesCustom 1 and 2,", Nodes_custom1, "and", Nodes_custom2

        print "x0", x[0]
     
        #determine the distance between two nodes - I assume this function is correct
        def dist(x,y):   
            return np.sqrt(np.sum((x-y)**2))

        def newdist(x1,x2,y1,y2):
            return np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
        
        count = 0
        while count < len(Nodes_custom2):
            a = np.array((Nodes_custom1[count],Nodes_custom2[count]))#, Nodes_custom[count,2])) (not needed since we will use 2-D now)
            print "a", a
            
            #in the above line, "count" refers to the node number and the 0 or 1 refers to the x or y coordinate
            #G_dirAnt.add_node(count) #we are adding a node whose sequential number is "count" #,pos=(Nodes_custom[count,0],Nodes_custom[count,1])) #always will add node 0 
            print "count value", count
            dist_b_count=0
            #dist_b_count is basically the 2nd node, separately counted
            while dist_b_count < len(Nodes_custom2):
                b = np.array((Nodes_custom1[dist_b_count],Nodes_custom2[dist_b_count]))#,Nodes_custom[dist_b_count,2]))
#                 print "a", a
#                 print "b", b
#                 print "dist_b_count", dist_b_count
                dist_a_b = dist(a,b)
                dist_a_b_new = newdist(Nodes_custom1[count],Nodes_custom1[dist_b_count],Nodes_custom2[count],Nodes_custom2[dist_b_count])
                print "dist_ab:", dist_a_b, "dist_a_b_new:", dist_a_b_new
             
                #my code edits from 1-15-17 below
                xdist = Nodes_custom1[dist_b_count] - Nodes_custom1[count]
                ydist = Nodes_custom2[dist_b_count] - Nodes_custom2[count]
                print "xdist", xdist, "ydist", ydist
                #check the above
                #print "xdist", xdist
                #print "ydist", ydist
                #if xdist == 0 and ydist == 0:
                    #transmissiondistance = 0
                if xdist == 0 and ydist != 0:
                    if Nodes_custom2[dist_b_count] > Nodes_custom2[count]:
                        anglefound = pi/2
                    elif Nodes_custom2[dist_b_count] < Nodes_custom2[count]:
                        anglefound = (3/2) * pi
                elif ydist == 0 and xdist != 0:
                    if Nodes_custom1[dist_b_count] > Nodes_custom1[count]:
                        anglefound = 0
                    elif Nodes_custom1[dist_b_count] < Nodes_custom1[count]:
                        anglefound = pi
                elif xdist == 0 and ydist == 0:
                        anglefound = 0.000001 #placeholder value to deal with node compared to itself (i.e., same x and y coordinates)
                else:
                    if ydist > 0 and xdist > 0:
                        anglefound = np.arctan((ydist / xdist))
                    elif ydist > 0 and xdist < 0:
                        anglefound = (np.arctan((ydist / xdist))) + pi
                    elif ydist < 0 and xdist < 0:
                        anglefound = (np.arctan((ydist / xdist))) + pi
                    elif ydist < 0 and xdist > 0: 
                        anglefound = (np.arctan((ydist / xdist))) + (2*pi)
                    else:
                        print "Error"
                #print "anglefound in radians:", anglefound, ",", "anglefound in degrees", int(math.degrees(anglefound))
                #print "dist_a_b",dist_a_b, "dist_a_b_new", dist_a_b_new
                #if dist_a_b <= Transmitter[Nodes_custom[count,3]-1,1]:#verifying the transmitter is transmitting far enough
                #if Nodes_custom[count,2]>=Nodes_custom[dist_b_count,2]-1: #and Nodes_custom[count,2]<=Nodes_custom[dist_b_count,2]+1:
                #print 'ok' #verifying that along the z axis there is line of sight with the antenna within some tolerance
                if anglefound == 0.000001: #0.000001 is a placeholder value to deal with having to compare each node with itself
                        transmissiondistance = 0
                        print "anglefound is 0000001"
                else: #int(degreenumber[int(anglefound-1)]) == int(math.degrees(anglefound)): #math.degrees(math.int(anglefound)):
                                #print "counter", counter
                        #nonlocal(distancenumber)
                        transmissiondistance = distancenumber[int(math.degrees(anglefound))-1]
                
                print "transmission distance is ", transmissiondistance, "and dist_a_b_new is ", dist_a_b_new
                #print "transdist", transmissiondistance
                if transmissiondistance == 0:
                        #dist_b_count = dist_b_count + 1
                        print "same node", count, "=", dist_b_count
                elif float(transmissiondistance) >= dist_a_b_new: 
                        #print "test printing of transmissiondistance and dist_a_b", transmissiondistance, ",", dist_a_b_new
                        self.nodeOnlyGraph.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
                        self.nodeOnlyGraph.add_node(self.ids[dist_b_count], coor = coorsTuples[dist_b_count], distToClosest = self.distToClosest[dist_b_count], tcurr = self.tcurr[dist_b_count], trec = self.trec[dist_b_count], battCap = self.battCap[dist_b_count])
                        self.nodeOnlyGraph.add_edge(count, dist_b_count, dist = dist_a_b_new, capacity = 1, number = count)
                        
                        #using G_dirAnt for plotting
                        G_dirAnt.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
                        
                        G_dirAnt.add_edge(count, dist_b_count)#,weight = 1) #weight = power consumption
                        #self.nodeOnlyGraph.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])
                        #dist_b_count = dist_b_count + 1
                        #print "success with edge between", count, "and ", dist_b_count
                        #print "selfids", self.ids[count]
                        #print "count", count
                        #comment
                        
                elif float(transmissiondistance) < dist_a_b_new:
                        #dist_b_count = dist_b_count + 1
                        print "failure with edge between", count, "and ", dist_b_count
                dist_b_count = dist_b_count + 1
            count = count + 1
        
        #print "lentempBEFORE", len(self.nodeOnlyGraph.edges())
        
        for node in self.nodeOnlyGraph.edges():
            if self.nodeOnlyGraph.has_edge(node[1],node[0]) == False:
                self.nodeOnlyGraph.remove_edge(node[0],node[1])

        #print "lentempAFTER", len(self.nodeOnlyGraph.edges())
        #print "nodeOnlygraph" ,self.nodeOnlyGraph.nodes(data=True)
        #commented out the above line on 5-29-17
        
        #print "the tcurr values are ", self.nodeOnlyGraph.nodes(self.tcurr)[1]
        #tcurr = self.nodeOnlyGraph.nodes(self.tcurr[0])[0]
        #print "the tcurr values now are ", tcurr
        
        print"The following line begins the data from the new code: "
        #for m in range(0,o):
        #if 1 == 1:
        numberofitems = list(self.nodeOnlyGraph.nodes(self.tcurr[0]))
        print "The number of nodes for this file is", len(numberofitems), "."
        m = 0 #added
        #dict = {} #added
        
        while m < len(numberofitems): 
            #print "The current node number for which the following info applies is node", m, ":"
            tcurr1new = list(self.nodeOnlyGraph.nodes(self.tcurr[m]))
            #print "The values of all nodes (nodes 0 to m) are all included in this list/ dictionary: ", tcurr1new
            #print "original tcurr1[0] is ", self.nodeOnlyGraph.nodes(self.tcurr[m])[m]
            #print "the original list with a dictionary inside is ", tcurr1new
            #print "The list/ dictionary only including the info for node", m, "is: ", tcurr1new[m]
            ##tcurr1new.remove([m][0])
            listx = list(tcurr1new[m])
            #print "tuple has now become list listx, which is the following:", listx
            #listx.pop(0)
            del listx[0]
            #print "new list listx after removal of node # is", listx
            tcurr1new = tuple(listx)
            #del tcurr1new[m][0]
            #print "The dictionary within a tuple for node", m,  "after removal of the node # is", str(tcurr1new)
            #print "tcurr1new[0] is: ", tcurr1new[0]
            
            tcurr1new = str(tcurr1new[0]).strip('[]')
            #NOTE TO SELF: confirm thqt this should be [0] and not [m]
            
            tcurr1new = ast.literal_eval(tcurr1new)
            #print "The final dictionary containing tcurr, trec, etc., is: ", tcurr1new, "."
            print "Therefore, the tcurr value for node", m, "is", tcurr1new['tcurr'], "."
            #wrong: print "#ofkeys", len(tcurr1new)
            #m = m + 1
        
            #tcurr1new = list(self.nodeOnlyGraph.nodes(self.tcurr[0])[0])
            #print "the original tcurr1 is ", tcurr1new
            #tcurr1new.remove(0)
            #print "tcurr1new is ", str(tcurr1new)
            #tcurr1new=str(tcurr1new[0]).strip('[]')
            #tcurr1new = ast.literal_eval(tcurr1new)
            #print "and the final tcurr value is ", tcurr1new['tcurr']
            
            #print "G_dirAnt nodes", G_dirAnt.nodes(data=True)
            #print "G_dirAnt edges", G_dirAnt.edges()
            #commented out the previous 2 lines on 5-29-17
            #sys.exit()
        
            print "NEW CODE BEGINS HERE"
            
            
            if tcurr1new['tcurr'] < 50:
                with open('/Users/wbl62/Desktop/directed-code/Transmitter_directed_med.csv', 'rU') as f:
                    reader1 = csv.reader(f)
                    mycsvlist1 = list(reader1)
                    degreenumber = [x[0] for x in mycsvlist1]
                    distancenumber = [x[1] for x in mycsvlist1]
                    print "degreenumber:", degreenumber, 
                    print "distancenumber:", distancenumber[180]
            else:
                with open('/Users/wbl62/Desktop/directed-code/Transmitter_directed_high.csv', 'rU') as f:
                    reader1 = csv.reader(f)
                    mycsvlist1 = list(reader1)
                    degreenumber = [x[0] for x in mycsvlist1]
                    distancenumber = [x[1] for x in mycsvlist1]
                    print "degreenumber:", degreenumber, 
                    print "distancenumber:", distancenumber[180] #example only to use [180]
            
            
            G_dirAnt = nx.Graph()
            
            d = etree.parse(open(path))
                    #print "d", d
                    #facility information
            self.ids = [int(i) for i in d.xpath('//node/id[1]/text()')]
            self.numNodes = len(self.ids)
            self.coors = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
            self.distToClosest = [ast.literal_eval(i) for i in d.xpath('//node/distToClosest[1]/text()')]
            self.tcurr = [ast.literal_eval(i) for i in d.xpath('//node/tcurr[1]/text()')]
                #tcurr = [ast.literal_eval(i) for i in d.xpath('//node/tcurr[1]/text()')]
                #print "initial tcurr value is ", self.tcurr
                #print "tcurr1 value is ", tcurr1
                #print "distToClosest", self.distToClosest
                #sys.exit()
            self.trec = [ast.literal_eval(i) for i in d.xpath('//node/trec[1]/text()')]
            self.battCap = [ast.literal_eval(i) for i in d.xpath('//node/battCap[1]/text()')]
            coorsTuples = [(i[0], i[1]) for i in self.coors]
                #print "coors", self.coors
                #commented out: print "coorsTuples are as follows", coorsTuples
            
                #wrong: tcurrTuples = [i for i in self.tcurr]
                #wrong: print "tcurrTuples are", tcurrTuples
                #test comment
            self.commodities = {}
            
                #Start directed antenna processing to determine the nodes and edges in the graph based on the antenna pattern for each node
            coordinate1 = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
                    #for above: related syntax from wnj_dataset is coors = [ast.literal_eval(i) for i in d.xpath('//node/coor[1]/text()')]
            print "c1", coordinate1
    
            List_of_lists = map(list, coordinate1)
            print "list of lists", List_of_lists
                    #Nodes_custom = ()
                    #Nodes_custom[0] = coordinate1
                    #Nodes_custom = None
            Nodes_custom1 = [x[0] for x in coordinate1]
                    #x coordinates of nodes
            Nodes_custom2 = [x[1] for x in coordinate1]
                    #y coordinates of nodes
            print "Nodes_custom 1 and 2 values are,", Nodes_custom1, "and", Nodes_custom2
    
            print "x0", x[0]
         
                #determine the distance between two nodes - I assume this function is correct
                #REPETITIVE BELOW - COMMENTED OUT
                #def dist(x,y):   
                    #return np.sqrt(np.sum((x-y)**2))
    
                #def newdist(x1,x2,y1,y2):
                    #return np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
            
            count = 0
            while count < len(Nodes_custom2):
                a = np.array((Nodes_custom1[count],Nodes_custom2[count]))#, Nodes_custom[count,2])) (not needed since we will use 2-D now)
                print "a", a
                
                #in the above line, "count" refers to the node number and the 0 or 1 refers to the x or y coordinate
                #G_dirAnt.add_node(count) #we are adding a node whose sequential number is "count" #,pos=(Nodes_custom[count,0],Nodes_custom[count,1])) #always will add node 0 
                print "count value", count
                dist_b_count=0
                #dist_b_count is basically the 2nd node, separately counted
                while dist_b_count < len(Nodes_custom2):
                    b = np.array((Nodes_custom1[dist_b_count],Nodes_custom2[dist_b_count]))#,Nodes_custom[dist_b_count,2]))
    #                 print "a", a
    #                 print "b", b
    #                 print "dist_b_count", dist_b_count
                    dist_a_b = dist(a,b)
                    dist_a_b_new = newdist(Nodes_custom1[count],Nodes_custom1[dist_b_count],Nodes_custom2[count],Nodes_custom2[dist_b_count])
                    print "dist_ab:", dist_a_b, "dist_a_b_new:", dist_a_b_new
                 
                    #my code edits from 1-15-17 below
                    xdist = Nodes_custom1[dist_b_count] - Nodes_custom1[count]
                    ydist = Nodes_custom2[dist_b_count] - Nodes_custom2[count]
                    print "xdist", xdist, "ydist", ydist
                    #check the above
                    #print "xdist", xdist
                    #print "ydist", ydist
                    #if xdist == 0 and ydist == 0:
                        #transmissiondistance = 0
                    if xdist == 0 and ydist != 0:
                        if Nodes_custom2[dist_b_count] > Nodes_custom2[count]:
                            anglefound = pi/2
                        elif Nodes_custom2[dist_b_count] < Nodes_custom2[count]:
                            anglefound = (3/2) * pi
                    elif ydist == 0 and xdist != 0:
                        if Nodes_custom1[dist_b_count] > Nodes_custom1[count]:
                            anglefound = 0
                        elif Nodes_custom1[dist_b_count] < Nodes_custom1[count]:
                            anglefound = pi
                    elif xdist == 0 and ydist == 0:
                            anglefound = 0.000001 #placeholder value to deal with node compared to itself (i.e., same x and y coordinates)
                    else:
                        if ydist > 0 and xdist > 0:
                            anglefound = np.arctan((ydist / xdist))
                        elif ydist > 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + pi
                        elif ydist < 0 and xdist < 0:
                            anglefound = (np.arctan((ydist / xdist))) + pi
                        elif ydist < 0 and xdist > 0: 
                            anglefound = (np.arctan((ydist / xdist))) + (2*pi)
                        else:
                            print "Error"
                    print "anglefound in radians:", anglefound, ",", "anglefound in degrees", int(math.degrees(anglefound))
                    #print "dist_a_b",dist_a_b, "dist_a_b_new", dist_a_b_new
                    #if dist_a_b <= Transmitter[Nodes_custom[count,3]-1,1]:#verifying the transmitter is transmitting far enough
                    #if Nodes_custom[count,2]>=Nodes_custom[dist_b_count,2]-1: #and Nodes_custom[count,2]<=Nodes_custom[dist_b_count,2]+1:
                    #print 'ok' #verifying that along the z axis there is line of sight with the antenna within some tolerance
                    if anglefound == 0.000001: #0.000001 is a placeholder value to deal with having to compare each node with itself
                            transmissiondistance = 0
                            print "anglefound is 0000001"
                    else: #int(degreenumber[int(anglefound-1)]) == int(math.degrees(anglefound)): #math.degrees(math.int(anglefound)):
                                    #print "counter", counter
                            transmissiondistance = distancenumber[int(math.degrees(anglefound))-1]
                    
                    print "transmission distance is ", transmissiondistance, "and dist_a_b_new is ", dist_a_b_new
                    #print "transdist", transmissiondistance
                    if transmissiondistance == 0:
                            #dist_b_count = dist_b_count + 1
                            print "same node", count, "=", dist_b_count
                    elif float(transmissiondistance) >= dist_a_b_new: 
                            #print "test printing of transmissiondistance and dist_a_b", transmissiondistance, ",", dist_a_b_new
                            self.nodeOnlyGraph.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
                            self.nodeOnlyGraph.add_node(self.ids[dist_b_count], coor = coorsTuples[dist_b_count], distToClosest = self.distToClosest[dist_b_count], tcurr = self.tcurr[dist_b_count], trec = self.trec[dist_b_count], battCap = self.battCap[dist_b_count])
                            self.nodeOnlyGraph.add_edge(count, dist_b_count, dist = dist_a_b_new, capacity = 1, number = count)
                            
                            #using G_dirAnt for plotting
                            G_dirAnt.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])#,pos=(Nodes_custom[dist_b_count,0],Nodes_custom[dist_b_count,1]))
                            
                            G_dirAnt.add_edge(count, dist_b_count)#,weight = 1) #weight = power consumption
                            #self.nodeOnlyGraph.add_node(self.ids[count], coor = coorsTuples[count], distToClosest = self.distToClosest[count], tcurr = self.tcurr[count], trec = self.trec[count], battCap = self.battCap[count])
                            #dist_b_count = dist_b_count + 1
                            #print "success with edge between", count, "and ", dist_b_count
                            #print "selfids", self.ids[count]
                            #print "count", count
                            #comment
                            
                    elif float(transmissiondistance) < dist_a_b_new:
                            #dist_b_count = dist_b_count + 1
                            print "failure with edge between", count, "and ", dist_b_count
                    dist_b_count = dist_b_count + 1
                count = count + 1
            m = m + 1
        
        #print "lentempBEFORE", len(self.nodeOnlyGraph.edges())
        
        for node in self.nodeOnlyGraph.edges():
            if self.nodeOnlyGraph.has_edge(node[1],node[0]) == False:
                self.nodeOnlyGraph.remove_edge(node[0],node[1])
        
        print "NEW CODE ENDS HERE"
        
        #for index in range(0,len(self.ids)):#G_dirAnt.nodes():
         #   self.nodeOnlyGraph.add_node(self.ids[index], coor = coorsTuples[index], distToClosest = self.distToClosest[index], tcurr = self.tcurr[index], trec = self.trec[index], battCap = self.battCap[index])
            #self.nodeOnlyGraph.add_node(index)
        for commod in d.xpath('//odPair'):
            id = int(commod.xpath('./id/text()')[0])
            self.commodities[id] = {}
            origin = int(commod.xpath('./origin/text()')[0])
            destination = int(commod.xpath('./destination/text()')[0])
            demand = float(commod.xpath('./demand/text()')[0])
            self.commodities[id]['odPair'] = (origin, destination)
            self.commodities[id]['demand'] = demand 
            
        #print "nodesonlygraphedges", self.nodeOnlyGraph.edges(data=True)
        #print "SORTEDnodesonlygraphedges", self.nodeOnlyGraph.edges(data=True)
        #print "G_dirAnt nodes", G_dirAnt.nodes(data=True)
        #commented out the previous 3 lines on 5-29-17
        
        #sys.exit()
#         nx.draw(G_dirAnt,with_labels=True)
#         #just drawing the graph with the new nodes & edges
#         nx.write_gml(G_dirAnt, "/Users/jhuff/Desktop/directed-code/transmitter_2_results.txt")
#         print "gml_written"
#         plt.show()
