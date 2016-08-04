'''
Created on Sep 30, 2013

@author: hmedal
'''

import sqlite3 as lite
import sys
import datetime
import os
    
def getBranchAndCutTableCreateString(tableName):
    string =         'CREATE TABLE ' + tableName + '''(date timestamp, 
                     inst tinytext,
                     clusterName tinytext,
                     flowVarsType tinytext,
                     datasetName tinytext,
                     datasetIndex tinyint(100),
                     numCommod tinyint(100),
                     commodLocType tinytext,
                     commodDemandType tinytext,
                     numJamLocs tinyint(100),
                     commRange real,
                     infRange real,
                     multRadiosPerNode tinyint(100),
                     numChannels tinyint(100),
                     jamBudget tinyint(100),
                     jamRange real,
                     maxNumHops tinyint(100),
                     algType tinytext,
                     runTime real,
                     lb real,
                     ub real,
                     objValue real)'''
    return string

def getComputationResultsTableCreateString(tableName):
    string =         'CREATE TABLE ' + tableName + '''(date timestamp, 
                     inst tinytext,
                     clusterName tinytext,
                     flowVarsType tinytext,
                     datasetName tinytext,
                     datasetIndex tinyint(100),
                     numCommod tinyint(100),
                     commodLocType tinytext,
                     commodDemandType tinytext,
                     numJamLocs tinyint(100),
                     commRange real,
                     infRange real,
                     multRadiosPerNode tinyint(100),
                     numChannels tinyint(100),
                     jamBudget tinyint(100),
                     jamRange real,
                     maxNumHops tinyint(100),
                     algType tinytext,
                     runTime real,
                     numNodesExplored tinyint(100),
                     numMaxWtIndSetProbsSolved tinyint(100),
                     lb real,
                     ub real,
                     objValue real)'''
    return string
    
def getHeuristicResultsTableCreateString():
    tableName = 'NoInterferenceHeuristic'
    string =         'CREATE TABLE ' + tableName + '''(date timestamp, 
                     inst tinytext,
                     clusterName tinytext,
                     flowVarsType tinytext,
                     datasetName tinytext,
                     datasetIndex tinyint(100),
                     numCommod tinyint(100),
                     commodLocType tinytext,
                     commodDemandType tinytext,
                     numJamLocs tinyint(100),
                     commRange real,
                     infRange real,
                     multRadiosPerNode tinyint(100),
                     numChannels tinyint(100),
                     jamBudget tinyint(100),
                     jamRange real,
                     maxNumHops tinyint(100),
                     algType tinytext,
                     runTime real,
                     heuristicSoln tinytext,
                     bestSoln tinytext,
                     heuristicThroughput real,
                     actualThroughput real,
                     bestThroughput real);'''
    return string
    
def createTable_BranchAnCut(databaseName, tableName):
    con = None
    try:
        con = lite.connect(databaseName)
        
        c = con.cursor()
        
        c.execute('drop table if exists ' + tableName)
        
        # Create table
        c.execute(getBranchAndCutTableCreateString(tableName))

        con.commit()
    
    except lite.Error, e:
        
        if con:
            con.rollback()
            
        print "Error %s:" % e.args[0]
        sys.exit(1)
        
    finally:
        if con:
            con.close()
            
def printResultsToDB(databaseName, tableName, infrastructureInfo, modelInfo, dataSetInfo, instanceInfo, algParams, algOutput, solnOutput):
    con = None
    date = datetime.datetime.now()
    combinedList = [date] + infrastructureInfo + modelInfo + dataSetInfo + instanceInfo + algParams + algOutput + solnOutput
    print combinedList
    valuesStringQuestionMark = ""
    for index in range(len(combinedList) - 1):
        valuesStringQuestionMark += "?" + ", "
    valuesStringQuestionMark += "?"
    #print " VALUES(" + valuesStringQuestionMark + ")"
    print "databaseName", databaseName
    try:
        con = lite.connect(databaseName)
        con.execute("INSERT INTO " + tableName +" VALUES(" + valuesStringQuestionMark + ")", combinedList)
        con.commit()
    
    except lite.Error, e:
        
        if con:
            con.rollback()
            
        print "Error %s:" % e.args
        sys.exit(1)
        
    finally:
        
        if con:
            con.close()

def testDataEntry_BranchAndCut(databaseName, tableName):
    None
    
def createAllTables(databaseName):
    createTable_BranchAnCut(databaseName, "BranchAndCut")
    
def getStringsToCreateAllTables():
    print getBranchAndCutTableCreateString("BranchAndCut")
    
if __name__ == '__main__':
    databaseName = "/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/expr_output/wnj-results_local.db"
    
    print getHeuristicResultsTableCreateString()
    #createAllTables(databaseName)
    #getStringsToCreateAllTables()
    print "finished"
