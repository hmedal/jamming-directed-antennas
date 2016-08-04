import xml.etree.cElementTree as ET
import os
import numpy as np
import itertools
import math
from src.wnj_proj.wnj_jamInf.mod import wnj_interference

#paths and stringsdatabaseName
dataPath = None
ending = None
exprFilePath = None
outputFilePath = None
exeStatement = None
runFilesPath = None
runFilesPathName = None
localExprFilePath = None
databaseName = None
modelType = None
g_inst = None
g_clusterName = None
#hazardsFileEnding = "_conditional"

#names and statements
#cluster = None
moduleStatement = None
#timeString = None
myNumThreads = 0

#constants
g_timeLimitHours = 0

FUZZ = 0.0001

exeMap = {'cormican' : "mod/wnj_interference.py", 'regular' : "mod/wnj_interference.py"}
abbrev = {'cormican': 'cor', 'noInterferenceHeuristic' : 'heur', 'berkeley': 'berk', 'b-and-c-reg' : 'no-cuts', 'branch-and-cut' : 'bandc', 'benders' : 'bend', 'benders-callback' : 'bend-c','benders-tr' : 'bend-tr'}
datasetsMap = {}
jobFilesList = []
computersInfoMap = {}

createForLocal = False

def createComputersInfoMap():
    global computersInfoMap 
    computersInfoMap = {'hmedal' : {}, 'uark': {}, 'msu' : {}}
    computersInfoMap['hmedal']['local-day'] = {}
    computersInfoMap['hmedal']['local-day']['numThreads'] = 4
    computersInfoMap['hmedal']['local-day']['programToExecuteWith'] = 'python'
    computersInfoMap['hmedal']['local-overnight'] = {}
    computersInfoMap['hmedal']['local-overnight']['numThreads'] = 8
    computersInfoMap['hmedal']['local-overnight']['programToExecuteWith'] = 'python'
    #timeStringServers = '2:00:00'
    def createUARKInfo():
        uarkSchedulerCommand = '/share/apps/torque/bin/qsub'
        uarkProgramToExecuteWith = '/home/hmedal/proj/wnj/scripts/gurobi.sh'
        
        computersInfoMap['uark']['star'] = {} # note: star is not currently working with Gurobi because of license issues
        computersInfoMap['uark']['star']['numThreads'] = 8
        computersInfoMap['uark']['star']['programToExecuteWith'] = uarkProgramToExecuteWith
        computersInfoMap['uark']['star']['queueName'] = 'short8core'
        #computersInfoMap['uark']['razor']['timeString'] = timeStringServers
        computersInfoMap['uark']['star']['schedulerCmd'] = uarkSchedulerCommand
        
        computersInfoMap['uark']['razor-12'] = {}
        computersInfoMap['uark']['razor-12']['numThreads'] = 12
        computersInfoMap['uark']['razor-12']['programToExecuteWith'] = uarkProgramToExecuteWith
        computersInfoMap['uark']['razor-12']['queueName'] = 'tiny12core'
        #computersInfoMap['uark']['razor']['timeString'] = timeStringServers
        computersInfoMap['uark']['razor-12']['schedulerCmd'] = uarkSchedulerCommand
        
        computersInfoMap['uark']['razor-16'] = {}
        computersInfoMap['uark']['razor-16']['numThreads'] = 16
        computersInfoMap['uark']['razor-16']['programToExecuteWith'] = uarkProgramToExecuteWith
        computersInfoMap['uark']['razor-16']['queueName'] = 'tiny16core'
        #computersInfoMap['uark']['razor']['timeString'] = timeStringServers
        computersInfoMap['uark']['razor-16']['schedulerCmd'] = uarkSchedulerCommand
        
        computersInfoMap['uark']['debug'] = {}
        computersInfoMap['uark']['debug']['numThreads'] = 12
        computersInfoMap['uark']['debug']['programToExecuteWith'] = uarkProgramToExecuteWith
        computersInfoMap['uark']['debug']['queueName'] = 'debug12core'
        #computersInfoMap['uark']['debug']['timeString'] = timeStringServers
        computersInfoMap['uark']['debug']['schedulerCmd'] = uarkSchedulerCommand
    def createMSUInfo():
        msuSchedulerCommand = 'qsub'
        msuProgramToExecuteWith = 'python'
        computersInfoMap['msu']['raptor'] = {}
        computersInfoMap['msu']['raptor']['numThreads'] = 4
        computersInfoMap['msu']['raptor']['programToExecuteWith'] = msuProgramToExecuteWith
        computersInfoMap['msu']['raptor']['queueName'] = 'q128p48h'
        #computersInfoMap['msu']['raptor']['timeString'] = timeStringServers
        computersInfoMap['msu']['raptor']['schedulerCmd'] = msuSchedulerCommand
        
        computersInfoMap['msu']['talon'] = {}
        computersInfoMap['msu']['talon']['numThreads'] = 12
        computersInfoMap['msu']['talon']['programToExecuteWith'] = msuProgramToExecuteWith
        computersInfoMap['msu']['talon']['queueName'] = 'q192p48h'
        #computersInfoMap['msu']['talon']['timeString'] = timeStringServers
        computersInfoMap['msu']['talon']['schedulerCmd'] = msuSchedulerCommand
        
        computersInfoMap['msu']['shadow'] = {}
        computersInfoMap['msu']['shadow']['numThreads'] = 12
        computersInfoMap['msu']['shadow']['programToExecuteWith'] = msuProgramToExecuteWith
        computersInfoMap['msu']['shadow']['queueName'] = 'q192p48h'
        #computersInfoMap['msu']['talon']['timeString'] = timeStringServers
        computersInfoMap['msu']['shadow']['schedulerCmd'] = msuSchedulerCommand
    createUARKInfo()
    createMSUInfo()
    
def createDatasetsMap():
    global datasetsMap
    #create grid datasets
    def createGridMaps():
        for gridSize in range (3, 11):
            name = 'grid-' + str(gridSize) + 'x' + str(gridSize)
            distToClosest = 1 / (gridSize - 1.0)
            datasetsMap[name] = {}
            datasetsMap[name]['numNodes'] = gridSize**2
            datasetsMap[name]['commodLocType'] = 'corners'
            datasetsMap[name]['commRangeLevels'] = {}
            datasetsMap[name]['commRangeLevels'][1] = round(distToClosest + FUZZ, 4) # can communicate with N-S-E-W neighbors
            datasetsMap[name]['commRangeLevels'][2] = round(distToClosest * math.sqrt(2) + FUZZ, 4) # can also communicate with NE-NW-SE-SW neighbors
            datasetsMap[name]['commRangeLevels'][3] = round(2 * (distToClosest + FUZZ), 4) #can communicate with nodes two links removed
            datasetsMap[name]['jamRangeLevels'] = {}
            datasetsMap[name]['jamRangeLevels'][1] = 0.0 + FUZZ # can only jam node if located right on top of it
            datasetsMap[name]['jamRangeLevels'][2] = round(distToClosest / 2.0 + FUZZ, 4) # can jam two nodes if placed between them
            datasetsMap[name]['jamRangeLevels'][3] = round(distToClosest * 0.5 * math.sqrt(2) + FUZZ, 4) # can jam four nodes if placed in the middle of them
            datasetsMap[name]['jamRangeLevels'][4] = round(distToClosest + FUZZ, 4) # can jam N-S-E-W nodes
            datasetsMap[name]['jamRangeLevels'][5] = round(distToClosest * math.sqrt(2) + FUZZ, 4) # can also jam NE-NW-SE-SW nodes
    def createMITMap():
        name = 'mit'
        datasetsMap[name] = {}
        datasetsMap[name]['numNodes'] = 92
        datasetsMap[name]['baselineCommRange'] = 0.3572
        datasetsMap[name]['commodLocType'] = 'rand'
        datasetsMap[name]['jamRangeLevels'] = {}
        gridSizeOfComparisonDataset = int(round(math.sqrt(datasetsMap[name]['numNodes']), 0))
        nameOfComparisonDataset = 'grid-' + str(gridSizeOfComparisonDataset) + 'x' + str(gridSizeOfComparisonDataset)
        datasetsMap[name]['jamRangeLevels'][1] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][1]
        datasetsMap[name]['jamRangeLevels'][2] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][2]
        datasetsMap[name]['jamRangeLevels'][3] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][3]
        datasetsMap[name]['jamRangeLevels'][4] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][4]
    def createBerkeleyMap():
        name = 'berkeley'
        datasetsMap[name] = {}
        datasetsMap[name]['numNodes'] = 54
        datasetsMap[name]['baselineCommRange'] = 0.1667
        datasetsMap[name]['commodLocType'] = 'rand'
        datasetsMap[name]['jamRangeLevels'] = {}
        gridSizeOfComparisonDataset = int(round(math.sqrt(datasetsMap[name]['numNodes']), 0))
        nameOfComparisonDataset = 'grid-' + str(gridSizeOfComparisonDataset) + 'x' + str(gridSizeOfComparisonDataset)
        datasetsMap[name]['jamRangeLevels'][1] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][1]
        datasetsMap[name]['jamRangeLevels'][2] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][2]
        datasetsMap[name]['jamRangeLevels'][3] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][3]
        datasetsMap[name]['jamRangeLevels'][4] = datasetsMap[nameOfComparisonDataset]['jamRangeLevels'][4]
    
    createGridMaps()
    createMITMap()
    createBerkeleyMap()
    #print "datasetsMap", datasetsMap
    
# def setQueueInfo(debug):
#     global numThreadsHPC
#     global cluster, timeString
#     if(debug):
#         cluster = 'debug12core' # there is no debug for 16 core
#         timeString = '30:00'
#         numThreadsHPC = 12
#     else:
#         cluster = 'tiny16core'
#         timeString = '2:00:00'
        
def setPathsAndNumThreads(inst, clusterName, exeMapKey):
    global dataPath, ending, exprFilePath, outputFilePath, exeStatement, runFilesPath, localExprFilePath, runFilesPathName, modelType
    global moduleStatement, myNumThreads, databaseName
    global createForLocal
    global g_inst, g_clusterName
    g_inst = inst
    g_clusterName = clusterName
    modelType = exeMapKey
    pathEnding = inst + '/' + clusterName
    localExprFilePath = '/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/exprFiles/' + pathEnding
    runFilesPath = '/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/runFiles/' + pathEnding
    runFilesPathName = '/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/runFiles/' + pathEnding
    myNumThreads = computersInfoMap[inst][clusterName]['numThreads']
    if inst == 'hmedal':
        print exeMapKey, exeMap, "exeMap[exeMapKey]", exeMap[exeMapKey]
        createForLocal = True
        moduleStatement = ''
        exeStatement = '/home/hmedal/Documents/2_msu/research_manager/code/ide/eclipse_new3/WNJ3/src/wnj_proj/wnj_jamInf/' + exeMap[exeMapKey]
        dataPath = '/home/hmedal/Documents/2_msu/research_manager/data/'
        ending = '_local'
        exprFilePath = localExprFilePath
        outputFilePath = '/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/outputFiles'
        databaseName = "/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/expr_output/wnj-results_local.db"
    else:
        createForLocal = False
        exeStatement = '~/code/src/python/wnj/src/wnj_proj/wnj_jamInf/' + exeMap[exeMapKey]
        dataPath = '/home/hmedal/data/'
        ending = ''
        exprFilePath = '/home/hmedal/proj/wnj/exprFiles/' + clusterName
        outputFilePath = '/home/hmedal/proj/wnj/outputFiles'
        runFilesPathName = '/home/hmedal/proj/wnj/runFiles/' + clusterName
        databaseName = "/home/hmedal/proj/wnj/outputFiles/wnj.db"
        if(inst == 'uark'):
            moduleStatement = "module load gcc/4.6.3 mkl/13.1.0 python/2.7.5 gurobi/5.5.0"
        elif(inst == 'msu'):
            moduleStatement = ''
        else:
            raise Exception('invalid computer name')

def getDatasetPart(baseline, name, datasetIndex, numCommod, commodLocType, commodDemandType):
    if(baseline):
        return 'base'
    else:
        return name + '_index-' + str(datasetIndex) + "_commods-" + str(numCommod) \
        + "_cLoc-" + str(commodLocType) + "_dem-" + str(commodDemandType)
        
def getInstancePart(baseline, numJamLocs, commRange, infRange, multRadiosPerNode, numChannels, jamBudget, jamRange, maxNumHopsMult):
    if(baseline):
        return 'base'
    else:
        return 'jamLocs-' + str(numJamLocs) + '_commRg-' + str(round(commRange, 2)) + '_infRange-' + str(round(infRange, 2)) + \
            '_multRadio-' + str(int(multRadiosPerNode)) + '_chan-' + str(numChannels) + '_jamB-' + str(jamBudget) + '_jamRng-' + \
            str(jamRange) + '_hop-' + str(maxNumHopsMult)
    
def getInstancePartMinor(baseline):
    if(baseline):
        return ''
    else:
        return ''

def getAlgorithmPart(baseline):
    if(baseline):
        return ''
    else:
        return ''
    
def getModelPart(baseline, jamVarsType, flowVarsType):
    if(baseline):
        return ''
    else:
        return str(jamVarsType) + '_' + str(flowVarsType)

def createExprAndJobFiles(algType, paramsTuple, baselineFlagsMap, modelName, numLevels=3, isLastJob = False):
    global g_timeLimitHours
    #print "createExprAndJobFiles isLastJob", isLastJob
    modelBaseline = baselineFlagsMap['dataset']
    datasetBaseline = baselineFlagsMap['dataset']
    instanceBaseline = baselineFlagsMap['instance']
    instanceBaselineMinor = baselineFlagsMap['instanceMinor']
    algorithmBaseline = baselineFlagsMap['algorithm']
    index = 0
    datasetName = paramsTuple[index]
    index += 1
    datasetIndex = paramsTuple[index]
    index += 1
    numCommod = paramsTuple[index]
    index += 1
    commodDemandType = paramsTuple[index]
    index += 1
    numJamLocs = paramsTuple[index]
    index += 1
    infRangeMult = paramsTuple[index]
    index += 1
    commRange = paramsTuple[index]
    index += 1
    multRadiosPerNode = paramsTuple[index]
    index += 1
    numChannels = paramsTuple[index]
    index += 1
    jamBudget = paramsTuple[index]
    index += 1
    jamRange = paramsTuple[index]
    index += 1
    maxNumHopsMult = paramsTuple[index]
    index += 1
    jamVarsType = paramsTuple[index]
    index += 1
    flowVarsType = paramsTuple[index]
    index += 1
    timeLimitHours = paramsTuple[index]
    index += 1
    # derived values
    numNodes = datasetsMap[datasetName]['numNodes']
    #commRange = commRangeMult * (datasetsMap[datasetName]['baselineCommRange'])
    infRange = infRangeMult * commRange
    g_timeLimitHours = timeLimitHours
    
    root = ET.Element("experimentData")
    dataset = ET.SubElement(root, "dataset")
    commodLocType = datasetsMap[datasetName]['commodLocType']
    pathName = dataPath + 'wirelessNet/' + datasetName + '_index-' + str(datasetIndex) \
        + "_commods-" + str(numCommod) \
        + "_commodLoc-" + str(commodLocType) + "_commodDemand-" + str(commodDemandType) + ".xml"
        
    ET.SubElement(dataset, "path").text = pathName
    ET.SubElement(dataset, "name").text = datasetName
    ET.SubElement(dataset, "numNodes").text = str(numNodes)
    ET.SubElement(dataset, "datasetIndex").text = str(datasetIndex)
    ET.SubElement(dataset, "numCommod").text = str(numCommod)
    ET.SubElement(dataset, "commodLocType").text = str(commodLocType)
    ET.SubElement(dataset, "commodDemandType").text = str(commodDemandType)
    
    instance = ET.SubElement(root, "instance")
    ET.SubElement(instance, "numJamLocs").text = str(numJamLocs)
    #ET.SubElement(instance, "commRange").text = str(commRange)
    ET.SubElement(instance, "commRange").text = str(commRange)
    #ET.SubElement(instance, "infRange").text = str(infRange)
    ET.SubElement(instance, "infRange").text = str(infRange)
    ET.SubElement(instance, "multRadiosPerNode").text = str(int(multRadiosPerNode))
    ET.SubElement(instance, "numChannels").text = str(numChannels)
    ET.SubElement(instance, "jamBudget").text = str(jamBudget)
    ET.SubElement(instance, "jamRange").text = str(jamRange)
    ET.SubElement(instance, "maxNumHopsMult").text = str(maxNumHopsMult)
    
    algorithm = ET.SubElement(root, "algorithm")
    ET.SubElement(algorithm, "timeLimit").text = str(timeLimitHours * 3600)
    ET.SubElement(algorithm, "numThreads").text = str(myNumThreads)
    
    otherSubElement = ET.SubElement(root, "other")
    ET.SubElement(otherSubElement, "databaseName").text = databaseName
    ET.SubElement(otherSubElement, "inst").text = g_inst
    ET.SubElement(otherSubElement, "clusterName").text = g_clusterName
    
    modelSubElement = ET.SubElement(root, "model")
    ET.SubElement(modelSubElement, "modelType").text = modelType
    ET.SubElement(modelSubElement, "jamVarsType").text = jamVarsType
    ET.SubElement(modelSubElement, "flowVarsType").text = flowVarsType
    
    tree = ET.ElementTree(root)
    sep = "_"
    headString = "WnjInf_" +abbrev[modelName] + sep + abbrev[algType] + sep
    headStringExpr = "WnjInf_" +abbrev[modelName] + sep
    middleString = getModelPart(modelBaseline, jamVarsType, flowVarsType)  + sep + \
        getDatasetPart(datasetBaseline, datasetName, datasetIndex, numCommod, commodLocType, commodDemandType)  + sep + \
        getInstancePart(instanceBaseline, numJamLocs, commRange, infRange, multRadiosPerNode, numChannels, jamBudget, jamRange, maxNumHopsMult) + sep + \
        getInstancePartMinor(instanceBaselineMinor) + sep + getAlgorithmPart(algorithmBaseline)
    exprFileName = exprFilePath + '/' + headStringExpr + middleString + ending + '.xml'
    outputFileName = outputFilePath + '/' + headString  + middleString + ending + '.log'
    #print "exprFileName", exprFileName
    #print "outputFileName", outputFileName
    exprFile = localExprFilePath + '/' + headStringExpr  + middleString + ending + '.xml'
    tree.write(exprFile)
    createJobFile(algType, exprFileName, outputFileName, headString + middleString + "_Job", runFilesPath + '/' + headString + middleString + '_Job.pbs', isLastJob)
    jobFilesList.append(runFilesPathName + '/' + headString + middleString + '_Job.pbs')

def createJobFile(algType, exprFile, outputFileName, exprName, pbsOutputFile, isLastJob = False):
    if g_inst == 'hmedal':
        createJobFile_Local(algType, exprFile, outputFileName, pbsOutputFile)
        os.system("chmod a+x " + pbsOutputFile)
    elif g_inst == 'uark' or g_inst == 'msu':
        #print "IsLastJob", isLastJob
        createJobFile_PBSServer(algType, exprFile, outputFileName, exprName,  pbsOutputFile, isLastJob)
    else:
        raise Exception('invalid computer name')

def createJobFile_Local(algType, exprFile, outputFileName, pbsOutputFile):
    global jobFilesList
    f = open(pbsOutputFile, 'w')
    myStr = '#!/bin/bash\n' 
    myStr += 'export PYTHONPATH="$PYTHONPATH:/home/hmedal/Documents/2_msu/research_manager/code/ide/eclipse_new3/WNJ3/src/wnj_proj/wnj_jamInf/mod:/home/hmedal/Documents/2_msu/research_manager/code/ide/eclipse_new3/WNJ3:/home/hmedal/Documents/2_msu/research_manager/code/ide/eclipse_new3/WNJ3/src"\n'
    myStr += computersInfoMap[g_inst][g_clusterName]['programToExecuteWith'] + ' ' + exeStatement + " --algType " + algType + " --exprfile " + exprFile + " > " + outputFileName
    f.write(myStr)
    
def createJobFile_PBSServer(algType, exprFile, outputFileName, exprName, pbsOutputFile, isLastJob = False):
    #print "exprFile again", exprFile
    #print "outputFileName again", outputFileName
    global jobFilesList
    f = open(pbsOutputFile, 'w')
    myStr = "#PBS -N " + exprName + "\n"
    myStr += "#PBS -q " + computersInfoMap[g_inst][g_clusterName]['queueName'] + "\n"
    myStr += "\n"
    myStr += "#PBS -j oe\n"
    myStr += "\n"
    myStr += "#PBS -M hugh.medal@msstate.edu\n" # send me email when job aborts (with an error)
    if isLastJob:
        myStr += "#PBS -m ae\n"
    else:
        myStr += "#PBS -m a\n"
    myStr += "#PBS -o " + exprName +".$PBS_JOBID\n"
    myStr += "#PBS -l nodes=1:ppn=" + str(myNumThreads) + "\n"
    myStr += "#PBS -l walltime=" + str(g_timeLimitHours + 1) +":00:00" + "\n"
    myStr += "\n"
    myStr += "cd $PBS_O_WORKDIR\n"
    myStr += moduleStatement + "\n"
    myStr += computersInfoMap[g_inst][g_clusterName]['programToExecuteWith'] + ' ' + exeStatement + " --algType " + algType + " --exprfile " + exprFile + " > " + outputFileName
    if isLastJob:
        print myStr
    f.write(myStr)

def createRunScript(scriptFile, exprFiles):
    if(createForLocal):
        createRunScript_Local(scriptFile, exprFiles)
        os.system("chmod a+x " + scriptFile)
    else:
        createRunScript_Server(scriptFile, exprFiles)

def createRunScript_Local(scriptFile, exprFiles):
    #print "createRunScript"
    f = open(scriptFile, 'w')
    myStr = "#!/bin/bash\n"
    for exprFile in exprFiles:
        print exprFile
        myStr += exprFile + "\n"
    f.write(myStr)
    #print "scriptFile", scriptFile
    
def createRunScript_Server(scriptFile, exprFiles):
    #print "createRunScript"
    f = open(scriptFile, 'w')
    myStr = "#!/bin/sh\n"
    myStr += ". ~/.bashrc"
    myStr += "\n"
    for exprFile in exprFiles:
        print computersInfoMap[g_inst][g_clusterName]['schedulerCmd'] + " " + exprFile
        myStr += computersInfoMap[g_inst][g_clusterName]['schedulerCmd'] + " " + exprFile + "\n"
    f.write(myStr)
    #print "scriptFile", scriptFile

def createFilesForJob(algType, tuple, flagsTuple, modelName, isTheLastJob = False):
    global jobFilesList
    #print "createFilesForJob isTheLastJob", isTheLastJob
    createExprAndJobFiles(algType, tuple, flagsTuple, modelName, isLastJob = isTheLastJob)
    
def createFilesForJob_Batch(algType, arrays, flagsTuple, modelName, scenario):
    cartProd = itertools.product(*arrays)
    lengthCtr = 1
    for array in arrays:
        lengthCtr *= len(array)
    lengthOfIterator = lengthCtr
    counter = 0
    for array in cartProd:
        counter += 1
        createFilesForJob(algType, array, flagsTuple, modelName, counter == lengthOfIterator)
    createRunScript(runFilesPath + '/scripts/' + modelName + '_' + algType + '_' + scenario + '_script.sh', jobFilesList)
    
def getModelParamsArray(param, variation):
    baselineDataset = 'grid-7x7'
    if(variation == 'example'):
        myDataset = 'grid-5x5'
        if(param == 'datasetName'):
            return [myDataset]
        elif(param == 'datasetIndex'):
            return [1]
        elif(param == 'numCommod'):
            return [8,16]
        elif(param == 'commodDemandType'):
            return ['rand02']
        elif(param == 'numJamLocs'):
            return [9, 16, 25, 36, 49, 64, 81]
        elif(param == 'commRange'):
            return [datasetsMap[myDataset]['commRangeLevels'][1]]
        #elif(param == 'commRangeMult'):
        #    return [1, 1.25, 1.5]
        elif(param == 'infRangeMult'):
            return [1]
        elif(param == 'multRadiosPerNode'):
            return [False]
        elif(param == 'numChannels'):
            return [1]
        elif(param == 'jamBudget'):
            return [1, 2]
        elif(param == 'jamRange'):
            print "jamRanges", datasetsMap[myDataset]['jamRangeLevels'][1], datasetsMap[myDataset]['jamRangeLevels'][2]
            return [datasetsMap[myDataset]['jamRangeLevels'][1], datasetsMap[myDataset]['jamRangeLevels'][2]]
        elif(param == 'maxNumHopsMult'):
            return [2.0]
        elif(param == 'jamVarsType'):
            return ['jam-vars']
        elif(param == 'flowVarsType'):
            return ['arc-based']
    if(variation == 'testing'):
        if(param == 'datasetName'):
            return ['grid-5x5']
        elif(param == 'datasetIndex'):
            return [1]
        elif(param == 'numCommod'):
            return [16]
        elif(param == 'commodDemandType'):
            return ['unif5']
        elif(param == 'numJamLocs'):
            return [49]
        elif(param == 'commRange'):
            return [datasetsMap[baselineDataset]['commRangeLevels'][1]]
        #elif(param == 'commRangeMult'):
        #    return [1, 1.25, 1.5]
        elif(param == 'infRangeMult'):
            return [1.75]
        elif(param == 'multRadiosPerNode'):
            return [False]
        elif(param == 'numChannels'):
            return [1]
        elif(param == 'jamBudget'):
            return [1]
        elif(param == 'jamRange'):
            return [datasetsMap['grid-5x5']['jamRangeLevels'][1]]
        elif(param == 'maxNumHopsMult'):
            return [2.0]
        elif(param == 'jamVarsType'):
            return ['jam-vars']
        elif(param == 'flowVarsType'):
            return ['path-based', 'arc-based']
    elif(variation == 'all'):
        if(param == 'datasetName'):
            return ['grid-7x7', 'berkeley']
    elif(variation == 'makeup'):
        if(param == 'datasetName'):
            return ['berkeley']
        elif(param == 'timeLimitHours'):
            return [12]
    elif(variation == 'makeup_density_of_nodes_and_ranges'):
        if(param == 'datasetName'):
            return ['grid-8x8', 'grid-9x9', 'berkeley']
    elif(variation == 'base'):
        if(param == 'datasetName'):
            return [baselineDataset]
        elif(param == 'datasetIndex'):
            return [1]
        elif(param == 'numCommod'):
            return [16]
        elif(param == 'commodDemandType'):
            return ['rand02']
        elif(param == 'numJamLocs'):
            return [49]
        elif(param == 'commRange'):
            return [datasetsMap[baselineDataset]['commRangeLevels'][1]]
        #elif(param == 'commRangeMult'):
        #    return [1.5]
        elif(param == 'infRangeMult'):
            return [1.75]
        elif(param == 'multRadiosPerNode'):
            return [False]
        elif(param == 'numChannels'):
            return [1]
        elif(param == 'jamBudget'):
            return [2]
        elif(param == 'jamRange'):
            return [datasetsMap[baselineDataset]['jamRangeLevels'][2]]
        elif(param == 'maxNumHopsMult'):
            return [2.0]
        elif(param == 'jamVarsType'):
            return ['jam-vars']
        elif(param == 'flowVarsType'):
            return ['arc-based']
        elif(param == 'timeLimitHours'):
            return [5]
    elif(variation == 'latency'):
        if(param == 'datasetName'):
            return ['grid-7x7', 'berkeley']
        elif(param == 'jamBudget'):
            return [1, 2, 3]
        elif(param == 'maxNumHopsMult'):
            return [1.0, 1.25, 2]
        elif(param == 'flowVarsType'):
            return ['path-based']
        elif(param == 'timeLimitHours'):
            return [5]
    elif(variation == 'runtime'):# 9/11/14: removed multiple radios per node (can add back if algorithm speeds up)
        if(param == 'datasetName'):
            return ['grid-7x7', 'grid-8x8', 'grid-9x9', 'berkeley']
        elif(param == 'numChannels'):
            return [1, 2]
        elif(param == 'jamBudget'):
            return [1, 3]
        elif(param == 'flowVarsType'):
            return ['arc-based']
        elif(param == 'timeLimitHours'):
            return [2]
    elif(variation == 'density_of_nodes_and_ranges'):
        if(param == 'datasetName'):
            return ['grid-7x7', 'grid-8x8', 'grid-9x9', 'berkeley']
        elif(param == 'commRange'):
            return [datasetsMap[baselineDataset]['commRangeLevels'][1], datasetsMap[baselineDataset]['commRangeLevels'][2]]
        elif(param == 'infRangeMult'):
            return [0, 1, 1.75]
        elif(param == 'timeLimitHours'):
            return [5]
    elif(variation == 'num_of_channels'):# no longer have multiple radios per node
        if(param == 'datasetName'):
            return ['grid-7x7', 'berkeley']
        elif(param == 'numChannels'):
            return [1, 2, 3]
        elif(param == 'timeLimitHours'):
            return [10]
    elif(variation == 'noInterferenceHeuristic'):# no longer have multiple radios per node
        if(param == 'numCommod'):
            return [1,16]
        elif(param == 'jamBudget'):
            return [1, 3]
        elif(param == 'numChannels'):
            return [1, 2]
    #elif(variation == 'radio_range'):
    #    if(param == 'commRange'):
    #        return [datasetsMap[baselineDataset]['commRangeLevels'][1], 
    #                datasetsMap[baselineDataset]['commRangeLevels'][2], datasetsMap[baselineDataset]['commRangeLevels'][3]]
    #elif(variation == 'inf_range'):
    #    if(param == 'infRangeMult'):
    #        return [0, 1, 1.5]
    elif(variation == 'jam_locs'):
        if(param == 'datasetName'):
            return ['grid-7x7', 'berkeley']
        elif(param == 'numJamLocs'):
            return [9, 16, 25, 36, 64]
        elif(param == 'timeLimitHours'):
            return [5]
    elif(variation == 'heuristic-testing'):
        if(param == 'infRangeMult'):
            return [0]
        elif(param == 'datasetName'):
            return ['berkeley']
    elif(variation == 'num_jam_devices_and_jam_range'):
        if(param == 'datasetName'):
            return ['grid-7x7', 'berkeley']
        elif(param == 'jamBudget'):
            return [1, 2, 3, 4, 5]
        elif(param == 'jamRange'):
            def createList():
                myList = []
                for index in [1, 2, 3, 4, 5]:
                    myList.append(datasetsMap[baselineDataset]['jamRangeLevels'][index])
                return myList
            return createList()
        elif(param == 'timeLimitHours'):
            return [5]
    
def createFiles(modelName, algType, scenario):
    global jobFilesList
    datasetNameArray = getModelParamsArray('datasetName', 'base')
    datasetIndexArray = getModelParamsArray('datasetIndex', 'base')
    numCommodArray = getModelParamsArray('numCommod', 'base')
    commodDemandTypeArray = getModelParamsArray('commodDemandType', 'base')
    #instance
    numJamLocsArray = getModelParamsArray('numJamLocs', 'base')
    infRangeMultsArray = getModelParamsArray('infRangeMult', 'base')
    commRangeArray = getModelParamsArray('commRange', 'base')
    multRadiosPerNodeArray = getModelParamsArray('multRadiosPerNode', 'base')
    numChannelsArray = getModelParamsArray('numChannels', 'base')
    jamBudgetArray = getModelParamsArray('jamBudget', 'base')
    jamRangeArray = getModelParamsArray('jamRange', 'base')
    maxNumHopsMultArray = getModelParamsArray('maxNumHopsMult', 'base')
    jamVarsTypeArray = getModelParamsArray('jamVarsType', 'base')
    flowVarsTypeArray = getModelParamsArray('flowVarsType', 'base')
    timeLimitArray = getModelParamsArray('timeLimitHours', 'base')
    if(scenario == 'example'):
        #dataset
        datasetNameArray = getModelParamsArray('datasetName', 'example')
        datasetIndexArray = getModelParamsArray('datasetIndex', 'example')
        numCommodArray = getModelParamsArray('numCommod', 'example')
        commodDemandTypeArray = getModelParamsArray('commodDemandType', 'example')
        #instance
        numJamLocsArray = getModelParamsArray('numJamLocs', 'example')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'example')
        commRangeArray = getModelParamsArray('commRange', 'example')
        multRadiosPerNodeArray = getModelParamsArray('multRadiosPerNode', 'example')
        numChannelsArray = getModelParamsArray('numChannels', 'example')
        jamBudgetArray = getModelParamsArray('jamBudget', 'example')
        jamRangeArray = getModelParamsArray('jamRange', 'example')
        maxNumHopsMultArray = getModelParamsArray('maxNumHopsMult', 'example')
        jamVarsTypeArray = getModelParamsArray('jamVarsType', 'example')
        flowVarsTypeArray = getModelParamsArray('flowVarsType', 'example')
    elif(scenario == 'testing'):
        #dataset
        datasetNameArray = getModelParamsArray('datasetName', 'testing')
        datasetIndexArray = getModelParamsArray('datasetIndex', 'testing')
        numCommodArray = getModelParamsArray('numCommod', 'testing')
        commodDemandTypeArray = getModelParamsArray('commodDemandType', 'testing')
        #instance
        numJamLocsArray = getModelParamsArray('numJamLocs', 'testing')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'testing')
        commRangeArray = getModelParamsArray('commRange', 'testing')
        multRadiosPerNodeArray = getModelParamsArray('multRadiosPerNode', 'testing')
        numChannelsArray = getModelParamsArray('numChannels', 'testing')
        jamBudgetArray = getModelParamsArray('jamBudget', 'testing')
        jamRangeArray = getModelParamsArray('jamRange', 'testing')
        maxNumHopsMultArray = getModelParamsArray('maxNumHopsMult', 'testing')
        jamVarsTypeArray = getModelParamsArray('jamVarsType', 'testing')
        flowVarsTypeArray = getModelParamsArray('flowVarsType', 'testing')
    elif(scenario == 'base'):
        None
    elif(scenario == 'latency'):
        datasetNameArray = getModelParamsArray('datasetName', 'all')
        jamBudgetArray = getModelParamsArray('jamBudget', 'latency')
        flowVarsTypeArray = getModelParamsArray('flowVarsType', 'latency')
        maxNumHopsMultArray = getModelParamsArray('maxNumHopsMult', 'latency')
    elif(scenario == 'runtime'):
        datasetNameArray = getModelParamsArray('datasetName', 'runtime')
        numChannelsArray = getModelParamsArray('numChannels', 'runtime')
        jamBudgetArray = getModelParamsArray('jamBudget', 'runtime')
        flowVarsTypeArray = getModelParamsArray('flowVarsType', 'runtime')
        timeLimitArray = getModelParamsArray('timeLimitHours', 'runtime')
    elif(scenario == 'density_of_nodes_and_ranges'):
        datasetNameArray = getModelParamsArray('datasetName', 'density_of_nodes_and_ranges')
        commRangeArray = getModelParamsArray('commRange', 'density_of_nodes_and_ranges')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'density_of_nodes_and_ranges')
    elif(scenario == 'density_of_nodes_and_ranges_makeup'):
        datasetNameArray = getModelParamsArray('datasetName', 'makeup_density_of_nodes_and_ranges')
        commRangeArray = getModelParamsArray('commRange', 'density_of_nodes_and_ranges')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'density_of_nodes_and_ranges')
        timeLimitArray = getModelParamsArray('timeLimitHours', 'makeup')
    elif(scenario == 'num_of_channels'):
        datasetNameArray = getModelParamsArray('datasetName', 'all')
        numChannelsArray = getModelParamsArray('numChannels', 'num_of_channels')
    elif(scenario == 'jam_locs'):
        datasetNameArray = getModelParamsArray('datasetName', 'all')
        numJamLocsArray = getModelParamsArray('numJamLocs', 'jam_locs')
    elif(scenario == 'noInterferenceHeuristic'):
        datasetNameArray = getModelParamsArray('datasetName', 'num_jam_devices_and_jam_range')
        numChannelsArray = getModelParamsArray('numChannels', 'noInterferenceHeuristic')
        numCommodArray = getModelParamsArray('numCommod', 'noInterferenceHeuristic')
        jamBudgetArray = getModelParamsArray('jamBudget', 'noInterferenceHeuristic')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'density_of_nodes_and_ranges')
    elif(scenario == 'noInterferenceHeuristic-testing'):
        datasetNameArray = getModelParamsArray('datasetName', 'heuristic-testing')
        numChannelsArray = getModelParamsArray('numChannels', 'noInterferenceHeuristic')
        numCommodArray = getModelParamsArray('numCommod', 'noInterferenceHeuristic')
        jamBudgetArray = getModelParamsArray('jamBudget', 'noInterferenceHeuristic')
        infRangeMultsArray = getModelParamsArray('infRangeMult', 'heuristic-testing')
    elif(scenario == 'num_jam_devices_and_jam_range'):
        datasetNameArray = getModelParamsArray('datasetName', 'num_jam_devices_and_jam_range')
        jamBudgetArray = getModelParamsArray('jamBudget', 'num_jam_devices_and_jam_range')
        jamRangeArray = getModelParamsArray('jamRange', 'num_jam_devices_and_jam_range')
    else:
        raise Exception("scenario not found")
    arrays = [datasetNameArray, datasetIndexArray, numCommodArray, commodDemandTypeArray, 
              numJamLocsArray, infRangeMultsArray, commRangeArray, multRadiosPerNodeArray, numChannelsArray, jamBudgetArray, \
              jamRangeArray, maxNumHopsMultArray, 
              jamVarsTypeArray, flowVarsTypeArray, timeLimitArray]
    print "arrays", arrays
    baselineFlagsMap = {'model' : False, 'dataset' : False, 'instance' : False, 'instanceMinor' : True, 'algorithm' : True}
    createFilesForJob_Batch(algType, arrays, baselineFlagsMap, modelName, scenario)
    jobFilesList = []

def clearOld(inst, clusterName, flag):
    if flag:
        projectDir = "/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/"
        runDir = projectDir + 'runFiles/' + inst + '/' + clusterName + "/"
        exprDir = projectDir + 'exprFiles/' + inst + '/' + clusterName + "/"
        bashCommand1 = "for i in " + runDir + "*.pbs; do rm -f $i; done"
        bashCommand2 = "rm " + runDir + "scripts/*.sh"
        bashCommand3 = "for i in " + exprDir + "*.xml; do rm -f $i; done"
        os.system(bashCommand1 + '; '+ bashCommand2 + ';' + bashCommand3)
    
def archiveLocalOutput(flag, clusterName):
    if flag:
        projectDir = "/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/"
        outputDir = projectDir + 'outputFiles/hmedal/' + clusterName + '/'
        bashCommand = "for i in " + outputDir + '*.log; do mv -f $i ' + outputDir + 'archive; done'
        print "bash command: ", bashCommand
        os.system(bashCommand)

if __name__ == "__main__":
    createDatasetsMap()
    createComputersInfoMap()  
    
    clearOldFlag = True
    archiveLocalOutputFlag = True
    debug = False
    
    print "starting"
    #myModelAndAlgNames = [[['cormican-path-based-jam-vars'], ['branch-and-cut']], 
    #                      [['cormican-path-based-jam-and-interdict-vars'], ['branch-and-cut']]]
    cormicanBendersClassic = [['cormican', 'benders'], ['cormican', 'benders-ki-tr'], ['cormican', 'benders-ki-pareto'], ['cormican', 'benders-pareto-tr'], ['cormican', 'benders-ki-pareto-tr']]
    runTime2Set = [[['cormican'], ['benders']], [['cormican'], ['benders-tr']], [['cormican'],['branch-and-cut']], [['cormican'],['b-and-c-reg']],[['cormican'],['benders-callback']]]
    runTime3Set = [[['cormican'],['b-and-c-reg']]]
    runTime4Set = [[['cormican'],['branch-and-cut']]]
    
    experimentsMap = {1 : {}, 2 : {}, 3 : {}, 4 : {}, 5 : {}, 6 : {}, 7 : {}, 8 : {}, 9 : {}, 10: {}}
    index = 1
    experimentsMap[index]['myModelAndAlgNames'] = cormicanBendersClassic
    experimentsMap[index]['scenarios'] = ['example']
    experimentsMap[index]['inst'] = 'hmedal'
    experimentsMap[index]['clusterName'] = 'local-day'
    experimentsMap[index]['runScriptPath'] = ''
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['benders']]]
    experimentsMap[index]['scenarios'] = ['jam_locs']
    experimentsMap[index]['inst'] = 'hmedal'
    experimentsMap[index]['clusterName'] = 'local-overnight'
    experimentsMap[index]['runScriptPath'] = ''
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['branch-and-cut']]]
    experimentsMap[index]['scenarios'] = ['latency']
    experimentsMap[index]['inst'] = 'uark'
    experimentsMap[index]['clusterName'] = 'star'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m sr -i uark -c star -s razor.uark.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['branch-and-cut']]]
    experimentsMap[index]['scenarios'] = ['num_of_channels']
    experimentsMap[index]['inst'] = 'uark'
    experimentsMap[index]['clusterName'] = 'razor-12'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m sr -i uark -c razor-12 -s razor.uark.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['benders']], [['cormican'], ['benders-tr']], [['cormican'], ['benders-callback']]]
    experimentsMap[index]['scenarios'] = ['runtime']
    experimentsMap[index]['inst'] = 'uark'
    experimentsMap[index]['clusterName'] = 'razor-16'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m sr -i uark -c razor-16 -s razor.uark.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['branch-and-cut']]]
    experimentsMap[index]['scenarios'] = ['num_jam_devices_and_jam_range']
    experimentsMap[index]['inst'] = 'msu'
    experimentsMap[index]['clusterName'] = 'raptor'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m sr -i msu -c raptor -s titan.hpc.msstate.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['benders']], [['cormican'], ['benders-tr']], [['cormican'], ['benders-callback']]]
    experimentsMap[index]['scenarios'] = ['jam_locs']
    experimentsMap[index]['inst'] = 'msu'
    experimentsMap[index]['clusterName'] = 'talon'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m s -i msu -c talon -s titan.hpc.msstate.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['noInterferenceHeuristic']]]
    experimentsMap[index]['scenarios'] = ['noInterferenceHeuristic']
    experimentsMap[index]['inst'] = 'uark'
    experimentsMap[index]['clusterName'] = 'razor-16'
    experimentsMap[index]['runScriptPath'] = '~/scripts/crWnj.sh -m sr -i uark -c razor-16 -s razor.uark.edu'
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = [[['cormican'], ['noInterferenceHeuristic']]]
    experimentsMap[index]['scenarios'] = ['noInterferenceHeuristic-testing']
    experimentsMap[index]['inst'] = 'hmedal'
    experimentsMap[index]['clusterName'] = 'local-day'
    experimentsMap[index]['runScriptPath'] = ''
    
    index += 1
    experimentsMap[index]['myModelAndAlgNames'] = runTime3Set
    experimentsMap[index]['scenarios'] = ['runtime']
    experimentsMap[index]['inst'] = 'hmedal'
    experimentsMap[index]['clusterName'] = 'local-day'
    experimentsMap[index]['runScriptPath'] = ''
    
    exprLogPath = '/home/hmedal/Documents/2_msu/1_MSU_Projects/Papers/PAPER_JammingSpatialInterference/exprLog.log'
    experimentNumsToRun = [5]
    
    if wnj_interference.debug is True or wnj_interference.writeToFile is True:
        raise Exception("cannot submit jobs when debug or writeToFile are on")
    for experimentNum in experimentNumsToRun:
        if experimentsMap[experimentNum]['inst'] == 'hmedal':
            archiveLocalOutput(archiveLocalOutputFlag, experimentsMap[experimentNum]['clusterName'])
        clearOld(experimentsMap[experimentNum]['inst'], experimentsMap[experimentNum]['clusterName'], clearOldFlag)
        cartProd = itertools.product(experimentsMap[experimentNum]['myModelAndAlgNames'], experimentsMap[experimentNum]['scenarios'])
        # create run files
        for myModelAndAlgName, scenario in cartProd:
            os.system('echo "$(date) ' + experimentsMap[experimentNum]['inst'] + ' ' + experimentsMap[experimentNum]['clusterName'] + ' ' + scenario + '" >> ' + exprLogPath)
            print "myModelAndAlgName", myModelAndAlgName
            myModelName = myModelAndAlgName[0][0]
            myAlgType = myModelAndAlgName[1][0]
            print "myModelName", myModelName
            print "myAlgType", myAlgType
            print myAlgType, myModelName, scenario
            setPathsAndNumThreads(inst = experimentsMap[experimentNum]['inst'], clusterName = experimentsMap[experimentNum]['clusterName'], 
                                  exeMapKey = myModelName)
            createFiles(myModelName, myAlgType, scenario)
            print "files created"
        os.system(experimentsMap[experimentNum]['runScriptPath'])