'''
Created on Jun 14, 2014

@author: hmedal
'''

import numpy as np
from scipy import stats
import gurobipy
import logging
import time
from multiprocessing import Pool

writeToFile = False

def getAllocLevelsVectorFromBinaryVector(binaryVector, instance):
    return [sum([int(a*b) for a,b in zip(list,range(instance.numAllocLevels))]) for list in binaryVector]
    
def getProbabilityFromAllocationComplex(allocation, fac):
    return 1 - np.exp(-allocation)

def getProbabilityFromAllocation(allocation, numAllocLevels):
    return min(1.0, allocation/(numAllocLevels + 0.0))

def OLD_RESIDUAL_PROB_IF_UNPROTECTED_getProbabilityFromAllocation(allocation, numAllocLevels):
    return min(1.0,(allocation + 1.0)/numAllocLevels)

def getProbabilityFromAllocationAndHazardLevel(allocation, numAllocLevels, hazardLevel, numHazardLevels):
    return getProbabilityFromAllocation(allocation, numAllocLevels)**(hazardLevel/(numHazardLevels + 0.0))

def getRandomIndependentBinomials(allocationVector, instance, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    pArray = [myGetProbabilityFromAllocation(allocationVector[fac], instance.numAllocLevels) for fac in range(instance.numFacs)]
    #print "pArray", pArray
    #print "instance.numCapLevels", instance.numCapLevels
    return [np.random.binomial(instance.numCapLevels - 1 ,p) for p in pArray]

def getRandomIndependentBinomials_withHazards(allocationVector, instance, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    hazardLevelVector = instance.getRandomHazardLevelsVector()
    #print "hazardLevelVector", hazardLevelVector
    pArray = [getProbabilityFromAllocationAndHazardLevel(allocationVector[fac], instance.numAllocLevels, hazardLevelVector[fac], instance.numHazardLevels) 
              for fac in range(instance.numFacs)]
    #print "pArray", pArray
    #print "instance.numCapLevels", instance.numCapLevels
    #return [np.random.binomial(instance.numCapLevels - 1 ,p) for p in pArray]
    return [[hazardLevelVector[fac], np.random.binomial(instance.numCapLevels - 1 , pArray[fac])] for fac in range(instance.numFacs)]


def getCapLevelProbabilityForAllocationLevel(fac, capLevel, allocLevel, instance, myNumCapLevels, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    p = myGetProbabilityFromAllocation(allocLevel, instance.numAllocLevels)
    
    prob = 0.0
    if(p == 0.0):
        if(capLevel == 0):
            prob = 1
        else:
            prob = 0
    elif(p == 1.0):
        if(capLevel == myNumCapLevels -1):
            prob = 1
        else:
            prob = 0
    else:
        prob = stats.binom.pmf(capLevel, myNumCapLevels - 1, p)
    return prob

def getCapLevelProbabilityForAllocationLevelAndHazardLevel(fac, capLevel, allocLevel, hazardLevel, instance, 
                                                           myNumCapLevels, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    p = getProbabilityFromAllocationAndHazardLevel(allocLevel, instance.numAllocLevels, hazardLevel, instance.numHazardLevels)
    #print "fac:", fac, "p:", p
    prob = 0.0
    if(p == 0.0):
        if(capLevel == 0):
            prob = 1
        else:
            prob = 0
    elif(p == 1.0):
        if(capLevel == myNumCapLevels -1):
            prob = 1
        else:
            prob = 0
    else:
        prob = stats.binom.pmf(capLevel, myNumCapLevels - 1, p)
    return prob

def createFacilityProbabilitiesForAllocationLevels(instance):
    facilityProbabilitiesForAllocationLevels = [[[getCapLevelProbabilityForAllocationLevel(j,l,k, instance, instance.numCapLevels) for k in range(instance.numAllocLevels)] for l in range(instance.numCapLevels)] for j in range(instance.numFacs)]
    return facilityProbabilitiesForAllocationLevels

def calculateFacilityStateProbabilities(scen, instance, scenariosSet, facilityProbabilitiesForAllocationLevels):
        return [[sum(facilityProbabilitiesForAllocationLevels[j][l][k] * scenariosSet[scen][j][l] for l in range(instance.numCapLevels)) 
                                               for k in range(instance.numAllocLevels)] for j in range(instance.numFacs)]
        
def getCapLevelProbabilityForAllocationLevel_AltNumAllocLevels(fac, capLevel, allocLevel, myNumAllocLevels, myNumCapLevels, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    p = myGetProbabilityFromAllocation(allocLevel, myNumAllocLevels)
    prob = 0.0
    if(p == 0.0):
        if(capLevel == 0):
            prob = 1
        else:
            prob = 0
    elif(p == 1.0):
        if(capLevel == myNumCapLevels -1):
            prob = 1
        else:
            prob = 0
    else:
        prob = stats.binom.pmf(capLevel, myNumCapLevels - 1, p)
    return prob

def getCapLevelProbabilityForAllocationLevelAndHazardLevel_AltNumAllocLevels(fac, capLevel, allocLevel, hazardLevel, myNumAllocLevels, myNumCapLevels, myNumHazardLevels,
                                                                             myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    p = getProbabilityFromAllocationAndHazardLevel(allocLevel, myNumAllocLevels, hazardLevel, myNumHazardLevels)
    prob = 0.0
    if(p == 0.0):
        if(capLevel == 0):
            prob = 1
        else:
            prob = 0
    elif(p == 1.0):
        if(capLevel == myNumCapLevels -1):
            prob = 1
        else:
            prob = 0
    else:
        prob = stats.binom.pmf(capLevel, myNumCapLevels - 1, p)
    return prob

def getProbabilityOfCapLevelScenConsideringHazards(capLevelsVector, allocLevelsVector, instance, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForHazardScenarios = [getProbabilityOfCapLevelScenario_withHazards(capLevelsVector,allocLevelsVector,hazardLevelsVector, instance, myGetProbabilityFromAllocation)
                                       for hazardLevelsVector in instance.scenarioHazardLevels]
    meanValue = sum(a*b for a,b in zip(probabilitiesForHazardScenarios, instance.scenarioProbs))
    return meanValue

def getProbabilityOfCapLevelScenConsideringHazards_AltNumCapLevels(capLevelsVector, allocLevelsVector, instance, myNumAllocLevels, myNumCapLevels, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForHazardScenarios = [getProbabilityOfCapLevelScenario_AltNumCapLevels_withHazards(capLevelsVector,allocLevelsVector,hazardLevelsVector, instance, myNumAllocLevels, myNumCapLevels, myGetProbabilityFromAllocation)
                                       for hazardLevelsVector in instance.scenarioHazardLevels]
    meanValue = sum(a*b for a,b in zip(probabilitiesForHazardScenarios, instance.scenarioProbs))
    return meanValue

def getProbabilityOfCapLevelScenario(capLevelsVector, allocLevelsVector, instance, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForFacs = [getCapLevelProbabilityForAllocationLevel(fac,capLevelsVector[fac], allocLevelsVector[fac], instance, instance.numCapLevels, myGetProbabilityFromAllocation) 
                            for fac in range(len(allocLevelsVector))]
    #print "capLevelsVector", capLevelsVector, "probabilitiesForFacs", probabilitiesForFacs, "product: ", np.prod(probabilitiesForFacs)
    return np.prod(probabilitiesForFacs)

def getProbabilityOfCapLevelScenario_withHazards(capLevelsVector, allocLevelsVector, hazardLevelsVector, instance, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForFacs = [getCapLevelProbabilityForAllocationLevelAndHazardLevel(fac,capLevelsVector[fac], allocLevelsVector[fac], hazardLevelsVector[fac], instance, instance.numCapLevels, myGetProbabilityFromAllocation) 
                            for fac in range(len(allocLevelsVector))]
    #print "capLevelsVector", capLevelsVector, "hazardLevelsVector", hazardLevelsVector, "probabilitiesForFacs", probabilitiesForFacs, "product: ", np.prod(probabilitiesForFacs)
    return np.prod(probabilitiesForFacs)

def getProbabilityOfCapLevelScenario_AltNumCapLevels(capLevelsVector, allocLevelsVector, instance, myNumAllocLevels, myNumCapLevels, myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForFacs = [getCapLevelProbabilityForAllocationLevel_AltNumAllocLevels(fac,capLevelsVector[fac], allocLevelsVector[fac], myNumAllocLevels, myNumCapLevels, myGetProbabilityFromAllocation) 
                            for fac in range(len(allocLevelsVector))]
    #print "capLevelsVector", capLevelsVector, "probabilitiesForFacs", probabilitiesForFacs, "product: ", np.prod(probabilitiesForFacs)
    return np.prod(probabilitiesForFacs)

#is the actual probability (not the conditional probability)
def getProbabilityOfCapLevelScenario_AltNumCapLevels_withHazards(capLevelsVector, allocLevelsVector, hazardLevelsVector, instance, myNumAllocLevels, myNumCapLevels,
                                                                 myGetProbabilityFromAllocation = getProbabilityFromAllocation):
    probabilitiesForFacs = [getCapLevelProbabilityForAllocationLevelAndHazardLevel_AltNumAllocLevels(fac,capLevelsVector[fac], allocLevelsVector[fac], hazardLevelsVector[fac],
                                                                                       myNumAllocLevels, myNumCapLevels, instance.numHazardLevels, myGetProbabilityFromAllocation) 
                            for fac in range(len(allocLevelsVector))]
    #print "capLevelsVector", capLevelsVector, "probabilitiesForFacs", probabilitiesForFacs, "product: ", np.prod(probabilitiesForFacs)
    return np.prod(probabilitiesForFacs)

def getConstraintExpressionsForBendersCut(coefficientDualTermsAllScenarios, constantDualTerms, iteration, numBunches, allocVarsArg, thetaVarsArg, instance):
    expressions = []
    numScenarios = len(constantDualTerms)
    bunchSize = numScenarios/numBunches
    remainder = numScenarios%numBunches
    
    for bunch in range(numBunches):
        rangeForBunch = range(bunch * bunchSize, (bunch + 1) * bunchSize)
        coefficientDualTermsSumForBunch = [[sum([coefficientDualTermsAllScenarios[scenarioIndex][j][k] for scenarioIndex in rangeForBunch]) for k in range(instance.numAllocLevels)] 
                                           for j in range(instance.numFacs)]
        constantDualsSumForBunch = sum(constantDualTerms[bunch * bunchSize:(bunch + 1) * bunchSize])
        expressions.append(sum([allocVarsArg[j][k] * coefficientDualTermsSumForBunch[j][k] for k in range(instance.numAllocLevels) for j in range(instance.numFacs)]) + constantDualsSumForBunch 
                                        >= thetaVarsArg[bunch])
    if(remainder > 0):
        bunch = numBunches
        rangeForLastBunch = range(numBunches * bunchSize, numScenarios)
        coefficientDualTermsSumForBunch = [[sum([coefficientDualTermsAllScenarios[scenarioIndex][j][k] for scenarioIndex in rangeForLastBunch]) for k in range(instance.numAllocLevels)] 
                                           for j in range(instance.numFacs)]
        constantDualsSumForBunch = sum(constantDualTerms[numBunches * bunchSize : numScenarios])
        expressions.append(sum([allocVarsArg[j][k] * coefficientDualTermsSumForBunch[j][k] for k in range(instance.numAllocLevels) for j in range(instance.numFacs)]) + constantDualsSumForBunch 
                                        >= thetaVarsArg[bunch])
    print "expressions", expressions
    return expressions

def generateRandomScenariosFromAllocVector(allocationSoln, count, instance):
    #sum([a*b for a,b in zip(allocationSoln,range(instance.numAllocLevels))])
    allocLevelsSet = [sum([int(a*b) for a,b in zip(list,range(instance.numAllocLevels))]) for list in allocationSoln]
    probabilities = [getProbabilityFromAllocation(allocLevelsSet[fac], instance.numAllocLevels) for fac in allocLevelsSet]
    return [[np.random.binomial(instance.numCapLevels - 1, probabilities[fac]) for fac in range(instance.numFacs)] for i in range(count)]

class FirstStageModel(object):
    
    gurobiModel = None
    allocVars = None
    thetaVars = None
    costOfAlloc = None
    
    def __init__(self, instance, costOfAlloc, timeLimit, numScenarios, initialUpperBound, numBunches, showOutput = False):
        #bunchSize = numScenarios/numBunches
        #print "numBunches", numBunches
        remainder = numScenarios%numBunches
        #print "numBunches", numBunches
        #print "remainder", remainder
        self.gurobiModel = gurobipy.Model("master problem single cut")
        try:
            self.allocVars = [[self.gurobiModel.addVar(0,1,vtype = gurobipy.GRB.BINARY,name="y_"+str(j)+","+str(k)) for k in range(instance.numAllocLevels)] for j in range(instance.numFacs)]
            self.thetaVars = [self.gurobiModel.addVar(0, initialUpperBound,vtype = gurobipy.GRB.CONTINUOUS,name="theta_"+str(scenNumber)) for scenNumber in range(numBunches)]
            if(remainder > 0):
                self.thetaVars.append(self.gurobiModel.addVar(0, initialUpperBound,vtype = gurobipy.GRB.CONTINUOUS,name="theta_"+str(numBunches)))
            
            self.gurobiModel.update()
            #print "thetaVars", thetaVars
            self.gurobiModel.setObjective(sum(self.thetaVars), gurobipy.GRB.MAXIMIZE)
    
            self.gurobiModel.update()
            # alloc levels constraint
            [self.gurobiModel.addConstr(sum([self.allocVars[j][k] for k in range(instance.numAllocLevels)]) == 1, "one alloc level_" + str(j)) for j in range(instance.numFacs)]
            self.gurobiModel.update()
            # instance.budget constraint
            self.gurobiModel.addConstr(sum([costOfAlloc[j][k] * self.allocVars[j][k] for k in range(instance.numAllocLevels) for j in range(instance.numFacs)]) <= instance.budget, "allocation instance.budget")
            self.gurobiModel.update()
            self.gurobiModel.setParam('OutputFlag', showOutput)
            #print "time_limit", time_limit
            self.gurobiModel.setParam('TimeLimit', timeLimit)
        except gurobipy.GurobiError as e:
            print 'createProbabilityChainMasterProblemModel_SingleCut_GRB: Gurobi Error reported ' + e
    
class SecondStageProblem(object):
    '''
    classdocs
    '''
    gurobiModel = None
    capacityConstraints = None
    instance = None
    assignVars = None
    
    def setInstance(self, instance):
        self.instance = instance
        
    def createModelGurobi(self):
        global gurobiModel,capacityConstraints, assignVars
        print "createModelGurobi"
        gurobiModel = gurobipy.Model("myLP")
        #print "GUROBI MODEL CREATED"
        try:
            # Create variables)
            assignVars = [[gurobiModel.addVar(0,1,vtype=gurobipy.GRB.CONTINUOUS,name="x_"+str(i)+","+str(j)) for j in range(self.instance.numFacs+1)] for i in range(self.instance.numDemPts)]
            # Integrate new variables
            gurobiModel.update()
            # Set objective
            gurobiModel.setObjective(sum([self.instance.demPtWts[i]*self.instance.pairsUtilityMatrix[i][j]*assignVars[i][j] for j in range(self.instance.numFacs+1) for i in range(self.instance.numDemPts)]), gurobipy.GRB.MAXIMIZE)
            gurobiModel.update()
            capacityConstraints = [gurobiModel.addConstr(sum([self.instance.demPtWts[i]*assignVars[i][j] for i in range(self.instance.numDemPts)]) <= self.instance.capacities[j], "capacity_"+str(j)) for j in range(self.instance.numFacs)]
            for i in range(self.instance.numDemPts):
                gurobiModel.addConstr(sum([assignVars[i]  [j] for j in range(self.instance.numFacs+1)]) == 1, "demand_met"+str(i))
            gurobiModel.update()
            gurobiModel.setParam('OutputFlag', False )
            
        except gurobipy.GurobiError as e:
            print 'createTransportationModelGurobi: Gurobi Error reported' + str(e)
            #logging.error('Error reported')
            
    def resetRHSCapacities(self, capacityLevels, myNumCapLevels):
        #print "resetRHSCapacities"
        for fac in range(self.instance.numFacs):
            capacityConstraints[fac].setAttr("rhs", (capacityLevels[fac]/float(myNumCapLevels - 1)) * float(self.instance.capacities[fac]))
        gurobiModel.update()
        
    def computeSecondStageUtility(self, capacityLevels, myNumCapLevels):
        #print "computeSecondStageUtility", capacityLevels
        self.resetRHSCapacities(capacityLevels, myNumCapLevels)
        #print "after reset"
        gurobiModel.optimize()
        #print "capLevels: ", [capacityConstraints[fac].getAttr("rhs") for fac in range(self.instance.numFacs)]
        #print "capLevels", capacityLevels, "objVal", gurobiModel.objVal
        if(writeToFile):
            gurobiModel.write("/home/hmedal/Documents/Temp/secondStage_" + str(capacityLevels) + str(gurobiModel.objVal) + "_" +"_testOutput.lp")
        return gurobiModel.objVal
    
    def computeSecondStageUtilityAndSolution(self, capacityLevels, myNumCapLevels):
        #print "computeSecondStageUtility", capacityLevels
        self.resetRHSCapacities(capacityLevels, myNumCapLevels)
        #print "after reset"
        gurobiModel.optimize()
        #print "capLevels: ", [capacityConstraints[fac].getAttr("rhs") for fac in range(self.instance.numFacs)]
        #print "capLevels", capacityLevels, "objVal", gurobiModel.objVal
        if(writeToFile):
            gurobiModel.write("/home/hmedal/Documents/Temp/secondStage_" + str(capacityLevels) + str(gurobiModel.objVal) + "_" +"_testOutput.lp")
        assignSoln = [[assignVars[i][j].X for j in range(self.instance.numFacs+1)] for i in range(self.instance.numDemPts)]
        #return assignSoln, 
        return assignSoln, gurobiModel.objVal
    
    def computeAverageSecondStageUtility_Sampling(self, allocationVector, numSamples, numThreads, instance):
        localPoolVar = Pool(processes = int(numThreads))
        samples = [getRandomIndependentBinomials(allocationVector, instance.numAllocLevels, instance.numCapLevels, instance.numFacs) for i in range(numSamples)]
        meanValue = np.mean(localPoolVar.map(self.computeSecondStageUtility, samples))
        localPoolVar.terminate()
        return meanValue
    
    def computeAverageSecondStageUtility_SamplesGiven(self, allocationVector, samples, instance):
        localPoolVar = Pool()
        objVals = localPoolVar.map(self.computeSecondStageUtility, samples)
        scenarioProbs = [getProbabilityOfCapLevelScenario(sample, allocationVector, instance) for sample in samples]
        meanValue = sum([a*b for a,b in zip(objVals, scenarioProbs)])
        localPoolVar.terminate()
        return meanValue