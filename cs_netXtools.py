#!/usr/bin/python

# A collection of additional calculations that I find useful for networkX objects
import random as rand
import networkx as nx
import math
import numpy as np
import json
from cs_dictTools import populateDictOfLists as popDoL

# Takes a networkX graph object and calculates the size of the largest connected component.
def size_lcc(myGraph):
    nodeList_cc = sorted(nx.connected_components(myGraph))
    # give the sizes of the lists
    size_cc = [len(nodeList) for nodeList in nodeList_cc]
    # find index of the largest one (should be first one...)
    i_largest = size_cc.index(max(size_cc))
    # return the number of nodes in the largest connected component
    cc_largest = len(nodeList_cc[i_largest])

    return cc_largest

# Takes a networkX graph object and returns subgraph containing largest connected component.
def subgraph_lcc(myGraph):
    nodeList_cc = sorted(nx.connected_components(myGraph))
    # give the sizes of the lists
    size_cc = [len(nodeList) for nodeList in nodeList_cc]
    # find index of the largest one (should be first one...)
    i_largest = size_cc.index(max(size_cc))
    # return the subgraph of largest connected component
    newGraph = myGraph.subgraph(nodeList_cc[i_largest])

    return newGraph

# Takes a networkX graph and subgraph and calculates the significance of the lcc in the subgraph
def sig_lccOfSubgraph(fullgraph, nodeList, numBins=10, numSimulations=10000):
    lccSizes_preserveDeg = []
    tempGraph = fullgraph.subgraph(nodeList)
    nodesInGraph = tempGraph.nodes()
    targetDegreeDist = getDict_degree2numNodes(fullgraph, nodeList=nodesInGraph, numBins=numBins)
    lcc_subgraph = size_lcc(fullgraph.subgraph(nodeList))
    for x in range(numSimulations):
        currSimGraph = randSubgraph_preserveDegree(fullgraph, targetDegreeDist, binFlag='yes')
        # calculate lcc for random network subgraph
        cc_largest = size_lcc(currSimGraph)
        lccSizes_preserveDeg.append(cc_largest)
    pVal = (sum(i >= lcc_subgraph for i in lccSizes_preserveDeg))/float(numSimulations)
    zScore = (lcc_subgraph - np.mean(lccSizes_preserveDeg)) / np.std(lccSizes_preserveDeg)

    retStruct = {}
    retStruct['lcc'] = lcc_subgraph
    retStruct['pVal'] = pVal
    retStruct['zScore'] = zScore
    retStruct['simResults'] = lccSizes_preserveDeg

    return retStruct

# Takes a networkX graph and subgraph and calculates the significance of the diameter in the subgraph
def sig_diameter(fullgraph, nodeList, numBins=10, numSimulations=10000):
    diameter_preserveDeg = []
    currSubgraph = fullgraph.subgraph(nodeList)
    currSubgraph_lcc = subgraph_lcc(currSubgraph)
    nodesInGraph = currSubgraph_lcc.nodes()
    targetDegreeDist = getDict_degree2numNodes(fullgraph, nodeList=nodesInGraph, numBins=numBins)
    diameter_subgraph = nx.diameter(currSubgraph_lcc)
    for x in range(numSimulations):
        currSimGraph = randSubgraph_preserveDegree(fullgraph, targetDegreeDist, binFlag='yes')
        # calculate distance for lcc of random network subgraph
        currDiameter = nx.diameter(subgraph_lcc(currSimGraph))
        diameter_preserveDeg.append(currDiameter)
    pVal = (sum(i >= diameter_subgraph for i in diameter_preserveDeg))/float(numSimulations)
    zScore = (diameter_subgraph - np.mean(diameter_preserveDeg)) / np.std(diameter_preserveDeg)

    retStruct = {}
    retStruct['diameter'] = diameter_subgraph
    retStruct['pVal'] = pVal
    retStruct['zScore'] = zScore
    retStruct['simResults'] = diameter_preserveDeg

    return retStruct


# Takes a networkX graph object and returns dictionary with the degree distribution of the graph.  The keys are the degrees and the values are the nodes.
def getDict_degree2nodeNames(myGraph, **keyword_parameters):
    deg2nodeNames = {}
    if ('nodeList' in keyword_parameters):
        # remove entries from nodeList that are not in graph
        nodeList = nx.subgraph(myGraph, keyword_parameters['nodeList']).nodes()
        for currNode in nodeList:
            #make dictionary where key is degree and value is the list of nodes
            currDeg = myGraph.degree(currNode)
            popDoL(deg2nodeNames, currDeg, currNode)
    else:
        for currNode in myGraph.nodes():
            # make dictionary where key is degree, and value is list of nodes
            currDeg = myGraph.degree(currNode)
            popDoL(deg2nodeNames, currDeg, currNode)
    # at this point, we have generated deg2nodeNames dictionary either for all nodes, or for the subset of nodes defined by optional argument nodeList

    if ('numBins' in keyword_parameters):
        degDist = deg2nodeNames
        numBins = keyword_parameters['numBins']
        # calculate desired bin width
        binWidth = (math.log10(max(degDist.keys())) - math.log10(min(degDist.keys()))) / numBins

        lb = 10 ** (binWidth * np.array(range(numBins)))
        ub = 10 ** (binWidth * (np.array(range(numBins))+1))

        # rewrite degDist dictionary updated for binning
        degDist_binned = {}
        for binNum in range(numBins):
            # initiate dictionary space
            binKey = ub[binNum] #round(ub[binNum],3)
            degDist_binned[binKey] = []
            # function to evaluate whether key is in bin
            def meetCondition(element):
                return bool(element > lb[binNum] and element < ub[binNum])
            # return keys that are in the bin
            keysInBin = filter(meetCondition, degDist.keys())
            summedNumNodes = 0
            for currKey in keysInBin:
                degDist_binned[binKey].extend(degDist[currKey])
        return degDist_binned

    elif ('binLimit' in keyword_parameters):
        degDist = deg2nodeNames
        ub = sorted(keyword_parameters['binLimit'])
        lb = [-1]
        lb.extend(ub[:-1])

        # rewrite degDist dictionary updated for binning
        degDist_binned = {}

        for binNum in range(len(ub)):
            # initiate dictionary space
            binKey = ub[binNum]
            degDist_binned[binKey] = []
            # function to evaluate whether key is in bin
            def meetCondition(element):
                return bool(element > lb[binNum] and element < ub[binNum])
            # return keys that are in the bin
            keysInBin = filter(meetCondition, degDist.keys())
            summedNumNodes = 0
            for currKey in keysInBin:
                degDist_binned[binKey].extend(degDist[currKey])
        return degDist_binned

    else:
        return deg2nodeNames

# Takes networkX graph object and returns dictionary with the degree distribution of the graph.  The keys are the degrees and the values are the number of nodes with that degree.
def getDict_degree2numNodes(myGraph, **keyword_parameters):
    # in this function, we just run getDict_degree2nodeNames and then count the number of nodes inside.
    if ('nodeList' in keyword_parameters):
        if ('numBins' in keyword_parameters):     # nodeList & numBins
            deg2nodeNames = getDict_degree2nodeNames(myGraph,                   nodeList=keyword_parameters['nodeList'], numBins=keyword_parameters['numBins'])
        elif ('binLimit' in keyword_parameters):    # nodeList & binLimits
            deg2nodeNames = getDict_degree2nodeNames(myGraph,                   nodeList=keyword_parameters['nodeList'], binLimit=keyword_parameters['binLimit'])
        else:   # only nodeList
            deg2nodeNames = getDict_degree2nodeNames(myGraph,               nodeList=keyword_parameters['nodeList'])
    else:   # no nodeList
        if ('numBins' in keyword_parameters):     # nodeList & numBins
            deg2nodeNames = getDict_degree2nodeNames(myGraph,                   numBins=keyword_parameters['numBins'])
        elif ('binLimit' in keyword_parameters):    # nodeList & binLimits
            deg2nodeNames = getDict_degree2nodeNames(myGraph, binLimit=keyword_parameters['binLimit'])
        else:
            deg2nodeNames = getDict_degree2nodeNames(myGraph)

    deg2numNodes = {}
    for currKey in deg2nodeNames:
        # make dictionary where key is degree, and value is list of nodes
        deg2numNodes[currKey] = len(deg2nodeNames[currKey])

    return deg2numNodes

# Takes networkX graph object and returns random subgraph with numNodes nodes.
def randSubgraph(ogGraph, numNodes):
    # random sampling of networks with size same as iuis list
    list_rand = rand.sample(list(ogGraph), numNodes)
    # generate subgraph of the randomly selected genes
    randSubgraph = ogGraph.subgraph(list_rand)

    return randSubgraph

# Takes networkX graph object and returns a subgraph with user defined degree distribution of nodes.
# ogGraph is the source graph.  degDist is as dictionary containing the degrees as keys and the number required of each degree as values.  (generate this with getDict_degree2numNodes(myGraph). don't forget to get the degree distribution on the OG GRAPH!!)
def randSubgraph_preserveDegree(ogGraph, degDist, binFlag):
    # 0print(degDist)
    if binFlag == 'yes':
        currBinLimits = sorted(degDist.keys())
        # print('we made it into preserveDegree yes binFlag')
        currBinDict = getDict_degree2nodeNames(ogGraph, binLimit=currBinLimits)
        currBinDictNums = getDict_degree2numNodes(ogGraph, binLimit=currBinLimits)
        # print('donor graph degree distribution')
        # print(currBinDictNums)
        # print('desired degree distribution')
        # print(degDist)
    elif binflag == 'no':
        currBinDict = getDict_degree2nodeNames(ogGraph)
    else:
        print('ERROR: binFlag must be yes or no')

    degMatchedRandNodeNames = []
    # for each value in the desired degree bin,
    for currBin in degDist:
        numNodesNeeded = degDist[currBin]
        candNodes = rand.sample(list(currBinDict[currBin]), numNodesNeeded)
        degMatchedRandNodeNames.extend(candNodes)
        # print(len(candNodes), ' | ', len(degMatchedRandNodeNames))

    newGraph = ogGraph.subgraph(degMatchedRandNodeNames)
    return newGraph

# Takes networkX graph and returns backbone of the graph
def graphBackbone(ogGraph, co_alpha, weightAttr):
    G_backbone = ogGraph.copy()
    # calculate sum of weights at every node
    sumWeights_node = {}
    for eaNode in G_backbone.nodes():
        currSum = 0
        for eaNeighbor in nx.neighbors(G_backbone, eaNode):
            weight = G_backbone[eaNode][eaNeighbor][weightAttr]
            currSum += weight
        sumWeights_node[eaNode] = currSum
    # check edges from both directions for sig fr uni, store new graph in G_result
    G_result = nx.Graph()
    for u,v in G_backbone.edges():
        w = G_backbone[u][v][weightAttr]
        pa = w/sumWeights_node[u]
        ka = nx.degree(G_backbone, u)
        alpha1 = pa + (1.-pa)**(ka-1)
        pb = w/sumWeights_node[v]
        kb = nx.degree(G_backbone, v)
        alpha2 = pb + (1.-pb)**(kb-1)

        if alpha1 < co_alpha or alpha2 < co_alpha:
            G_result.add_edge(u, v, weight=w)
    return G_result


# Takes networkX graph + two node lists and returns the distances between all the points.
def dist_btw2subgraphs(myGraph, nodes1, nodes2, moduleDef='allNodes', distDef='pairwise',returnType = 'average'):
    # if flag 'center' is indicated, the distances between the centers of the subgraphs will be calculated.  Definition of center requires that all the components in the subgraph are connected, hence call to subgraph_lcc.
    if moduleDef == 'center':
        leftNodes = nx.center(subgraph_lcc(myGraph.subgraph(nodes1)))
        rightNodes = nx.center(subgraph_lcc(myGraph.subgraph(nodes2)))
    elif moduleDef == 'allNodes':
        leftNodes = nodes1
        rightNodes = nodes2
    else:
        raise NameError('moduleDef must be allNodes or center')

    # calculate pairwise distances
    pairDist = {}
    for eaOrigin in leftNodes:
        pairDist[eaOrigin] = []
        for eaTarget in rightNodes:
            pairDist[eaOrigin].append(nx.shortest_path_length(myGraph, source=eaOrigin, target=eaTarget))
    # if closest flag on, then find the min in each dict of pairwise distances
    pooledDist = []
    if distDef == 'pairwise':
        pooledDist = [dist for sublist in pairDist.values() for dist in sublist]
    elif distDef == 'closest':
        pooledDist = [min(pairDist[node]) for node in pairDist]
    else:
        raise NameError('distDef must be pairwise or closest')

    if returnType == 'average':
        aveDist = sum(pooledDist)/float(len(pooledDist))
        return aveDist
    elif returnType == 'all':
        return pooledDist
    else:
        raise NameError('returnType must be average or all')


def getPosfrJson(fileName):
    posDict = {}
    ogPos = json.load(open(fileName))
    for eaNode in ogPos["elements"]["nodes"]:
        nodeName = str(eaNode["data"]["shared_name"])
        posTuple = (eaNode["position"]["x"], eaNode["position"]["y"])
        posDict[nodeName] = posTuple
    return posDict

def getPosfrJson3D(fileName):
    posDict = {}
    ogPos = json.load(open(fileName))
    for eaNode in ogPos["elements"]["nodes"]:
        nodeName = str(eaNode["data"]["shared_name"])
        posTuple = (eaNode["position"]["x"], eaNode["position"]["y"], eaNode["position"]["z"])
        posDict[nodeName] = posTuple
    return posDict
