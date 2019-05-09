import utilities
import copy
import itertools as it
from random import randrange, seed, choice

# build dictionaries containing the prefix and suffix values associated with
# each middle substring
def constructEdges(kLen, frags):
    edges = {}
    for ndx in range(0, len(frags)):
        key = frags[ndx][0:kLen - 1]
        val = frags[ndx][1:kLen]
        utilities.appendToDict(edges, key, val)
    return edges

def constructGraph(edges):
    graph = [(k, v) for k, v in edges.items()]
    return graph

def printEdges(path):
    for ndx in range(0, len(path) - 1):
        print(str(path[ndx]) + ' -> ' + str(path[ndx+1])) 

def pickRandomChild(children):
    maxKid = len(children) - 1
    if maxKid == 0:
        return 0
    return randrange(maxKid)

def buildCycle(graph, nextNode):
    path = []
    adjacentNodes = workingCopy[nextNode]
    #print('Working copy: {0}'.format(workingCopy))

    # Create the initial cycle
    while len(workingCopy[nextNode]) > 0:
        #print('-----------------------')
        currentNode = nextNode
        path.append(currentNode)
        adjacentNodes = workingCopy[currentNode]
        #print('Current node: {0}'.format(currentNode))
        randomChild = pickRandomChild(adjacentNodes)  
        nextNode = adjacentNodes[randomChild]
        #print('Next node: {0}'.format(nextNode))
        #path.append(nextNode)
        del adjacentNodes[randomChild]
        #print('Path: {0}'.format(path))
        #print('Adjacent nodes: {0}'.format(adjacentNodes)) 
        #print('Working copy: {0}'.format(workingCopy))
    # once we run out of nodes to append to the main cycle we should be back at the 
    # starting node
    # path.append(nextNode)
    return path


# find epicycles
# walk the path to find a nodes that have unexplored edges
def findEpicycles(workingCopy, edgeCount, path):
    print('Edge count: {0}'.format(edgeCount))
    for ndx in range(0, edgeCount+1):
        # if a node in the cycle has an unexplored edge...
        print('ndx = {0}, path len = {1}'.format(ndx, len(path)))
        print('path: {0}'.format(path))
        node = path[ndx]
        if len(workingCopy[node]) > 0:
            # build a cycle starting from node
            subPath = buildCycle(workingCopy, node)
            path = path[0:ndx] + subPath + path[ndx+1:]
    return path

def buildStringFromPath(path):
    stringArr = []
    stringArr.extend(path[0])
    for base in path[1::]:
        stringArr.append(base[-1])
    print(stringArr)
    return ''.join(stringArr)

def printEdges(dictionary):
    for key in dictionary.keys():
        print('{0} -> {1}'.format(key, dictionary[key]))

kmers = []
with open('kUniversal.txt', 'r') as infile:
    kLen = int(infile.readline().strip())
    alphabet = ['0', '1']
    lex = it.product(alphabet, repeat=kLen)
    # note: iterating over lex destroys it
    #for item in lex:
        #print(item)
    for kmer in lex:
        kmers.append(kmer) 
    
edges = constructEdges(kLen, kmers)
graph = constructGraph(edges)
print(graph[0])


edgeCount = 2**kLen - kLen

with open('kUniversal.results.txt', 'w') as outfile:
    #strPath = [str(item) for item in path]
    #output = '->'.join(strPath)
    #outfile.write(output)
    pass
