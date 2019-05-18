import utilities
import copy
import itertools as it
from random import randrange, seed, choice

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
    #print('Edge count: {0}'.format(edgeCount))
    for ndx in range(0, edgeCount+1):
        # if a node in the cycle has an unexplored edge...
        #print('ndx = {0}, path len = {1}'.format(ndx, len(path)))
        #print('path: {0}'.format(path))
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

edges = {}
with open('kUniversal.txt', 'r') as infile:
    kLen = int(infile.readline().strip())
    alphabet = ['0', '1']
    #edges = {}
    degrees = {}
    lex = it.product(alphabet, repeat=kLen)
    # note: iterating over lex destroys it
    #for item in lex:
        #print(item)
    for kmer in lex:
        key = kmer[0:kLen-1]
        val = kmer[1:kLen]
        if (val != key):
            utilities.appendToDict(edges, key, val)
        #utilities.appendToDict(edges, key, val)
        if val not in edges:
            edges[val] = []
        utilities.decrementDict(degrees, key)
        utilities.incrementDict(degrees, val)
    #print('edges: {0}'.format(edges))
    
    #printEdges(edges)

#print('edges: {0}'.format(edges))
workingCopy = copy.deepcopy(edges)
#print('workingCopy: {0}'.format(workingCopy))
nextNode = choice(list(workingCopy.keys()))
adjacentNodes = workingCopy[nextNode]

# create the initial cycle:
path = buildCycle(workingCopy, nextNode)

#edgeCount = 0
#for item in edges.items():
#    edgeCount += len(item[1])

edgeCount = 2**kLen - kLen - 1

#print('path: {0}'.format(path))
#print('Edges in graph: {0}, edges in path: {1}'.format(edgeCount, len(path)-1))
path = findEpicycles(workingCopy, edgeCount, path)

#print(adjacentNodes)     
#print("Working copy: ", workingCopy)
#printEdges(path)
#print('Edges in graph: {0}, edges in path: {1}'.format(edgeCount, len(path)-1))
print('Path: {0}'.format(path))
output = buildStringFromPath(path)

print(output)
with open('kUniversal.results.txt', 'w') as outfile:
    #strPath = [str(item) for item in path]
    #output = '->'.join(strPath)
    outfile.write(output)
