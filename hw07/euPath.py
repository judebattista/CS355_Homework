import utilities
import copy
import re
from random import randrange, seed, choice

def printEdges(path):
    for ndx in range(0, len(path) - 1):
        print(str(path[ndx]) + ' -> ' + str(path[ndx+1])) 

def pickRandomChild(children):
    maxKid = len(children) - 1
    if maxKid == 0:
        return 0
    return randrange(maxKid)

def findEndNodes(degrees):
    startNode = 0
    endNode = 0
    for node in degrees.keys():
        if degrees[node] == -1:
            startNode = node
        if degrees[node] == 1:
            endNode = node
    return startNode, endNode


with open('euPath.txt', 'r') as infile:
    edges = {}
    revEdges = {}
    degrees = {}
    for line in infile:
        splitLine = re.split('->', line.strip())
        key = int(splitLine[0].strip()) 
        vals = re.split(',', splitLine[1].strip())
        #print(splitLine, key, vals)
        for val in vals:
            intVal = int(val)
            if intVal not in edges:
                edges[intVal] = []
            utilities.decrementDict(degrees, key)
            utilities.appendToDict(edges, key, intVal)
            utilities.appendToDict(revEdges, intVal, key)
            utilities.incrementDict(degrees, intVal)
    print('Degrees: {0}'.format(degrees))

workingCopy = copy.deepcopy(edges)
print('Edges: {0}'.format(workingCopy))

startNode, endNode = findEndNodes(degrees)
print('Start node: {0}, end node: {1}'.format(startNode, endNode))
path = []
stack = []
currentNode = startNode
neighbors = workingCopy[currentNode]
while len(neighbors) > 0 or len(stack) > 0:
    #neighbors = workingCopy[currentNode]
    if len(neighbors) > 0: # add to stack, pick neighbor, remove edge, neighbor = current
        stack.append(currentNode)
        childNdx = pickRandomChild(neighbors)
        nextNode = neighbors[childNdx]
        print('Removing {0} from {1}'.format(nextNode, neighbors))
        workingCopy[currentNode].remove(nextNode)
        currentNode = nextNode
    else:   # No neighbors: add to path, remove last item from stack and set as current node
        path.append(currentNode)
        currentNode = stack.pop(-1)
    neighbors = workingCopy[currentNode]
    print('Path: {0}'.format(path))
    print('Stack: {0}'.format(stack))
    print('Current node: {0}, neighbors: {1}'.format(currentNode, neighbors))

path.append(currentNode)
path = path[::-1]
print('Final path: {0}'.format(path))


with open('euPath.results.txt', 'w') as outfile:
    strPath = [str(item) for item in path]
    output = '->'.join(strPath)
    outfile.write(output)
