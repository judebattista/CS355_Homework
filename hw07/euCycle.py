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
    path.append(nextNode)
    return path

with open('euCycle.txt', 'r') as infile:
    edges = {}
    for line in infile:
        splitLine = re.split('->', line.strip())
        key = splitLine[0].strip() 
        vals = re.split(',', splitLine[1].strip())
        #print(splitLine, key, vals)
        for val in vals:
            utilities.appendToDict(edges, int(key), int(val))
    print(edges)

workingCopy = copy.deepcopy(edges)
nextNode = choice(list(workingCopy.keys()))
adjacentNodes = workingCopy[nextNode]

# create the initial cycle:
path = buildCycle(workingCopy, nextNode)

edgeCount = 0
for item in edges.items():
    edgeCount += len(item[1])

# find epicycles
# walk the path to find a node that has a 
for ndx in range(0, edgeCount+1):
    # if a node in the cycle has an unexplored edge...
    print('ndx = {0}, path len = {1}'.format(ndx, len(path)))
    node = path[ndx]
    if len(workingCopy[node]) > 0:
        # build a cycle starting from node
        subPath = buildCycle(workingCopy, node)
        path = path[0:ndx] + subPath + path[ndx+1:]
        

#print(adjacentNodes)     
#print("Working copy: ", workingCopy)
#print(path)
#printEdges(path)
print('Edges in graph: {0}, edges in path: {1}'.format(edgeCount, len(path)-1))
with open('euCycle.results.txt', 'w') as outfile:
    strPath = [str(item) for item in path]
    output = '->'.join(strPath)
    outfile.write(output)
