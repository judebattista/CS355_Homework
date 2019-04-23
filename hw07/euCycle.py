import utilities
import copy
import re
from random import randrange, seed

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

keys = edges.keys()
maxKey = max(keys)
minKey = min(keys)
keyRange = maxKey - minKey
path = []
startNode = randrange(minKey, maxKey + 1)
path.append(startNode)
nextNode = (startNode + 1) % keyRange
workingCopy = copy.deepcopy(edges)
currentNode = startNode
print(workingCopy)
print(startNode)
adjacentNodes = workingCopy[startNode]
print(adjacentNodes)
while len(adjacentNodes) > 0:
    nextNode = adjacentNodes[0]
    path.append(nextNode)
    print(path)
    del adjacentNodes[0]
    adjacentNodes = workingCopy[nextNode]
print(adjacentNodes)     
print(workingCopy)
print(path)
with open('euCycle.results.txt', 'w') as outfile:
    pass
