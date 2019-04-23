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
startNode = randrange(minKey, maxKey + 1)
nextNode = startNode + 1
workingCopy = copy.deepcopy(edges)
print(workingCopy)
adjacentNodes = workingCopy[startNode]
if len(adjacentNodes) > 0:
    nextNode = adjacentNodes[0]
    del adjacentNodes[0]
    

with open('euCycle.results.txt', 'w') as outfile:
    pass
