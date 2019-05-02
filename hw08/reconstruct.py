import utilities
import copy
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
    #startNode = choice(degrees.keys())
    #endNode = choice(degrees.keys()) 
    # if we don't care about getting random end points we can just use the first key
    startNode = next(iter(degrees))
    endNode = next(iter(degrees))
    
    for node in degrees.keys():
        if degrees[node] == -1:
            startNode = node
        if degrees[node] == 1:
            endNode = node
    return startNode, endNode

def genPath(bases, degrees):
    workingCopy = copy.deepcopy(bases)
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
            #print('Removing {0} from {1}'.format(nextNode, neighbors))
            workingCopy[currentNode].remove(nextNode)
            currentNode = nextNode
        else:   # No neighbors: add to path, remove last item from stack and set as current node
            path.append(currentNode)
            currentNode = stack.pop(-1)
        #print('Working copy: {0}'.format(workingCopy))
        #print('Path: {0}'.format(path))
        #print('Stack: {0}'.format(stack))
        #print('Current node: {0}, neighbors: {1}'.format(currentNode, neighbors))
        neighbors = workingCopy[currentNode]

    path.append(currentNode)
    return path[::-1]

def buildStringFromPath(path):
    stringArr = []
    stringArr.append(path[0])
    for base in path[1::]:
        stringArr.append(base[-1])
    return ''.join(stringArr)


with open('reconstruct.txt', 'r') as infile:
    klen = int(infile.readline())
    kmers = []
    for line in infile:
        kmers.append(line.strip())

bases = {}
degrees = {}

for kmer in kmers:
    key = kmer[0:klen-1]
    val = kmer[1:klen]
    utilities.appendToDict(bases, key, val)
    if val not in bases:
        bases[val] = []
    utilities.decrementDict(degrees, key)
    utilities.incrementDict(degrees, val)


path = genPath(bases, degrees)
print(path)
ntide = buildStringFromPath(path)
print(ntide)

with open('reconstruct.results.txt', 'w') as outfile:
    outfile.write(ntide)

