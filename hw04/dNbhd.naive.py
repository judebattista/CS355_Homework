import utilities
import itertools

def getIndiceTuples(kmerLen, distance):
    return list(itertools.combinations(range(kmerLen), distance))

def mutSlices(listOfSlices, alphabet, tail):
    fragCount = len(listOfSlices)
    itMutations = itertools.product(listOfSlices, alphabet)
    mutations = [''.join(list(tup)) for tup in itMutations]
    splitMutes = [mutations[ndx * 4:ndx*4 + 4] for ndx in range(fragCount)]
    cartesian = itertools.product(*splitMutes)
    return [''.join(list(tup)) + tail for tup in cartesian]

with open('dNbhd.naive.txt', 'r') as infile:
    root = infile.readline().strip()
    print(root)
    distance = int(infile.readline().strip())
'''
algorithm:
    Find all the indices where we need to do a replacement
    Split each string at the indices, discarding the character at each index
    Make four new strings from each fragment and append A, C, G, and T to each
    Append each fragment to construct a cousin string
    
'''

# The letters used as replacements in our string
alphabet = 'ACGT'
kmerLen = len(root)
indices = getIndiceTuples(kmerLen, distance)
cousinFrags = []
lastNdx = 0
cousins = []
# break apart the pattern string at each index in the set of indices
for indexSet in indices:
    for ndx in indexSet:
        kmerSlice = root[lastNdx:ndx]
        cousinFrags.append(kmerSlice)
        lastNdx = ndx + 1
    # everything left over after the final split
    tail = (root[lastNdx :])
    # Append every letter of the alphabet to an instance of every fragment
    print(cousinFrags)
    #print(cousinFragMutes)
    cousins.append(mutSlices(cousinFrags, alphabet, tail))

with open('dNbhd.naive.results.txt', 'w') as outfile:
    utilities.writeListToFileOnNewlines(outfile, cousins)
#cousin = ''.join(cousinFragMutes)
#print(cousin)
