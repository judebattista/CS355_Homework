import utilities
import itertools

# Given a set of fragments, a length k, and a distance d
# Find every kmer of length k within a distance d of every substring of length k in the first fragment
# Remove it from the list if it is not within distance d of a substring in every fragment in the list

subfragDict = {}
fragments = []
nTides = 'ACGT' 

# for a given kmer, calculate all the cousins of degree dist using the specified alphabet
# store the results in a dictionary
def genCousinsDremoved(kmer, alphabet, dist, dictionary):
    # Generate all the possible locations for replacement
    for positions in itertools.combinations(range(len(kmer)), dist):
        #generate a every possible set of values to fill in each set of locations
        # ... kind of. We will not include the last letter in the alphabet
        for replacements in itertools.product(range(len(alphabet) - 1), repeat=dist):
            # work with a list rather than a string for efficiency
            cousin = list(kmer)
            for pos, rep in zip(positions, replacements):
                # if the letter is the same as the replacement, replace it with the 
                # 'spare' letter from the alphabet
                if cousin[pos] == alphabet[rep]:
                    cousin[pos] = alphabet[-1]
                else:
                    cousin[pos] = alphabet[rep]
            dictionary[''.join(cousin)] = 0

# check to see if a pattern is in hamming distance of any substring in fragment
# fragment and pattern should be lists of chars, distance an int
def patternNearFrag(fragment, pattern, distance):
    length = len(pattern)
    #pattern = list(pattern)
    subfrags = [fragment[foo:foo+length] for foo in range(0, len(fragment) - length + 1)] 
    #matchingFrags = [1 for x in subfrags if utilities.hammingDistance(pattern, x) <= distance]
    for subfrag in subfrags:
        #utilities.hammingDistance(pattern, subfrag) 
        if utilities.hammingDistance(pattern, subfrag) <= distance:
            return True
    #print('Using fragment %s and pattern %s, the minimum Hamming distance was %d' % (fragment, pattern, bar))
    return False

# for given fragment, look at each substring of length kLen
# Then create all the possible cousins within distance degrees of each substring
def genSubfragCousins(fragment, kLen, distance):
    subfrags = [fragment[foo:foo+kLen] for foo in range(0, len(fragment) - kLen + 1)] 
    for subfrag in subfrags:
        for d in range(0, distance+1):
            genCousinsDremoved(subfrag, nTides, d, subfragDict)

# for a given fragment, check each cousin to see if the cousin is within degree distance 
# from any substring of the fragment with the appropriate length
def filterCousinsByFrag(fragment, cousins, dist):
    filters = [patternNearFrag(fragment, cousin, dist) for cousin in cousins]
    cousins = list(itertools.compress(cousins, filters))
    return cousins

with open('motifEnum.txt', 'r') as infile:
    intParams = utilities.readIntListFromFile(infile)
    kLen = intParams[0]
    distance = intParams[1]
    for line in infile:
        fragments.append(list(line.strip()))

# generate all the cousins of the subfragments of fragment 0
genSubfragCousins(fragments[0], kLen, distance)
cousins =list(subfragDict.keys())
for fragment in fragments[1:]:
    cousins = filterCousinsByFrag(fragment, cousins, distance)

with open('motifEnum.results.txt', 'w') as outfile:
    utilities.writeListToFile(outfile, cousins)
