import itertools
import utilities


nTides = ['A','C','G','T']




with open('medianString.txt', 'r') as infile:
    kLen = int(infile.readline().strip())
    fragments = []
    for line in infile:
        fragments.append(list(line.strip()))

allKmers = list(itertools.product(nTides, repeat=kLen))
zeros = itertools.cycle([0])
kmerDict = {k:v for k, v in zip(allKmers, zeros)}
bestKmer = None
minDist = kLen
for fragment in fragments:
    fragMinDist = kLen
    subfrags = [fragment[foo:foo+kLen] for foo in range(0, len(fragment) - kLen + 1)]
    for kmer in allKmers:
        kmerMinDist = kLen
        #kmerDict[kmer] += utilities.hammingDistance(kmer, subfrag)
        for subfrag in subfrags:
            currentDist = utilities.hammingDistance(kmer, subfrag)
            if currentDist < kmerMinDist:
                kmerMinDist = currentDist
        kmerDict[kmer] += kmerMinDist
medianString = ''.join(min(kmerDict, key=kmerDict.get)) 

with open('medianString.results.txt', 'w') as outfile:
    outfile.write(medianString)
