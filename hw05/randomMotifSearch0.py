import itertools
from random import randrange, seed
import utilities

letterToNum = {'A':0, 'C':1, 'G': 2, 'T':3}
alphabet = ['A','C','G','T']

# Score a set of strings
def scoreMotifs(motifs, alphabet):
    pass

# Generate a profile
def buildProfile(kLen, motifs, alphabet):
    pass

# Generate motifs based on a profile
def buildMotifsFromProfile(kLen, profile, alphabet):
    pass

def distillMotifs(kLen, motifs, dna, alphabet):
    bestScore = kLen * len(motifs)
    bestMotifs = []
    while True:
        profile = buildProfile(kLen, motifs, alphabet)
        newMotifs = buildMotifsFromProfile(kLen, profile, dna, alphabet)
        newScore = scoreMotifs(newMotifs, alphabet)
        if newScore < bestScore:
            bestScore = newScore
            bestMotifs = newMotifs
        else:
            return bestScore, bestMotifs

def selectRandomMotifs(kLen, dna):
    maxNdx = len(dna) - kLen
    randomMotifs = []
    for frag in dna:
        ndx = randrange(maxNdx)
        randomMotifs.append(frag[ndx:ndx+kLen])
    return randomMotifs

def run(iterations, alphabet):
    with open('randomMotifSearch.txt', 'r') as infile:
        intParams = utilities.readIntListFromFile(infile)
        kLen = intParams[0]
        numFrags = intParams[1]
        dna = []
        for line in infile:
            dna.append(list(line.strip()))
    bestScore = kLen * numFrags # worst possible score
    bestMotifs = []

    # pick random kmers as seed
    # distill the kmers to find a locally best set of kmers
    # check to see if they are better than our current set
    for ndx in range(0, iterations):
        if ndx % 1000 == 0:
            seed()
        seedMotifs = selectRandomMotifs(kLen, dna)
        localBestScore, localBestMotifs = distillMotifs(kLen, seedMotifs, dna, alphabet)
        if localBestScore < bestScore:
            bestScore = localBestScore
            bestMotifs = localBestMotifs

    with open('randomMotifSearch.results.txt', 'w') as outfile:
        bestMotifStrings = [''.join(motifs) for motif in bestMotifs]
        utilities.writeListToFineOnNewlines(outfile, bestMotifStrings)
