import itertools
from random import randrange
import utilities

letterToNum = {'A':0, 'C':1, 'G': 2, 'T':3}
alphabet = ['A','C','G','T']

def fragmentProbability(fragment, matrix):
    kLen = len(fragment)
    probability = 1
    for foo in range(0, kLen):
        nTide = letterToNum[fragment[foo]]
        probability *= matrix[nTide][foo]
    #print('Probability of fragment %s is %f.' % (fragment, probability))
    return probability

def findMostProbable(text, kLen, matrix):
    fragments = [text[foo:foo+kLen] for foo in range(0, len(text) - kLen + 1)]
    probabilities = {''.join(fragment):fragmentProbability(fragment, matrix) for fragment in fragments}
    #print('----------')
    return list(max(probabilities, key=probabilities.get))
    
# profile will have len(alphabet) rows and kLen cols
def buildProfile(kLen, motifs, alphabet):
    # The denominator for every fraction in our profile should be the number of rows in our set of motifs + 1 for every letter of our alphabet.
    denom = len(motifs) + len(alphabet)
    profile = []
    # initialize profile to ones
    for ndx in range(0, len(alphabet)):
        profile.append(list(itertools.repeat(1, kLen)))
    # tally the number of characters at each position
    for col in range(0, kLen):
        for row in motifs:
            nTide = row[col]
            profile[letterToNum[nTide]][col] += 1
    # divide each tally by the number of rows in the motifs to get a percentage
    percentProfile = list(map(lambda x: list(map(lambda y: y / denom, x)), profile))
    return percentProfile

# Scores set of motifs based on the hamming distance from a consensus string
def scoreMotifs(motifs, profile, alphabet):
    # find the most popular string -> highest value in each column of profile
    rows = len(profile)
    cols = len(profile[0])
    consensus = []
    for col in range(0, cols):
        maxProb = 0
        popular = None
        for row in range(0, rows):
            if profile[row][col] > maxProb:
                maxProb = profile[row][col]
                popular = alphabet[row]
        consensus.append(popular)
    score = 0
    for row in motifs:
        score += utilities.hammingDistance(row, consensus)
    return score


# Score a set of motifs based on the probability of the consensus string
def scoreMotifsPercentile(motifs, profile, alphabet):
    # find the consensus string
    rows = len(profile)
    cols = len(profile[0])
    consensus = []
    score = 1.0
    for col in range(0, cols):
        maxProb = 0
        popular = None
        for row in range(0, rows):
            if profile[row][col] > maxProb:
                maxProb = profile[row][col]
                popular = alphabet[row]
        consensus.append(popular)
        score = score * profile[letterToNum[popular]][col]
    return score

def chooseRandomMotifs(dna, kLen):
    # get a random starting index into each string
    # We can start a maximum of kLen characters from the end
    maxNdx = len(dna[0]) - kLen
    motifs = []
    #indices = []
    # for each fragment in DNA
    for fragment in dna:
        # select a random kmer of length kLen
        startNdx = randrange(maxNdx)
        motifs.append(fragment[startNdx : startNdx+kLen])
        #indices.append(startNdx)
    #print(indices)
    return motifs

def distillMotifs(kLen, motifs, alphabet):
    bestMotifs = motifs
    bestScore = kLen
    while True:
        profile = buildProfile(kLen, motifs, alphabet)
        motifs = [findMostProbable(motif, kLen, profile) for motif in motifs]
        score = scoreMotifs(motifs, profile, alphabet)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
        else:
            return bestScore, bestMotifs
             

def run(iterations, alphabet):
    with open('randomMotifSearch.txt', 'r') as infile:
        intParams = utilities.readIntListFromFile(infile)
        kLen = intParams[0]
        numFrags = intParams[1]
        dna = []
        for line in infile:
            dna.append(list(line.strip()))
    #print(dna)

    # form a motif from random kmers in each fragment
    bestMotifs = chooseRandomMotifs(dna, kLen)
    initialProfile = buildProfile(kLen, bestMotifs, alphabet)
    #bestScore = scoreMotifsPercentile(bestMotifs, initialProfile, alphabet)
    bestScore = scoreMotifs(bestMotifs, initialProfile, alphabet)
    #bestScore = kLen

    #print(bestScore)

    # pick random kmers from each fragment iterations times
    # Then distill the kmers to their 'best' motifs
    # if the 'best' motifs beat the best score, use them as our new baseline
    for ndx in range(0, iterations):
        seedMotifs = chooseRandomMotifs(dna, kLen)
        localBestScore, localBestMotifs = distillMotifs(kLen, seedMotifs, alphabet)
        if localBestScore < bestScore:
            bestScore = localBestScore
            bestMotifs = localBestMotifs

    with open('randomMotifSearch.results.txt', 'w') as outfile:
        bestMotifStrings = [''.join(item) for item in bestMotifs]
        utilities.writeListToFileOnNewlines(outfile, bestMotifStrings)

# Run the program
run(10000, alphabet)
