import itertools
from random import randrange, seed
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

def findMostProbable(matrix, alphabet):
    rows = len(matrix)
    cols = len(matrix[0])
    newFrag = []
    for col in range(cols):
        maxColValue = 0
        maxRow = 0
        # in a given column, find out which row has the highest probability
        for row in range(rows):
            if matrix[row][col] > maxColValue:
                maxColValue = matrix[row][col]
                maxRow = row
        # when we get the highest probability, append the corresponding letter to our kmer
        newFrag.append(alphabet[maxRow])
    return newFrag

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

# Note: This scoring system is bad! Instead, just score the motifs like we scored
# a set of strings!
# However, this system does still produce the right results
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

def gibbs(kLen, motifs, alphabet):
    numFrags = len(motifs)
    fragToReplace = randrange(numFrags)
    newMotifs = motifs
    del newMotifs[fragToReplace]
    profile = buildProfile(kLen, newMotifs, alphabet)
    newFrag = findMostProbable(profile, alphabet)
    newMotifs.insert(fragToReplace, newFrag)
    score = scoreMotifs(newMotifs, profile, alphabet)
    return newMotifs, score
             

def run(alphabet):
    with open('gibbs.txt', 'r') as infile:
        intParams = utilities.readIntListFromFile(infile)
        kLen = intParams[0]
        numFrags = intParams[1]
        iterations = intParams[2]
        dna = []
        for line in infile:
            dna.append(list(line.strip()))
    #print(kLen, numFrags, iterations) 
    #print(dna)
    
    # form a motif from random kmers in each fragment
    #bestMotifs = chooseRandomMotifs(dna, kLen)
    #initialProfile = buildProfile(kLen, bestMotifs, alphabet)
    #bestScore = scoreMotifsPercentile(bestMotifs, initialProfile, alphabet)
    #bestScore = scoreMotifs(bestMotifs, initialProfile, alphabet)
    bestScore = kLen * numFrags
    bestMotifs = []
    #print(bestScore)

    # pick random kmers from each fragment iterations times
    # Then distill the kmers to their 'best' motifs
    # if the 'best' motifs beat the best score, use them as our new baseline
    for ndx in range(0, iterations):
        if ndx % 1000 == 0:
            seed()
        seedMotifs = chooseRandomMotifs(dna, kLen)
        newMotifs, newScore = gibbs(kLen, seedMotifs, alphabet)
        if newScore < bestScore:
            bestScore = newScore
            bestMotifs = newMotifs

    with open('gibbs.results.txt', 'w') as outfile:
        bestMotifStrings = [''.join(item) for item in bestMotifs]
        utilities.writeListToFileOnNewlines(outfile, bestMotifStrings)

# Run the program
run(alphabet)
