import itertools
from random import randrange, seed, uniform
import time
import utilities


letterToNum = {'A':0, 'C':1, 'G': 2, 'T':3}
alphabet = ['A','C','G','T']

def scoreKmerWithProfile(kmer, profile, alphabet):
    score = 0
    for ndx in range(0, len(kmer)):
        letter = kmer[ndx]
        letterRow = letterToNum[letter]
        letterScore = profile[letterRow][ndx]
        score += letterScore
    return score

def findMostProbable(testKmers, profile, alphabet):
    bestKmer = ''
    bestScore = -1
    for kmer in testKmers:
        score = scoreKmerWithProfile(kmer, profile, alphabet)
        if score > bestScore:
            bestScore = score
            bestKmer = kmer
    return bestKmer


def pickKmer(testKmers, profile, alphabet):
    # running total of the score
    totalScore = 0
    kmerCount = len(testKmers)
    # we are interested in tracking the total score, not individual scores
    runningScores = []
    for kmer in testKmers:
        totalScore += scoreKmerWithProfile(kmer, profile, alphabet)
        runningScores.append(totalScore)
    #print(runningScores)
    roll = uniform(0, runningScores[kmerCount - 1])
    #print(roll)
    kmerScore = runningScores[0]
    ndx = 0
    while kmerScore < roll:
        ndx += 1
        kmerScore = runningScores[ndx]
    #print(ndx, testKmers)
    return testKmers[ndx - 1]

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

# Score motifs based on difference from most common string
# low scores are good.
# make counts column-major
def scoreMotifs(motifs, alphabet):
    numChars = len(alphabet)
    numMotifs = len(motifs)
    lenMotif = len(motifs[0])
    counts = []
    # set counts to zero
    for ndx in range(0, lenMotif):
        counts.append([0]* numChars)
    # count each character in motifs
    for row in range(0, numMotifs):
        for col in range(0, lenMotif):
            nTide = motifs[row][col]
            countsRow = letterToNum[nTide]
            counts[col][countsRow] += 1
    #print(motifs)
    #print(counts)
    score = 0
    for col in range(0, lenMotif):
        colMax = max(counts[col])
        colScore = numMotifs - colMax
        score += colScore
    #for motif in motifs:
    #    print(motif)
    #print(score)
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

def gibbs(kLen, motifs, alphabet, iterations):
    numFrags = len(motifs)
    fragLen = len(motifs[0])
    bestMotifs = []
    
    #choose random kmers from each motif as a starting point
    seedMotifs = chooseRandomMotifs(motifs, kLen)
    newMotifs = seedMotifs
    bestScore = scoreMotifs(seedMotifs, alphabet)
    startingScore = bestScore
    
    # iterate over our set of seed motifs
    #if startingScore < int(numFrags / 2) * fragLen:
    if True:
        for ndx in range(0, iterations):
            ndxToReplace = randrange(numFrags)
            fragToReplace = motifs[ndxToReplace]

            testKmers = [fragToReplace[ndx : ndx+kLen] for ndx in range(0, fragLen - kLen + 1)]

            del newMotifs[ndxToReplace]
            profile = buildProfile(kLen, newMotifs, alphabet)
            newFrag = pickKmer(testKmers, profile, alphabet)
            newMotifs.insert(ndxToReplace, newFrag)
            score = scoreMotifs(newMotifs, alphabet)
            if score < bestScore:
                bestScore = score
                bestMotifs = newMotifs
    return bestMotifs, bestScore, startingScore
             

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
    
    bestScore = kLen * numFrags
    bestMotifs = []
    startingScoreForBest = bestScore 
    #print(bestScore)
    
    iterationsPerStartingPoint = 5 
    iterations = 50000
    #iterations *= 10
    #iterations *= 50
    startTime = time.clock()
    for ndx in range(0, iterations):
        #print('----Run %d----' % ndx)
        if ndx % 1000 == 0:
            seed()
        newMotifs, newScore, startingScore = gibbs(kLen, dna, alphabet, iterationsPerStartingPoint)
        if newScore < bestScore:
            bestScore = newScore
            bestMotifs = newMotifs
            startingScoreForBest = startingScore
    print('Best score: ', bestScore, '. Starting score: ', startingScoreForBest)
    endTime = time.clock()
    print('Time elapsed: ', endTime - startTime)
    with open('gibbs.results.txt', 'w') as outfile:
        bestMotifStrings = [''.join(item) for item in bestMotifs]
        utilities.writeListToFileOnNewlines(outfile, bestMotifStrings)
    

# Run the program
run(alphabet)
