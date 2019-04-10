import itertools
from random import randrange, seed, uniform
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
    # we are interested in tracking the total score, no individual scores
    runningScores = []
    for kmer in testKmers:
        totalScore += scoreKmerWithProfile(kmer, profile, alphabet)
        runningScores.append(totalScore)
    print(runningScores)
    roll = uniform(0, runningScores[kmerCount - 1])
    print(roll)
    kmerScore = runningScores[0]
    ndx = 0
    while kmerScore < roll:
        ndx += 1
        kmerScore = runningScores[ndx]
    print(ndx, testKmers)
    return testKmers[ndx]

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
def scoreMotifsiWithProfile(motifs, profile, alphabet):
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
# make counts column-major
def scoreMotifs(motifs, alphabet):
    numChars = len(alphabet)
    numMotifs = len(motifs)
    lenMotif = len(motifs[0])
    counts = []
    # set counts to zero
    for ndx in range(0, lenMotif):
        counts.append([0]* numCharsOB)
    # count each character in motifs
    for row in range(0, numMotifs)
        for col in range(0, lenMotif):
            nTide = motifs[row][col]
            countsRow = letterToNum[nTide]
            counts[col][countsRow] += 1
    print(motifs)
    print(counts)
    score = 0
    for col in range(0, numMotifs):
        colMax = max(count[col])
        colScore = numMotifs - colMax
        score += colScore



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
    ndxToReplace = randrange(numFrags)
    fragToReplace = motifs[ndxToReplace]

    #choose random kmers from each motif
    seedMotifs = chooseRandomMotifs(motifs, kLen)
    newMotifs = seedMotifs
    testKmers = [fragToReplace[ndx : ndx+kLen] for ndx in range(0, len(fragToReplace) - kLen + 1)]

    del newMotifs[ndxToReplace]
    profile = buildProfile(kLen, newMotifs, alphabet)
    newFrag = pickKmer(testKmers, profile, alphabet)
    newMotifs.insert(ndxToReplace, newFrag)
    scoreProfile = buildProfile(kLen, newMotifs, alphabet)
    score = scoreMotifs(newMotifs, scoreProfile, alphabet)
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
    bestScore = kLen * numFrags
    bestMotifs = []
    #print(bestScore)
    
    iterations = 2000
    # pick random kmers from each fragment iterations times
    # Then distill the kmers to their 'best' motifs
    # if the 'best' motifs beat the best score, use them as our new baseline
    for ndx in range(0, iterations):
        if ndx % 1000 == 0:
            seed()
        newMotifs, newScore = gibbs(kLen, dna, alphabet)
        if newScore < bestScore:
            bestScore = newScore
            bestMotifs = newMotifs

    with open('gibbs.results.txt', 'w') as outfile:
        bestMotifStrings = [''.join(item) for item in bestMotifs]
        utilities.writeListToFileOnNewlines(outfile, bestMotifStrings)

# Run the program
run(alphabet)
