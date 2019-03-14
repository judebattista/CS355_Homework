import itertools
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
    motifRows = len(motifs)
    profile = []
    # initialize profile to zeros
    for ndx in range(0, len(alphabet)):
        profile.append(list(itertools.repeat(0, kLen)))
    # tally the number of characters at each position
    for col in range(0, kLen):
        for row in motifs:
            nTide = row[col]
            profile[letterToNum[nTide]][col] += 1
    # divide each tally by the number of rows in the motifs to get a percentage
    percentProfile = list(map(lambda x: list(map(lambda y: y / motifRows, x)), profile))
    #for row in profile:
    #    for col in row:
    #        col /= float(motifRows)
    return percentProfile

def buildProfile1(kLen, motif, alphabet):
    return False

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

with open('greedyMotifSearch.txt', 'r') as infile:
    intParams = utilities.readIntListFromFile(infile)
    kLen = intParams[0]
    numFrags = intParams[1]
    dna = []
    for line in infile:
        dna.append(list(line.strip()))
print(dna)

# form a motif from the first kmers in each fragment
bestMotifs = []
for ndx in range(0, numFrags):
    bestMotifs.append(dna[ndx][0:kLen])
print('Best Motifs: ', list(map(lambda x: ''.join(x), bestMotifs)))
profile = buildProfile(kLen, bestMotifs, alphabet)
bestScore = scoreMotifs(bestMotifs, profile, alphabet)
print(bestScore)

#build a new set of motifs based around every kmer in the first fragment
for foo in range(0,len(dna[0]) - kLen + 1):
    motifs = []
    # start the motifs with a kmer from the first dna fragment
    motifs.append(dna[0][foo:foo+kLen])
    # for every fragment after the first one
    for bar in range(1,numFrags): 
        # build a profile using the motif
        profile = buildProfile(kLen, motifs, alphabet)
        nextMotif = findMostProbable(dna[bar], kLen, profile)
        motifs.append(nextMotif)
    print('Profile: ', profile)
    print('Motifs: ', list(map(lambda x: ''.join(x), motifs)))
    #newScore = scoreMotifsPercentile(motifs, profile, alphabet) 
    newScore = scoreMotifs(motifs, profile, alphabet) 
    print('New score: ', newScore)
    # hamming distance score needs to be less than the best
    if newScore < bestScore:
    # probability score needs to be greater than the best
    #if newScore > bestScore:
        bestScore = newScore
        bestMotifs = motifs
    print('Best Motifs: ', list(map(lambda x: ''.join(x), bestMotifs)))
    print('Best score: ', bestScore)

with open('greedyMotifSearch.results.txt', 'w') as outfile:
    bestMotifStrings = [''.join(item) for item in bestMotifs]
    utilities.writeListToFileOnNewlines(outfile, bestMotifStrings)
