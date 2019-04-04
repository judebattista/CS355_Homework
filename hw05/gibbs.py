import utilities
from random import randrange, seed

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

def run():
    with open('gibbs.txt', 'r') as infile:
        intParams = utilities.readIntListFromFile(infile)
        kLen = intParams[0]
        numFrags = intParams[1]
        n = intParams[2]
        dna = []
        for line in infile:
            dna.append(list(line.strip()))

        # randomly select a kmer from each string
        seedMotifs = chooseRandomMotifs(dna, kLen) 

run()
