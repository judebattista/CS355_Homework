from itertools import product
import utilities

# One way to generate all possible kmers of length length
def allKmers(length):
    # Generate all possible strings of length using the character A, T, C, and G
    return map(''.join, product('ATCG', repeat=length))

# Generate a dictionary of all kmers of length length that occur in text
def allKmersInString(text, length, defaultValue):
    kmers = {}
    for ndx in range(0, len(text) - patLen + 1):
        kmer = text[ndx: ndx + patLen]
        kmers[kmer] = defaultValue

# generate a list of all kmers of length length that occur in text
def listAllKmersInString(text, length):
    kmers = []
    for ndx in range(0, len(text) - length + 1):
        kmers.append(text[ndx: ndx + length]) 
    return kmers

# compute the Hamming distance between two strings
def hammingDistance(text0, text1):
    return sum([ntide0 != ntide1 for ntide0, ntide1 in zip(text0, text1)])

# For a given pattern, find out how often it matches a substring of text within 
# a Hamming distance of threshold
def generateFreqs(dictionary, text, pattern, threshold):
    patLen = len(pattern)
    for ndx in range(0, len(text) - patLen + 1):
        testKmer = text[ndx: ndx+patLen]
        if hammingDistance(testKmer, pattern) <= threshold:
            utilities.incrementDict(dictionary, pattern)
    return dictionary

with open('frequentWithMismatches.txt', 'r') as infile:
    text = infile.readline().strip()
    kmerLength, maxDistance = map(lambda x: int(x), infile.readline().strip().split(' '))
    kmers = allKmers(kmerLength)
    freqs = {}
    for kmer in kmers:
        generateFreqs(freqs, text, kmer, maxDistance)
    maxFreq = max(freqs.values())
    mostFrequentKmers = []
    for key, value in freqs.items():
        if value == maxFreq:
            mostFrequentKmers.append(key)
with open('frequentWithMismatchesResults.txt', 'w') as outfile:
    utilities.writeListToFile(outfile, mostFrequentKmers)
