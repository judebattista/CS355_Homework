# function to count the kmers in a string
# this version should run in |text| * k time
# text: the string to be searched
# k: the length of kmers to search for
# kmerCounts: a dictionary containing each kmer and the number of times it appears in the text
def findKmers(text, k):
    # Array to hold a kmer
    kmer = []
    # dictionary to hold all the kmers and their counts
    kmerCounts = {}
    for ndx in range(0, len(text)- k + 1):
        kmer = []
        kmer.append(text[ndx:ndx+k])
        #print(kmer)    # uncomment to debug
        stringKmer = ''.join(kmer)
        if stringKmer in kmerCounts:
            kmerCounts[stringKmer] += 1
        else:
            kmerCounts[stringKmer] = 1
    return kmerCounts

# TODO: Find a |text| time version
# Eli Min pointed me to this article on searching strings using tries on Medium.com
# https://medium.com/basecs/trying-to-understand-tries-3ec6bede0014
# Assuming a four-letter alphabet, this seems like a potentially ideal solution
# Conclusion: Works quite nicely for storing kmers!
#   Unfortunately, not well suited for sorting kmers by frequency
#   Probably better off using a dictionary
#   Fun experiment though!
class trieNode:
    def __init__(self):
        self.value = None
        self.count = 0
        self.children = {'A':None, 'C':None, 'G':None, 'T':None}

    def dump(self):
        print(self.count)
        for key in self.children:
            print(key)
            if self.children[key] is not None:
                self.children[key].dump()

#Naive approach to adding kmers to the trie
#Should really add the current letter to k tries at the same time
#TODO: add parallel walk
def addKmerToTrie(root, kmer):
    workingNode = root
    for ntide in kmer:
        #walk the trie, checking each letter of the kmer
        #if the letter doesn't exist, create it
        if workingNode.children[ntide] is None:
            workingNode.children[ntide] = trieNode()
        #if the letter does exist, move to it
        else:
            workingNode = workingNode.children[ntide]
    #After the for loop, we should be at the leaf so increment the count
    workingNode.count += 1

# Any way to do this with a map?
# Could also use recursion. More fun than a for loop at least.
def buildTrie(text, k):
    root = trieNode()
    numberOfKmers = len(text) - k + 1;
    for ndx in range(0, numberOfKmers):
        addKmerToTrie(root, text[ndx:ndx+k])
    return root
    
# Function to print the kmer counts prettily
# kmerCounts: dictionary with kmer strings as keys and integer counts as values
# Used W3Schools.com for dictionary functions
# https://www.w3schools.com/python/python_dictionaries.asp
def printKmerCounts(kmerCounts):
    for key, value in kmerCounts.items():
        print(key + ': ' + str(value))

# Used python.org's sorted documentation to figure out how to sort a dictionary by value
# https://docs.python.org/3/library/functions.html#sorted
# Doesn't work quite right: currently returns a list of key value pairs
# Could build a dictionary based off the list, but that seems horribly inefficient.
# Will leave it as a list of kvps unless a compelling reason to do otherwise crops up
def sortKmerCountsByCount(kmerCounts):
    kvps = kmerCounts.items()
    kvps = sorted(kvps, key=lambda kvp: kvp[1], reverse=True)
    return kvps

# Function to print the sorted kmer counts prettily
# kmerCounts: list of kvps with kmer strings as keys and integer counts as values
def printSortedKmerCounts(kmerCounts):
    for kvp in kmerCounts:
        print(kvp[0] + ': ' + str(kvp[1]))


def writeMostCommonKmers(kmerCounts, filename):
    targetCount = kmerCounts[0][1]
    ndx = 0
    with open(filename, 'w') as outfile:
        while kmerCounts[ndx][1] == targetCount:
            outfile.write(kmerCounts[ndx][0] + '\n')
            ndx += 1

with open('frequentWords.txt') as infile:
    text = infile.readline().strip()
    k = int(infile.readline().strip())
    kmerCounts = sortKmerCountsByCount(findKmers(text, k))
    # Currently no good way to sort results by frequency, so not using tries
    #root = buildTrie(text, k)
    #root.dump()

filename = 'frequentWordsResults.txt'
writeMostCommonKmers(kmerCounts, filename)
