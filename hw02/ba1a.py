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

sample = "This is a test. This is only a test"
kmerCounts = findKmers(sample, 4)
sortedKmerCounts = sortKmerCountsByCount(kmerCounts)
printSortedKmerCounts(sortedKmerCounts)
