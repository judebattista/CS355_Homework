kmers = {}

# adds a kmer occurrence to a dictionary
# kmer: the kmer to add
# index: the index where the kmer was found
# dictionary: uses kmers as keys, matching each with a list of indices at which the kmer has been found
def addKmer(kmer, index, dictionary):
    if kmer in dictionary:
        dictionary[kmer].append(index)
    else:
        dictionary[kmer] = [index]
    return dictionary

'''
#test for addKmer. Really need to get some TDD going...
addKmer('test', kmers, 0);
addKmer('test', kmers, 1);
print(kmers)
'''

# finds all kmers within text, storing them in kmerDictionary
#   using the kmer as the key and the list of indices at which it was found as the value
# text: the text to search
# k: the length of kmer we are searching for
# kmerDictionary: Storage with kmers as keys and the list of indices where the key can be found in text as the value 
def findKmers(text, k, kmerDictionary):
    # Array to hold a kmer
    kmer = []
    for ndx in range(0, len(text)- k + 1):
        kmer = []
        kmer.append(text[ndx:ndx+k])
        #print(kmer)    # uncomment to debug
        stringKmer = ''.join(kmer)
        addKmer(stringKmer, ndx, kmerDictionary)
    return kmerDictionary

# Checks to see if the values ever occur at least freq times in span distance of each other returns true if any clump is found
# values: a list of indices at which the item was found
# span: the distance in which the occurences must be found to count
# freq: how many times the item must occur within span to count
def clumped(values, span, freq):
    result = False;
    for ndx in range(0, len(values) - freq + 1):
        print('Checking ' + str(values[ndx + freq - 1]) + ' and ' + str(values[ndx]) + ' against ' + str(span) + ' ndx=' + str(ndx) + ' freq=' + str(freq))
        if values[ndx + freq - 1] - values[ndx] <= span:
            result = True;
            return result
    return result


# check kmerDictionary for any clumps, defined by freq occurrences within span distance
# kmerDictionary: Storage with kmers as keys and the list of indices where the key can be found in text as the value 
# span: the distance in which the occurences must be found to count
# freq: how many times the item must occur within span to count
def findClumps(kmerDictionary, span, freq):
    clumps = []
    for k, v in kmerDictionary.items():
        # only check kmers with at least freq occurrences
        # and yeah, I mostly typed that to make the pun
        if len(v) >= freq and clumped(v, span, freq):
            clumps.append(k)
    return clumps

with open('patternClumps.txt', 'r') as infile:
    data = infile.readline().strip()
    kmerLen, span, freq = infile.readline().strip().split(' ')
    #print(str(kmerLen) + ', ' + str(span) + ', ' + str(freq))
    kmers = {}
    kmers = findKmers(data, int(kmerLen), kmers)
    #print(kmers)
    # based on the failure of my first submission, I'm guessing all occurrences of the kmer must be completely contained within span
    #   not simply start within span
    # No problem, we'll just shorten the span by the k value
    clumps = findClumps(kmers, int(span) - int(kmerLen), int(freq))
    print(clumps)

with open('patternClumpsResult.txt', 'w') as outfile:
    #map(lambda x: outfile.write(x + ' '), clumps)
    # so that didn't work. Let's try using a join...
    outfile.write(' '.join(clumps))
