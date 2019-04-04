import utilities

def splitIntoKmers(fragment, kLen):
    fragLen = len(fragment)
    kmers = []
    for ndx in range(0, fragLen - kLen + 1):
        kmers.append(fragment[ndx : ndx + kLen])
    return sorted([''.join(kmer) for kmer in kmers])

def run():
    with open('kmerComposition.txt', 'r') as infile:
        kLen = int(infile.readline().strip())
        fragment = list(infile.readline().strip())

    #print(pattern)
    #print(dna)
    kmers = splitIntoKmers(fragment, kLen)
    print(kmers)
    with open('kmerComposition.results.txt', 'w') as outfile:
        utilities.writeListToFileOnNewlines(outfile, kmers)
run()
