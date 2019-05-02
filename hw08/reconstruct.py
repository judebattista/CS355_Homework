import utilities

with open('reconstruct.txt', 'r') as infile:
    klen = int(infile.readline())
    kmers = []
    for line in infile:
        kmers.append(line.strip())

prefixes = {}
suffixes = {}
kDict = {}

for kmer in kmers:
    utilities.appendToDict(prefixes, kmer[1:klen], kmer[0])
    utilities.appendToDict(suffixes, kmer[0:klen-1], kmer[klen-1])
    utilities.incrementDict(kDict, kmer)
print(suffixes)
print(prefixes)

# select a kmer
# peel off its suffix to create a base.
# look up its base in suffixes to grab the trailing character of a matching kmer.
# create the matching kmer by appending the character to the end of base.
# create a new base and repeat
# This should resemble an Eularian path solution. Guess I'd better get that done.