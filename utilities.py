def appendToDict(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]

def incrementDict(dictionary, key):
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = 1

def writeListToFile(outfile, data):
    outfile.write(' '.join(map(lambda value: str(value), data)))

def writeListToFileOnNewlines(outfile, data):
    for item in data:
        outfile.write(item.strip())
        outfile.write('\n')

def readIntListFromFile(infile):
    return list(map(lambda x: int(x), infile.readline().strip().split(' ')))

# compute the Hamming distance between two strings (or lists. Lists are better.)
def hammingDistance(text0, text1):
    return sum([ntide0 != ntide1 for ntide0, ntide1 in zip(text0, text1)])

def genCousinsOfD(kmer, alphabet, distance):
    # Generate all the possible locations for replacement
    for positions in combinations(range(len(kmer)), distance):
        #generate a every possible set of values to fill in each set of locations
        # ... kind of. We will not include the last letter in the alphabet
        for replacements in product(range(len(alphabet) - 1), repeat=distance):
            # work with a list rather than a string for efficiency
            cousin = list(kmer)
            for pos, rep in zip(positions, replacements):
                # if the letter is the same as the replacement, replace it with the 
                # 'spare' letter from the alphabet
                if cousin[pos] == alphabet[rep]:
                    cousin[pos] = alphabet[-1]
                else:
                    cousin[pos] = alphabet[rep]
            yield ''.join(cousin)


