import utilities
from itertools import combinations, product, chain

# After talking to Dr. Jones and investigating itertools, I opted for the 
# Hamming ball/circle approach as explained by Dr. Jones and sourced from:
# https://codereview.stackexchange.com/questions/88912/create-a-list-of-all-strings-within-hamming-distance-of-a-reference-string-with
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

with open('dNbhd.txt', 'r') as infile:
    root = infile.readline().strip()
    distance = int(infile.readline().strip())
    alphabet = 'ACGT'

cousins = chain.from_iterable(genCousinsOfD(root, alphabet, dist) for dist in range(distance + 1))

with open('dNbhd.results.txt', 'w') as outfile:
    utilities.writeListToFileOnNewlines(outfile, cousins)
     
