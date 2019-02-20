# Calculate the Hamming distance between sequences 0 and 1
# seq0: string of length n
# seq1: string of length n
# note: will still work with strings of uneven length but the distance will be
#   inflated by the difference in lengths
# returns an integer representing the Hamming distance between the strings
def hamDist(seq0, seq1):
    return sum([ntide0 != ntide1 for ntide0, ntide1 in zip(seq0, seq1)])

# Appends value to a list at key in dictionary
def appendToDict(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]

with open('approximateOccurrences.txt', 'r') as infile:
    pattern = infile.readline().strip()
    text = infile.readline().strip()
    maxDist = int(infile.readline().strip())
pLen = len(pattern)
tLen = len(text)
#print(str(hamDist('abcde', 'abcdf')))
#print(str(hamDist('abcde', 'abcde')))
# Create a dictionary that uses the Hamming index between pattern and every equal length
#   substring of text as keys and a list of the indices of the corresponding substrings
#   as the value for the entry
distances = {}
for ndx in range(0, tLen - pLen + 1):
    #print(pattern)
    #print(text[ndx:ndx+pLen])
    appendToDict(distances, (hamDist(pattern, text[ndx:ndx+pLen])), ndx)

# Collapse the lists of indices at the valid keys into a single list
indices = []
[indices.extend(value) for (key, value) in distances.items() if key <= maxDist]
#print(indices)

with open('approximateOccurencesResults.txt', 'w') as outfile:
    outfile.write(' '.join(map(lambda num: str(num), sorted(indices))))
