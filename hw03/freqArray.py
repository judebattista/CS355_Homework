from itertools import product
import utilities


def allKmers(length):
    # Generate all possible strings of length using the character A, T, C, and G
    return map(''.join, product('ACGT', repeat=length))

def generateDictFromKeys(listOfKeys, defaultValue):
    dictionary = {}
    for key in listOfKeys:
        dictionary[key] = defaultValue
    return dictionary

# Find the frequency of every substring of length length in text
def generateFreqs(dictionary, text, length):
    for ndx in range(0, len(text) - length + 1):
        utilities.incrementDict(dictionary, text[ndx: ndx+length])
    return dictionary

with open('freqArray.txt', 'r') as infile:
    text = infile.readline().strip()
    length = utilities.readIntListFromFile(infile)[0]
    freqs = generateDictFromKeys(allKmers(length), 0)
    generateFreqs(freqs, text, length)
    frequencies = list(freqs.values())

with open('freqArrayResults.txt', 'w') as outfile:
    utilities.writeListToFile(outfile, frequencies)
