# function to find the reverse complement of a DNA strand
# pattern: A string representation of DNA consisting of the chars 'A', 'C, 'G', and 'T'
# complements: a dictionary translating a nucleotide to its complement
# revComp: The pattern's complement, complete with 5' to 3' reversal
# We may need to store complement in pattern in order to comply with the book's instructions
# However, this makes for much less readable code, so we'll stick with complement as our return value for now
def revComp(pattern, complements):
    comp = []
    for ntide in pattern[len(pattern) - 1 : : -1]:
        comp.append(complements[ntide])
    revComp = ''.join(comp)
    return revComp

def simpleReverseComp(pattern, complements):
    revComp = ''
    for ntide in reversed(pattern):
        revComp += complements[ntide]
    return revComp
        
def mapReverseComp(pattern, complements):
    mapComp = list(map(lambda ntide: complements[ntide], reversed(pattern)))
    revComp = ''.join(mapComp)
    return revComp


# function to print the pattern, its reversed string, and its reversed complement
# pattern: A string representation of DNA consisting of the chars 'A', 'C, 'G', and 'T'
# complements: a dictionary translating a nucleotide to its complement
def printComplementDetails(pattern, complements):
    print('Source:      ' + pattern)
    print('Rev source:  ' + pattern[len(pattern)-1::-1])
    print('Complement:  ' + revComp(pattern, complements))
    print('Simple Comp: ' + simpleReverseComp(pattern, complements))
    print('Map Comp:    ' + mapReverseComp(pattern, complements))

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
source = "ACCGGGTTTT"
printComplementDetails(source, complements)

with open('reverseComplement.txt', 'r') as infile:
    text = infile.readline().strip()
    #print(len(text))
    revComp = mapReverseComp(text, complements)
    #print(len(revComp))
with open('reverseComplementResults.txt', 'w') as outfile:
    outfile.write(revComp)

