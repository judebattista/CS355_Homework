# function to find the reverse complement of a DNA strand
# pattern: A string representation of DNA consisting of the chars 'A', 'C, 'G', and 'T'
# complement: The pattern's complement, complete with 5' to 3' reversal
# We may need to store complement in pattern in order to comply with the book's instructions
# However, this makes for much less readable code, so we'll stick with complement as our return value for now
def revComp(pattern):
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    comp = []
    for ntide in pattern[len(pattern) - 1 : : -1]:
        comp.append(complements[ntide])
    complement = ''.join(comp)
    return complement

# function to print the pattern, its reversed string, and its reversed complement
# pattern: A string representation of DNA consisting of the chars 'A', 'C, 'G', and 'T'
def printComplementDetails(pattern):
    print('Source:      ' + source)
    print('Rev source:  ' + source[len(source)-1::-1])
    print('Complement:  ' + revComp(source))

source = "ACCGGGTTTT"
printComplementDetails(source)    
