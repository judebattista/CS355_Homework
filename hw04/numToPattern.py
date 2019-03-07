translate = ['A', 'C', 'G', 'T']

def peelRemainder(value, radix, dna):
    if value > 0:
        dna = peelRemainder(int(value / radix), radix, translate[value % radix] + dna)
    return dna 
        

with open('numToPattern.txt', 'r') as infile:
    value = int(infile.readline().strip())
    length = int(infile.readline().strip())
    dna = peelRemainder(value, 4, '').rjust(length, translate[0])
    print(dna)

with open('numToPattern.results.txt', 'w') as outfile:
    outfile.write(dna) 
