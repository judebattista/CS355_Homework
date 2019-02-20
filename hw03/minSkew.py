class Accumulator:
    def __init__(self):
        self.value = 0
        self.skewDict = {}
    def inc(self, step = 1):
        self.value += step
    def dec(self, step = 1):
        self.value -= step
    #def addSkewToDict(self, ndx, skew):

    def tallySkew(self, ndx, char):
        char = char.upper()
        if char == 'C':
            self.dec()
        elif char == 'G':
            self.inc()
        if self.value in self.skewDict:
            self.skewDict[self.value].append(ndx)
        else:
            self.skewDict[self.value] = [ndx]

#sample = '0' + 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
with open('minSkew.txt', 'r') as infile:
    sample = '0' + infile.readline()
    accum = Accumulator()
    skews = [accum.tallySkew(ndx, char) for (ndx, char) in enumerate(sample)]
    indices = accum.skewDict[min(accum.skewDict.keys())]

with open('minSkewResults.txt', 'w') as outfile:
    outfile.write((' '.join(map(lambda number: str(number), indices))))

