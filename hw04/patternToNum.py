import utilities
from functools import reduce

class radixAccum:
    translate = {'A':0, 'C':1, 'G':2, 'T':3}
    
    def __init__(self, radix):
        self.radix = radix
        self.columnValue = 1
    
    def calculate(self, letter):
        number = self.translate[letter] * self.columnValue
        self.columnValue *= self.radix
        return int(number)                

def patternToNum(pattern, accum):
    return reduce(lambda x, y: int(x) + accum.calculate(y), pattern)

with open('patternToNum.txt', 'r') as infile:
    pattern = '0' + infile.readline().strip().upper()[::-1]
    accum = radixAccum(4)
    value = patternToNum(pattern, accum)
    print(value)

with open('patternToNum.results.txt','w') as outfile:
    outfile.write(str(value))
