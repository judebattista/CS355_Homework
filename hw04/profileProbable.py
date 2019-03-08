import utilities

letterToNum = {'A':0, 'C':1, 'G': 2, 'T':3}

def fragmentProbability(fragment, matrix):
    kLen = len(fragment)
    probability = 0
    for foo in range(0, kLen):
        nTide = letterToNum[fragment[foo]]
        probability += matrix[nTide][foo]
    return probability

def findMostProbable(text, kLen, matrix):
    fragments = [dna[foo:foo+kLen] for foo in range(0, len(dna) - kLen + 1)]
    probabilities = {''.join(fragment):fragmentProbability(fragment, matrix) for fragment in fragments}
    return max(probabilities, key=probabilities.get)


with open('profileProbable.txt', 'r') as infile:
    dna = list(infile.readline().strip())
    kLen = int(infile.readline().strip())
    matrix = []
    for line in infile:
        matrix.append(list(map(lambda num: float(num), line.strip().split())))

mostProbable = findMostProbable(dna, kLen, matrix)

with open('profileProbable.results.txt', 'w') as outfile:
    outfile.write(mostProbable)
