letterToNum = {'A':0, 'C':1, 'G':2, 'T':3}

def scoreMotifs(motifs, alphabet):
     numChars = len(alphabet)
     numMotifs = len(motifs)
     lenMotif = len(motifs[0])
     counts = []
     # set counts to zero
     for ndx in range(0, lenMotif):
         counts.append([0]* numChars)
     # count each character in motifs
     for row in range(0, numMotifs):
         for col in range(0, lenMotif):
             nTide = motifs[row][col]
             countsRow = letterToNum[nTide]
             counts[col][countsRow] += 1
     #print(motifs)
     #print(counts)
     score = 0
     for col in range(0, lenMotif):
         colMax = max(counts[col])
         colScore = numMotifs - colMax
         score += colScore
     #for motif in motifs:
     #    print(motif)
     #print(score)
     return score

with open('test.txt') as infile:
    frags = []
    for line in infile:
        frags.append(line.strip())

score = scoreMotifs(frags, ['A','C','G','T'])
print(score)
