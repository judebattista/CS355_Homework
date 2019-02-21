with open('hamming.txt', 'r') as infile:
    seq0 = infile.readline()
    seq1 = infile.readline()
    discrepancies = [ntide0 != ntide1 for ntide0, ntide1 in zip(seq0, seq1)]
    #print(discrepancies)
    hdist = sum(discrepancies)
    # So Python 3 got rid of reduce in favor of explicit accumulation loops
    #print(hdist)

with open('hammingResults.txt', 'w') as outfile:
    outfile.write(str(hdist))
