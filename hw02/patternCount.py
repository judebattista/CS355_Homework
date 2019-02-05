# Count the number of times that pattern appears in text
# Let's do it with thoroughly unnecessary recursion. Because we can. Hooray!
# text: a string to search for the pattern
# pattern: the pattern to search for

def PatternCount(text, pattern):
    ndx = text.find(pattern)
    count = 0
    if ndx > -1:
        count += 1 + PatternCount(text[ndx+1:], pattern)
    return count

# make sure the test data has overlapping occurences of pattern
# otherwise we could just use string.count()

#Run vs test data sets
'''
text = 'GCGCG'
pattern = 'GCG'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)

text = 'AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT'
pattern = 'AAA'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)

text = 'AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'
pattern = 'TTT'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)

text = 'GGACTTACTGACGTACG'
pattern = 'ACT'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)


text = 'ATCCGATCCCATGCCCATG'
pattern = 'CC'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)

text = 'CTGTTTTTGATCCATGATATGTTATCTCTCCGTCATCAGAAGAACAGTGACGGATCGCCCTCTCTCTTGGTCAGGCGACCGTTTGCCATAATGCCCATGCTTTCCAGCCAGCTCTCAAACTCCGGTGACTCGCGCAGGTTGAGTA'
pattern = 'CTC'
print(pattern  + ' occurs ' + str(PatternCount(text, pattern)) + ' times in ' + text)
'''

count = 0
with open('patternCount.txt', 'r') as filein:
    text = filein.readline().strip()
    print(text)
    pattern = filein.readline().strip()
    print(pattern)
    count = PatternCount(text, pattern)
    print(count)
with open('patternCountResults.txt', 'w') as fileout:
    fileout.write(str(count))

