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

# So, I misread the problem, turns out we need to keep the indices of matches, not just count them. 
# Find the indices where the pattern appears in text
# Let's do it with thoroughly unnecessary recursion. Because I've already committed to this bit back when I thought it would be a simple counting gig. Yay?
# text: a string to search for the pattern
# pattern: the pattern to search for
# Well, that got uglier than I had hoped, need to include some state information
# indices: the list of indices where the pattern is found
# offset: the offset into the total text where the current search takes place
def PatternIndices(text, pattern, indices, offset):
    ndx = text.find(pattern)
    if ndx > -1:
        indices.append(ndx + offset)
        indices = PatternIndices(text[ndx+1:], pattern, indices, offset+ndx+1)
    return indices


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
    pattern = filein.readline().strip()
    #print('Searching for ' + pattern)
    text = filein.readline().strip()
    #print('Searching in ' + text)
    indices = []
    indices = PatternIndices(text, pattern, indices, 0)
    count = PatternCount(text, pattern)
    print('Counted: ' + str(count))
    print(indices)
    print('Found ' + str(len(indices)))
with open('patternCountResults.txt', 'w') as fileout:
    for ndx in indices:
        fileout.write(str(ndx) + ' ')
