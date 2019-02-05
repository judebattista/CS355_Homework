import os

#open the data file
with open('rosalind_ini6.txt', 'r') as infile:
    #get the first line out of the file and remove whitespace
    s = infile.readline()
s = s.strip()
print(s)
#create a new dictionary
d = {}
#create a list of all the words in the sentence
wordList = s.split(' ')
#for every word in the sentence, make it a key in the dictionary with an initial value of 1
#if the key already exists, add one to it.
for word in wordList:
    d[word] = d.get(word, 0) + 1
print(d)

#whoops. Didn't need to output it as a dictionary, just as KVPs
'''
with open('exercise06.txt', 'w') as outfile: 
    outfile.write('{')
    for key in d:
        outfile.write(key + ': ' + str(d[key]) + ', ')
    #remove the trailing ', '
    # the 2 says seek from the end of the file
    # we can't use a non-zero seek from the end unless we open the file in binary mode
    # so let's go to the end of the file
    outfile.seek(0, os.SEEK_END) #os.SEEK_END = 2
    # ... then go backwards from our current position by using a whence value of 1 instead of 2
    outfile.seek(outfile.tell() -2, os.SEEK_SET) #os.SEEK_SET = 1
    outfile.truncate()
    outfile.write('}')
'''
with open('exercise06.txt', 'w') as outfile: 
    for key in d:
        outfile.write(key + ' ' + str(d[key]) + '\n')
