infile = open('rosalind_ini5.txt', 'r')
#for line in infile:
#   print(line)
outfile = open('output.txt', 'w')
for line in infile:
    nextline = infile.readline()
    outfile.write(nextline)
    print('Read: ', line)
    print('Wrote: ', nextline)

infile.close()
outfile.close()

