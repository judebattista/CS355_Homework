import utilities

def run():
    with open('reconstructFromPath.txt', 'r') as infile:
        frag = []
        for line in infile:
            kmer = list(line.strip().split())
            frag.append(kmer[-1])

    #print(pattern)
    #print(dna)
    print(frag)
    with open('reconstructFromPath.results.txt', 'w') as outfile:
        outfile.write(''.join(frag))

run()
