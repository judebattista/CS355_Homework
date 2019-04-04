import utilities

def minDistanceFromString(pattern, fragment):
    fragLen = len(fragment)
    patLen = len(pattern) 
    minDist = patLen
    for ndx in range(0, fragLen - patLen + 1):
        dist = utilities.hammingDistance(fragment[ndx : ndx + patLen], pattern)
        if dist < minDist:
            minDist = dist
    return minDist

def run():
    with open('distanceBetween.txt', 'r') as infile:
        pattern = list(infile.readline().strip())
        dna = []
        for line in infile:
            stringFrags = line.strip().split()
            listFrags = [list(frag) for frag in stringFrags]
            dna.extend(listFrags)

    #print(pattern)
    #print(dna)
    dist = 0
    for fragment in dna:
        dist += minDistanceFromString(pattern, fragment)
    print(dist)
    with open('distanceBetween.results.txt', 'w') as outfile:
        outfile.write(str(dist))

run()
