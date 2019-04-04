import utilities

# build dictionaries containing the prefix and suffix values associated with
# each middle substring
def constructFixes(frags, prefixes, suffixes):
    kLen = len(frags[0])
    numFrags = len(frags)
    for ndx in range(0, numFrags):
        frag = frags[ndx]
        utilities.appendToDict(prefixes, ''.join(frag[1:]), frag[0])
        utilities.appendToDict(suffixes, ''.join(frag[0:kLen-1]), frag[kLen - 1])

# using the dictionaries, generate a list of pairs representing the graph
def constructGraph(prefixes, suffixes):
    graph = []
    for key in prefixes:
        prefs = prefixes[key]
        if key in suffixes.keys():
            sufs = suffixes[key]
            for foo in prefs:
                for bar in sufs:
                    pair = (foo + key, key + bar)
                    graph.append(pair)
    return graph

# print the graph in rosalind format
def printGraph(graph):
    for pair in graph:
        print(pair[0], ' -> ', pair[1])

def run():
    with open('overlapGraph.txt', 'r') as infile:
        frags = []
        strFrags = []
        for line in infile:
            frags.append(list(line.strip()))
            strFrags.append(line.strip())
    #print(frags)
    #print(strFrags)
    prefixes = {}
    suffixes = {}
    constructFixes(frags, prefixes, suffixes)
    print(prefixes)
    print(suffixes)
    graph = constructGraph(prefixes, suffixes)
    graph = sorted(graph, key=lambda pair: pair[0])
    printGraph(graph) 
    with open('overlapGraph.results.txt', 'w') as outfile:
        for pair in graph:
            outfile.write(pair[0] + ' -> ' + pair[1] + '\n')

run()
