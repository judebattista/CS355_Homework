import utilities

# build dictionaries containing the prefix and suffix values associated with
# each middle substring
def constructEdges(kLen, frags, edges):
    for ndx in range(0, len(frags)):
        key = frags[ndx][0:kLen - 1]
        val = frags[ndx][1:kLen]
        utilities.appendToDict(edges, key, val)

def constructGraph(edges):
    graph = [(k, v) for k, v in edges.items()]
    return graph

# print the graph in rosalind format
def printGraph(graph):
    print(graph)
    for pair in graph:
        print(pair[0], ' -> ', ','.join(sorted(pair[1])))

def run():
    with open('kmerDeBruijn.txt', 'r') as infile:
        frags = []
        for line in infile:
            frags.append(infile.readline().strip())
    kLen = len(frags[0])
    print(kLen, frags)
    edges = {}
    constructEdges(kLen, frags, edges)
    print(edges)
    graph = constructGraph(edges)
    graph = sorted(graph, key=lambda pair: pair[0])
    printGraph(graph) 
    with open('kmerDeBruijn.results.txt', 'w') as outfile:
        for pair in graph:
            valueString = ','.join(sorted(pair[1]))
            outfile.write(pair[0] + ' -> ' + valueString + '\n')

run()
