import utilities

# build dictionaries containing the prefix and suffix values associated with
# each middle substring
def constructEdges(kLen, text, edges):
    k = kLen - 1
    for ndx in range(0, len(text) - kLen + 1):
        key = text[ndx : ndx + kLen - 1]
        val = text[ndx + 1 : ndx + kLen]
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
    with open('stringDeBruijn.txt', 'r') as infile:
        kLen = int(infile.readline().strip())
        text = infile.readline().strip()
    print(kLen, text)
    edges = {}
    constructEdges(kLen, text, edges)
    print(edges)
    graph = constructGraph(edges)
    graph = sorted(graph, key=lambda pair: pair[0])
    printGraph(graph) 
    with open('stringDeBruijn.results.txt', 'w') as outfile:
        for pair in graph:
            valueString = ','.join(sorted(pair[1]))
            outfile.write(pair[0] + ' -> ' + valueString + '\n')

run()
