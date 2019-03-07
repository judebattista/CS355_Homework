def appendToDict(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]

def incrementDict(dictionary, key):
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = 1

def writeListToFile(outfile, data):
    outfile.write(' '.join(map(lambda value: str(value), data)))

def writeListToFileOnNewlines(outfile, data):
    for item in data:
        outfile.write(item.strip())
        outfile.write('\n')

def readIntListFromFile(infile):
    return list(map(lambda x: int(x), infile.readline().strip().split(' ')))

