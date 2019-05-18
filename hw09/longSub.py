def printMatrix(matrix):
    for row in matrix:
        print(row)

def evalNode(row, col, downVals, rightVals):
    upCost = -1
    leftCost = -1
    if row >  0:
        upCost = downVals[row - 1][col]
    if col > 0:
        leftCost = rightVals[row][col - 1]
    return upCost, leftCost

#start at node row, col
def buildPaths(row, col, downVals, rightVals, totals):
    total = totals[row][col]
    upCost, leftCost = evalNode(row, col, downVals, rightVals)
    if row > 0:
        nextUp = row - 1
        nextUpTotal = totals[nextUp][col]
        proposedUp = totals[row][col] + upCost
        if proposedUp > nextUpTotal:
            totals[nextUp][col] = proposedUp
        totals = buildPaths(nextUp, col, downVals, rightVals, totals)
    if col > 0:
        nextLeft = col - 1
        nextLeftTotal = totals[row][nextLeft]
        proposedLeft = totals[row][col] + leftCost
        if proposedLeft > nextLeftTotal:
            totals[row][nextLeft] = proposedLeft
        totals = buildPaths(row, nextLeft, downVals, rightVals, totals)
    return totals

with open('longPath.txt', 'r') as infile:
    dims = infile.readline().strip().split()
    rows = int(dims[0])
    cols = int(dims[1])
    downVals = []
    rightVals = []
    for ndx in range(0, rows):
        row = infile.readline().strip()
        row = row.split()
        row = map(int, row)
        downVals.append(row)
    # eat the - separator between matrices
    infile.readline()
    for ndx in range(0, rows+1):
        row = infile.readline().strip()
        row = row.split()
        row = map(int, row)
        rightVals.append(row)

printMatrix(downVals)
printMatrix(rightVals)

totals = [[0 for col in range(cols+1)] for row in range(rows+1)]
printMatrix(totals)

totals = buildPaths(rows, cols, downVals, rightVals, totals)
printMatrix(totals)

with open('longPath.results.txt', 'w') as outfile:
    outfile.write(str(totals[0][0]))


