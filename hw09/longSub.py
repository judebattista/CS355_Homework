import sys

def printMatrix(matrix):
    for row in matrix:
        print(row)

def evalNode(row, col, rowStr, colStr):
    upCost = -1
    leftCost = -1
    diagCost = -1
    if row >  0:
        upCost = 0 
    if col > 0:
        leftCost = 0
    if col > 0 and row > 0:
        if rowStr[row-1] == colStr[col-1]:
            diagCost = 1
            #print('comparing {0} with {1} to get {2}'.format(rowStr[row-1], colStr[col-1], diagCost))
        else:
            diagCost == 0
    return upCost, leftCost, diagCost

#start at node row, col
def buildPaths(row, col, totals, rowStr, colStr):
    upCost, leftCost, diagCost = evalNode(row, col, rowStr, colStr)
    if col > 0 and row > 0:
        proposedDiag = totals[row][col] + diagCost
        #print(proposedDiag)
        if totals[row-1][col-1] < proposedDiag:
            totals[row-1][col-1] = proposedDiag
        totals = buildPaths(row-1, col-1, totals, rowStr, colStr)
    if row > 0:
        if totals[row-1][col] < totals[row][col]:
            totals[row-1][col] = totals[row][col]
        totals = buildPaths(row-1, col, totals, rowStr, colStr)
    if col > 0:
        if totals[row][col-1] < totals[row][col]:
            totals[row][col-1] = totals[row][col]
        totals = buildPaths(row, col-1, totals, rowStr, colStr)
    return totals

def reconstructPath00(path, row, col, matrix, rowStr, colStr):
    #if rowStr[row] == colStr[col]:
    #    path.append(rowStr[row])
    if row <= 0 and col <= 0:
        print('Reached origin')
        return path
    nbhd = []
    if row > 0 and col > 0:
        # add diagonal
        nbhd.append((matrix[row-1][col-1], row-1, col-1))
    if row > 0:
        nbhd.append((matrix[row-1][col], row-1, col))
        # add up
    if col > 0:
        nbhd.append((matrix[row][col-1], row, col-1))
        # add left
    target = max(nbhd)
    print(target)
    newRow = target[1]
    newCol = target[2]

    return reconstructPath(path, newRow, newCol, matrix, rowStr, colStr)


def reconstructPath01(path, row, col, matrix, rowStr, colStr):
    if row <= 0 or col <= 0:
        return path 
    nbhd = []
    nbhd.append((matrix[row-1][col-1], row-1, col-1))
    nbhd.append((matrix[row][col-1], row, col-1))
    nbhd.append((matrix[row-1][col], row-1, col))
    
    if matrix[row-1][col-1] == matrix[row][col] + 1:
        #path.append(rowStr[row])
        pass 
    target = max(nbhd)
    newRow = target[1]
    newCol = target[2]
    path.append((newRow, newCol))
    return reconstructPath(path, newRow, newCol, matrix, rowStr, colStr)

def buildBtm(rowStr, colStr):
    rows = len(rowStr)
    cols = len(colStr)
    path = [[0 for col in range(cols)] for row in range(rows)]
    btm = [[0 for col in range(cols)] for row in range(rows)]
    for row in range(1, rows):
        for col in range(1, cols):
            up = (path[row-1][col], 'down')
            left = (path[row][col-1], 'right')
            if rowStr[row] == colStr[col]:
                diag = (path[row-1][col-1] +1, 'diag')
            else:
                diag = (-1, 'diag')
            target = max([up, left, diag])
            path[row][col] = target[0]
            btm[row][col] = target[1]
    return btm

def buildLcs(btm, rowStr, i, j, lcs):
    print('row: {0}, col: {1}, val: {2}'.format(i, j, btm[i][j]))
    if i == 0 or j == 0:
        print('Returning LCS: {0}'.format(lcs))
        return lcs
    if btm[i][j] == 'down':
        buildLcs(btm, rowStr, i-1, j, lcs)
    if btm[i][j] == 'right':
        buildLcs(btm, rowStr, i, j-1, lcs)
    if btm[i][j] == 'diag':
        print('Appending {0} to lcs: {1}'.format(rowStr[i], lcs))
        lcs.append(rowStr[i])
        buildLcs(btm, rowStr, i-1, j-1, lcs)

def buildStr(path, rowStr, colStr):
    pathStr = []
    for row, col in path:
        if rowStr[row] == colStr[col]:
            pathStr.append(rowStr[row])
    return pathStr
    
with open('longStr.txt', 'r') as infile:
    rowStr = list(infile.readline().strip())
    colStr = list(infile.readline().strip())
    print('{0} by {1}'.format(rowStr, colStr))
    rows = len(rowStr) - 1
    cols = len(colStr) - 1
    downVals = []
    rightVals = []

#totals = [[0 for col in range(cols+1)] for row in range(rows+1)]
#printMatrix(totals)

sys.setrecursionlimit(2000)
btm = buildBtm(rowStr, colStr)
printMatrix(btm)
lcs = []
buildLcs(btm, rowStr, rows, cols, lcs)
if rowStr[0] == colStr[0]:
    lcs.append(rowStr[0])
strLcs = ''.join(lcs[::-1])
print(strLcs)
with open('longStr.results.txt', 'w') as outfile:
    outfile.write(strLcs)


