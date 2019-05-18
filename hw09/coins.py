with open('coins.txt', 'r') as infile:
    total = int(infile.readline().strip())
    coinStr = infile.readline().strip().split(',')
    coins = map(int, coinStr)

print(coins)
minCoins = [total] * (total+1)
minCoins[0] = 0
coinTypes = len(coins)

for foo in range(0, total+1):
    for bar in range(0, coinTypes):
        nextNdx = foo + coins[bar]
        newCount = minCoins[foo] + 1
        print('foo: {0}, nextNdx: {1}, newCount: {2}'.format(foo, nextNdx, newCount))
        if nextNdx < len(minCoins) and minCoins[nextNdx] > newCount:
            print('Current value: {0}, current count: {1}, next value: {2}, nexti count: {3}, proposed next count: {4}'.format(foo, minCoins[foo], nextNdx, minCoins[nextNdx], newCount))
            minCoins[nextNdx] = newCount

print(minCoins)
with open('coins.results.txt', 'w') as outfile:
    outfile.write(str(minCoins[total]))


