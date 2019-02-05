a = 9
if a < 10:
    print('The number is less than 10')          
    print('Also, do another thing in the if block')
else: 
    print('The number is greater than or equal to 10')
    print('And do something else in the else block')
print('Always printed')

greetings = 1
while greetings <= 3:
    print('Hello! ' * greetings)
    greetings += 1
print('While loop terminated')

names = ['Alpha', 'Bravo', 'Charlie']
for name in names:
    print('Hello, ' + name)
print('for loop terminated')

foo = 10
for bar in range(foo):
    print(bar)
print('Numeric range for loop terminated')

#interestingly, these just output the string 'range(5, 12)'
print(range(5, 12))
rangeOfInts = range(5, 12)
print(rangeOfInts)
#I assume in python 2, this behaved differently. Check with Dr. Jones or Alyssa.

print('Come with me if you want to live.')

'''
Assumptions: 
    a < b < 10000
    sum all the odd integers on [a,b]
'''
'''
Need to handle four cases:
    O - O
    O - E
    E - E
    E - O
'''

a = 4141 
b = 8814

# this fails on O - O
sum = 0
start = a // 2
end = b // 2
for foo in range(start, end):
    sum += foo * 2 + 1
print(sum)

# this also fails on O - O
delta = b - a
kRange = (delta + 1) // 2 #effectively takes ceiling(delta / 2). I hope.
sum = 0
for foo in range(start, start + kRange):
    sum += foo * 2 + 1
print(sum)

# Pretty sure this works on all four cases
sum = 0
start = a // 2          #floor of a / 2
end = (b + 1) // 2      #ceiling of b / 2
delta = end - start
sum = 0
for foo in range(start, end):
    sum += foo * 2 + 1
print (sum)

