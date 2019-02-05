list_name = [0, 1, 2, 'three']
print(list_name[2:4])
print(list_name[2:])
print(list_name[:2])

quote = 'I swear by my pretty floral bonnet, I will end you.';
print(quote[:14])
print(quote[14:])

s = 'xgK1Ci2L5u6UZ5bPelusiosoW8qgqwc3Bu89vEWWZlhPa8c3NAdmNV5gmelanuroidesIxQQJxCO68KlXahLrQ2hMRn5yUbj4szGrQ7K5mpCMErJxM5tAbYgg1JZcWRSZvfJkd1l9QzOE7NXW4mNg6b2LMGKkfylHwev2.' 
a = 15
b = 22
c = 56
d = 67
answer = s[a : b+1] + ' ' + s[c : d+1]
print(answer)
