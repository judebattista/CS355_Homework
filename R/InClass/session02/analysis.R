bacteria <- readLines("bacteria.csv")
shortNames <- bacteria[lapply(bacteria, nchar) < 10]
print(shortNames)

frag_counts = c(1000, 5000, 49000, NA, 70000, 500, NA)
frag_names = c('f1','f2','f3','f4','f5','f6','f7')
frag_names [frag_counts > 1000]
frag_names [frag_counts > 1000 & !is.na(frag_counts)]

print('d - density, returns the height of the distribution')
print('p - returns the cumulative density function')
print('q - return the inverse density function (quartiles)')
print('r - return randomly generated numbers')
print('Use these in conjunction with distribution names, e.g. dnorm')
print('distributions: t, binom, chisq')

for (i in -10:10) {
    # defaults to norm(0, 1)
    print(dnorm(i))
}

for (i in -10:10) {
    cat(i, ' ', dnorm(i, mean=4, sd=2), '\n')
}

# lower.tail is a binary switch for ltpn vs utpn
for (i in -10:10) {
    cat(i, ' ', pnorm(i, lower.tail=FALSE), '\n')
}

for (i in -10:10) {
    cat(i, ' ', pnorm(i, lower.tail=TRUE), '\n')
}

# can use seq to generate non-integer sequences
for (i in seq(0.0, 1.0, 0.1)) {
    cat(i, ' ', qnorm(i, lower.tail=FALSE, mean=4, sd=4), '\n')
}


