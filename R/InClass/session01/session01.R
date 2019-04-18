dna = read.csv("AncestryDNA.csv", head=TRUE, sep=",")
names <- c("alice", "bob", "cara", "don")
print("names:")
print(names)
class <- factor(c(1, 1, 2, 4))
print("class")
print(class)
ages <- c(20, 20, 21, 22)
print("ages:")
print(ages)

print("class is a factor. Class summary:")
print(summary(class))

print("ages is a numeric list. Ages summary:")
print(summary(ages))

print("frames store multiple objects of different types")
students = data.frame(name = names, year = class, age = ages)
summary(students)

print("tables are enhanced versions of frames")
grades <- factor(c("A","A","B","B","C","A","D","B","C","A"))
results <- table(grades)
results

print("We can query a table")
print(grades == "A")



