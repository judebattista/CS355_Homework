heisenberg <- read.csv(file="simple.csv", head=TRUE, sep=",")
#Summary data! How cool is that?
summaryData <- summary(heisenberg)
print("dir() = ls")
dir()
print("getwd = pwd")
getwd()

print("Contents of the simple.csv file:")
print(heisenberg)
print("names(var) gets the column names for the csv")
print("You can reference columns by varname$colname")
cols = names(heisenberg)
heisenberg$cols[1]      # this does not work, even though cols[1] = trial
heisenberg$trial        # this does work

print("ls() gets you all the variable names in your workspace")
ls()

