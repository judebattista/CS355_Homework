cars <- read.csv('auto_mpg.csv', header=TRUE)
summary(cars)

print(cars)
# Yeah, that's too much data
# Let's get the column names
names(cars)

# attach will bring in all the columns of the table cars into separate vectors:
attach(cars)
# now we have vectors named after each column in cars!
#print(mpg)

# Let's create a linear model that compares weight to mpg using our newly created vectors
model0 <- lm(mpg ~ weight)
summary(model0)

# we can create plots directly in R
# make the plot area:
plot(0.0, 0.0, ylim=c(0.0, 50.0), xlim = c(0.0, 5500), xlab="weight", ylab="mpg")

# add a title
title(main="Mileage Characteristics of Common Automobiles")

# add the points
points(weight, mpg, col='red', pch=10)

# can we make a better model?
model1 <- lm( mpg ~ weight + cylinders + displacement + horsepower + model_year )
summary(model1)

# we can use a model and a data set to make predictions
predicted_mpg <- predict(model1, cars)

plot(0.0, 0.0, ylim = c(0.0,50.0), xlim = c(0.0,5500), xlab="weight", ylab="mpg")
title(main="Predicted Mileage Characteristics of Common Autos")
points(weight, mpg, col = "red", pch=10)
# Add predicted points in another color
points( weight, predicted_mpg, col = "blue" )

# We can do a Principle Components Analysis (PCA)
# princomp_fit <- princomp(cars, cor=TRUE)
# throws an error: we have NA values in our data

# princomp_fit <- princomp(na.omit(cars), cor=TRUE)
# still have an error. Oh hey, we have a text column for car name. That's... not numeric

data <- na.omit(cars[1:8])
princomp_fit <- princomp(data, cor=TRUE)
summary(princomp_fit)
loadings(princomp_fit)
