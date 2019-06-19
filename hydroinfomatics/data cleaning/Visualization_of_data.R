# This R script contains some examples regarding barplot, boxplot, simple.scatter plot and some lattice plots

## boxplot

# The data we are going to use is the data set "PlantGrowth", which is included in R base
data("PlantGrowth") # There are 3 groups, one control and two treatments. For each group,weight are recorded.
PlantGrowth # The data is generated in this way,by recording a weight and group for each plant.

# We can attach the data set to our search path, such that we can directly refer to its variable
attach(PlantGrowth)

# We want to re-organize the data, say, with  variables recorded corresponding to the levels of the factor
# this can be done using the unstack()
unstack(PlantGrowth)
# Then we can do a box plot of these three groups
library(ggplot2)
ggplot(PlantGrowth,aes(group,weight))+geom_boxplot()
boxplot(unstack(PlantGrowth))
# This can also be done by the model formula notation that R uses allows this to be done in a systematic manner
boxplot(weight~group)

## barplot
# First, let us consider a univariate categorical data
# suppose, a group of 25 people are surveyed as to their beer-drinking preference. The categories were (1) Domestic
# can, (2) Domestic bottle, (3) Microbrew and (4) import, the raw data is 
# 3 4 1 1 3 4 3 3 1 3 2 1 2 1 2 3 2 3 1 1 1 1 4 3 1
beer<-c(3,4,1,1,3,4,3,3,1,3,2,1,2,1,2,3,2,3,1,1,1,1,4,3,1)
# when using the barplot(), we need to summarize the data first
table(beer) # table() build a contingency table of the counts at each combination of factor levels
barplot(table(beer))
# we can also plot the proportion for each of the group
barplot(table(beer)/length(beer))
# Then, we consider the bivariate categorical data
# suppose a student survey is done to evaluate if students who smoke study less, the data recorded is
# Person  Smokes   Amount of studying
#  1       Y       less than 5 hours
#  2       N       5-10 hours
#  3       N       5-10 hours
#  4       Y       more than 10 hours
#  5       N       more than 10 hours
#  6       Y       less than 5 hours
#  7       Y       5-10 hours
#  8       Y       less than 5 hours
#  9       N       more than 5 hours
#  10      Y       5-10 hours

# We can handle this in R by creating two vectors to hold our data
smokes = c("Y","N","N","Y","N","Y","Y","Y","N","Y")
amount = c(1,2,2,3,3,1,2,1,3,2)
table(smokes,amount)
barplot(table(smokes,amount))
barplot(table(amount,smokes))
smokes=factor(smokes) # for names
barplot(table(smokes,amount),beside = TRUE,legend.text = TRUE) # put beside not stacked. The default one is stacked
barplot(table(amount,smokes),main="table(amount,smokes)",beside = TRUE,legend.text = c('less than 5','5-10','more than 10'))

# Violinplots and density plots, compared with boxplot
# The functions simple.violinplot and simple.densityplot can be used in place
# of side-by-side boxplots to compare different distributions.
# A violinplot is very similar to a boxplot, only the box is replaced by a density which is given a mirror image for
# clarity. A densityplot plots several densities on the same scale. Multiple histograms would look really awful,
# but multiple densities are manageable.

par(mfrow=c(1,3)) # 3 graphs per page
data(InsectSprays) # load in the data
InsectSprays
boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
install.packages("UsingR")
library(UsingR)
simple.violinplot(count ~ spray, data = InsectSprays, col = "lightgray")
simple.densityplot(count ~ spray, data = InsectSprays)

# simple scatterplot and paired scatterplot
# we use the data set "emissions" which contains data on the GDP and CO2 emissions for several European countries and the United states for the year 1999
data(emissions)
emissions
attach(emissions)
simple.scatterplot(perCapita,CO2) # Shows scatterplot of x vs y with histograms of each on sides of graph

pairs(emissions) # A matrix of scatterplots is produced. 
#The pairs command will produce scatterplots for each possible pair


# lattice package
# we have introduced making beautiful plots using ggplot2. There are also some other packages available, such as the package "lattice"

install.packages('lattice')
library(lattice)
# for the box plot we have created, it can also be created using bwplot() in lattice
bwplot( ~ weight | group , data =PlantGrowth )
# histogram
histogram(~count | spray, data = InsectSprays)
# density plot
densityplot(~count | spray, data = InsectSprays)
# scatter plot matrix
splom(emissions)
# barplot
barchart(amount~smokes,groups=smokes)

