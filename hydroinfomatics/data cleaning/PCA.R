# This R script exemplifies the method of PCA.
# The data set we used is Edgar Anderson's Iris Data, which gives the 
# measurements in centimeters of the variables sepal length and width
# and petal length and width, respectively, for 50 flowers from each of 3 species of iris.
# The species are Iris setosa, versicolor and virginica
# Please download "iris_data.csv" from IVLE and import it to working directory
my_data<-read.table('iris_data.csv',sep=',',header = TRUE)
# We can perform a principal component analysis with the princomp(dataset).
# To be noted, the dataset should contain numeric variables only. If there are
# non-numeric variables in your data set, you must exclude them with bracket notation or with the subset() function
my_data<-subset(my_data,select = c('Sepal.Length','Sepal.Width','Petal.Length','Petal.Width'))
# or just use my_data<-my_data[-5]

my_model<-princomp(my_data)

# we can use the summary() function to view the proportion of the total variance explained by each component

summary(my_model)

# From the output we can see that 92.4% of the variation in the dataset is explained by the first component alone, and 97.8% is explained by the first two components.

# To use the loadings of each component, we use this command
my_model$loadings

# Similarly, we can view the scores for each of the observations as shown
my_model$scores

# We can also plot the "scree plot", which plots the variances against the number of the principal component

screeplot(my_model)



























# reference: http://www.instantr.com/2012/12/18/performing-a-principal-component-analysis-in-r/