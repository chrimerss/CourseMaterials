# This script exemplifies the methods to normalize the data, including Min-Max Normalization, Z-score Normalisation and Decimal-scaling Normalisation

# We use the data "OwascoInlet_data.csv" from IVLE

my_data<-read.table('OwascoInlet_data.csv',sep=' ',header = TRUE)

# To make it simple, we just use a subset of this data, i.e, data during the period from 2009-05-01 to 2009-05-31
my_data<-my_data[10:40,]


# Now we want to normalize the data for each of its columns. To attain this, we can create our own function, and apply it to each of the column of the data.frame using lapply()
# Since the class of my_data$date is a factor, we want to change it to numeric such that our defined functio can work on the whole data.frame
my_data$date<-as.numeric(my_data$date)
# you can also use the subset() to exclude the first column. See "PCA.R"
# Or just use my_data[-1] to exclude the first column

# Min-Max normalization

min_max <- function(x){
y<-(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x, na.rm=TRUE)) 
return (y)
} # normalize the data to the range [0,1]

my_data_normed1<-as.data.frame(lapply(my_data,min_max))

# Z-score normalization

Z_score<-function(x){
  y<-(x - mean(x, na.rm=TRUE))/(sd(x,na.rm = TRUE)) 
return(y)
}
my_data_normed2<-as.data.frame(lapply(my_data,Z_score))

# Decimal normalization
# we need to define a function to round up a number to the nearest power of 10
roundUp <- function(x) 10^ceiling(log10(x)) # type ?ceiling to see what does it mean

Decimal_norm<-function(x){
y<-x/roundUp(max(x)) # say, the maximum number of x is 87, we round it up to the nearest power of 10, which is 100. 
return(y)            # then, for decimal normalization, all the values of x should be divided by 100
}
my_data_normed3<-as.data.frame(lapply(my_data,Decimal_norm))

