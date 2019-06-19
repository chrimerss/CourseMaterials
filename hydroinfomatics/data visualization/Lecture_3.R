#xinli@u.nus.edu

## This R scripts includes some of the basic R commands to extract the descriptive statistics of data, as well as some R statements for data visualization;
# To run this script, we need a lot of additional packages. However, some of the packages are not available from the default global (CDN) CRAN mirror, please set the CRAN mirror as below
## Click Tools/Global Options/Packages/; Then Click "Change" button on the right-hand side of the default "Global (CDN) - RStudio" text within "CRAN mirror" category, you will find a list of CRAN mirrors, choose "USA(CA1)[https]-University of California, Berkeley,CA (the 6th one from the bottom), then click "Ok", and "Ok"

## Please login to the IVLE, and download the dataset with the name "OwascoInlet_data.csv" (Files/HYDROINFORMATICS Files/Data Sets/OwascoInlet_data.csv) 

## Put the downloaded data set into your working directory

# At the bottom right panel of your RStudio, Click "More/Go To Working Directory", put the downloaded data set into this directory

# Read this data, store the data in a data frame named "my_data"

# Sometimes, if you run the script on your laptop, the window for plotting is too small, and RStudio may not function well. Please adjust the size of the bottom right window of RStudio, to make the plotting window larger 

# Try to select all the content in this script (Ctrl+all),and then press"Ctrl+Enter", see whether it runs fluently (if all the lines run well, you can see the last plot is scatter-histogram plot)
# First we install the packages we need to run this script. I would illustrate why we need this package in the corresponding context
install.packages(c('modeest','moments','psych','UsingR','ggplot2',
                   +'plotly','scales','ggthemes','ggExtra'))


my_data<-read.table('OwascoInlet_data.csv',header=TRUE,sep=' ') # Since this data are delimited/seperated by white space, sep should be ' '. If the data are seperated by ',', sep should be ','.

# Now we want to see the structure of this data set
str(my_data) # you can see this data set contains 888 observations of 6 variables, date, P_mm (Precipitation with unit of milimeter),Streamflow_m3s (streamflow with unit of cubic meter per sec); baseflow_m3s (baseflow with unit of cubic per sec),Tmax_C and Tmin_C (max and min daily temperature with unit of Celsius)

#### Descriptive Statistics

## Central Tendency (mean, median, mode)

# We just use Precipitation as an illustrative example, since my_data is a data.frame, the "$" can be used to refer to some of the variables (columns), Precipitation is corresponding to the variable "P_mm"; we assign it to a variable called Pre

Pre<-my_data$P_mm; Pre  # the precipitation data are stored in the variable "Pre", you can type the variable name "Pre", to print the values of Pre on the console

# mean

mean(Pre)

# trimmed mean (mean is sensitive to outliers, in order to tackle this issue, we can use trimmed mean)

mean(Pre, trim=1/10) # trim 1/10 off top and bottom

# median

median(Pre)

# mode
# In R, the function mode() return the type of an object. If you type mode(Pre), you will find it returns "numeric"
# However, there is a package called "modeest", which can calulate the mode of univariate unimodal data
# first install the package "modeest"
#install.packages('modeest') # This has already been done in a batch from the beginning. Same for the other packages
# load this package into working space
library('modeest')
# use the function mlv to calculate the mode. Have no idea about mlv, type ?mlv
?mlv  #it is a continuous variable
# read through the description, you may figure out an easy command to compute the mode of Pre
mlv(Pre,method='mfv') # mfv, which returns the most frequent values in a given mumerical vector
#for discrete variable


## Measures of variation (range, interquartile range, variance, standard deviation)

# range
range(Pre)

# 1st and 3rd quartile
quantile(Pre, probs=c(0.25,0.75))

# interquartile range
IQR(Pre)

# variance
var(Pre)

# standard deviation

sd(Pre)  

# We can use the function fivenum() to return Tukey's five number summary (min, lower-hinge, median, upper-hinge, maximum)

fivenum(Pre)



## High-order statistics (skewness and kurtosis)

# high-order moments like skewness and kurtosis can be calculated using the package called "moments"

#install.packages('moments')
library(moments)

skewness(Pre) 
kurtosis(Pre) 

# Examples above are for vector data (Precipitation), for a data frame which contains different variables, we can use the summary() function to extract some of the descriptive statistics

summary(my_data)


# You may recall from Lecture_1 that we can also use the function describe() in the package called "psych" to obtain a more complete set of descriptive statistics

#install.packages('psych')
library(psych)
describe(my_data)


#### Data Visualization

## Data visualization can be handled by different graphical systems, base graphics, lattice and ggplot2, etc. 
## We firstly introduce how to make basic plots (scatter plots, histograms, box plots, time series) using base graphics in R; Further, we show how to make decent and more beautiful plots using the package ggplot2


# There is a column "date" in my_data, this is in conflict with the date() in R, we rename it to date_1
names(my_data)[1]<-c('date_1')
# Before we start, we can use attach() function to attach our data (in this case,"my_data") to the R search path (In this way, objects in my_data can be accessed by simply giving their names); We can use the function detach() to detach it.
attach(my_data)
## scatter plot
# scatter plots are used to display the relationship between two continous variables, say, the relationship between Tmax_C and Tmin_C
plot(Tmax_C,Tmin_C,type='p',main='scatter plot of Tmax and Tmin',xlab='Tmax (°C)',ylab='Tmin (°C)') # type ?plot for the detailed argument of plot()

# histogram

hist(Pre,main='Histogram of Precipitation')

# box plot

boxplot(Pre,horizontal=TRUE,main='boxplot of Precipitation')

# we can even plot histogram and box plot together, using the function simple.hist.and.boxplot in a package called "UsingR"

#install.packages('UsingR')
library(UsingR)
simple.hist.and.boxplot(Pre,main='Histogram and box plot of Precipitation')

# Time series 

plot(as.Date(date_1),Streamflow_m3s,'o',main='Time series of streamflow')
# ISOdatetime() to create a standard time series
#or use strptime()

# ggplot2 has gained its fame in the R community in recent few years because of its versatility, clear and consistent interface,and beautiful output
#install.packages('ggplot2') 
library(ggplot2)

# scatter plot
p <- ggplot(my_data, aes(Tmax_C, Tmin_C),main='hao') # This serves a initialization of a ggplot object. Type ?ggplot for detailed description of the argument of ggplot()
p + geom_point()+ggtitle('scatter plot of Tmax and Tmin')+xlab('Tmax (°C)')+ylab('Tmin (°C)') # The point geom is used to create scatterplots

# histogram

p <- ggplot(my_data, aes(P_mm))
p + geom_histogram()+ggtitle('Histogram and box plot of Precipitation')+xlab('Precipitation (mm)')

# boxplot
my_string<-rep('Owasco',length(Pre)) 
station<-factor(my_string) # create a factor variable
my_data$station<-station 

p <- ggplot(my_data, aes(x=station,y=P_mm)) # we create a factor variable called "station" in order to make it easier to produce a boxplot in ggplot2
p + geom_boxplot()

# time series
# we need to install a package called "plotly" to facilitate plotting of time series
#install.packages('plotly')
#install.packages("scales")
library(ggplot2)
library(plotly)
library(scales)
ggplot(my_data, aes(x = as.Date(date_1), y =Streamflow_m3s )) +
  geom_point() +
  geom_line() +
  scale_x_date(breaks=date_breaks("6 months"),labels=date_format("%b %y"))+ # scale_x_date is aiming for positioning the dates in a more readable manner
  ggtitle('Time series of streamflow')+xlab('Date')+ylab('Streamflow (m^3/s)')


#### Tufte in R; please read the text in http://motioninsocial.com/tufte/ first before proceed to the following script

## Here I introduce a bit more about the art of visualization, the visualization practice developed by Edward Tufte.
# The idea behind "Tufte in R" is to produce a graph with minimal data-ink (total ink used to print the graph). This is important, because it illustrates the most elemental principle.

# Here I produced several plots for our data set

## Minimal Line Plot
# Say, using my_data,we want to see how the daily maximum temperature changes from 2009-04-22 to 2009-04-30 (The first 9 elements of Tmax_C), how to plot the minimum line plot using ggplot2,according to the instruction in http://motioninsocial.com/tufte/?
# First, we install a package called "ggthemes"
#install.packages('ggthemes')
library('ggthemes')
ggplot(my_data[1:9,], aes(x=as.Date(date_1[1:9]),y=Tmax_C[1:9])) + geom_line() + geom_point(size=3) + theme_tufte(base_size = 15) +
  theme(axis.title=element_blank())+scale_y_continuous(breaks=seq(5, 35, 5),label=sprintf('%s °C',seq(5,35,5)))

## Dot-dash plot in ggplot2
# say, we want to plot the scatter plot of Tmax and Tmin, to plot a dot-dash plot, we can use the following command:
ggplot(my_data, aes(Tmax_C, Tmin_C)) + geom_point() + geom_rug() + theme_tufte(ticks=F) + 
  xlab("Tmax (°C)") + ylab("Tmin (°C)") 

## Marginal histogram scatterplot
# To plot marginal histogram scatterplot, we need to use the "ggMarginal" function in package "ggExtra"
#install.packages('ggExtra')
library(ggExtra)
p <- ggplot(my_data, aes(Tmax_C, Tmin_C)) + geom_point() + theme_tufte(ticks=F)+xlab("Tmax (°C)") + ylab("Tmin (°C)") 
ggMarginal(p, type = "histogram", fill="transparent")










  









