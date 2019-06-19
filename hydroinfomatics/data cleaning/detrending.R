# This script exemplifies how to detrend the global annual mean temperature anomaly time series using temporal differentiation method

# The global annual mean temperature anomaly data can be downloaded from http://data.okfn.org/data/core/global-temp

# Import the "Global_annual_mean_temperature_annomaly.csv" data into working directory.

my_data<-read.table("Global_annual_mean_temperature_annomaly.csv",sep=',',header = TRUE)

# The original data come from two souces, global component of Climate at a Glance (GCAG) and GISS Surface Temperature (GISTEMP)
# We seperate the data into two data sets, and sort them in ascending order of "Year".
attach(my_data) # attach() function helps to load "my_data" into the working space, such that you can refer to its variable (Source,Year, Mean) directly
GCAG_data<-my_data[Source=='GCAG',]
GIS_data<-my_data[Source=='GISTEMP',]
GCAG_data_new<-GCAG_data[order(GCAG_data$Year),] # order a data.frame by one of the column, using order(). The default is in ascending order
GIS_data_new<-GIS_data[order(GIS_data$Year),]
par(mfrow=c(2,1)) # Display figures in two rows, one column
plot(GCAG_data_new$Year,GCAG_data_new$Mean,type='o',xlab='year',ylab='Global Annual Mean Temp Anomaly',main='GCAG')
plot(GIS_data_new$Year,GIS_data_new$Mean,type='o',xlab='year',ylab='Global Annual Mean Temp Anomaly',main='GIS')

# Now we want to detrend the time series using temporal differentiation method diff()

GCAG_detrend<-diff(GCAG_data_new$Mean) # You can find  that the dimension of GCAG_detrend is (n-1), since the differentiation operates like y2-y1,y3-y2,...y(n)-y(n-1)
GIS_detrend<-diff(GIS_data_new$Mean)

par(mfrow=c(2,1))
# To be noted, Year: 2:n
plot(GCAG_data_new$Year[2:length(GCAG_data_new$Year)],GCAG_detrend,type='o',xlab='year',ylab='Detrended Global Annual Mean Temp Anomaly',main='GCAG: Detrended')
plot(GIS_data_new$Year[2:length(GIS_data_new$Year)],GIS_detrend,type='o',xlab='year',ylab=' Detrended Global Annual Mean Temp Anomaly',main='GIS: Detrended')
