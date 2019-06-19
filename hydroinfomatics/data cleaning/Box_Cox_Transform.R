# This R script exemplifies the method of Box-Cox Transformation
# To perform this transformation, we need to install a package called "forecast"
# We use the Streamflow data from "OwascoInlet_data.csv" from IVLE
install.packages('forecast')
library(forecast)
my_data<-read.table('OwascoInlet_data.csv',sep=' ',header = TRUE)
Streamflow<-my_data$Streamflow_m3s
# We can plot the histogram of precipitation. You can see that the distribtion of precipitaiton is not normal
hist(Streamflow)
#normalization
log_function<-function(x){
  y<-log10(x)
}
stream_new<-log_function(Streamflow)
hist(stream_log)
stream_log<-unlist(stream_log)
class(stream_log)
stream_log<-lapply(Streamflow,log10)
# Now we perform the Box-Cox transformation
# Details of Box-Cox transformation
# f(x;lambda)=(x^lambda - 1)/lambda if lambda is not equal to 0. For lambda=0,f(x;0)=log(x)

# Hence, the first step is to estimate the parameter, i.e. lambda, we use the function BoxCox.lambda() from package "forecast"
lambda<-BoxCox.lambda(Streamflow, method=c("loglik"), lower=-3, upper=3);lambda
# Then, we perform the Box-Cox transformation using the BoxCox() function from package "forecast"

Streamflow_transformed<-BoxCox(Streamflow, lambda);Streamflow_transformed
# see whether the transformed data is normal
hist(Streamflow_transformed)
# Maybe we can put these two histograms together for comparison
par(mfrow=c(1,2))
hist(Streamflow)
hist(Streamflow_transformed)
# We can also check the normality of data using Q-Q plot. Type ?qqnorm for details.
par(mfrow=c(1,2))
qqnorm(Streamflow)
qqline(Streamflow,col=2) # col=2 means the red color
qqnorm(Streamflow_transformed)
qqline(Streamflow_transformed,col=2)
