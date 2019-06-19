# This R script exemplifies the smoothing method using LOESS
# We use the 'OwascoInlet_data.csv' data from IVLE
my_data<-read.table('OwascoInlet_data.csv',sep=' ',header = TRUE)
Tmax<-my_data$Tmax_C
Tmin<-my_data$Tmin_C
# plot the scatterplot
plot(Tmax,Tmin,type='p',main='scatter plot of Tmax and Tmin',xlab='Tmax (°C)',ylab='Tmin (°C)')

# smooth with Loess using scatter.smooth() (Plot and add a smooth curve computed by loess to a scatter plot.)
scatter.smooth(Tmax, Tmin, span = 0.5, degree = 1,
             family = c("symmetric", "gaussian"), evaluation = 50,col=2)
# span: smoothness parameter for loess
# degree: degree of local polynomial used
# family: if "gaussian", fitting by least-squares
# evaluation: number of points at which to evaluate the smooth curve

# if we want to perform a quadratic fit using loess
scatter.smooth(Tmax, Tmin, span = 0.5, degree = 2,
               family = c("symmetric", "gaussian"), evaluation = 50,col=2)

