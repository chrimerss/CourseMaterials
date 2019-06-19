# Simple Linear regression.
# Suppose we know that the maximum heart rate of a person is related to age, now we want to build a simple linear regression model to study the relationship between these two variables.
# And we have a sample of (age, maximum heart rate) as below:
age<-c(18,23,25,35,65,54,34,56,72,19,23,42,18,39,37)
max_heart_rate<-c(202,186,187,180,156,169,174,172,153,199,193,174,198,183,178)

# First, we plot the scatterplot
plot(age,max_heart_rate)
# the simple regression model can be built using lm(y~x). y is the dependent variable, x is the independent variable or predictor
lm(max_heart_rate~age)
# the regression line can be plotted using abline(lm(y~x))
abline(lm(max_heart_rate~age))
install.packages('UsingR')
library(UsingR)
# This can also be done using simple.lm() function
lm.result<-simple.lm(age,max_heart_rate) # the result of simple.lm(x,y) is of class lm, so the plot and summary commands adapt themselves to that
# if you want to plot confidence interval, set show.ci=TRUE
simple.lm(x=age,y=max_heart_rate,show.ci=TRUE) 
# summary of the results
summary(lm.result) # This give you a whole summary of the statistics, regression coef, standard error, p-value, etc

# Or you want to extract some of them
# coef
coef(lm.result)
# Residuals.
lm.res=resid(lm.result)
summary(lm.res)
# if you want to get all the coefficients
coefficients(summary(lm.result))
# Check with model validity. 
# The validity of the model can be checked graphically via eda. The assumption on the errors being i.i.d. normal
# random variables translates into the residuals being normally distributed. They are not independent as they add to
# 0 and their variance is not uniform, but they should show no serial correlations.

# We can test for normality with eda tricks: histograms, boxplots and normal plots. We can test for correlations
# by looking if there are trends in the data. This can be done with plots of the residuals vs. time and order. We can
# test the assumption that the errors have the same variance with plots of residuals vs. time order and fitted values.
# The plot command will do these tests for us if we give it the result of the regression

par(mfrow=c(2,2)) # plot 4 separate graphs placed in 2*2
plot(lm.result)

#Residuals vs. fitted: This plots the fitted (by) values against the residuals. Look for spread around the line y = 0
# and no obvious trend.
# Normal qqplot: The residuals are normal if this graph falls close to a straight line.
# Scale-Location: This plot shows the square root of the standardized residuals. The tallest points, are the largest
#residuals.
#Cook's distance: This plot identifies points which have a lot of influence
# in the regression line.

# If you want to get the fitted values
fitted(lm.result) # yi(fitted)=b0+b1*xi

# to use model to get predicted values. Say if we input new predictor value c(50,60,70,80,90), we can use the predict() to get the corresponding prediction
new<-data.frame(age=c(50,60,70,20,60))
predict(lm(max_heart_rate~age),new)



# say, if we want to add both confidence interval and prediction interval; see http://www2.stat.duke.edu/~tjl13/s101/slides/unit6lec3H.pdf for descripition of prediction interval and confidence interval for simple linear regression

pred.w.plim <-predict(lm(max_heart_rate~age),new,interval="prediction") # add prediction interval
pred.w.clim <-predict(lm(max_heart_rate~age),new,interval="confidence") # add confidence interval

# matplot is used to plot columns of matrices
matplot(new$age, cbind(pred.w.clim, pred.w.plim[,-1]), col=1:5,lty = c(1,2,2,3,3),type = "l", ylab = "predicted y") # lty to set line type, see http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
legend("topright", inset=.05, legend=c("model_fit", "lower confidence limit","upper confidence limit","lower prediction limit","upper prediction limit"), col=1:5,lty = c(1,2,2,3,3),horiz=FALSE)


# Ref: https://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf