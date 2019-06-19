# Multiple Linear Regression
# The data set we used is "enrollmentForecast.csv" which can be downloaded from IVLE.
# This data set contains information used to estimate undergraduate enrollment at the University of New Mexico (Office of Institutional Research, 1990).
# The variables in this data set
# Roll: fall enrollment
# UNEM: unemployment rate
# HGRAD: number of spring high school graduates
# INC: per capita income 
my_data<-read.table('enrollmentForecast.csv',sep=',',header=TRUE)

# we want to see the the pair scatter plot 
plot(my_data) # when the plot() function operates on a data.frame, it gives the pair scatter plot.
# this can also be done by using
pairs(my_data)
 
# the lm() function can be used to create a multiple linear regression model, there are two parameters in lm(), formula and data
# formula: YVAR~XVAR1+XVAR2+XVAR3+...XVARi where YVAR is the dependent, or predicted variable and XVARs are the independent or predictors variable
# data, the variable that contains the data set.

# It is recommended that you save a newly created linear model into a variable, by doing so, the model can be used in subsequent calculations and analyses without having to retype the entire lm() function each time.

# In this particular case, we are using the unemployment rate (UNEM) and number of spring high school graduates (HGRAD) to predict the fall enrollment (ROLL)

# Creating a linear model with two predictors
twoPredictorModel <- lm(ROLL ~ UNEM + HGRAD,data=my_data)
# what is the model
twoPredictorModel # as we can see, the model is: Fall Enrollment = -8255.8 + 698.2 * Unemployment Rate + 0.9 * Number of Spring High School Graduate
# we can use summary() function to display a wealth of important information of a linear model
summary(twoPredictorModel)
# The summary(OBJECT) function has provided us with t-test, F-test, R-squared, residual, and significance values. All of this data can be used to answer important questions related to our models.

# We want to plot the results of the regression, to see whether our model is valid (i.i.d assumption of residuals)
par(mfrow=c(2,2))
plot(twoPredictorModel)

# Creating a linear model with three predictors
# Simply, one can just continue to add variables to the FORMULA argument until all of them are accounted for
# say, if we want to add another variable, per capita income INC
threePredictorModel <- lm(ROLL ~ UNEM + HGRAD + INC, data=my_data)
# what is the model
threePredictorModel
# use summary() function to display information about this model
summary(threePredictorModel)
# From the results, we can see, compared with the two-predictor model, performance of three-predictor model is better (larger R^2 and adjusted R^2, which means the three-predictor model can explain more variance than the two-predictor model)
# We want to plot the results of the regression, to see whether our model is valid (i.i.d assumption of residuals)
par(mfrow=c(2,2))
plot(threePredictorModel)


# Model prediction, confidence and prediciton intervals
# We often use our regression models to estimate the mean response or predict
# future values of the response variable for certain values of the response
# variables. The function predict() can be used to make both confidence intervals
# for the mean response and prediction intervals. To make confidence intervals for
# the mean response use the option interval="confidence". To make a prediction
# interval use the option interval="prediction". By default this makes 95%
# confidence and prediction intervals. If you instead want to make a 99%
# confidence or prediction interval use the option level=0.99

# say if we want to predict the mean fall enrollment with UNEM=7.8, HGRAD=19500 and INC=3000
new=data.frame(UNEM=7.8,HGRAD=19500,INC=3000)

# model prediction
model_prediction<-predict(threePredictorModel,new);model_prediction

# this is the same as  

model_prediction<-predict(lm(ROLL ~ UNEM + HGRAD + INC, data=my_data),new);model_prediction

# confidence interval
model_prediction_confidence<-predict(threePredictorModel,new,interval = "confidence");model_prediction_confidence
# prediction interval
model_prediction_predicinterval<-predict(threePredictorModel,new,interval = "prediction");model_prediction_predicinterval


# Finally, we want to see how our model fit the obervation
model_fit<-threePredictorModel$fitted.values
observation<-my_data$ROLL
model_obs<-data.frame(observation,model_fit)

# matplot is used to plot columns of matrices
matplot(my_data$YEAR, model_obs, lty = c(1,2),col=1:2,type = "l", ylab = "Roll")

legend("bottomright", inset=.05, legend=c("observation", "model_fit"), col=1:2,lty=c(1,2), horiz=TRUE)










# Reference 1: https://www.r-bloggers.com/r-tutorial-series-multiple-linear-regression/
# Reference 2: https://www.r-bloggers.com/wp-content/uploads/2009/12/multiLinRegExample1.txt
# Reference 3: http://www.stat.columbia.edu/~martin/W2024/R6.pdf