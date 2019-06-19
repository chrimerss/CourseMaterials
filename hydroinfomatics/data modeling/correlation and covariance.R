# covariance and correlation.
# Suppose we know that the maximum heart rate of a person is correlated with age.
# And we have a sample of (age, maximum heart rate) as below, how do we calculate the covariance and correlation?
# There are three types of covariance and correlation, pearson r correlation (linear relationship), kendall rank correlation, and spearman rank correlation.
# For the detailed introduction of these three measures, please see http://www.statisticssolutions.com/correlation-pearson-kendall-spearman/

age<-c(18,23,25,35,65,54,34,56,72,19,23,42,18,39,37)
max_heart_rate<-c(202,186,187,180,156,169,174,172,153,199,193,174,198,183,178)

# First, we plot a scatterplot to see what is the relationship like

plot(age,max_heart_rate,type='p')

# Now we can calculate the three measures of covariance and correlation

# pearson r cov and cor
cov(age,max_heart_rate,method='pearson')
cor(age,max_heart_rate,method='pearson')
# Kendall rank cov and cor
cov(age,max_heart_rate,method='kendall')
cor(age,max_heart_rate,method='kendall')
# Spearman rank cov and cor\
cov(age,max_heart_rate,method='spearman')
cor(age,max_heart_rate,method='spearman')

# If you want to calculate the confidence interval or p-value, use cor.test()
cor.test(age,max_heart_rate,method='pearson')
cor.test(age,max_heart_rate,method='kendall')
cor.test(age,max_heart_rate,method='spearman')


# the cor() and cov() can be applied directly to a data.frame, which would generate a correlation or covariance matrix










# Reference: https://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf