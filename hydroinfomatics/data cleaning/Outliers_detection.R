# This R script exemplifies some the methods to detect outliers

# Box-and-whisker method

# For more or less unimodal and symmetrically distributed data, Tukey's box-and-whisker
# method for outlier detection is often appropriate. In this method, an observation is an outlier
# when it is larger than the so-called ``whiskers'' of the set of observations. The upper whisker is
# computed by adding 1.5 times the interquartile range to the third quartile and rounding to the
# nearest lower observation. The lower whisker is computed likewise.


# we can use the function boxplot.stats() to list the outliers
x<-c(1:10, 20, 30)
boxplot.stats(x)$out

# Here, 20 and 30 are detected as outliers since they are above the upper whisker of the
# observations in x. The factor 1.5 used to compute the whisker is to an extent arbitrary and it can
# be altered by setting the coef option of boxplot.stats. A higher coefficient means a higher
# outlier detection limit (so for the same dataset, generally less upper or lower outliers will be detected).
boxplot.stats(x, coef = 2)$out

# Hiridoglou and Berthelot's method
# This method takes skewness into account. It is appropriate for finding outliers on both sides of the distribution.
# Moreover, because of the different behaviour for small and large x-values, it is
# appropriate for skewed (long-tailed) distributions.

hboutlier <- function(x,r){
  x <- x[is.finite(x)]
  stopifnot(
    length(x) > 0
    , all(x>0)
  )
  xref <- median(x)
  if (xref <= sqrt(.Machine$double.eps))
    warning("Reference value close to zero: results may be inaccurate")
  pmax(x/xref, xref/x) > r
} # This function returns a logical value for each of the elements in x.
  # "True" means the corresponding element is an outlier

z<-rexp(10,1) # generate 10 numbers from exponential distribution
hboutlier(z,5) # here the reference value r is set to be 5.

# reference: https://cran.r-project.org/doc/contrib/de_Jonge+van_der_Loo-Introduction_to_data_cleaning_with_R.pdf
