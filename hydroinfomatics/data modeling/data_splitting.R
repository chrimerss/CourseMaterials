# data splitting
# Content below is from https://ragrawal.wordpress.com/2012/01/14/dividing-data-into-training-and-testing-dataset-in-r/

# During machine learning one often needs to divide the two different data sets, namely training and testing datasets.
# While you can't directly use the "sample" command in R, there is a simple workaround for this. 
# Essentially, use the "sample" command to randomly select certain index number and then use the selected index numbers to divide the dataset into training and testing dataset. 
# Below is the sample code for doing this. In the code below I use 20% of the data for testing and rest of the 80% for training.

# By default R comes with few datasets. 
data = mtcars
dim(data)  # 32 11

#Sample Indexes
indexes = sample(1:nrow(data), size=round(0.2*nrow(data)))

# Split data
test = data[indexes,]
dim(test)  # 6 11
train = data[-indexes,]
dim(train) # 26 11





# Sometimes, we need to divide the data frames into three groups, trainning data, validation data, and testing data

# say, we need to divide the data into 0.6 trainning data, 0.2 validation data, 0.2 testing data
nr<-dim(data)[1]
# shuffle the data
my_data<-data[sample.int(nr),] # if we just want to select the first 60% as training data, the subsequent 20% and 20% as validation and testing data,respectively, then we do not need to do the shuffling.
train_index=round(nr*0.6); # training index
val_index=round(nr*(0.6+0.2)); # validation index
# the remaining is testing index

train_data<-my_data[1:train_index,]
val_data<-my_data[(train_index+1):val_index,]
test_data<-my_data[(val_index+1):nr,]
