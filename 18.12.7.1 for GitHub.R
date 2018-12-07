library(devtools)
library(ggord)
library(labeling)
library(ggplot2)

library(e1071)
library(MASS)
library(sampling)
library(readxl)
library(caret)
library(dplyr)

library(mlr)
library(selectr)
library(FSelector)
library(kernlab)
library(kknn)
library(LiblineaR)

#Count accuracy
count_accuracy <- function(result,data_test){  
  n <- length(result)  
  count_correct <- 0  
  i <- 1  
  for (i in 1:n){  
    if (result[i]==data_test[i,1]){  
      count_correct = count_correct +1  
    }  
  }  
  print(count_correct/n)  
}  
#Read dataset and data split
trainingSetnormal<-read.csv("path", header=TRUE)
index<-createDataPartition(trainingSetnormal$Name, time=1, p=0.625, list=F)
trainingSet<-trainingSetnormal[index, ]
testSet <-trainingSetnormal[-index, ]

#feature selection for PC3/nucleotides
#Feature number vs accuracy plot
#LDA
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sfbs", max.features =1, beta=0)
rdesc = makeResampleDesc("Bootstrap", predict = 'test', iters =100, stratify = TRUE)
lrn = makeLearner("classif.lda")
sfeatspoldaPC1 = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatspoldaPC1)
#SVM
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sfbs", max.features =1)
rdesc = makeResampleDesc("Bootstrap",  predict = 'test', iters =100,stratify = TRUE)
lrn = makeLearner("classif.svm")
sfeatsposvmlinearPC1 = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatsposvmlinearPC1)
#KNN
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sfbs", max.features =1)
rdesc = makeResampleDesc("Bootstrap",  predict = 'test', iters =100, stratify = TRUE)
lrn = makeLearner("classif.kknn",k=3)
sfeatsknnPC1 = selectFeatures(learner = lrn, classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatsknnPC1)
#Forward search
#LDA
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sffs", max.features = 5)
rdesc = makeResampleDesc("Bootstrap", predict = 'test', iters =100, stratify = TRUE)
lrn = makeLearner("classif.lda")
sfeatspoldaPC1 = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatspoldaPC1)
#SVM
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sffs", max.features = 5)
rdesc = makeResampleDesc("Bootstrap",  predict = 'test', iters =100,stratify = TRUE)
lrn = makeLearner("classif.svm")
sfeatsposvmlinearPC1 = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatsposvmlinearPC1)
#KNN
classif.task = makeClassifTask(data = trainingSet, target = "Name")
ctrl = makeFeatSelControlSequential(method = "sffs", max.features = 5)
rdesc = makeResampleDesc("Bootstrap",  predict = 'test', iters =100, stratify = TRUE)
lrn = makeLearner("classif.kknn",k=3)
sfeatsknnPC1 = selectFeatures(learner = lrn, classif.task, resampling = rdesc, control = ctrl, show.info = FALSE)
analyzeFeatSelResult(sfeatsknnPC1)

#Test set accuracy
#LDA
trainingSet<-trainTransformed[index, ]
testSet <-trainTransformed[-index, ]
trainingSet<-trainingSet[,c('Name',sfeatspoldaPCfs$x)]
testSet <-testSet[,c('Name',sfeatspoldaPCfs$x)]
ions.lda=lda(formula = Name ~ ., data=trainingSet)
ions.ldapredict=predict(ions.lda,testSet[,-1],decision.values=T)
ions.ldapredict
table(testSet$Name,ions.ldapredict$class)
count_accuracy(ions.ldapredict$class,testSet)
#LDA plot
p<-ggord(ions.lda,trainingSet$Name)
p
#KNN
trainingSet<-trainTransformed[index, ]
testSet <-trainTransformed[-index, ]

trainingSet<-trainingSet[,c('Name',sfeatsknnPCfs$x)]
testSet <-testSet[,c('Name',sfeatsknnPCfs$x)]
data.fit <- knn3(factor(Name) ~ ., trainingSet, k = 3)
data.fit
ions.ldapredict=predict(data.fit,testSet[,-1],type = "class")
ions.ldapredict
table(testSet$Name,ions.ldapredict)
count_accuracy(ions.ldapredict,testSet)
#SVM
trainingSet<-trainTransformed[index, ]
testSet <-trainTransformed[-index, ]
trainingSet<-trainingSet[,c('Name',sfeatsposvmlinearPCfs$x)]
testSet <-testSet[,c('Name',sfeatsposvmlinearPCfs$x)]
ions.svm=svm(Name ~ .,data=trainingSet, kernel="linear")
summary(ions.svm)
ions.svmpredict=predict(ions.svm,testSet[,-1],decision.values=T)
ions.svmpredict
table(testSet$Name,ions.svmpredict)
count_accuracy(ions.svmpredict,testSet)

#Repeated cross-validation
accuracyldaall=0
accuracyknnall=0
accuracysvmall=0
for(i in 1:30){
  index<-createDataPartition(trainingSetnormal$Name, time=1, p=0.625, list=F)
  trainingSet<-trainingSetnormal[index, c('Name', sfeatspoldaPCfs$x)]
  testSet <-trainingSetnormal[-index, c('Name', sfeatspoldaPCfs$x)]
  #LDA
  ions.lda=lda(formula = Name ~ ., data=trainingSet)
  ions.ldapredict=predict(ions.lda,testSet[,-1],decision.values=T)
  ions.ldapredict
  AccuracyLDA = count_accuracy(ions.ldapredict$class,testSet)
  accuracyldaall = accuracyldaall+AccuracyLDA
  #KNN
  trainingSet<-trainTransformed[index, c('Name', sfeatsknnPCfs$x)]
  testSet <-trainTransformed[-index, c('Name', sfeatsknnPCfs$x)]
  data.fit <- knn3(factor(Name) ~ ., trainingSet, k = 3)
  ions.knnpredict=predict(data.fit,testSet[,-1],type = "class")
  Accuracyknn = count_accuracy(ions.knnpredict,testSet)
  accuracyknnall = accuracyknnall+Accuracyknn
  #SVM  
  trainingSet<-trainTransformed[index, c('Name', sfeatsposvmlinearPCfs$x)]
  testSet <-trainTransformed[-index, c('Name', sfeatsposvmlinearPCfs$x)]
  ions.svm=svm(Name ~ .,data=trainingSet, kernel="linear")
  ions.svmpredict=predict(ions.svm,testSet[,-1],decision.values=T)
  Accuracysvm = count_accuracy(ions.ldapredict$class,testSet)
  accuracysvmall = accuracysvmall+Accuracysvm
}
accuracyldaall/30
accuracyknnall/30
accuracysvmall/30
