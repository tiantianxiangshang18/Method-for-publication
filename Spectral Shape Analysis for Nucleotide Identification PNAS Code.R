#Code for 'Spectral Shape Analysis for Nucleotide Identification'
# www.pnas.org/cgi/doi/10.1073/pnas.1820713116
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
trainingSetnormal<-read.csv("C:/Users/Yun/OneDrive/Schanze group one drive/18.9.30 paper writing/paper writing/origin data/fluorescence/18.6.7 final/all normal.csv", header=TRUE)
index<-createDataPartition(trainingSetnormal$Name, time=1, p=0.625, list=F)
trainingSet<-trainingSetnormal[index, ]
testSet <-trainingSetnormal[-index, ]

#Feature selection
#Feature number vs accuracy plot
#Sequential Backward Floating Selection (SBFS)

SBFS <- function(dataset, algorithm){
  classif.task = makeClassifTask(data = dataset, target = "Name")
  ctrl = makeFeatSelControlSequential(method = "sfbs", max.features =1, beta=0)
  rdesc = makeResampleDesc("Bootstrap", predict = 'test', iters =100, stratify = TRUE)
  lrn = algorithm
  result = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
  return(result)
}

#LDA SBFS_LDA_result
SBFS_LDA_result = SBFS(dataset = trainingSet, algorithm = makeLearner('classif.lda'))
analyzeFeatSelResult(SBFS_LDA_result)

#SVM SBFS_SVM_result
SBFS_SVM_result = SBFS(dataset = trainingSet, algorithm = makeLearner('classif.svm'))
analyzeFeatSelResult(SBFS_SVM_result)

#KNN SBFS_KNN_result
SBFS_KNN_result = SBFS(dataset = trainingSet, algorithm = makeLearner('classif.kknn', k=3))
analyzeFeatSelResult(SBFS_KNN_result)

#Sequential Forward Floating Selection (SBFS)
SFFS <- function(dataset, algorithm){
  classif.task = makeClassifTask(data = dataset, target = "Name")
  ctrl = makeFeatSelControlSequential(method = "sffs", max.features = 5)
  rdesc = makeResampleDesc("Bootstrap", predict = 'test', iters =100, stratify = TRUE)
  lrn = algorithm
  result = selectFeatures(learner = lrn, task = classif.task, resampling = rdesc,control = ctrl, show.info = FALSE)
  return(result)
}

#LDA SFFS_LDA_result
SFFS_LDA_result = SFFS(dataset = trainingSet, algorithm = makeLearner('classif.lda'))
analyzeFeatSelResult(SFFS_LDA_result)

#SVM SFFS_SVM_result
SFFS_SVM_result = SFFS(dataset = trainingSet, algorithm = makeLearner('classif.svm'))
analyzeFeatSelResult(SFFS_SVM_result)

#KNN SFFS_KNN_result
SFFS_KNN_result = SFFS(dataset = trainingSet, algorithm = makeLearner('classif.kknn', k=3))
analyzeFeatSelResult(SFFS_KNN_result)


#Testset accuracy
trainTransformed <- trainingSetnormal
#LDA
trainingSet<-trainTransformed[index, ]
testSet <-trainTransformed[-index, ]
trainingSet<-trainingSet[,c('Name',SFFS_LDA_result$x)]
testSet <-testSet[,c('Name',SFFS_LDA_result$x)]
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
trainingSet<-trainingSet[,c('Name',SFFS_SVM_result$x)]
testSet <-testSet[,c('Name',SFFS_SVM_result$x)]
data.fit <- knn3(factor(Name) ~ ., trainingSet, k = 3)
data.fit
ions.ldapredict=predict(data.fit,testSet[,-1],type = "class")
ions.ldapredict
table(testSet$Name,ions.ldapredict)
count_accuracy(ions.ldapredict,testSet)
#SVM
trainingSet<-trainTransformed[index, ]
testSet <-trainTransformed[-index, ]
trainingSet<-trainingSet[,c('Name',SFFS_KNN_result$x)]
testSet <-testSet[,c('Name',SFFS_KNN_result$x)]
ions.svm=svm(Name ~ .,data=trainingSet, kernel="linear")
summary(ions.svm)
ions.svmpredict=predict(ions.svm, testSet[,-1],decision.values=T)
ions.svmpredict
table(testSet$Name,ions.svmpredict)
count_accuracy(ions.svmpredict,testSet)

#Repeated cross-validation
accuracyldaall=0
accuracyknnall=0
accuracysvmall=0
for(i in 1:30){
  index<-createDataPartition(trainingSetnormal$Name, time=1, p=0.625, list=F)
  #LDA
  trainingSet<-trainingSetnormal[index, c('Name', SFFS_LDA_result$x)]
  testSet <-trainingSetnormal[-index, c('Name', SFFS_LDA_result$x)]
  ions.lda=lda(formula = Name ~ ., data=trainingSet)
  ions.ldapredict=predict(ions.lda,testSet[,-1],decision.values=T)
  AccuracyLDA = count_accuracy(ions.ldapredict$class,testSet)
  accuracyldaall = accuracyldaall+AccuracyLDA
  #KNN
  trainingSet<-trainTransformed[index, c('Name', SFFS_KNN_result$x)]
  testSet <-trainTransformed[-index, c('Name', SFFS_KNN_result$x)]
  data.fit <- knn3(factor(Name) ~ ., trainingSet, k = 3)
  ions.knnpredict=predict(data.fit,testSet[,-1],type = "class")
  Accuracyknn = count_accuracy(ions.knnpredict,testSet)
  accuracyknnall = accuracyknnall+Accuracyknn
  #SVM  
  trainingSet<-trainTransformed[index, c('Name', SFFS_SVM_result$x)]
  testSet <-trainTransformed[-index, c('Name', SFFS_SVM_result$x)]
  ions.svm=svm(Name ~ .,data=trainingSet, kernel="linear")
  ions.svmpredict=predict(ions.svm,testSet[,-1],decision.values=T)
  Accuracysvm = count_accuracy(ions.ldapredict$class,testSet)
  accuracysvmall = accuracysvmall+Accuracysvm
}
accuracyldaall/30
accuracyknnall/30
accuracysvmall/30
