#######set directory
setwd('D:\\QSAR workflow')
#######Load package
library(RWeka)
library(caret)
library(randomForest)
library(kernlab)
library(corrplot)
library(C50)
library(nnet)
library(e1071)
library(GA)
library(cvTools) 
library(Metrics)
library(MASS)
library(pls)
library(Interpol)
library(protr)
library(seqinr)
library(Peptides)
library(AUC)
library(ROCR)
library(prospectr)

x1 <-read.fasta('AFP.fasta', seqtype="AA", as.string = TRUE)
x2<- read.fasta('nonAFP.fasta', seqtype="AA", as.string = TRUE)

x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]
label = c(rep("AFP", length(x1)),rep("nonAFP", length(x2)))

AAC1 <- t(sapply(x1, extractAAC))
AAC2 <- t(sapply(x2, extractAAC))
AAC = rbind(AAC1,AAC2)

df = data.frame(scale(AAC),Class = label)
data = df[!duplicated(df), ]
data2 = data[,!sapply(data,function(x) any(is.na(x)))]
                         
################ Random
Pos = subset(data2, Class == 'AFP')
Neg = subset(data2, Class == 'nonAFP')

#######  Dividing Training and Testing sets on positive and negative classes
sample1 <- c(sample(1:nrow(Pos) ,300))
sample2 <- c(sample(1:nrow(Neg),300))
  train1  <- Pos[sample1,] ####Positive set for training
  train2  <- Neg[sample2,] ####Negative set for training
  test1 <-   Pos[-sample1,]    ####Positive set for testing
  test2 <-   Neg[-sample2,]    ####Negative set for testing 
  
  TR <- rbind(train1,train2) ####combining for internal set
  TS <- rbind(test1,test2)    ####combining for external set

################ Kennard-Stone algorithm
  
npos = round(nrow(Pos)*0.75)
nneg = round(nrow(Neg)*0.75)                                     
kenX <- kenStone(Pos[,-ncol(Pos)],k=npos,metric='euclid',pc=5)
kenY <- kenStone(Neg[,-ncol(Neg)],k=nneg,metric='euclid',pc=5)

internal = rbind(Pos[kenX$model,], Neg[kenX$model,])                     
external = rbind(Pos[kenY$test,], Neg[kenY$test,])         
                      
                      
