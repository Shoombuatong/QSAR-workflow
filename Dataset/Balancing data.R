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

x1 <-read.fasta('AFP.fasta', seqtype="AA", as.string = TRUE)
x2<- read.fasta('nonAFP.fasta', seqtype="AA", as.string = TRUE)

x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]
label = c(rep("AFP", length(x1)),rep("nonTHP", length(x2)))

AAC1 <- t(sapply(x1, extractAAC))
AAC2 <- t(sapply(x2, extractAAC))
AAC = rbind(AAC1,AAC2)

df = data.frame(scale(AAC),Class = label)
data = df[!duplicated(df), ]
internal = data[,!sapply(data,function(x) any(is.na(x)))]
                         
################ Random
Pos = read.csv("ALL-AFP.csv", header = TRUE) 
Neg = read.csv("ALL-nonAFP.csv", header = TRUE)

#######  Dividing Training and Testing sets on positive and negative classes
sample1 <- c(sample(1:nrow(Pos) ,300))
sample2 <- c(sample(1:nrow(Neg),300))
  train1  <- Pos[sample1,] ####Positive set for training
  train2  <- Neg[sample2,] ####Negative set for training
  test1 <-   Pos[-sample1,]    ####Positive set for testing
  test2 <-   Neg[-sample2,]    ####Negative set for testing 
  
  TR <- rbind(train1,train2) ####combining for internal set
  TS <- rbind(test1,test2)    ####combining for external set



           
