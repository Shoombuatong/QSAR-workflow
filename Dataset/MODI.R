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

x1 <- read.fasta('AFP.fasta', seqtype="AA", as.string = TRUE)
x2<- read.fasta('nonAFP.fasta', seqtype="AA", as.string = TRUE)

x1 <- x1[(sapply(x1, protcheck))]
x2 <- x2[(sapply(x2, protcheck))]
label = c(rep("AFP", length(x1)),rep("nonAFP", length(x2)))

AAC1 <- t(sapply(x1, extractAAC))
AAC2 <- t(sapply(x2, extractAAC))
AAC = rbind(AAC1,AAC2)

df = data.frame(scale(AAC),Class = label)
data = df[!duplicated(df), ]
internal = data[,!sapply(data,function(x) any(is.na(x)))]

############# MODI for classification task
Dat <- internal
n = ncol(Dat)
m= n-1
AA <- Dat[,1:m]
class = as.numeric(Dat[,ncol(Dat)])

d1 <- dist(AA, upper=TRUE, diag=TRUE, method = "euclidean")
nd1 <- scale(d1)
nd2 = ((nd1-min(nd1))/(max(nd1)-min(nd1)))
MOBI <- matrix(nrow = nrow(Dat), ncol = 1)

for (i in 1:nrow(Dat)){
MOBI[i,] <- Dat[order(nd2[i,]),][2,n]
}

result = data.frame(class,MOBI)
X <- subset(result, result[,1] == '1')
Y <- subset(result, result[,1] == '2')

MODIclass = (table(X)[1]/nrow(X)+ table(Y)[2]/nrow(Y))/2

##### MODI > 0.65 
############# MODI for classification task
feature <- read.csv("Model.csv", header = TRUE) 
label <- read.csv("Koc.csv", header = TRUE) 
Dat = data.frame(feature, label)
n = ncol(Dat)
m= n-1
AA <- Dat[,1:m]

d1 <- dist(AA, upper=TRUE, diag=TRUE, method = "euclidean")
nd1 <- scale(d1)
nd2 = ((nd1-min(nd1))/(max(nd1)-min(nd1)))
MOBI3knn <- matrix(nrow = nrow(Dat), ncol =3)
MOBI5knn <- matrix(nrow = nrow(Dat), ncol =5)
MOBI3knn_aver <- matrix(nrow = nrow(Dat), ncol =1)
MOBI5knn_aver <- matrix(nrow = nrow(Dat), ncol =1)

for (i in 1:nrow(Dat)){
MOBI3knn[i,] <- Dat[order(nd2[i,]),][1:3,n]
MOBI5knn[i,] <- Dat[order(nd2[i,]),][1:5,n]
}

for (i in 1:nrow(Dat)){
MOBI3knn_aver[i,] <- mean(MOBI3knn[i,])
MOBI5knn_aver[i,] <- mean(MOBI5knn[i,])
}

result = data.frame(MOBI3knn_aver,MOBI5knn_aver,Dat[,n])

MODI_q2_3knn = cor(result[1],result[3])*cor(result[1],result[3])
MODI_q2_5knn = cor(result[2],result[3])*cor(result[2],result[3])

MODI_q2 = data.frame(MODI_q2_3knn,MODI_q2_5knn)

##### MODI (3knn) > 0.46 and MODI (5knn) > 0.47

