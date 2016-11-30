library(R.matlab)
library(nnls)
library(glmnet)
setwd("~/Research/Fall 2016 circuit");

# dataMat<-readMat('forR_som_cell12loc1all.mat');
dataMat<-readMat('test.mat');
design1<-dataMat$design1;
design2<-dataMat$design2;
Y<-dataMat$Y;

########################
net.fit<-glmnet(cbind(design1[,2]*-1,design2),Y,lower=0,family='poisson');
betafit<-coef(net.fit,0);
betafit[2]<-betafit[2]*-1;

betafit

write(betafit[,1], file = "betafit.txt", sep = " ");
