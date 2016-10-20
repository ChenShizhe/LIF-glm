library(R.matlab)
library(nnls)
library(glmnet)
setwd("~/Research/Fall 2016 circuit");

dataMat<-readMat('forR.mat');
design1<-dataMat$design1;
design2<-dataMat$design2;
Y<-dataMat$Y;

# vlog <- function() {
#   ## link
#   linkfun <- function(y) log(exp(y)-1)
#   ## inverse link
#   linkinv <- function(eta)  log(exp(eta)+1)
#   ## derivative of invlink wrt eta
#   mu.eta <- function(eta) { 1/(exp(-eta) + 1) }
#   valideta <- function(eta) TRUE
#   link <- "log(exp(y)-1)"
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta, 
#                  name = link),
#             class = "link-glm")
# }
# vv <- vlog();
# 
# g<-glm(Y ~ ., family = poisson(link=vv),data=data.frame(cbind(design1[,2]*-1,design2)));
# 
# 
# vstep <- function() {
#   ## link
#   linkfun <- function(y){linkfun=matrix(0,nrow(y),ncol(y));linkfun[y>=1]=y[y>=1]-1;linkfun[y<1]=log(y[y<1]);}
#   ## inverse link
#   linkinv <- function(eta){
#     if (eta>=0){1+eta}
#     else{exp(eta)}
#   }
#   ## derivative of invlink wrt eta
#   mu.eta <- function(eta){
#     if (eta>=0){1}
#     else{log(eta)}
#   }
#   valideta <- function(eta) TRUE
#   link <- "log(exp(y)-1)"
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta, 
#                  name = link),
#             class = "link-glm")
# }
# vv <- vstep();
# 
# g<-glm(Y ~ ., family = poisson(link=vv),data=data.frame(cbind(design1[,2]*-1,design2)));


########################

# design=cbind(design1[,2],design2);
# 
# nnglm<-nnls(design,Y);
# coef.nn <- coef(nnglm)
# 
# 
# 
# irls =
#   function(A, b, family=poisson, maxit=25, tol=1e-08)
#   {
#     x = rep(0,ncol(A))
#     for(j in 1:maxit)
#     {
#       eta    = A %*% x
#       g      = family()$linkinv(eta)
#       gprime = family()$mu.eta(eta)
#       z      = eta + (b - g) / gprime
#       W      = as.vector(gprime^2 / family()$variance(g))
#       xold   = x
#       #x      = solve(crossprod(A,W*A), crossprod(A,W*z), tol=2*.Machine$double.eps)
#       x      = coef(nnls(sqrt(W)*A, sqrt(W)*b))
#       
#       if(sqrt(crossprod(x-xold)) < tol) break
#     }
#     list(coefficients=x,iterations=j)
#   }
# 
# g.fit<-glm(Y ~ ., family = poisson,,x=TRUE,data=data.frame(cbind(design1[,2]*-1,design2)));
# coef(g.fit)
# 
# irls.fit<-irls(g.fit$x,g$y,family=poisson);
# coef(irls.fit)



########################
net.fit<-glmnet(cbind(design1[,2]*-1,design2),Y,lower=0,family='poisson');
betafit<-coef(net.fit,0);

write(betafit[,1], file = "betafit.txt", sep = " ")