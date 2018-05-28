\name{nnet.coxph}
\alias{nnet.coxph}
\title{Feedforward neural network for right-censored survival data}
\description{Implements a simple (leaky) ReLU feedforward neural network using the partial likelihood at the output node during backprop to tune model weights thereby allowing for censoring}
\usage{
nnet.coxph(X, So, layers, strat = NULL, X.test=NULL, So.test=NULL, strat.test=NULL,
                       learning.rate=0.01, eta=1, num.epochs=1500, leaky=0,
                       weights = NULL, sig=1,
                       lam=0, keep.prob = NULL, adam=F, b1=0.9, b2=0.999, eps=1e-8, ss=1,
                       standardize=T, verbose=T, checkgrads=F, coxbench=F)
}
\arguments{
        \item{X}{A design matrix where each column is an observation (transpose of usual format in R)}
        \item{So}{A Surv() object}
        \item{layers}{A vector with the number of nodes in each layer including the terminal node which is one}
        \item{strat}{Index of different strata}
        \item{X.test}{Desgin matrix to estimate generalization accuracy during training}
        \item{So.test}{Surv() object to estimate generalization accuracy during training}
        \item{strat.test}{Strata index to estimate generalization accuracy during training}
        \item{learning.rate}{Gradient descent size}
        \item{eta}{Geometric decay rate for learning rate: rate <- eta * rate}
        \item{num.epochs}{Number of gradient steps (no SGD implementation)}
        \item{leaky}{Activation function is pmax(x,0) + pmin(x,0)*leaky)}
        \item{weights}{Values to initialize the weights at}
        \item{sig}{Standard deviation of noise for weight initialization}
        \item{lam}{Lambda penalty on L2-norm of the weights}
        \item{keep.prob}{Defaults to 1, otherwise implements dropout}
        \item{adam}{Boolean whether to use the Adam gradient descent optimizer}
        \item{b1/b2}{Hyperparameters for Adam optimizer}
        \item{eps}{Noise in denominator for Adam optimizer}
        \item{ss}{Seeding number}
        \item{standardize}{Boolean whether to standardize the features}
        \item{verbose}{Boolean whether to print output during training}
        \item{checkgrads}{Boolean whether to check backprop gradients match numerically derived estimates}
        \item{coxbench}{Boolean whether to print loss function for classical model}
}
\value{Returns a list with: weighted (weights), terminal function (aterminal), predicted risk scores (ahat), leaky parameter (leaky), column means and SDs (vecMu and vecSd) and training rate (cost.df).}
\author{Erik Drysdale}

\examples{
# library(stringr)
# library(glmnet)
# library(BoutrosLab.statistics.survival)
# # Generate survival data on a "hazard ball"
# set.seed(1)
# n.train <- 1000
# n.test <- 500
# x1 <- runif(n.train+n.test,-1,1)
# x2 <- sample(sign(x1))*sqrt(1 - x1^2)
# X <- cbind(x1,x2)
# eta0 <- 1 + 0.5*X[,1]^2 + 2*X[,2]^2
# t0 <- -log(runif(n.train+n.test))/exp(eta0)
# c0 <- rexp(n.train+n.test,rate=5)
# tobs <- ifelse(c0 < t0, c0, t0)
# delta <- ifelse(c0 < t0, 0, 1)
# So <- Surv(tobs,delta)
# So <- So[order(tobs)]
# X <- X[order(tobs),]
# # Training/test
# X.train <- X[1:n.train,]
# X.test <- X[-(1:n.train),]
# So.train <- So[1:n.train]
# So.test <- So[-(1:n.train)]
# # Fit glmnet
# mdl.train.glmnet <- cv.glmnet(x=X.train,y=as.matrix(So.train),family='cox')
# eta.glmnet <- as.vector(predict.cv.glmnet(mdl.train.glmnet,newx=X.test,s='lambda.min'))
# # Fit coxph
# mdl.coxph <- coxph(So.train ~ .,data=data.frame(X.train))
# eta.coxph <- as.vector(predict(mdl.coxph,newdata=data.frame(X.test)))
# # Fit neural network
# mdl.snet <- nnet.coxph(X=t(X.train),So=So.train,learning.rate = 0.01,
#                        layers = c(5,3,1),verbose = FALSE,
#                        num.epochs = 1500,adam = F,lam = 0)
# eta.snet <- as.vector(predict.snet(t(X.test),mdl.snet))
# # Compare accuracy
# acc.glmnet <- survConcordance(So.test ~ eta.glmnet)$conc
# acc.coxph <- survConcordance(So.test ~ eta.coxph)$conc
# acc.snet <- survConcordance(So.test ~ eta.snet)$conc
# cat(sprintf('Test C-index for: \n glmnet: %0.2f\n Cox-PH: %0.2f\n Neural-Net: %0.2f',
#               acc.glmnet,acc.coxph,acc.snet))
}
