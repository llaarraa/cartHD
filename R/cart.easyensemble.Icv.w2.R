#' CART for binary outcomes using easy ensemble with internal cross-validation and prediction for new data
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.data.test new data; samples by rows, variables by columns. Must have the same number of variables as my.data.train, the variables must be in the same order.
#' @param y.test class membreship in test set
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param nfolds.Icv number of folds to be used in the internal cross-validation (used to assess the error within boosting)
#' @param num.boost number of boosting iterations
#' @param num.ee number of easy ensemble iterations
#' @return list list containing the predictions on new data (pred), the true class membership of new data (y.test), the overall accuracy (accuracy) and the boosting score (y.score) 
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(20*1000), ncol=1000)
#' y=rbinom(20, size=1, prob=0.7)
#' my.data.test=matrix(rnorm(20*1000), ncol=1000)
#' y.test=rbinom(20, size=1, prob=0.7)
#' pred.mds=cart.easyensemble.Icv(my.data.train, y, my.data.test, y.test, nfolds=2)

######### function that performs regular easy ensemble (boosting), with internal CV used to evaluate the prediction error on the training set 
####### the output is the prediction on an independent test set


cart.easyensemble.Icv.w2 <-
function(my.data.train, y,  my.data.test, y.test, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			nfolds.Icv=5, 
			num.boost=10, num.ee=9){




			
######### function that performs regular easy ensemble (boosting), with internal CV used to evaluate the prediction error on the training set 


####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
####max.depth: maximum number of levels of the tree
#####n.folds: number of folds to be used in CV
#num.boost: number of boosting iterations

####### my.data.test: test data, will be used during the boosting step
####### y.test: class membreship in test set
###num.ee: number of easy ensemble iterations


################ function that performs easy ensemble
############ 13/11 and performs boosting, 
############ different from the original CART/CV function because it re-uses some quantities calculated for the original cart




######################### beginning of the function #####################################



########################## step 0:  initializations, calculate quantities needed for cart and EE ###############################

						#### saving the number of samples and number of variables
						#number of samples
						N=num.samples=nrow(my.data.train)


						############## added for Easy Ensamble ####################

						#### sample size from each class
						#class with y=1
						N1=sum(y==1)
						#class with y=0
						N0=N-N1

						#which is the majority class, 1 or 0? should it give a warning if the sample size is exactly the same??
						maj.class=ifelse(N1>N0, 1, 0)

						#########which samples belong to the majority class
						which.maj.class=which(y==maj.class)

						#########how many from the majority class have to be removed
						#

						#### list with the sample IDs to remove from the majority class to obtain a balanced training set
						##### each element of the list (with num.ee elements) contains the IDs to REMOVE! to be used with the function cart.for.ee.internal
						
						#10/6/2013: added for debugging
						#set.seed(1)
	
						my.folds=lapply(1:num.ee, function(i) sample(which.maj.class, abs(N1-N0)))
						my.folds=lapply(my.folds, sort)	
						
						############## end added for Easy Ensamble ####################


						###################### quantities needed for cart ##################################


						#number of variables
						num.var=ncol(my.data.train)


						#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
						my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var)

						#sorting the training data using the calculated indexes
						#sorted training data set, each variable is sorted
						my.tr.s=matrix(unlist(lapply(1:num.var, function(i) my.data.train[,i][my.ind.s[,i]])), ncol=num.var)







						#weights in matrix format
						my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
						#class membership in matrix format
						my.y.s<-matrix(y[my.ind.s], ncol=num.var)


						#weighted relative frequencies, for class 1 and class 2
						n.l<-apply(my.w.s,2,cumsum)
						n.r<-1-n.l

						##### utility matrices: calculated once and reused later
						my.y.0=my.y.s==0
						my.w.times.y.0=my.w.s*my.y.0
						my.w.times.y.s=my.w.s*my.y.s
						## derives the ranks - to backtrasform the sorted data into the original data, uses a function from library(matrixStats)
						#my.ranks=t(colRanks(my.data.train))
						my.ranks=apply(my.data.train, 2, rank, ties="first")

						#defining the prior if it was not specified by the user, using the empirical prior; 
						#the empirical prior is equal to the WEIGHTED proporion of cases in each class
						#order for my.prior: class 1, class 0
						if(is.null(my.prior)) my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

						############## to include the prior in the calculation of the gini gain
						N.l.w=sum(my.weights[y==1])
						N.r.w=1-N.l.w







				########################## step 1: fit cart on original data ##########################		
				########### necessary??
				results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
														 
				########################## end fit cart on original data ##########################		
		

				############################# step 2: perform cross-validation ############################

				############### obtain balanced folds
				#the output is a list contaning the sample ids to be used in each fold

				#my.folds=.balanced.folds(y, nfolds)
				########## add an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list
				#### deleted, not necessary for regular boostin
				#if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

				################ obtain cross-validated prediction using nfolds-CV
				##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
				#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]
				#### modified to obtain regular boosting, does not use CV!
				########### necessary??
				prediction.cv=my.predict.cart(results.cart, my.data.train)



#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))



 
				#number of samples in test set
				 num.samples.test=nrow(my.data.test)






######################### beginning of the easy enemble step ##########################

#### performs Easy Ensemble with internal cross-validation to evaluate the error in the boosting step, my.folds is the list containing the IDs of the samples to be removed at each of the num.ee steps

			
		#	debug(cart.for.ee.cv.internal)
		#debug(.cart.for.cv.internal.with.arguments)

my.scores=matrix(unlist(lapply(my.folds, .cart.for.ee.cv.internal.w2.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, nfolds.Icv)), ncol=num.samples.test, byrow=T)

#derive the global score for the test samples, combining the 
my.scores.final=apply(my.scores, 2, sum)

y.pred=ifelse(my.scores.final>0, 1, 0)


##### added 27/2/2013: checks if there are any scores equal to exactly 0
	#are there any ties?  if there are ties assign the class at random
	
	#10/6/2013: added for debugging
	#set.seed(1)
	
	num.ties=sum(my.scores.final==0, na.rm=T)
	if(num.ties>0)  y.pred[which(my.scores.final==0)]=sample(c(0,1), num.ties, replace=T)




return(list(pred=y.pred, class.test=y.test, accuracy.test.all=mean(y.pred==y.test), accuracy.test.class1=mean(y.pred[y.test==1]==y.test[y.test==1]), accuracy.test.class0=mean(y.pred[y.test==0]==y.test[y.test==0]), y.score=my.scores.final))

}
