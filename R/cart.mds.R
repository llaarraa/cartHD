#' CART for binary outcomes using multiple downsizing with prediction on new data
#'
#' Fit of CART on multiple downsized (class-balanced) training sets and prediction on new data based on majority voting. 
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.data.test new data; samples by rows, variables by columns. Must have the same number of variables as my.data.train, the variables must be in the same order.
#' @param y.test class membreship in test set
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param num.mds number of downsized data sets to be used 
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list containing the predictions on new data and some accuracy measures
#' \item{prediction.test}{class predictions for new samples}
#' \item{class.test}{true class membreship for new samples}
#' \item{score.test.class1}{scores assigned to each new sample, evaluated as the average number of times that the test sample was assigned to class 1 using downsized training sets}
#' \item{accuracy.test.all}{overall predictive accuracy}
#' \item{accuracy.test.class1}{predictive accuracy for class 1}
#' \item{accuracy.test.class0}{predictive accuracy for class 0}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @note The downsized training sets are obtained randomly removing some of the majority class samples as to obtain class balanced training sets. The class membership of the test set samples is obtained using majority voting: the samples are classified in the class to which they were assigned most frequently by the CART models fitted using the downsized training sets. The scores (scores.test) are the proportion of times that a sample was classified in class 1. The class membership is assigned randomly in case of ties (score=0.50, possible outcome if the number of downsized training sets (num.mds) is an even number). 
#' @export
#' @seealso \code{\link{cart.mds.cv}}
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(30*1000), ncol=1000)
#' y=rbinom(30, size=1, prob=0.7)
#' my.data.test=matrix(rnorm(10*1000), ncol=1000)
#' y.test=rbinom(10, size=1, prob=0.5)
#' pred.mds=cart.mds(my.data.train, y, my.data.test, y.test)



			
######### function that performs multiple downsizing with CART
#return(list(pred=y.pred, class.test=y.test, error=mean(y.pred==y.test)))

cart.mds <-
function(my.data.train, y,  my.data.test, y.test, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			#nfolds=10, 
			#num.boost=10, 
			num.mds=9, check.data=TRUE, verbose=TRUE){


			
######### function that performs multiple downsizing with CART


####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
####max.depth: maximum number of levels of the tree
#####n.folds - not relevant here: number of folds to be used in CV
#num.boost  - not relevant here: number of boosting iterations

####### my.data.test: test data, will be used during the boosting step
####### y.test: class membreship in test set
###num.mds: number of downsized datasets to combine


################ function that performs multiple downsizing for cart
############ 
############ 



######################### beginning of the function #####################################



####################### step 00: check data #############################
			
			#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
				.check.test.data(my.data.test, y.test, my.data.train)
				
				}


############################## step 0: calculates some quantities needed for CART estimation


			#### saving the number of samples and number of variables
			#number of samples
			N=num.samples=nrow(my.data.train)


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



		#number of samples in test set
		 num.samples.test=nrow(my.data.test)




########################## fit cart on original data ##########################		
########### necessary??
#results.cart=cart.internal()
#results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
########################## fit cart on original data ##########################		
		

############################# perform cross-validation ############################

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
#prediction.cv=my.predict.cart(results.cart, my.data.train)



#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))


 #maxerror <- min(1 - max(summary(y))/sum(summary(y)), 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance! 
 #maxerror <- min(1 - max(table(y))/num.samples, 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance!   

 
 






######################### beginning of the multiple downsizing step ##########################


############## added for multiple downsizing ####################

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
my.folds=lapply(1:num.mds, function(i) sample(which.maj.class, abs(N1-N0)))
#9/5/2013: sort the indexes within the folds
my.folds=lapply(my.folds, sort)
############## end added for multiple downsizing ####################



#the re-ordering of the ouput is not necessary, as the predictions are obtained on the independent test sets

	#if the MDS progress must be printed on screen, an index is initialized
			if(verbose==TRUE) {cat("Performing multiple downsizing, using", num.mds , "downsized training sets. \n")
								my.env=environment()
								#.internal.index=1						
								.internal.index=1						
									}



#my.predictions.new=(matrix(unlist(lapply(my.folds, .cart.for.mds.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, my.data.test, verbose)), ncol=num.samples.test, byrow=T))
my.predictions.new=(matrix(unlist(lapply(my.folds, .cart.for.mds.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, my.data.test, my.env, verbose)), ncol=num.samples.test, byrow=T))

#derive the global score for the test samples, combining the predictions from the downsized analyses
my.scores.final=apply(my.predictions.new, 2, sum)

y.pred=ifelse(my.scores.final>num.mds/2, 1, 0)


##### added 27/2/2013: checks if there are any scores equal to exactly 0
	#are there any ties?  if there are ties assign the class at random
	num.ties=sum(my.scores.final==num.mds/2, na.rm=T)
	if(num.ties>0)  y.pred[which(my.scores.final==num.mds/2)]=sample(c(0,1), num.ties, replace=T)


return(list(prediction.test=y.pred, class.test=y.test, score.test.class1=my.scores.final/num.mds, accuracy.test.all=mean(y.pred==y.test), accuracy.test.class1=mean(y.pred[y.test==1]==y.test[y.test==1]), accuracy.test.class0=mean(y.pred[y.test==0]==y.test[y.test==0]) ))

}
