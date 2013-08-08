#' CART for binary outcomes using boosting with internal cross-validation and prediction of the class membership for new data
#'
#' Function that fits CART using boosting with internal cross-validation for error estimation (within boosting) and produces cross-validated class predictions for the traning data
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param nfolds number of folds to be used for cross-validation
#' @param num.boost number of boosting iterations
#' @param nfolds.Icv number of folds to be used for internal cross-validation (within boosting for the assessment of predictive accuracy)
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list with the cross-validated predictions on training data and additional information (see below)
#' \item{prediction.cv}{cross-validated class predictions for training samples obtained using boosting}
#' \item{class.train}{true class membreship for training set samples}
#' \item{score.cv.class1}{cross-validated scores assigned to each training sample, evaluated using boosting - samples with positive scores are classified in Class 1.}
#' \item{accuracy.cv.all}{overall cross-validated predictive accuracy}
#' \item{accuracy.cv.class1}{cross-validated predictive accuracy for class 1}
#' \item{accuracy.cv.class0}{cross-validated predictive accuracy for class 0}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples 
#' set.seed(2)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000)
#' y=rbinom(10, size=1, prob=0.5)
#' pred.boosting.Icv.cv=cart.boosting.Icv.cv(my.data.train, y, nfolds=2)


			
######### function that performs boosting with CART and produces cross-validated predictions on the training set - it does not use new data!
#return(list(pred=y.pred, class.train=y, accuracy=mean(y.pred==y)))

cart.boosting.Icv.cv <-
function(my.data.train, y,  
			#my.data.test, y.test, 
			my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			nfolds=5, 
			num.boost=10,#, 
			nfolds.Icv=5,
			check.data=TRUE, verbose=TRUE){


			
######### function that performs multiple downsizing with CART

####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
####max.depth: maximum number of levels of the tree
#####n.folds -  number of folds to be used in CV
#num.boost  - number of boosting iterations

####### my.data.test: not relevant here: test data, will be used during the boosting step
####### y.test: not relevant here: class membreship in test set

###num.mds: not relevant here: number of downsized datasets to combine


################ function that performs multiple downsizing for cart
############ 
############ 



######################### beginning of the function #####################################


############################## step 00: check the data

		#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
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



			## derives the ranks - to backtrasform the sorted data into the original data
			#my.ranks=t(colRanks(my.data.train))

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

 
 
#number of samples in test set
 #num.samples.test=nrow(my.data.test)






######################### step 1: beginning of the cross-validation step ##########################

#18/3/2013: the IDs of the samples to leave out are reordered within the folds, necessary to obtain a correct order of the predictions
	
#10/6/2013: for debugging
#set.seed(1)	
	if(nfolds<num.samples) {my.folds.cv=.balanced.folds(y, nfolds) 
							my.folds.cv=lapply(my.folds.cv, sort)}	else {my.folds.cv=as.list(1:num.samples)
																	nfolds=num.samples
																	}




#my.predictions.new=matrix(unlist(lapply(my.folds, .cart.for.mds.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, my.data.test)), ncol=num.samples.test, byrow=T)

if(verbose==TRUE) {cat("Performing boosting with cross-validation, using", nfolds, "folds \n")
					my.env=environment()
					.internal.index=1
				}

#changed 7/5/2013: the new version returns the cross-validated classification scores and CV-class membership to the main function
#y.pred=unlist(lapply(my.folds.cv, .cart.for.cv.boosting.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ind.s, num.boost, verbose))[order(unlist(my.folds.cv))]
#return(list(pred.cv=prediction.new.data.boosting, score.cv=diff.alpha))

#8/5/2013: defined the new function that uses internal CV for the evaluation of the error within boosting, nfolds.Icv is an additional argument compared to the non-Icv function
y.pred=lapply(my.folds.cv, .cart.for.cv.boosting.Icv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ind.s, num.boost, nfolds.Icv, verbose, my.env)

score.cv=unlist(lapply(y.pred, function(x) x$score.cv))[order(unlist(my.folds.cv))]
y.pred=unlist(lapply(y.pred, function(x) x$pred.cv))[order(unlist(my.folds.cv))]


#.cart.for.cv.mds.internal.with.arguments=function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.mds){



#derive the global score for the test samples, combining the 
#my.scores.final=apply(my.predictions.new, 2, sum)

#y.pred=ifelse(my.scores.final>num.mds/2, 1, 0)


##### added 27/2/2013: checks if there are any scores equal to exactly 0
	#are there any ties?  if there are ties assign the class at random
#	num.ties=sum(my.scores.final==num.mds/2, na.rm=T)
#	if(num.ties>0)  y.pred[which(my.scores.final==num.mds/2)]=sample(c(0,1), num.ties, replace=T)




return(list(prediction.cv=y.pred, class.train=y, score.cv.class1=score.cv,
		accuracy.cv.all=mean(y.pred==y), accuracy.cv.class1=mean(y.pred[y==1]==y[y==1]), accuracy.cv.class0=mean(y.pred[y==0]==y[y==0])))

}
