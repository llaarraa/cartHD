#' CART for binary outcomes using easy ensemble with cross-validated prediction 
#'
#' Fit of CART using easy ensemble - a method that uses downsized class balanced training sets and boosting - with class cross-validated prediction for training data.  
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param nfolds number of folds to be used for cross-validation 
#' @param num.boost number of boosting iterations
#' @param num.ee number of easy ensemble iterations
#' @param replace.ee logical indicator (TRUE/FALSE) indicating if the downsized training sets should be obtained with resampling with replacement (TRUE) or not (FALSE). Note that only resampling without replacement is currently implemented.
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list with the cross-validated predictions on training data and additional information (see below)
#' \item{prediction.cv}{cross-validated class predictions for training samples}
#' \item{class.train}{true class membreship for training set samples}
#' \item{score.cv.class1}{cross-validated scores assigned to each training sample, evaluated using boosting and easy ensemble - samples with positive scores are classified in Class 1.}
#' \item{accuracy.cv.all}{overall cross-validated predictive accuracy}
#' \item{accuracy.cv.class1}{cross-validated predictive accuracy for class 1}
#' \item{accuracy.cv.class0}{cross-validated predictive accuracy for class 0}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(100*1000), ncol=1000)
#' y=rbinom(100, size=1, prob=0.7)
#' pred.ee=cart.easyensemble.cv(my.data.train, y, max.depth=1, nfolds=2, num.boost=5, num.ee=2)



######### function that performs regular boosting on multiple down-sized data sets, with no internal CV for the evaluation of the error on the training set
########## the output is the cross-validated prediction on the samples from the test set
#return(list(pred=y.pred, class.test=y.test, accuracy=mean(y.pred==y.test),y.score=my.scores.final))

cart.easyensemble.cv <-
function(my.data.train, y, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			nfolds=10, 
			num.boost=10, num.ee=10, replace.ee=FALSE, check.data=TRUE, verbose=TRUE){

######### function that performs regular boosting, with CV

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

##### replace.ee: TRUE or FALSE; if FALSE the sampling for easy ensemble is performed without replacement, while it is performed with replacement if set to TRUE 
######### WARNING: Does not work correctly if replace.ee=TRUE. it does not use bootstrap samples, bug that needs to be fixed.

################ function that performs easy ensemble
############ 13/11 and performs boosting, 
############ different from the original CART/CV function because it re-uses some quantities calculated for the original cart


################ uses internal functions: 
#.cart.internal.with.arguments
#my.predict.cart
#.cart.for.ee.internal.with.arguments

##################### end of function declarations #############################################################



######################### beginning of the function #####################################

		####################### step 00 : check data ############################
			
			#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
						
				}




################## step 0: calculate some quantities needed to fit cart (multiple times)

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



			########################## step 1:  fit cart on original data ##########################		
			########### necessary?? is it used?? 
			results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		




			############################# step 2: obtain predictions on the training set (not cross-validated!) ############################
			prediction.cv=my.predict.cart(results.cart, my.data.train)
			########################## ############################################################## ##########################		


#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))


 #maxerror <- min(1 - max(summary(y))/sum(summary(y)), 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance! 
 #maxerror <- min(1 - max(table(y))/num.samples, 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance!   

 
 
			#number of samples in test set
			#num.samples.test=nrow(my.data.test)


 
 
 ################## step 3: perform Cross-validation  with embedded Easy Ensemble step
 
	#obtained folds used for CV
 	if(nfolds<num.samples) {my.folds.cv=.balanced.folds(y, nfolds) 
							my.folds.cv=lapply(my.folds.cv, sort)}	else my.folds.cv=as.list(1:num.samples)

	if(verbose==TRUE) {cat("Performing cross-validation on multiple downsized training sets, using", nfolds ,"folds \n")
						my.env=environment()
						#.internal.index<<- 1
						.internal.index=1
						}
 
		

	#y.pred=unlist(lapply(my.folds.cv, .cart.for.cv.ee.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.ee, num.boost, verbose))[order(unlist(my.folds.cv))]
	#updated 7/5/2013: returns not only the class predictions but also the scores used to generate them
	y.pred=lapply(my.folds.cv, .cart.for.cv.ee.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.ee, num.boost, verbose, my.env)
	
	scores.cv=unlist(lapply(y.pred, function(x) x$scores.final))[order(unlist(my.folds.cv))]
	
	y.pred=unlist(lapply(y.pred, function(x) x$y.pred))[order(unlist(my.folds.cv))]
		
	if(verbose==TRUE) rm(.internal.index)	

	return(list(prediction.cv=y.pred,  class.train=y,  score.cv.class1=scores.cv, 
	accuracy.cv.all=mean(y.pred==y),  accuracy.cv.class1=mean(y.pred[y==1]==y[y==1]), accuracy.cv.class0=mean(y.pred[y==0]==y[y==0])))
}
