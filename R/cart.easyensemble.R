#' CART for binary outcomes using easy ensemble with prediction for new data
#'
#' Fit of CART using easy ensemble - a method that uses downsized class balanced training sets and boosting - with class prediction for new data.  
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.data.test new data; samples by rows, variables by columns. Must have the same number of variables as my.data.train, the variables must be in the same order.
#' @param y.test class membreship in test set
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param num.boost number of boosting iterations
#' @param num.ee number of easy ensemble iterations
#' @param replace.ee logical indicator (TRUE/FALSE) indicating if the downsized training sets should be obtained with resampling with replacement (TRUE) or not (FALSE). Note that only resampling without replacement is currently implemented.
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list containing the predictions on new data (pred.test), the true class membership of new data (y.test), the overall accuracy (accuracy.all), the class specific predictive accuracies (accuracy.class1 and accuracy.class0 for classes y=1 and y=0, respectively) and the boosting score (y.score) 
#' \item{prediction.test}{class predictions for new samples}
#' \item{class.test}{true class membreship for new samples}
#' \item{score.test.class1}{scores assigned to each new sample, evaluated using boosting and easy ensemble - FIX THIS >>>> }
#' \item{accuracy.test.all}{overall predictive accuracy}
#' \item{accuracy.test.class1}{predictive accuracy for class 1}
#' \item{accuracy.test.class0}{predictive accuracy for class 0}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(100*1000), ncol=1000)
#' y=rbinom(100, size=1, prob=0.7)
#' my.data.test=matrix(rnorm(10*1000), ncol=1000)
#' y.test=rbinom(10, size=1, prob=0.7)
#' pred.ee=cart.easyensemble(my.data.train, y, my.data.test, y.test, max.depth=1)



######### function that performs regular boosting on multiple down-sized data sets, with no internal CV for the evaluation of the error on the training set
########## the output is the prediction on an indepdendent test set
#return(list(pred=y.pred, class.test=y.test, accuracy=mean(y.pred==y.test),y.score=my.scores.final))

cart.easyensemble <-
function(my.data.train, y,  my.data.test, y.test=NULL, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			#nfolds=10, 
			num.boost=9, num.ee=9, replace.ee=FALSE, check.data=TRUE, verbose=TRUE){



######### function that performs regular easyensemble, with no CV

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



####################### step 00: check data #############################
			
			#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
				.check.test.data(my.data.test, y.test, my.data.train)
				#added 24/4: stop EE if the classes are balanced
				if(sum(y==1)==sum(y==0)) stop("Easy ensemble is a method for class-imbalanced data sets. In your data set the number of samples in both classes is equal. Use cart() or cart.boosting() functions to derive the classification rule.") 
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
			if(verbose==TRUE) cat("Fitting CART on the complete data\n")
		
			results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		



			############################# step 2: obtain predictions on the training set (not cross-validated) ############################
			prediction.cv=my.predict.cart(results.cart, my.data.train)
			########################## ############################################################## ##########################		


#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))


 #maxerror <- min(1 - max(summary(y))/sum(summary(y)), 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance! 
 #maxerror <- min(1 - max(table(y))/num.samples, 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance!   

 
 
			#number of samples in test set
			 num.samples.test=nrow(my.data.test)


 
 
 
 ############## step 3: perform Easy Ensamble ####################

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
				####### WARNING: needs to be coded differently to be able to use bootstrap samples, does not work at the moment
				my.folds=lapply(1:num.ee, function(i) sample(which.maj.class, abs(N1-N0), replace=replace.ee))
				######## 18/3/2013: reordered the folds to obtain a correct order of the predictions
				my.folds=lapply(my.folds, sort)

				############## end added for Easy Ensamble ####################

				 



				######################### beginning of the easy enemble step ##########################
				#### performs Easy Ensemble, my.folds is the list containing the IDs of the samples to be removed at each of the num.ee steps, each element of the list is ordered

				#the reordering of the prediction output is not necessary here, the predictions are obtained on the test set
				
				
			if(verbose==TRUE) {cat("Performing easy ensemble using", num.ee, "downsized data sets \n")
					my.env=environment()
					.internal.index=1
					}

				my.scores=matrix(unlist(lapply(my.folds, .cart.for.ee.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, verbose, my.env)), ncol=num.samples.test, byrow=T)

				#my.scores=matrix(unlist(lapply(my.folds, .cart.for.ee.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s)), ncol=num.samples.test, byrow=T)[order(unlist(my.folds)),]

				#.cart.for.ee.internal.with.arguments=function(my.fold, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost){

				#derive the global score for the test samples, combining the diff.alpha obtained from single down-sized data sets (downsizing+boosting)
				my.scores.final=apply(my.scores, 2, sum)

				y.pred=ifelse(my.scores.final>0, 1, 0)


				##### added 27/2/2013: checks if there are any scores equal to exactly 0
					#are there any ties?  if there are ties assign the class at random
					num.ties=sum(my.scores.final==0, na.rm=T)
					if(num.ties>0)  y.pred[which(my.scores.final==0)]=sample(c(0,1), num.ties, replace=T)


		if(verbose==T) rm(.internal.index)					


	return(list(prediction.test=y.pred, class.test=y.test, score.test.class1=my.scores.final, accuracy.test.all=mean(y.pred==y.test),  accuracy.test.class1=mean(y.pred[y.test==1]==y.test[y.test==1]), accuracy.test.class0=mean(y.pred[y.test==0]==y.test[y.test==0])))
}
