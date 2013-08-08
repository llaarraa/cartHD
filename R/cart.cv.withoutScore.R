#' CART for binary outcomes with cross-validated class predictions	
#' REMOVE FROM THE PACKAGE
#' Old function for cross-validation, it does not produce the cross-validated class prediction scores.
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param nfolds number of folds to be used in k-fold cross-validation
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list containing the CART fit obtained on the original data and the cross-validated predictions
#' \item{Results.cart}{CART object containing the fit obtained on the complete data}
#' \item{prediction.cv}{cross-validated class predictions}
#' \item{accuracy.cv.all}{overall cross-validated predictive accuracy}
#' \item{accuracy.cv.class1}{cross-validated predictive accuracy for class 1}
#' \item{accuracy.cv.class0}{cross-validated predictive accuracy for class 0}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000)
#' y=rbinom(10, size=1, prob=0.5)
#' my.cart.fit=cart.cv(my.data.train, y, nfolds=2)


cart.cv.without.score<-
function(my.data.train, y, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  max.depth=5, nfolds=10, check.data=TRUE, verbose=TRUE){


################ function that performs k-fold cross validation after estimating CART


####### arguments ###########

####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
####max.depth: maximum number of levels of the tree


############## uses internal functions:
#.my.stump.forCART.oneStep.stepOne
#.my.stump.forCART.oneStep.stepNotOne
#.cart.internal.with.arguments
#.balanced.folds
#.cart.for.cv.internal.with.arguments

################################ beginning of the function ########################################

########### step 00: check the input data for correctness 
			#run checks on data only if the functions are not used for simulations		

		
		if(verbose==TRUE) cat("Checking data\n")

		if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
				}




########### step 0: calculate quantities needed for the cart fit


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
				#set.seed(1)
				my.ranks=apply(my.data.train, 2, rank, ties="first")

				
				#weighted relative frequencies, for class 1 and class 2

				n.l<-apply(my.w.s,2,cumsum)
				n.r<-1-n.l

				##### utility matrices: calculated once and reused later
				my.y.0=my.y.s==0
				my.w.times.y.0=my.w.s*my.y.0
				my.w.times.y.s=my.w.s*my.y.s
	
				#defining the prior if it was not specified by the user, using the empirical prior; 
				#the empirical prior is equal to the WEIGHTED proporion of cases in each class
				#order for my.prior: class 1, class 0
				if(is.null(my.prior)) my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

				############## to include the prior in the calculation of the gini gain
				N.l.w=sum(my.weights[y==1])
				N.r.w=1-N.l.w




########### step 1: obtain cart on the complete data set ###############


		if(verbose==TRUE) cat("Fitting CART on complete data\n")
		results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)

	

#results.cart=cart(my.data.train, y, my.weights, my.prior, min.samples, min.gain,  max.depth) 


############################# step 2: begin cross-validation ############################


############### obtain balanced folds
#the output is a list contaning the sample ids to be used in each fold

########## add an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list

		if(verbose==TRUE) cat("Performing cross-validation\n")


		if(nfolds<num.samples) {my.folds=.balanced.folds(y, nfolds) 
								my.folds=lapply(my.folds, sort)} else my.folds=as.list(1:num.samples)

		################ obtain cross-validated prediction using nfolds-CV
		##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
		#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]
		prediction.cv=unlist(lapply(my.folds, .cart.for.cv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds))]



		#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
		#return(list(results.cart=results.cart, prediction=prediction.cv))


			
		names(results.cart)=unlist(lapply(1:length(results.cart), function(k) .convert.node(k, results.cart[[k]]$depth )))

		#for the output: add the vector with the class membership of the training set and define the class of the object (class: cart)
		Results.cart=list(tree=results.cart, y=y)
			
			
			#define the class CART for the outout object	
			class(Results.cart)="CART"							




return(list(results.cart=Results.cart, prediction.cv=prediction.cv, accuracy.cv.all=mean(prediction.cv==y), accuracy.cv.class1=mean(prediction.cv[y==1]==y[y==1]), accuracy.cv.class0=mean(prediction.cv[y==0]==y[y==0]) ))




}
