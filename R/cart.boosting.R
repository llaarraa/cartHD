#' CART with boosting for binary outcomes with class prediction for new data
#'
#' Fit of CART for binary outcomes using boosting, the prediction of class membership is obtained for new data. 
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.data.test new data; samples by rows, variables by columns. Must have the same number of variables as my.data.train, the variables must be in the same order.
#' @param y.test class membreship in test set
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param num.boost number of boosting steps
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @param verbose logical indicator, if TRUE the progress of the computations is printed out in the console
#' @return list containing the predicted class membership for new samples and additional information (see below) 
##### \item{results.cart}{CART object with CART fit obtained on the original data}
##### \item{prediction.cv}{resubstitution class predictions on the training set - NB: obtained without CV, not with boosting!}
#' \item{prediction.test.boosting}{class predictions for new samples obtained with boosting}
#' \item{class.test}{true class membreship for new samples - if provided by the user}
#' \item{accuracy.test.all}{overall predictive accuracy for new data - obtained using boosting}
#' \item{accuracy.test.class1}{predictive accuracy for class 1}
#' \item{accuracy.test.class0}{predictive accuracy for class 0}
#' \item{beta}{boosting weight, defined at the t-th boosting iteration as 1/2 log(PE_t/(1-PE_t)), where PE_t is the weighted prediction error of the t-th boosting iteration}
#' \item{weights.after.boosting}{weights for each training set sample updated by boosting}
#' \item{error.boosting}{proportion of misclassified samples at each step of boosting}
#' \item{score.test.class1}{classification score for test set samples evaluated by boosting, samples with positive values are classified in class 1}
#' \item{pred.new.i}{predicted class membreship for new samples at each boosting iteration, rows are samples, columns are boosting iterations}
#' \item{iterations.boosting}{number of boosting iterations effectively performed}
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples
#' set.seed(1) 
#' my.data.train=matrix(rnorm(10*1000), ncol=1000)
#' y=rbinom(10, size=1, prob=0.5)
#' my.data.test=matrix(rnorm(10*1000), ncol=1000)
#' y.test=rbinom(10, size=1, prob=0.7)
#' pred.boosting=cart.boosting(my.data.train, y, my.data.test, y.test)



############### function that fits cart and performs classic ADAboost, it does not use internal cross-validation to evaluate the error on the training set, the output is the prediction on a test set
############## modifies cart.boosting3 to save time needed to process the non-performed boosting iterations
#return(list(results.cart=results.cart.original, prediction=prediction.cv, prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#weights.after.boosting=my.weights, error.boosting=error.boosting, prob=diff.alpha, pred.new.i=prediction.new[,1:m], m=m))


cart.boosting <-
function(my.data.train, y, my.data.test, y.test=NULL, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  
			max.depth=5, 
			#nfolds=10, 
			num.boost=10, check.data=TRUE, verbose=TRUE ){


######### function that performs regular boosting, with no CV

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

################ function that performs k-fold cross validation after estimating CART
############ 13/11 and performs boosting, 
############ different from the original CART/CV function because it re-uses some quantities calculated for the original cart

########### 20/3/2013: simplified the calculation of the error for each boosting step, it is calculated only for the effectively used boosting steps

################# uses internal functions
#.cart.internal.with.arguments
#my.predict.cart



######################### beginning of the function #####################################

			####################### step 00 : check data ############################
			
			#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
				.check.test.data(my.data.test, y.test, my.data.train)
				
				}




##################### step 0:  calculate some quantities needed later to fit cart (multiple times, and on subsets) ################################


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



			########################## step 1: fit cart on original data ###############################################################	
			#results.cart=cart.internal()
			
			if(verbose==TRUE) cat("Fitting CART on complete data\n")
			results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)


			######################### step 2: obtain  predictions on the training data, predictions are NOT cross-validated ##################
			#### modified to obtain regular boosting, does not use CV!

			prediction.cv=my.predict.cart(results.cart, my.data.train)




 #maxerror <- min(1 - max(summary(y))/sum(summary(y)), 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance! 
 #maxerror <- min(1 - max(table(y))/num.samples, 0.5)  ##tu poveš kakšna je napaka naivnega klasifikatroja - 0.5 èe balance ali celo manj èe ni balance!   

 
 
#number of samples in test set
 num.samples.test=nrow(my.data.test)
 
 #matrix storing the predictions for new data, the m-th columns of the matrices store the predictions obtained at the m-th boosting iterations
 #prediction.new: prediction obtained from the CART fitted at the m-th iteration
 #prediction.boosting: boosting class prediction at the m-th iteration
 prediction.new=prediction.boosting=matrix(NA, ncol=num.boost, nrow=num.samples.test)
 
 #save the name of the environment, will be used later to assign the predictions to the prediction.boosting matrix, within another function (lapply)
 my.env=environment()
 
#alpha from boosting
	alpha=numeric(num.boost)

#saving the cart fit from the original data
	results.cart.original=results.cart

	#results.cart.boosting=results.cart
 
 #initializing the indexes for the boosting

#cat("Boosting iteration #: ", 1, "\n") 



################################################## step 3: boosting ###################################


if(verbose==TRUE) cat("Fitting CART using boosting\n")

 m=1
 my.stop.boosting=FALSE
 prediction.cv.boosting=prediction.cv
 
 #saving the original weights
 my.weights.original=my.weights
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 #saving the prediction obtained at the 1st boosting step
#prediction.boosting[,1]=my.predict.cart(results.cart.original, my.data.test)
																
 
############################# begin boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

	

	if(verbose==TRUE) cat("Boosting iteration #: ", m, "\n")

#cat("Boosting iteration #: ", m, "\n") 



	############ naive error of the classifier
	#indicator variable, indicating if the CV-prediction was correct or not, logical indicator, 1 if wrong prediction, 0 otherwise
	incorrect.class.cv=prediction.cv.boosting!=y

	#calculates weighted error from the CV-estimate
	#weighted.error=sum(my.weights.boosting[incorrect.class.cv])
	weighted.error=sum(my.weights[incorrect.class.cv])
	
	#stop also error=0, setting a large alpha=10
	if(weighted.error==0 ) {alpha[m]=10 
							#stop the boosting if it achieves 0 error
							my.stop.boosting=TRUE} else alpha[m]=ifelse(weighted.error>=0.5, 0.001, log((1-weighted.error)/weighted.error))
							########## calculates alpha in the usual way if the error is not 0 and is not larger than 0.05
								#calculate alpha, 										
								#alpha, or 2*beta from Izenman, Adaboost algorith from page 513, 	
								#alpha is positive if the error is small (<0.5), negative if it is large (>0.5)
								#alpha[m]=log((1-weighted.error)/weighted.error), set to a small constant if the error is large
							
	#stop also if error=1, setting a small alpha=0.001
	if(weighted.error==1 ) {#alpha[m]=0.001 - already included in the previous command
							#stop the boosting if it achieves  error
							my.stop.boosting=TRUE} 
	

		

		
		
		
	#saving the prediction for new samples for the m-th step	
	prediction.new[,m]=my.predict.cart(results.cart, my.data.test)

		
	#go on to the next step of boosting
	if(!my.stop.boosting){
				
						#cat("Estimating\n")										
						#calculate alpha, 										
						#alpha, or 2*beta from Izenman, Adaboost algorith from page 513, 
						#alpha is positive if the error is small (<0.5), negative if it is large (>0.5)
						#alpha[m]=log((1-weighted.error)/weighted.error), set to a small constant if the error is large
						
												
						#update the weights
						#my.weights.boosting=ifelse(incorrect.class.cv, my.weights.boosting*(1-weighted.error)/weighted.error, my.weights.boosting)
						my.weights=ifelse(incorrect.class.cv, my.weights*(1-weighted.error)/weighted.error, my.weights)
						#renormalize the weights
						#my.weights.boosting=my.weights.boosting/sum(my.weights.boosting)
						my.weights=my.weights/sum(my.weights)
																	
											
					
					
					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
							
							#new weighted relative frequencies, for class 1 and class 2
							n.l<-apply(my.w.s,2,cumsum)
							n.r<-1-n.l

							##### utility matrices: calculated once and reused later
							my.w.times.y.0=my.w.s*my.y.0
							my.w.times.y.s=my.w.s*my.y.s
							
							#update not needed, my.ranks=t(colRanks(my.data.train)), #my.y.0=my.y.s==0,#my.y.s<-matrix(y[my.ind.s], ncol=num.var)
							
							#defining the prior if it was not specified by the user, using the empirical prior; 
							#the empirical prior is equal to the WEIGHTED proporion of cases in each class #order for my.prior: class 1, class 0
							
							#if(is.null(my.prior)) my.prior=c(sum(my.weights.boosting[y==1]), sum(my.weights.boosting[y==0]))
							if(is.null(my.prior)) my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

							############## to include the prior in the calculation of the gini gain
							#N.l.w=sum(my.weights.boosting[y==1])
							N.l.w=sum(my.weights[y==1])
							N.r.w=1-N.l.w

							############### end update of the matrices used in cart#####################											
							
						############### fit the cart model and obtain CV class prediction #####################											
							
							#using name results.cart because it is the default argument of cart.for.cv.internal
							#results.cart=cart.internal()
							results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)

	
							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)

						############### end fit the cart model and obtain CV class prediction #####################											
																	
			}#end if(!my.stop.boosting), updated the weights and alpha
										
			
											
			#increase the counter before exiting
			m=m+1
									
	
		#saving the prediction for new samples for the next step	
		#prediction.new[,m]=my.predict.cart(results.cart.original, my.data.test)

	
	
	
}#end boosting, while, stops if it achieves a sufficient number of iterations or if a classifier does not make any erros

################# decreasing the counter of the number of performed iterations
m=m-1

#reduce the boosting objects if less iterations were made
#if(m<num.boosting){
#	prediction.new=prediction.new[,1:(m-1)]
#	alpha=alpha[1:(m-1)]
#	}
	
	
	########################## prediction of class membership for new samples using boosting ########################
	
	#calculates the total sum of alphas, times 2
	sum.alpha.half=sum(alpha, na.rm=T)/2
	#calculates a matrix, values are 0 if the prediction was 0 and alpha_m for the predicted in class 1 in boosting iteration m
	mat.prediction.new.y1.times.alpha=prediction.new*rep(alpha, each=num.samples.test)
	#derive the class membership of new samples based on boosting, if the sum of the alphas for the predicted in class 1 is greater than the sum of the alphas for the predicted in class 0 assign to class 1
	#the difference between the sum of alpha[y=1] and sum of alpha[y=0] can be calculated as 2*sum(alpha[y=1])-sum(alpha)
	#diff.alpha=rowSums(mat.prediction.new.y1.times.alpha)-sum.alpha.half
	#problema missing
	diff.alpha=rowSums(mat.prediction.new.y1.times.alpha, na.rm=T)-sum.alpha.half
	prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)
	
	#are there any ties?  if there are ties assign the class at random
	num.ties=sum(diff.alpha==0, na.rm=T)
	if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	
	
#calculates the sum of alpha for those predicted 1, molti
#	prediction.new.data.boosting=unlist(lapply(c(1:num.samples.test), function(i) sum(prediction.new[i,]*alpha, na.rm=T)))
#prediction.new.data.boosting=ifelse(prediction.new.data.boosting>all.alpha-prediction.new.data.boosting, 1, 0)





########### if the class membership of the test set is provided, calculate the error rate for each number of boosting iterations

		if(!is.null(y.test)) {

				 if(m==1) {error.boosting=mean(y.test!=prediction.new[,1]) 
							accuracy.test.all=1-error.boosting		
							accuracy.test.class1=mean(prediction.new[y.test==1,1]==y.test[y.test==1])
							accuracy.test.class0=mean(prediction.new[y.test==0,1]==y.test[y.test==0])
							}	else{

									cumsum.alpha.half=cumsum(alpha)/2
					
									#######boosting error on the test set, one value for each possible value of boosting iterations, for example: the second value would be the average error obtained using 2 boosting iterations
									error.boosting=c(#first step 
											mean(y.test!=prediction.new[,1]),
											#next steps
											#unlist(lapply(2:num.boost, function(m) {
											#calculute the prediction only for the effective boosting iterations
											unlist(lapply(2:m, function(mm) {

																diff.alpha=rowSums(mat.prediction.new.y1.times.alpha[,1:mm])-cumsum.alpha.half[mm]
																prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)

																#are there any ties?  if there are ties assign the class at random
																num.ties=sum(diff.alpha==0, na.rm=T)
																if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
																
																#saving the prediction obtained at the mm-th boosting step
																#does not work within a function: prediction.boosting[,mm]=prediction.new.data.boosting  and #assign("prediction.boosting[,mm]", prediction.new.data.boosting, env=parent.frame() )														
																#uses the environment of the external function
																my.env$prediction.boosting[,mm]=prediction.new.data.boosting
																
																#evaluates the error rate at the mm-th step
																mean(prediction.new.data.boosting!=y.test)

																															
															}#end function(mm) 
												)#end lapply
											)#end unlist
										)#end error.boosting
								}#end else m==1

								#evaluate the final overall and class specific predictive accuracies
								accuracy.test.all=1-error.boosting[m]
								accuracy.test.class1=mean(prediction.new.data.boosting[y.test==1]==y.test[y.test==1])
								accuracy.test.class0=mean(prediction.new.data.boosting[y.test==0]==y.test[y.test==0])	

			} else {error.boosting<-NA
					accuracy.test.all=NA
					accuracy.test.class1=NA
					accuracy.test.class0=NA	}


	

		########## transform the CART fit in a CART object
		results.cart=vector("list",2)
		results.cart$tree=results.cart.original
		results.cart$y=y
		class(results.cart)="CART"


		#saving the result boosting prediction at step 1, same as regular prediction
		prediction.boosting[,1]=prediction.new[,1]		

		
		
		#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
		return(list(#results.cart=results.cart, #prediction.cv=prediction.cv, 
		prediction.test.boosting=prediction.new.data.boosting, class.test=y.test, 
		accuracy.test.all=accuracy.test.all, accuracy.test.class1=accuracy.test.class1, accuracy.test.class0=accuracy.test.class0,
		beta=alpha[1:m]/2, 
		#prediction.boosting=prediction.new, 
		#m: number of performed boosting iterations
		weights.after.boosting=my.weights, error.boosting=error.boosting,  
		#accuracy.test.all=,  accuracy.test.class1=,  accuracy.test.class0=1-error.boosting, 
		score.test.class1=diff.alpha,   prediction.test.boosting.by.iteration=prediction.boosting[,1:m], iterations.boosting=m))




}
