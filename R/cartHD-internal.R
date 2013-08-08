.balanced.folds <-
function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
         # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
         # nice we to get the ids in a list, split by class
###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
#cat(i)
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
### Now do a clever sort to mix up the NAs
   smallmat <- .permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }
 
 
.cart.for.cv.boosting.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ind.s, num.boost, verbose=FALSE, my.env=NULL){
#cart.for.cv.internal=function(my.fold){

	#################################################################################
	#################### step 1: reduce the data set ################################
	#################################################################################
	
	
	####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		
		#print the CV progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Fold: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Fold: ", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}
		
		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	my.data.test=my.data.train[-my.index.keep,]
	num.samples.test=nrow(my.data.test)
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			#results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		


	
	################### step 1: compute cart and non-cross validated prediction on the reduced training set
	
	#saving the cart fit from the original data
	results.cart.original=results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
	
	#obtaining non cross-validatated predictions 
	prediction.cv=my.predict.cart(results.cart, my.data.train)




	
	
	####################### NEXT: obtain prediction for the left out samples using boosting ##########################


#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
#return(list(prediction.cv, results.cart.cv))






	#################################################################################
	#################### step 2: Perform boosting on the reduced data set (training) and evaluate it on the left out sample (test) ################################
	#################################################################################
	
	
	
	
#number of samples in test set
 num.samples.test=nrow(my.data.test)
 
 #matrix storing the predictions for new data
 prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
#alpha from boosting
	alpha=numeric(num.boost)


	#results.cart.boosting=results.cart
 
 #initializing the indexes for the boosting

#cat("Boosting interation #: ", 1, "\n") 



################################################## step 3: boosting ###################################

	#matrix storing the predictions for new data
	prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
	#alpha from boosting
	alpha=numeric(num.boost)
	#results.cart.boosting=results.cart
 
	#initializing the indexes for the boosting

	#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.cv.boosting=prediction.cv
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

	


#cat("Boosting interation #: ", m, "\n") 



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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
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
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
						############### fit the cart model and obtain the non-CV class prediction on the training set#####################											
							
							#using name results.cart because it is the default argument of cart.for.cv.internal, 
							#also here we need to restrict the attention to the subset of data
							results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
								

							prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)

								

							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
						

							#prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)

							
							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
						


						################### removed, does not perform internal CV ################### 
						#	if(nfolds<num.samples) my.folds.boosting=.balanced.folds(y, nfolds) else my.folds.boosting=as.list(1:num.samples)
						#	my.folds.boosting=lapply(my.folds.boosting, sort)	
						#			prediction.cv.boosting=unlist(lapply(my.folds.boosting, .cart.for.cv.internal.with.arguments, N=N, my.ranks=my.ranks, num.var=num.var, my.tr.s=my.tr.s, my.y.s=my.y.s, my.w.s=my.w.s, my.y.0=my.y.0,my.data.train=my.data.train,  my.weights=my.weights, y=y, min.samples=min.samples, max.depth=max.depth, min.gain))[order(unlist(my.folds.boosting))] 
						################### ################### ################### ################### ################### 


						
							######## convert the selected indexes in the original indexes, problem: .balanced.folds does produce indexes that are in the range: 1 to num.samples (reduced by ee), must be re-transfromed in the original indexes to be used together with the old indeces.
	#						my.folds.internal.converted=lapply(my.folds.boosting, function(x) my.index.keep[x])
							
							#### add to the samples to remove for CV, those that really need to be cross-validated are reported first 
	#						my.folds.internal.all=lapply(my.folds.internal.converted, function(fold) c( fold, my.fold )) 
	
							
							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							####prediction.cv.boosting=unlist(lapply(my.folds.boosting, cart.for.cv.internal))[order(unlist(my.folds.boosting))]   
							
							###### updated 4/3/2013: the function cart.for.cv.internal was upadted to allow to arguments: the samples to leave out (all those that were not used in the easy-ensemble fold (my.fold[[i]]) and those that are left out because of the internal CV (my.folds.interal[[i]]), the evaluation is then made only for the samples left out for CV
	#						prediction.cv.boosting=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]])[1:length(my.folds.boosting[[i]])]))[order(unlist(my.folds.boosting))]


				
	

							
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
	

#return(prediction.new.data.boosting)
#changed 7/5/2013: returns cross-validated classification results and CV-classification scores to the main function
return(list(pred.cv=prediction.new.data.boosting, score.cv=diff.alpha))

}# end of cart.for.boosting.cv.internal.with.arguments






.cart.for.cv.boosting.Icv.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ind.s, num.boost, 
nfolds.Icv, #additional argument for Icv function
verbose=FALSE, my.env=NULL){
#cart.for.cv.internal=function(my.fold){

################ internal function that is used by cart.boosting.Icv.cv - boosting with internal CV and evaluation of the error using CV evaluated on training data

	#################################################################################
	#################### step 1: reduce the data set ################################
	#################################################################################
	
	
	####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		
		#print the CV progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Fold: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Fold: ", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}
		
		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	my.data.test=my.data.train[-my.index.keep,]
	num.samples.test=nrow(my.data.test)
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			#results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		


	
	################### step 1: compute cart and cross validated (!) prediction on the reduced training set
	
	#saving the cart fit from the original data
	results.cart.original=results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
	
	#previous version, used to obtain non-cross-validatated (!) predictions using test data
	#prediction.cv=my.predict.cart(results.cart, my.data.train) - used in the non-Icv version of the algorithm

	
	#folds for internal cross-validation
	#10/6/2013: for debugging
	#set.seed(1)
	if(nfolds.Icv<num.samples) {my.folds.Icv=.balanced.folds(y, nfolds.Icv) #
						my.folds.Icv=lapply(my.folds.Icv, sort)} else my.folds.Icv=as.list(1:N)

	
	#prediction.cv=unlist(lapply(my.folds.Icv, .cart.for.cv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
	prediction.cv=unlist(lapply(my.folds.Icv, .cart.for.cv.prediction.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]

	
	
	
	####################### NEXT: obtain prediction for the left out samples using boosting ##########################


		#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
		#return(list(prediction.cv, results.cart.cv))






	#################################################################################
	#################### step 2: Perform boosting on the reduced data set (training) and evaluate it on the left out sample (test) ################################
	#################################################################################
	
	
	
	
#number of samples in test set
 num.samples.test=nrow(my.data.test)
 
 #matrix storing the predictions for new data
 prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
#alpha from boosting
	alpha=numeric(num.boost)


	#results.cart.boosting=results.cart
 
 #initializing the indexes for the boosting

#cat("Boosting interation #: ", 1, "\n") 



################################################## step 3: boosting ###################################

	#matrix storing the predictions for new data
	prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
	#alpha from boosting
	alpha=numeric(num.boost)
	#results.cart.boosting=results.cart
 
	#initializing the indexes for the boosting

	#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.cv.boosting=prediction.cv
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

	


#cat("Boosting interation #: ", m, "\n") 



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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
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
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
						############### fit the cart model and obtain the CV class prediction on the training set#####################											
							
							#using name results.cart because it is the default argument of cart.for.cv.internal, 
							#also here we need to restrict the attention to the subset of data
							results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
								

							#version used in the non-Icv version of boosting
							#prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)
		
							#10/6/2013: added for debugging
							#important: to select different samples for CV at each boosting step
							# set.seed(1+m)
							#set.seed(1)
							if(nfolds.Icv<num.samples) {my.folds.Icv=.balanced.folds(y, nfolds.Icv) #
														my.folds.Icv=lapply(my.folds.Icv, sort)} else my.folds.Icv=as.list(1:N)
	
					
							#CV-class prediction
							#prediction.cv.boosting=unlist(lapply(my.folds.Icv, .cart.for.cv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
							prediction.cv.boosting=unlist(lapply(my.folds.Icv, .cart.for.cv.prediction.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
			
							
							
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
	#10/6/2013: for debugging
	#set.seed(1)
	num.ties=sum(diff.alpha==0, na.rm=T)
	if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	

#return(prediction.new.data.boosting)
#changed 7/5/2013: returns cross-validated classification results and CV-classification scores to the main function
return(list(pred.cv=prediction.new.data.boosting, score.cv=diff.alpha))

} #end .cart.for.boosting.cv.Icv.internal.with.arguments()



.cart.for.cv.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, verbose=FALSE, my.env=NULL){
	
#cart.for.cv.internal=function(my.fold){

		#print the CV progress on screen and increment the internal index
		if(verbose==TRUE) {	cat("Fold: ", my.env$.internal.index, "\n")
							#.internal.index<<-.internal.index+1
							assign(".internal.index", my.env$.internal.index+1, envir=my.env)
							}

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		





####################### prediction for left out samples ##########################


		#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
		#26/4/2013: added the output of the prediction scores, not only the class membership, the output is a matrix containing class prediction in the first row and the class 1 scores in the second row
		#warning: should be a vector in case of LOO-CV... deal with it correctly....
		prediction.cv=my.predict.cart.with.score(results.cart.cv, my.data.test)
		#return(list(prediction.cv, results.cart.cv))
		return(prediction.cv)
}#.end cart.for.cv.internal.with.arguments 









.cart.for.cv.prediction.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, verbose=FALSE, my.env=NULL){
#outputs only the prediction, not the scores

	
#cart.for.cv.internal=function(my.fold){

		#print the CV progress on screen and increment the internal index
		if(verbose==TRUE) {	cat("Fold: ", my.env$.internal.index, "\n")
							#.internal.index<<-.internal.index+1
							assign(".internal.index", my.env$.internal.index+1, envir=my.env)
							}

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		





####################### prediction for left out samples ##########################


		#prediction only, no scores
		prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
		#26/4/2013: added the output of the prediction scores, not only the class membership, the output is a matrix containing class prediction in the first row and the class 1 scores in the second row
		#warning: should be a vector in case of LOO-CV... deal with it correctly....
		#prediction and scores
		#prediction.cv=my.predict.cart.with.score(results.cart.cv, my.data.test)
		#return(list(prediction.cv, results.cart.cv))
		return(prediction.cv)
}#.end cart.for.cv.internal.with.arguments 



















.cart.for.cv.mds.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.mds, verbose=FALSE, my.env=NULL){
#cart.for.cv.internal=function(my.fold){

	#################################################################################
	#################### step 1: reduce the data set ################################
	#################################################################################
	
	#print the CV progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Fold: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Fold: ", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}

	
	####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#my.data.test=my.data.train[-my.index.keep,]
	#23/4/2013: fixed bug: keep a matrix format even if only one sample is included in the test set - as happens in LOO-CV
	my.data.test=matrix(my.data.train[-my.index.keep,], ncol=num.var, byrow=F)
	
	num.samples.test=nrow(my.data.test)
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			#results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		





####################### prediction for left out samples ##########################


#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
#return(list(prediction.cv, results.cart.cv))






	#################################################################################
	#################### step 2: Perform MDS ################################
	#################################################################################
	


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
		my.folds.mds=lapply(1:num.mds, function(i) sample(which.maj.class, abs(N1-N0)))
		#### 18/3/2013: added the sorting of the folds, the indexes within each fold are sorted
		my.folds.mds=lapply(my.folds.mds, sort)

		##### 18/3/2013: here it is not necessary to re-order the output because the test set is separate
		my.predictions.new=matrix(unlist(lapply(my.folds.mds, .cart.for.mds.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, my.data.test)), ncol=num.samples.test, byrow=T)
		#my.predictions.new=matrix(unlist(lapply(my.folds.mds, .cart.for.mds.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, my.data.test)), ncol=num.samples.test, byrow=T)[order(unlist(my.folds.mds),]

		#derive the global score for the test samples, combining the 
		my.scores.final=apply(my.predictions.new, 2, sum)

		
		## predicted class membership for the samples that were left out from cross-validation
		y.pred=ifelse(my.scores.final>num.mds/2, 1, 0)
		
		#are there any ties?  if there are ties assign the class at random
		num.ties=sum(my.scores.final==num.mds/2, na.rm=T)
		if(num.ties>0)  y.pred[which(my.scores.final==num.mds/2)]=sample(c(0,1), num.ties, replace=T)
	
		
		
		############## end added for multiple downsizing ####################



#return(y.pred)
#updated: 6/5/2013: the output contains also the score used for the classification
return(c(y.pred, my.scores.final/num.mds))
}# end of .cart.for.cv.mds.internal.with.arguments





.cart.for.cv.ee.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.ee, num.boost, verbose=FALSE, my.env=NULL){

#internal function that performs cross-validated class prediction using easy ensamble - without internal cross-validation for the evaluation of the prediction error within the boosting step



	#################################################################################
	#################### step 1: reduce the data set ################################
	#################################################################################
	
	if (verbose==TRUE) {	#cat("Fold: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Fold: ", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}
							
	####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#my.data.test=my.data.train[-my.index.keep,]
	#23/4/2013: fixed bug: keep a matrix format even if only one sample is included in the test set - as happens in LOO-CV
	my.data.test=matrix(my.data.train[-my.index.keep,], ncol=num.var, byrow=F)
	#class membership of the left out samples
	y.test=y[-my.index.keep]	
	
	num.samples.test=nrow(my.data.test)
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed for EE, re-evaluation of the order of the variables ###################################
	
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again
	#24/4/2013: check if it is needed or the original indexes are needed!!!!


	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			#results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		





####################### prediction for left out samples ##########################


#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
#return(list(prediction.cv, results.cart.cv))






	#################################################################################
	#################### step 2: Perform EE ################################
	#################################################################################
	


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
		my.folds.ee=lapply(1:num.ee, function(i) sample(which.maj.class, abs(N1-N0)))
		#### 18/3/2013: added the sorting of the folds, the indexes within each fold are sorted
		my.folds.ee=lapply(my.folds.ee, sort)

		##### 18/3/2013: here it is not necessary to re-order the output because the test set is separate
		#the output score is diff.alpha
		my.scores=matrix(unlist(lapply(my.folds.ee, .cart.for.ee.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s)), ncol=num.samples.test, byrow=T)
		#.cart.for.ee.interal.with.arguments=function(my.fold, 								    	N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s){

		
		#derive the global score for the test - left out - samples, combining the 
		my.scores.final=apply(my.scores, 2, sum)

		## predicted class membership for the samples that were left out from cross-validation
		y.pred=ifelse(my.scores.final>0, 1, 0)

		num.ties=sum(my.scores.final==0, na.rm=T)
		if(num.ties>0)  y.pred[which(my.scores.final==0)]=sample(c(0,1), num.ties, replace=T)
	
		
		
		## predicted class membership for the samples that were left out from cross-validation
		#### removed 7/5/2013: y.pred=ifelse(my.scores.final>num.ee/2, 1, 0) the classification rule is based on the sum of the diff.alpha only
		
		#are there any ties?  if there are ties assign the class at random
		#num.ties=sum(my.scores.final==num.ee/2, na.rm=T)
		#if(num.ties>0)  y.pred[which(my.scores.final==num.ee/2)]=sample(c(0,1), num.ties, replace=T)
	
		
		
		############## end added for multiple downsizing ####################



#return(y.pred)
#updated 7/5/2013: returns also the scores to the main function, not only the class prediction
return(list(y.pred=y.pred, scores.final=my.scores.final))
}# end of .cart.for.cv.ee.internal.with.arguments








.cart.for.ee.cv.internal.with.arguments <-
function(my.fold, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, nfolds){

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out
		
		#### the used data are balanced undersampled training sets, since we are using easy ensemble

		#########which samples to keep in the EE 	step
		my.index.keep=c(1:N)[-my.fold]

	
	
	################### step 0: resize the data#############################################
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)


	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#here the test 
	#my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################


		
		#if(length(table(y))==1) { #unnecessary check for EE! data are class-balanced
		#return(list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
		#classNode.l=NULL, classNode.r=NULL, 
		#		#QT.max=NULL, 
		#		which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=depth, parent=parent, parent.direction=parent.direction, 
		#		parent.history=parent.history, giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL,
		#		gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
		#		QT.rpart.max=NULL
		#		#, which.variables.used=NULL
		#		)#,change=change)
		#)#}



#step 1

#my.stump.mat1.Lara.duplicates.faster.internal<-function(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var,  
#		N, n.l, n.r, my.weights, my.index.keep=NULL, 
#		QT.max.old=NULL,  min.samples=2, min.gain=0.01, which.variables.used=NULL, dep
		
		
		
		
		
		########## perform cross-validated CART con the subset of the data, in this case on a balanced (easy ensemble) subset
		
		

		
		############################################################################
		############ step 1: estimate CART on the complete easy ensemble data
		############################################################################


		results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
		
												  
		


#saving the cart fit from the original data
	results.cart.original=results.cart


################################### step 2: obtain cross-validated class predictions #########################################
	
			#my.folds=.balanced.folds(y, nfolds)
			########## add an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list
			##### my.folds.internal: a list containing the IDs of the samples to be used in the 
			
			#10/6/2013: added for debugging
			#set.seed(1)
			
			if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
			######## convert the selected indexes in the original indexes, problem: .balanced.folds does produce indexes that are in the range: 1 to num.samples (reduced by ee), must be re-transfromed in the original indexes to be used together with the old indeces.
#			my.folds.internal.converted=lapply(my.folds.internal, function(x) my.index.keep[x])

			
			#### add to the samples to remove for CV, those that really need to be cross-validated are reported first 
#			my.folds.internal.all=lapply(my.folds.internal.converted, function(fold) c( fold, my.fold )) 
			

			################ obtain cross-validated prediction using nfolds-CV
			##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
			#prediction.cv=unlist(lapply(my.folds.internal, cart.for.cv.internal))[order(unlist(my.folds.internal))]
			###### updated 4/3/2013: the function cart.for.cv.internal was upadted to allow to arguments: the samples to leave out (all those that were not used in the easy-ensemble fold (my.fold[[i]]) and those that are left out because of the internal CV (my.folds.interal[[i]]), the evaluation is then made only for the samples left out for CV
			
#			prediction.cv=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]])[1:length(my.folds.internal[[i]])]))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(my.folds.internal, cart.for.cv.internal))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]], my.folds.internal[[i]] )))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))

			#prediction.cv=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
			prediction.cv=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]


#################################### step 3: boosting part ######################################################################



 
	#matrix storing the predictions for new data
	prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
	#alpha from boosting
	alpha=numeric(num.boost)
	#results.cart.boosting=results.cart
 
	#initializing the indexes for the boosting

	#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.cv.boosting=prediction.cv
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

	


#cat("Boosting interation #: ", m, "\n") 



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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
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
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
							#using name results.cart because it is the default argument of cart.for.cv.internal, 
							#also here we need to restrict the attention to the subset of data
							results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
																	  

							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
						

							#prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)

							
							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							#10/6/2013: added for debugging
							#set.seed(1+m)
							#set.seed(1)
							
							if(nfolds<num.samples) my.folds.boosting=.balanced.folds(y, nfolds) else my.folds.boosting=as.list(1:num.samples)
							my.folds.boosting=lapply(my.folds.boosting, sort)	
								
							######## convert the selected indexes in the original indexes, problem: .balanced.folds does produce indexes that are in the range: 1 to num.samples (reduced by ee), must be re-transfromed in the original indexes to be used together with the old indeces.
	#						my.folds.internal.converted=lapply(my.folds.boosting, function(x) my.index.keep[x])
							
							#### add to the samples to remove for CV, those that really need to be cross-validated are reported first 
	#						my.folds.internal.all=lapply(my.folds.internal.converted, function(fold) c( fold, my.fold )) 
	
							
							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							####prediction.cv.boosting=unlist(lapply(my.folds.boosting, cart.for.cv.internal))[order(unlist(my.folds.boosting))]   
							
							###### updated 4/3/2013: the function cart.for.cv.internal was upadted to allow to arguments: the samples to leave out (all those that were not used in the easy-ensemble fold (my.fold[[i]]) and those that are left out because of the internal CV (my.folds.interal[[i]]), the evaluation is then made only for the samples left out for CV
	#						prediction.cv.boosting=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]])[1:length(my.folds.boosting[[i]])]))[order(unlist(my.folds.boosting))]


					#prediction.cv.boosting=unlist(lapply(my.folds.boosting, .cart.for.cv.internal.with.arguments, N=N, my.ranks=my.ranks, num.var=num.var, my.tr.s=my.tr.s, my.y.s=my.y.s, my.w.s=my.w.s, my.y.0=my.y.0,my.data.train=my.data.train,  my.weights=my.weights, y=y, min.samples=min.samples, max.depth=max.depth, min.gain))[order(unlist(my.folds.boosting))] 
					prediction.cv.boosting=unlist(lapply(my.folds.boosting, .cart.for.cv.prediction.internal.with.arguments, N=N, my.ranks=my.ranks, num.var=num.var, my.tr.s=my.tr.s, my.y.s=my.y.s, my.w.s=my.w.s, my.y.0=my.y.0,my.data.train=my.data.train,  my.weights=my.weights, y=y, min.samples=min.samples, max.depth=max.depth, min.gain))[order(unlist(my.folds.boosting))] 

	

							
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
	
	#10/6/2013: added for debugging
	#set.seed(1)
	
	num.ties=sum(diff.alpha==0, na.rm=T)
	if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	
	
#calculates the sum of alpha for those predicted 1, molti
#	prediction.new.data.boosting=unlist(lapply(c(1:num.samples.test), function(i) sum(prediction.new[i,]*alpha, na.rm=T)))
#prediction.new.data.boosting=ifelse(prediction.new.data.boosting>all.alpha-prediction.new.data.boosting, 1, 0)





########### if the class membership of the test set is provided, calculate the error rate for each number of boosting iterations
############ removed, it is not used in the output
#if(!is.null(y.test)) {

# 	cumsum.alpha.half=cumsum(alpha)/2
	
#	if(m==1) error.boosting=mean(y.test!=prediction.new[,1]) else{
	
	#######boosting error on the test set, one value for each possible value of boosting iterations, for example: the second value would be the average error obtained using 2 boosting iterations
#	error.boosting=c(#first step 
#			mean(y.test!=prediction.new[,1]),
			#next steps
#			unlist(lapply(2:num.boost, function(mm) {
#
#								diff.alpha=rowSums(mat.prediction.new.y1.times.alpha[,1:mm])-cumsum.alpha.half[mm]
#								prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)
#
#								#are there any ties?  if there are ties assign the class at random
#								num.ties=sum(diff.alpha==0, na.rm=T)
#								if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
#								
#								mean(prediction.new.data.boosting!=y.test)
#
#
#							}#end function(m) 
#				)#end lapply
#			)#end unlist
#		)#end error.boosting
#}#end else m==1
#} else {error.boosting<-NA}


########### added for EE, calculates the score, defined as diff.alpha, a score is given for each sample 
return(score=diff.alpha)


#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
##return(list(results.cart=results.cart.original, prediction=prediction.cv, ##prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#prediction.boosting=prediction.new, 
#m: number of performed boosting iterations
##weights.after.boosting=my.weights, ##error.boosting=error.boosting,prob=diff.alpha,pred.new.i=prediction.new, m=m))







}


.cart.for.ee.internal.with.arguments <-
#browser()
function(my.fold, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, verbose=FALSE, my.env=NULL){

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#print the EE progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Easy ensemble step: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Easy ensemble step:", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}

		
		#########which samples to keep in the cross-validated step
		
		my.index.keep=c(1:N)[-my.fold]
	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)


	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#here the test 
	#my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			


		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################


	
	
			###################### step 1: evaluate cart############################### 
		
			########################## fit cart on original data ##########################		
			results.cart.ee=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		

	

			#prediction on the same data used to develop the predictor, the balanced training set, it is NOT cross-validated
			prediction.ee=my.predict.cart(results.cart.ee, my.data.train)

#return(list(prediction.cv, results.cart.cv))
#return(prediction.cv)





#################################### boosting part ###########################



 #matrix storing the predictions for new data
 prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
#alpha from boosting
	alpha=numeric(num.boost)

#saving the cart fit from the original data
	results.cart.original=results.cart.ee

	#results.cart.boosting=results.cart
 
 #initializing the indexes for the boosting

#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.ee.boosting=prediction.ee
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

		
#cat(m, " boosting \n")	


#cat("Boosting iteration #: ", m, "\n") 



	############ naive error of the classifier
	#indicator variable, indicating if the EE-prediction was correct or not, logical indicator, 1 if wrong prediction, 0 otherwise
	incorrect.class.ee=prediction.ee.boosting!=y

	#calculates weighted error from the EE-estimate
	#weighted.error=sum(my.weights.boosting[incorrect.class.ee])
	weighted.error=sum(my.weights[incorrect.class.ee])
	
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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
	prediction.new[,m]=my.predict.cart(results.cart.ee, my.data.test)

		
	#go on to the next step of boosting
	if(!my.stop.boosting){
				
						#cat("Estimating\n")										
						#calculate alpha, 										
						#alpha, or 2*beta from Izenman, Adaboost algorith from page 513, 
						#alpha is positive if the error is small (<0.5), negative if it is large (>0.5)
						#alpha[m]=log((1-weighted.error)/weighted.error), set to a small constant if the error is large
						
												
						#update the weights
						#my.weights.boosting=ifelse(incorrect.class.cv, my.weights.boosting*(1-weighted.error)/weighted.error, my.weights.boosting)
						my.weights=ifelse(incorrect.class.ee, my.weights*(1-weighted.error)/weighted.error, my.weights)
						#renormalize the weights
						#my.weights.boosting=my.weights.boosting/sum(my.weights.boosting)
						my.weights=my.weights/sum(my.weights)
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
						############### fit the cart model and obtain the (non cross validated) class prediction #####################											
							
							#using name results.cart because it is the default argument of cart.for.cv.internal
							#results.cart=cart.internal()
							results.cart.ee=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
															          

							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							prediction.ee.boosting=my.predict.cart(results.cart.ee, my.data.train)

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
	
	###prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)
	
	#are there any ties?  if there are ties assign the class at random
	###num.ties=sum(diff.alpha==0, na.rm=T)
	###if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	
	
#calculates the sum of alpha for those predicted 1, molti
#	prediction.new.data.boosting=unlist(lapply(c(1:num.samples.test), function(i) sum(prediction.new[i,]*alpha, na.rm=T)))
#prediction.new.data.boosting=ifelse(prediction.new.data.boosting>all.alpha-prediction.new.data.boosting, 1, 0)





########### if the class membership of the test set is provided, calculate the error rate for each number of boosting iterations
########## removed: not used in the output #####################
#if(!is.null(y.test)) {

 #	cumsum.alpha.half=cumsum(alpha)/2
	
	#######boosting error on the test set, one value for each possible value of boosting iterations, for example: the second value would be the average error obtained using 2 boosting iterations
# if(m==1) error.boosting=mean(y.test!=prediction.new[,1]) else{
#	error.boosting=c(#first step 
#			mean(y.test!=prediction.new[,1]),
#			#next steps
#			unlist(lapply(2:num.boost, function(m) {

#								diff.alpha=rowSums(mat.prediction.new.y1.times.alpha[,1:m])-cumsum.alpha.half[m]
#								prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)

								#are there any ties?  if there are ties assign the class at random
#								num.ties=sum(diff.alpha==0, na.rm=T)
#								if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
								
#								mean(prediction.new.data.boosting!=y.test)

 
#							}#end function(m) 
#				)#end lapply
#			)#end unlist
#		)#end error.boosting
#}# end else m==1
#} else {error.boosting<-NA}


########### added for EE, calculates the score, defined as diff.alpha, a score is given for each sample 
#return(list(score=diff.alpha,ff=mat.prediction.new.y1.times.alpha))
return(score=diff.alpha)


#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
##return(list(results.cart=results.cart.original, prediction=prediction.cv, ##prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#prediction.boosting=prediction.new, 
#m: number of performed boosting iterations
##weights.after.boosting=my.weights, ##error.boosting=error.boosting,prob=diff.alpha,pred.new.i=prediction.new, m=m))







}# end .cart.for.ee.internal.with.arguments



.cart.for.mds.internal.with.arguments <-
function(my.fold,  N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks,  my.data.test, my.env=NULL, verbose=FALSE){

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out, from 9/5/2013: the IDs are sorted within the folds!
		
		#### the used data are balanced undersampled training sets, since we are using easy ensemble

		#print the EE progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Easy ensemble step: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Multiple downsizing step:", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}

		
		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

		
		
		############## step 1: recalculate the quantities necessary for the estimation 
	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)


	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#here the test 
	#my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################


	
		
		############################################################################
		############ step 2: estimate CART on the complete downsized data
		############################################################################

		results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
												  




#saving the cart fit from the original data
	results.cart.original=results.cart


	
	############### step 2: obtain predictions for the new samples #####################
	
	prediction.new=my.predict.cart(results.cart.original, my.data.test)
	
	
	
	
########### returns the prediction for each sample of the test set 
return(prediction.new=prediction.new)


#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
##return(list(results.cart=results.cart.original, prediction=prediction.cv, ##prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#prediction.boosting=prediction.new, 
#m: number of performed boosting iterations
##weights.after.boosting=my.weights, ##error.boosting=error.boosting,prob=diff.alpha,pred.new.i=prediction.new, m=m))







}# end .cart.for.mds.internal.with.arguments


.cart.internal.with.arguments <-
function(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks){
#cart.internal.with.arguments=function(){

results.cart=vector("list")		

############# parent node, regular stump ######################

k=1

results.cart[[k]]=.my.stump.forCART.oneStep.stepOne(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, 
		N, n.l, n.r, my.weights,  #my.index.keep=1:N, 
		#QT.max.old=NULL, 
		min.samples=min.samples, min.gain=min.gain, #which.variables.used=NULL, depth=1, parent= 0, parent.direction=0, parent.history=1, 
		y=y, my.data.train=my.data.train, 
		#my.weigths.prob, 
		my.prior=my.prior, N.l.w=N.l.w, N.r.w=N.r.w)
		
		##### add predicted class of the node, in case it is a final node , it is based on the weighted frequency of class 1 vs class 0 only
		results.cart[[k]]$classNode.0=ifelse(sum(my.weights[y==1])>0.5, 1, 0)
		if(sum(my.weights[y==1])==1) results.cart[[k]]$classNode.0=sample(c(0, 1), 1) 
			
############################ next steps
#if(max.depth==1 | is.null(results.cart[[k]]$var))  return(results.cart) #return the list with the tree results


############################ next steps #####################################
##################### do this part only of depth was larger than 1, or a variable was selected in the first step

if(max.depth>1 & !is.null(results.cart[[k]]$var)){  #return(results.cart) #return the list with the tree results

for(i.depth in 2:max.depth){
		my.final.stop=0
	#	cat( "\n Depth: ", i.depth)
	for(i.horizontal in 1: (2^(i.depth-1)))	{
		k=k+1
	#	cat("H: ", i.horizontal, ",  ")
		#counter that shows how many nodes were not evaluated because the parent node was a final node
		#which is the parent node (k index) is calculated from the indexes i.depth and i.horizontal
		which.k.parent=2^(i.depth-2)+floor((i.horizontal-1)/2)
		which.direction.parent= i.horizontal %% 2 #1 if it considers the left node of the parent, n.l; 0 if it considers the right node of the parent: n.r, it essentially checks if a number is odd or even 

		
		#add to parent node the k of the left and right son				
		if(which.direction.parent==1) results.cart[[which.k.parent]]$k.son.l=k else results.cart[[which.k.parent]]$k.son.r=k
		
		#run stumps on i.depth level, i.horizontal element 

	
		#my.stop of parent node: 
		#TRUE: it was a final node
		#FALSE: further splitting can be done
		#NA: parent node was not split, stopped earlier
		if(results.cart[[which.k.parent]]$my.stop==T | is.na(results.cart[[which.k.parent]]$my.stop)) {
																		results.cart[[k]]=list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
																				classNode.l=NULL, classNode.r=NULL, 
																				#QT.max=QT.max, 
																				which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=i.depth, parent=which.k.parent, parent.direction=which.direction.parent, 
																				parent.history=c(results.cart[[which.k.parent]]$parent.history, k), 
																				#giniP=giniP, 
																				giniP=NULL,
																				p1Node.l.new=NULL, p1Node.r.new=NULL,
																				gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
																				QT.rpart.max=NULL
																				#, which.variables.used=NULL
																				)
																			my.final.stop=my.final.stop+1

																				}#end if

			  else {#perform stump on this node
								#derive which samples are going to be considered in this step, which.n.l [12] if going left, which.n.r [13] if going right
								My.index.keep=unlist(ifelse(which.direction.parent==1, results.cart[[which.k.parent]][5], results.cart[[which.k.parent]][6]))
								#p(A) of the parent node - probability of being in the parent node for future observations
								pA=ifelse(which.direction.parent==1, results.cart[[which.k.parent]]$pA.l,  results.cart[[which.k.parent]]$pA.r)
								
								#gini of the parent node
								giniP=ifelse(which.direction.parent==1, results.cart[[which.k.parent]]$gini.l, results.cart[[which.k.parent]]$gini.r)
		
							
		
							#calculate the next stump	
							results.cart[[k]]=.my.stump.forCART.oneStep.stepNotOne(my.tr.s, my.y.s, my.w.s, my.y.0, 
							#removed my.w.times.y.0, my.w.times.y.s,  ########
							num.var,
											N, 
											######## removed n.l,  n.r, ##############
											my.weights, my.index.keep=My.index.keep, 
											giniP=giniP,
											min.samples=min.samples, min.gain=min.gain,
											#which.variables.used=results.cart[[which.k.parent]]$which.variables.used, 
											depth=i.depth, parent=which.k.parent, parent.direction=which.direction.parent, 
											parent.history=c(results.cart[[which.k.parent]]$parent.history, k), y, 
											my.data.train=my.data.train,  
											#my.weigths.prob, 
											my.prior=my.prior, N.l.w=N.l.w, N.r.w=N.r.w, pA=pA, k=k, my.ranks=my.ranks) 
							
							
							#predicted class membership in case it is a final node
							if(which.direction.parent==1) results.cart[[k]]$classNode.0=results.cart[[which.k.parent]]$classNode.l else results.cart[[k]]$classNode.0=results.cart[[which.k.parent]]$classNode.r
		
							#check if there was a stopping at both son nodes, and the stopping indicator my.stop of the parent node needs to be changed, makes this check after running stump on the right node
							if(which.direction.parent==0 & is.na(results.cart[[k]]$my.stop) & is.na(results.cart[[k-1]]$my.stop))  results.cart[[which.k.parent]]$my.stop=TRUE
			
							#check if there was a stopping at one of the son nodes, but not both. in this case the tree is stopped in one direction only and the stopping indicator my.stop of the parent node needs to be changed
							#if: stop for left node but not for right node; else: stop for right node but not for left
							if(which.direction.parent==0 & !is.na(results.cart[[k]]$my.stop) & is.na(results.cart[[k-1]]$my.stop))  {results.cart[[k-1]]$var=1  #setting an aribitrary value for variable and threshold
																																	 results.cart[[k-1]]$thr=99999	
																																	 results.cart[[k-1]]$classNode.l=results.cart[[k-1]]$classNode.r=results.cart[[which.k.parent]]$classNode.l	
																																	 #added 6/5/2013: add also the probability of begin in class 1
																																	 results.cart[[k-1]]$p1Node.l.new=results.cart[[k-1]]$p1Node.r.new=results.cart[[which.k.parent]]$p1Node.l.new	
																														 																															 
																																	 results.cart[[k-1]]$my.stop=TRUE
																																	#else: stop for right node (k) but not for left node (k-1)
																																	} else {if(which.direction.parent==0 & is.na(results.cart[[k]]$my.stop) & !is.na(results.cart[[k-1]]$my.stop))   {
																																		results.cart[[k]]$var=1  #setting an aribitrary value for variable and threshold
																																		results.cart[[k]]$thr=99999	
																																		results.cart[[k]]$classNode.l=results.cart[[k]]$classNode.r=results.cart[[which.k.parent]]$classNode.r	
																																		#added 6/5/2013: add also the probability of begin in class 1
																																	 results.cart[[k]]$p1Node.l.new=results.cart[[k]]$p1Node.r.new=results.cart[[which.k.parent]]$p1Node.r.new	
																																	 
																																	 results.cart[[k]]$my.stop=TRUE
																																	
																																}#end else
																																}#end else2
			
			
			
							#count as stop also this case, where there was no sufficient gain in the node, and my.stop was set to NA
							if(is.na(results.cart[[k]]$my.stop)) my.final.stop=my.final.stop+1
							
							#check if the final depth was reached, in this case set stopping rule to TRUE, useful for prediction purposes
							#if(i.depth==max.depth & results.cart[[k]]$my.stop==FALSE) results.cart[[k]]$my.stop=TRUE#
							if(i.depth==max.depth & !is.na(results.cart[[k]]$my.stop) & results.cart[[k]]$my.stop==FALSE) results.cart[[k]]$my.stop=TRUE
							
							}# end else, perform stump on this node
									

		
	
		
}#end i.horizontal

if(my.final.stop==(2^(i.depth-1))) break
}#end i.depth


#end if(max.depth>1 & !is.null(results.cart[[k]]$var)) , performed steps after the first one
} else{#stopped at step 1
		results.cart[[1]]$k.son.l=2 
		results.cart[[1]]$k.son.r=3
		results.cart[[1]]$my.stop=TRUE

		} 

	return(results.cart)	
	}############# end cart.internal.with.arguments()
	
	
	
	
	
.my.stump.forCART.oneStep.stepNotOne <-
function(my.tr.s, my.y.s, my.w.s, my.y.0, #my.w.times.y.0, my.w.times.y.s,  
		num.var,  
		N, 
		#n.l, n.r, 
		my.weights, my.index.keep=NULL, 
		giniP=NULL,  min.samples=3, min.gain=0.01, 
		#which.variables.used=NULL, 
		depth=1, parent=0, 
		parent.direction=0, parent.history=NULL, y, my.data.train, 
		#my.weigths.prob, 
		my.prior, N.l.w, N.r.w, pA, k=NULL,  my.ranks=my.ranks){
#my.tr.s: sorted training set
#my.y.s: sorted class membership, in matrix format
#my.w.s: sorted weights, in matrix format
#my.y.0=my.y.s==0
#my.w.times.y.0=my.w.s*my.y.0
#my.w.times.y.s=my.w.s*my.y.s
#num.var: number of variables
#N: number of samples
#my.index.keep: index of the observations to keep in this node
#QT.max.old: gini impurity for future observations for the parent node
#min.samples: minimum number of samples in each class in order to procede with splitting
#min.gain: minimum gain obtained in order to procede with splitting
#depth: current depth of tree, 1 if at first step, etc
#parent: id of the parent node

#my.prior: class prior probabilities, vector: prior for class 1, prior for class0, must sum up to 1
#pA: P(A) of parent node, derived at a previous step
#pA.l and pA.r: matrices - removed now as arguments
#n.l.all: weighted relative frequencies

################## remove the observations not needed, version using the index to LEAVE OUT

#if(length(my.index.keep)<N) {


############### for prior incorporation ##################
#saving the original weights and y to 
my.weights.original=my.weights
y.original=y

#	my.w.times.y.s=my.w.s*my.y.s

 

########## resizing the objects, eliminating the samples that do not belong to this node


	y=y[my.index.keep]

		########## check if the class membership is the same for all the samples
		######## if all samples belong to the same class, stop

#		if(length(table(y))==1) {
#27/02/2013: changed to stop also if the number of sample is less than min.samples
		if(length(table(y))==1 | length(y)<min.samples) {

		return(list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
		classNode.l=NULL, classNode.r=NULL, 
				#QT.max=NULL, 
				which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=depth, parent=parent, parent.direction=parent.direction, 
				parent.history=parent.history, giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL,
				gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
				QT.rpart.max=NULL
				#, which.variables.used=NULL
				)#,change=change)
		)}


		
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)

	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	#which.index.rc=my.ranks[my.index.keep,]+rep(c(0:( num.var-1)), N)
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s

	
	#y=y[my.index.keep]
	#subsetting the weights
#	my.weights=my.weights[my.index.keep]
	## added
#	my.weights=my.weights/sum(my.weights)

	
##################### end of resizing ##############################
	
######## giniP, parent node Gini is retrieved from the previous nodes

## derive the indexes to retrieve the information of one of the observations
#pp1=sum(my.y.s[,1]*my.w.s[,1])
#pp2=sum(my.y.0[,1]*my.w.s[,1])
#giniP=pp1*(1-pp1)+pp2*(1-pp2)
#giniP=2*pp1*(1-pp1)	

#how many observations remain in the node
N.new=length(my.index.keep)

############ simplified calculation of the gini index ##############
	


################### simplified version for the calculation of the gini index ###################



 #p11=P(Y=0|Node=1), #p12=P(Y=1|Node=1) not calculated explicitly
 p11<-apply(my.w.times.y.0, 2, cumsum)/n.l
 
 #p22=P(Y=0|Node=0), ##p22=P(Y=1|Node=0) not calculated explicitly
 p21=(colSums(my.w.times.y.0)-p11*n.l)/n.r 
#p21.original=(colSums(my.w.times.y.0*my.norm.weights)-p11.original*n.l.original)/n.r.original 
 
 #gini impurity measure for each node
 gini.l=2*(p11*(1-p11))
 gini.r=2*(p21*(1-p21))

 ### QT evaluation using relative frequencies, previous version used this metric
#QT=(p11*(n.l-p11)+p12*(n.l-p12))/n.l+(p21*(n.r-p21)+p22*(n.r-p22))/n.r

 
 ###################################################
 

#N1.l.w.A=sum(my.weights.original[which.n.l][y.original[which.n.l]==1])
#N2.l.w.A=sum(my.weights.original[which.n.l][y.original[which.n.l]==0])

#N1.r.w.A=sum(my.weights.original[which.n.r][y.original[which.n.r]==1])
#N2.r.w.A=sum(my.weights.original[which.n.r][y.original[which.n.r]==0])

#####p(A) after the split, for each node
#pA.l=my.prior[1]*N1.l.w.A/N.l.w + my.prior[2]*N2.l.w.A/N.r.w 
#pA.r=my.prior[1]*N1.r.w.A/N.l.w + my.prior[2]*N2.r.w.A/N.r.w 

############# in comparison to step 1, the weights need to be set equal to their original values
pA.l=(my.prior[1]*(1-p11)*n.l/N.l.w + my.prior[2]*p11*n.l/N.r.w )*my.norm.weights
pA.r=(my.prior[1]*(1-p21)*n.r/N.l.w + my.prior[2]*p21*n.r/N.r.w )*my.norm.weights


#pA.l=my.prior[1]*(1-p11.original)*n.l.original/N.l.w + my.prior[2]*(p11.original)*n.l.original/N.r.w 
#pA.r=my.prior[1]*(1-p21.original)*n.r.original/N.l.w + my.prior[2]*p21.original*n.r.original/N.r.w 


#pA.l=my.prior[1]*(1-p11.original)/N.l.w + my.prior[2]*(p11.original)/N.r.w 
#pA.r=my.prior[1]*(1-p21.original)/N.l.w + my.prior[2]*p21.original/N.r.w 

#pA.l=my.prior[1]*(1-p11)/N.l.w + my.prior[2]*(p11)/N.r.w 
#pA.r=my.prior[1]*(1-p21)/N.l.w + my.prior[2]*p21/N.r.w 

QT.rpart=gini.l*pA.l+gini.r*pA.r


	
#### set to NA all the QT values that derive from duplicated variables values
#QT.rpart[apply(my.tr.s, 2, duplicated, fromLast=TRUE)]=NA
####modified 11/2/2013
##### setting the values of QT from duplicated variable values to very LARGE values, making sure that these values will not be selected (we select the minimum)
QT.rpart[apply(my.tr.s, 2, duplicated, fromLast=TRUE)]=10000


##########setting to extremely small value the QT of the variables already used, to make sure that they are not reselected,
## avoided NA bacause it would cause problems, having an entire column of NA

#unnecessary for step 1, removed from this version of the function the option of not considering the variables that had already been used in the same path of the tree
#if(!is.null(which.variables.used)) QT.rpart[,which.variables.used]=100000
#QT=QT[-N.new,]
#QT=QT[-N.new,]
#removing the last row

#QT.rpart=QT.rpart[-N.new,]
#Modified: 12/4/2013: if the number of samples within the node is only 2 with the previous step we obtain a vector instead of a matrix (only one row remains), retransform to a matrix
QT.rpart=matrix(QT.rpart[-N.new,], ncol=num.var)


#Modified: 11/2/2013: set to a very large number the values that are NaN, necessary because these values are generated when some variables produce a perfect separation, which in not uncommon when considering small subsets of the data
QT.rpart[is.na(QT.rpart)]=10000


 
#evaluate minimum gini (weighted with priors, as in rpart)
#QT.rpart.max=min(QT.rpart, na.rm=T)
#modified 11/2/2013: not necessary to take NA into account

############# added 27/2/2013: in case of variables producing exactly the same split of samples we would like to have the same gini value, this is not always the case beacuse of rounding, it can happen that very small (like e-17) differences arise, to avoid this we round QT.
QT.rpart=round(QT.rpart, 5)

QT.rpart.max=min(QT.rpart, na.rm=T)
	

#check if at level 1, otherwise compute the gain
#if(!is.null(giniP)){
	#check if Gini index is improved 
	gini.improvement=giniP*pA-QT.rpart.max
	
	
	#if(giniP- QT.rpart.max<min.gain)
	if(gini.improvement<min.gain)
		return(list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
		classNode.l=NULL, classNode.r=NULL, 
			#QT.max=QT.max, 
			which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=depth, parent=parent, parent.direction=parent.direction, 
			parent.history=parent.history, giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL, 
			gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
			QT.rpart.max=QT.rpart.max
			#, which.variables.used=NULL
			)#,change=change)
)




#gini.l=(p11*(1-p11)+p12*(1-p12))
#gini.r=(p21*(1-p21)+p22*(1-p22))
#QT=(-(n.l*gini.l+n.r*gini.r)/(n.l+n.r))[-N,]
#PA=prior[1]* n.rA/n.l + prior[2]* n.rA/n.r
#QT.improvement=giniP-(-(n.l*gini.l+n.r*gini.r)/(n.l+n.r))



#which(duplicated(my.tr.s[,1], fromLast=T))



#-((p11*(n.l-p11)+p12*(n.l-p12))/n.l+(p21*(n.r-p21)+p22*(n.r-p22))/n.r)/N



maxid<-apply(QT.rpart,2,which.min)


#####which.var<-which.max(apply(QT,2,max))
#modified for missing values



##################### re-check, to obtain the same results as the original cart, could be coded differently
 which.var<-which.min(apply(QT.rpart,2,min, na.rm=T))
#########################

 # which.var=tail(which(apply(QT.rpart,2,min, na.rm=T)==min(apply(QT.rpart,2,min, na.rm=T))), 1)



thr<-mean(c(my.tr.s[maxid[which.var],which.var],my.tr.s[maxid[which.var]+1,which.var]))

######## derive the indexes of the samples that go in the left/right nodes

#### derive the values of the 
#my.data.train.selected.var= my.tr.s[my.ranks[,which.var],which.var]


#if(!is.null(my.index.keep)){
			#id of the samples that go in node 1, just in this subset
		    which.n.l= which(my.data.train[,which.var]<thr & is.element(1:N, my.index.keep))
#			which.n.l= which(my.data.train.selected.var<thr & is.element(1:N, my.index.keep))
			
			#id of the samples that go in node 0, just in this subset
			which.n.r= which(my.data.train[,which.var]>=thr & is.element(1:N, my.index.keep))
#			} else {which.n.l= which(my.data.train[,which.var]<thr)
#					#id of the samples that go in node 0, just in this subset
#					which.n.r= which(my.data.train[,which.var]>=thr )
#					}


#predicted node for all the samples
#node.pred=ifelse(my.data.train[,which.var]<thr, 1, 0)
node.pred=ifelse(my.data.train[,which.var]<thr, 1, 0)[is.element(1:N, my.index.keep)]
#indicator of prediction in class 1
#idsn.l<-node.pred==1 
#idsn.l<-(node.pred==1 & is.element(1:N, my.index.keep))
 

 
#gini.improvement=giniP*pA-gini.l*pA.l-gini.r*pA.r
#gini.improvement=giniP*pA-QT.rpart.max

################## added to compute probabilities based on original class priors

#P(i=1|Node=1)=prior[1]*P(x in class 1| x in Node 1)/P(x in Node1); estimate: prior*n.l(Node1)/n.l /  (prior1*n.l(Node1)/n.l +  prior2*n.r(Node1)/n.r)
#P(i|A) in the description of rpart. 

#####modified Feb 2013: to correct a bug: situation none of the samples was assigned to one of the subnodes, these values were NAs
#p1Node.l.new=sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]/(sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]+sum(my.weights.original[which.n.l][y.original[which.n.l]==0])*my.prior[2])
if(length(which.n.l)>0) p1Node.l.new=sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]/(sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]+sum(my.weights.original[which.n.l][y.original[which.n.l]==0])*my.prior[2]) else p1Node.l.new=0

#P(i=1|Node=2)
#####modified Feb 2013: to correct a bug: same as for left node
#p1Node.r.new=sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]/(sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]+sum(my.weights.original[which.n.r][y.original[which.n.r]==0])*my.prior[2])
if(length(which.n.r)>0) p1Node.r.new=sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]/(sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]+sum(my.weights.original[which.n.r][y.original[which.n.r]==0])*my.prior[2]) else p1Node.r.new=0



############# determining the label of the class, using predictions for future observations, as in rpart


##node1: checking if there are any samples in this node, sum(node.pred) gives the number of samples classified in node1
##if(sum(node.pred)>0 & sum(node.pred)<N){
##### P(Y=1|Node=1)=P(Y=1|Node=2)
#if (p1Node.l.new==p1Node.r.new){ 
#cln.l<-sample(c(0,1),1) 
#cln.r<-1-cln.l
#} else { 
##### P(Y=1|Node=1)<P(Y=1|Node=2) -> Assign class 1 to node 2
#if (p1Node.l.new<p1Node.r.new) {cln.l<-0
#cln.r<-1} else {
#cln.l<-1
#cln.r<-0}
#}} else {if(sum(node.pred)==0) cln.l=NA 
#if(sum(node.pred==N)) cln.r=NA
##declare this a final node
##my.stop=TRUE - no, it could be a final node in one direction and not in another
#}





#####corrected, Rok (assign according to class prior)

#node1: checking if there are any samples in this node, sum(node.pred) gives the number of samples classified in node1


#### P(Y=1|Node=1)=prior.for.1
if (p1Node.l.new==my.prior[1]){ 
cln.l<-sample(c(0,1),1) 
} else { 
#### P(Y=1|Node=1)<prior.1 -> Assign class 1 to node 2
if (p1Node.l.new<my.prior[1]) cln.l<-0 else cln.l<-1 }

#### P(Y=1|Node=2)=prior.for.1
if (p1Node.r.new==my.prior[1]){ 
cln.r<-sample(c(0,1),1) 
} else { 
#### P(Y=1|Node=2)<prior.1 -> Assign class 1 to node 2
if (p1Node.r.new<my.prior[1]) cln.r<-0 else cln.r<-1 }



#### save the values of gini from each subnode


#index for the selected threshold
which.thr=maxid[which.var]

gini.l=gini.l[which.thr,which.var]
gini.r=gini.r[which.thr,which.var]

pA.l=pA.l[which.thr,which.var]
pA.r=pA.r[which.thr,which.var]

#system.time(res=my.stump())

#should this be a final node: set to TRUE if satisfying the stopping criteria
##### if not stopped earlier
#my.stop=ifelse(length(which.n.l)<min.samples | length(which.n.r)<min.samples, TRUE, FALSE) 

### 27/2/2013: changed to stop only if both nodes have less than min.samples
my.stop=ifelse(length(which.n.l)<min.samples & length(which.n.r)<min.samples, TRUE, FALSE) 



#if(!my.stop & !is.null(giniP)) my.stop=ifelse(QT.rpart.max-giniP<=min.gain, TRUE, FALSE)
#if(!my.stop) my.stop=ifelse(giniP-QT.rpart.max<=min.gain, TRUE, FALSE)

#if(which.var==398) i=a

list(var=which.var, thr=thr, #nNode1=n.l, nNode2=n.r, #p1Node1=p11,p2Node1=p12, p1Node2=p21,p2Node2=p22, 
	classNode.l=cln.l, classNode.r=cln.r, 
	#QT.max=QT.max, 
	which.n.l=which.n.l, which.n.r=which.n.r, my.stop=my.stop, depth=depth, parent=parent, parent.direction=parent.direction, 
	parent.history=parent.history, giniP=giniP, 
	p1Node.l.new=p1Node.l.new, p1Node.r.new=p1Node.r.new,  
	### quantities to check the improvement in gini
	gini.improvement=gini.improvement, gini.l=gini.l, gini.r=gini.r, pA=pA, pA.l=pA.l, pA.r=pA.r, QT.rpart.max=QT.rpart.max, #which.variables.used=c(which.variables.used, 
	which.var
	)#,change=change)

 
 
}
.my.stump.forCART.oneStep.stepOne <-
function(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var,  
		N, n.l, n.r, my.weights, 
		#my.index.keep=NULL, 
		#QT.max.old=NULL,  
		min.samples=min.samples, min.gain=min.gain, #which.variables.used=NULL, depth=1, parent=0, 
		#parent.direction=0, parent.history=NULL, 
		y, my.data.train, #my.weigths.prob, 
		my.prior, N.l.w, N.r.w){
#my.tr.s: sorted training set
#my.y.s: sorted class membership, in matrix format
#my.w.s: sorted weights, in matrix format
#my.y.0=my.y.s==0
#my.w.times.y.0=my.w.s*my.y.0
#my.w.times.y.s=my.w.s*my.y.s
#num.var: number of variables
#N: number of samples
#### removed: #my.index.keep: index of the observations to keep in this node
#### removed: #QT.max.old: QT at previous step
#min.samples: minimum number of samples in each class in order to procede with splitting
#min.gain: minimum gain obtained in order to procede with splitting
#### removed: #depth: current depth of tree, 1 if at first step, etc
#### removed: #parent: id of the parent node

#my.prior: class prior probabilities, vector: prior for class 1, prior for class0, must sum up to 1
#N.l.w: weighted relative frequency of class 1 samples in the complete dataset
#N.r.w: weighted relative frequency of class 0 samples in the complete dataset
#k: index 
#####removed: n1.all: weighted relative frequencies
################## remove the observations not needed, version using the index to LEAVE OUT

#if(length(my.index.keep)<N) {


############### for prior incorporation ##################
#saving the original weights and y 
my.weights.original=my.weights
y.original=y
#probability of being in this node
pA=1



## calculating the gini index for this node (all data)
pp1=sum(y*my.weights)
giniP=2*pp1*(1-pp1)	


############ older version, based on gini improvement calculated using relative frequencies ########################
#p11=P(Y=0|Node=1)
#p11<-apply(my.w.times.y.0, 2, cumsum)#/n1

#modified 19/10
#p12=P(Y=1|Node=1)
#p12<-apply(my.w.times.y.s,2,cumsum)#/n1


#p22=P(Y=1|Node=0)
#p22<-colSums(my.w.times.y.s)-p12


#p21=P(Y=0|Node=0)
#p21=colSums(my.w.times.y.0)-p11




### QT evaluation using relative frequencies, previous version used this metric
#QT=(p11*(n1-p11)+p12*(n1-p12))/n1+(p21*(n.r-p21)+p22*(n.r-p22))/n.r

################### simplified version ###################



 #p11=P(Y=0|Node=1), #p12=P(Y=1|Node=1) not calculated explicitly
 p11<-apply(my.w.times.y.0, 2, cumsum)/n.l
 
 #p21=P(Y=0|Node=0), ##p22=P(Y=1|Node=0) not calculated explicitly
 p21=(colSums(my.w.times.y.0)-p11*n.l)/n.r 
 
 #gini impurity measure for each node
 gini.l=2*(p11*(1-p11))
 gini.r=2*(p21*(1-p21))

 ###################################################
 

pA.l=my.prior[1]*(1-p11)*n.l/N.l.w + my.prior[2]*(p11)*n.l/N.r.w
pA.r=my.prior[1]*(1-p21)*n.r/N.l.w + my.prior[2]*p21*n.r/N.r.w 

###### weighted Gini, calculated as in rpart, using the predictions for future observations
QT.rpart=gini.l*pA.l+gini.r*pA.r
	
#### set to NA all the QT values that derive from duplicated variables values
#QT.rpart[apply(my.tr.s, 2, duplicated, fromLast=TRUE)]=NA

####modified 11/2/2013
##### setting the values of QT from duplicated variable values to very LARGE values, making sure that these values will not be selected (we select the minimum)
QT.rpart[apply(my.tr.s, 2, duplicated, fromLast=TRUE)]=10000

##########setting to extremely LARGE value the QT of the variables already used, to make sure that they are not reselected,
## avoided NA bacause it would cause problems, having an entire column of NA



#unnecessary for step 1
#if(!is.null(which.variables.used)) QT[,which.variables.used]=-100000
#QT=QT[-N.new,]
#potential problem with the following code: it becomes a vector if there are just two samples in the node
#QT.rpart=QT.rpart[-N,]

## Modified 22/4/2013: backtrasform to a matrix in case 
QT.rpart=matrix(QT.rpart[-N,], ncol=num.var)


#Modified: 11/2/2013: set to a very large number the values that are NaN, necessary because these values are generated when some variables produce a perfect separation, which in not uncommon when considering small subsets of the data
QT.rpart[is.na(QT.rpart)]=10000

############# added 27/2/2013: in case of variables producing exactly the same split of samples we would like to have the same gini value, this is not always the case beacuse of rounding, it can happen that very small (like e-17) differences arise, to avoid this we round QT.
QT.rpart=round(QT.rpart, 5)

 
#evaluate minimum gini (weighted with priors, as in rpart)
QT.rpart.max=min(QT.rpart, na.rm=T)
#changed 11/2/2013: not neccessary to take into account the min any more
#QT.rpart.max=min(QT.rpart)


##### stop if the gain is not sufficient
#### 18/2/2013: moved at a later step
#	if(giniP-QT.rpart.max<min.gain)
#		return(list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
#		classNode.l=NULL, classNode.r=NULL, 
#			#QT.max=QT.max, 
#			which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=1, parent=0, parent.direction=0, 
#			parent.history=1, giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL,
#			gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
#			QT.rpart.max=NULL
#			#,which.variables.used=NULL
#			)#,change=change)

						
#			)





#derive the id of the thresholds with the min ID for each variable
maxid<-apply(QT.rpart,2,which.min)


#####which.var<-which.max(apply(QT,2,max))
#modified for missing values#
#which.var<-which.max(apply(QT,2,max, na.rm=T))

#obtain the variable with min Gini
which.var<-which.min(apply(QT.rpart,2,min, na.rm=T))

#obtain the threshold value (as the mean value between the selected value and the successive value
thr<-mean(c(my.tr.s[maxid[which.var],which.var],my.tr.s[maxid[which.var]+1,which.var]))

####### derive the node membership for each sample
	which.n.l= which(my.data.train[,which.var]<thr)
	#id of the samples that go in node 0, just in this subset
	which.n.r= which(my.data.train[,which.var]>=thr )

	
#predicted node for all the samples
node.pred=ifelse(my.data.train[,which.var]<thr, 1, 0)#[is.element(1:N, my.index.keep)]

#gini.improvement=giniP*pA-gini.l*pA.l-gini.r*pA.r
gini.improvement=giniP*pA-QT.rpart.max

################## added to compute probabilities based on original class priors

#P(i=1|Node=1)=prior[1]*P(x in class 1| x in Node 1)/P(x in NodeLeft (1)); estimate: prior*n.l(Node1)/n.l /  (prior1*n.l(Node1)/n.l +  prior2*n.r(Node1)/n.r)
#P(i|A) in the description of rpart. 

#####modified Feb 2013: to correct a bug: situation none of the samples was assigned to one of the subnodes, these values were NAs

#p1Node.l.new=sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]/(sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]+sum(my.weights.original[which.n.l][y.original[which.n.l]==0])*my.prior[2])
if(length(which.n.l)>0) p1Node.l.new=sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]/(sum(my.weights.original[which.n.l][y.original[which.n.l]==1])*my.prior[1]+sum(my.weights.original[which.n.l][y.original[which.n.l]==0])*my.prior[2]) else p1Node.l.new=0

#P(i=1|Node=2(right))
#p1Node.r.new=sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]/(sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]+sum(my.weights.original[which.n.r][y.original[which.n.r]==0])*my.prior[2])
if(length(which.n.r)>0) p1Node.r.new=sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]/(sum(my.weights.original[which.n.r][y.original[which.n.r]==1])*my.prior[1]+sum(my.weights.original[which.n.r][y.original[which.n.r]==0])*my.prior[2]) else p1Node.r.new=0


############# determining the label of the class, using predictions for future observations, as in rpart


##node1: checking if there are any samples in this node, sum(node.pred) gives the number of samples classified in node1
#if(sum(node.pred)>0 & sum(node.pred)<N){
##### P(Y=1|Node=1)=P(Y=1|Node=2)
#if (p1Node.l.new==p1Node.r.new){ 
#cln.l<-sample(c(0,1),1) 
#cln.r<-1-cln.l
#} else { 
##### P(Y=1|Node=1)<P(Y=1|Node=2) -> Assign class 1 to node 2
#if (p1Node.l.new<p1Node.r.new) {cln.l<-0
#cln.r<-1} else {
#cln.l<-1
#cln.r<-0}
#}} else {if(sum(node.pred)==0) cln.l=NA 
#if(sum(node.pred==N)) cln.r=NA
##declare this a final node
##my.stop=TRUE - no, it could be a final node in one direction and not in another
#}


######### corrected,  using the threshold based on priors
#### P(Y=1|Node=1)=prior.for.1
if (p1Node.l.new==my.prior[1]){ 
cln.l<-sample(c(0,1),1) 
} else { 
#### P(Y=1|Node=1)<prior.1 -> Assign class 1 to node 2
if (p1Node.l.new<my.prior[1]) cln.l<-0 else cln.l<-1 }

#### P(Y=1|Node=2)=prior.for.1
if (p1Node.r.new==my.prior[1]){ 
cln.r<-sample(c(0,1),1) 
} else { 
#### P(Y=1|Node=2)<prior.1 -> Assign class 1 to node 2
if (p1Node.r.new<my.prior[1]) cln.r<-0 else cln.r<-1 }



#### save the values of gini from each subnode

#index for the selected threshold
which.thr=maxid[which.var]

gini.l=gini.l[which.thr,which.var]
gini.r=gini.r[which.thr,which.var]

pA.l=pA.l[which.thr,which.var]
pA.r=pA.r[which.thr,which.var]

#system.time(res=my.stump())

#should this be a final node: set to TRUE if satisfying the stopping criteria
##### if not stopped earlier, ok to set to true if it needs to be stopped, this node will still be used
#my.stop=ifelse(length(which.n.l)<min.samples | length(which.n.r)<min.samples , TRUE, FALSE) 

### 18/2/2013 moved here the stop due to small gain in Gini index
#my.stop=ifelse(length(which.n.l)<min.samples | length(which.n.r)<min.samples | giniP-QT.rpart.max<min.gain, TRUE, FALSE) 

### 27/2/2013: changed to stop only if both nodes have less than min.samples
my.stop=ifelse((length(which.n.l)<min.samples & length(which.n.r)<min.samples) | giniP-QT.rpart.max<min.gain, TRUE, FALSE) 

#if(!my.stop & !is.null(QT.max.old)) my.stop=ifelse(QT.rpart.max-giniP<=min.gain, TRUE, FALSE)
#if(!my.stop) my.stop=ifelse(giniP-QT.rpart.max<=min.gain, TRUE, FALSE)  # removed, checked earlier

list(var=which.var, thr=thr, #nNode1=n.l, nNode2=n.r, #p1Node1=p11,p2Node1=p12, p1Node2=p21,p2Node2=p22, 
	classNode.l=cln.l, classNode.r=cln.r, 
	#QT.max=QT.max, 
	which.n.l=which.n.l, which.n.r=which.n.r, my.stop=my.stop, 
	 depth=1, parent=0, parent.direction=0, parent.history=1, 
	giniP=giniP, 
	p1Node.l.new=p1Node.l.new, p1Node.r.new=p1Node.r.new,  
	### quantities to check the improvement in gini
	gini.improvement=gini.improvement, gini.l=gini.l, gini.r=gini.r, pA=pA, pA.l=pA.l, pA.r=pA.r, QT.rpart.max=QT.rpart.max 
	#, which.variables.used=which.var
	)#,change=change)


	
}



.permute.rows <-
function(x){
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}



.convert.node=function(k, my.depth){
######## function that converts the node index of a CART object in a decimal number: depth.position 
###### my.depth=$depth: depth of the node 
	#my.depth+(k-2^(my.depth-1)+1)*.1
	paste(my.depth, (k-2^(my.depth-1)+1), sep=".")
}




.convert.node.depth.unknown=function(k, max.depth=1000){
######## function that converts the node index of a CART object in a decimal number: depth.position 
###### to be used if the depth of the node is not known
	#my.depth+(k-2^(my.depth-1)+1)*.1

	my.depth=numeric(length(k))
	#calculates the depth of the node, 
	for(kk in 1:length(k))
	my.depth[kk]= which(k[kk]-(c(2^c(1:max.depth))[1:(max.depth)])<0)[1]

	paste(my.depth, (k-2^(my.depth-1)+1), sep=".")
}



.convert.node.depth.unknown.code0=function(k, max.depth=1000){
######## function that converts the node index of a CART object in a decimal number: depth.position 
###### to be used if the depth of the node is not known
	#my.depth+(k-2^(my.depth-1)+1)*.1
################# gives wrong results, do not use! bug fixed in function .convert.node.depth.unknown.code2
	
	my.depth=my.pos=numeric(length(k))
	my.names=rep(1, length(k))
	
	#calculates the depth of the node, 
	for(kk in 1:length(k)){
		my.depth[kk]= which(k[kk]-(c(2^c(1:max.depth))[1:(max.depth)])<0)[1]
		my.pos[kk]=k[kk]-2^(my.depth[kk]-1)+1
	
	if(my.depth[kk]==1) my.names[kk]=1 else{
		for(i in 1:(my.depth[kk]-1)){
		#my.names[kk]=paste(my.names[kk], ifelse( my.pos[kk]<= (2^(my.depth[kk]-1)/2^i), paste(i+1, "L", sep=""), paste(i+1, "R", sep="")), sep=".")
		my.pos.depth=k[kk]-2^(i-1)+1
		my.names[kk]=paste(my.names[kk], ifelse(k[kk] %% 2, paste(i+1, "R", sep=""), paste(i+1, "L", sep="")), sep=".")
		
		}}
	
	}
	
	
	
	my.names
}




.convert.node.depth.unknown.code2=function(k, max.depth=1000){
######## function that converts the node index of a CART object in a decimal number: depth.position 
###### to be used if the depth of the node is not known
	#my.depth+(k-2^(my.depth-1)+1)*.1
#updated 23/4/2013: fixed a bug
	
	my.depth=my.pos=numeric(length(k))
	my.names=rep(1, length(k))
	
	#calculates the depth of the node, 
	for(kk in 1:length(k)){
		my.depth[kk]= which(k[kk]-(c(2^c(1:max.depth))[1:(max.depth)])<0)[1]
		my.pos[kk]=k[kk]-2^(my.depth[kk]-1)+1
	

	if(my.depth[kk]==1) my.names[kk]=1 else{

	my.names.matrix=matrix(1, ncol=2^(my.depth[kk]-1), nrow=my.depth[kk])
	
		for(i in 2:my.depth[kk]){
		#my.names[kk]=paste(my.names[kk], ifelse( my.pos[kk]<= (2^(my.depth[kk]-1)/2^i), paste(i+1, "L", sep=""), paste(i+1, "R", sep="")), sep=".")
				my.names.matrix[i,]= rep(rep(c(paste(i, "L", sep=""), paste(i, "R", sep="")), each=2^(my.depth[kk]-i)), times= 2^(my.depth[kk]-(my.depth[kk]-i)-2))
		}
		my.names[kk]= paste(my.names.matrix[,my.pos[kk]], collapse=".")
	
		
		}
	
	}
	
	
	
	my.names
} # end function 



	.check.training.data=function(my.data.train, y, my.weights, my.prior, max.depth){
	################ function that checks the correctness of the training data
	#' @param my.data.train training data
	#' @param y class membership vector of the training data
	
			#missing values in the training data
			if (any(is.na(my.data.train))) { stop("There are some missing values in the training data. The handling of missing value is not implemented yet. Remove all the variables with missing values or impute the missing values.")
											 #stop()
											}
			#different number of samples in y and my.data.train 								
			if (dim(my.data.train)[1]!=length(y)) { stop("The number of samples included in the training data is not equal to the length of the class membership vector.")
											 #stop()
											}								

			#the number of classes is different than 2
			if (length(names(table(y)))!=2)  { stop("The function handles binary classification problems. Samples of two different classes should be included in the training data set.")
											 #stop()
											}								

			#coding of class membership, must be 0 or 1
			if (all(names(table(y))!=c("0","1")))  { stop("The class membership of the training data must be coded with 0 and 1.")
											 #stop()
												}

			#coding of class membership, must be 0 or 1
			if (class(y)!="numeric" & class(y)!="integer")  { stop("The class membership of the training data must be given  as a numeric or integer vector (not a factor or a character vector).")
											 #stop()
												}

												
			#any missing values in the y
			if (any(is.na(y)))  { stop("There are missing values in the class membership vector.")
											 #stop()
												}

			#weights not of the same length as y
			if (length(y)!=length(my.weights))  { stop("The length of the vector with weights is not equal to the number of samples in the training data set.")
											 #stop()
												}
												
			#the vector with prior not of length 2
			#check if priors were provided
			if (!is.null(my.prior)){
			#check if priors are a vector or lenght 2
				if (length(my.prior)!=2)  { stop("The class priors must be expressed as a vector with two elements.")
											 #stop()
										}#end if
												
				if (sum(my.prior)!=1)  { stop("The class priors must sum to 1.")
											 #stop()
										}#end if
								
								
			}#end !is.null(my.prior)	

			#check if the maximum depth of the tree is positive
			if (max.depth<=0) { stop("The maximum depth of the tree must be a positive integer")
											 #stop()
			
			
							}# end if(max.depth)	

			
			return()

	}#end .check.training.data


	
	.check.test.data=function(my.data.test, y.test, my.data.train=NULL){
	################ function that checks the correctness of the test data
	#' @param my.data.test test data
	#' @param y.test class membership vector of the test set
	#' @param my.data.train training data

			#missing values in the training data
			if (any(is.na(my.data.test))) {stop("There are some missing values in the test data. The handling of missing value is not implemented yet. Remove all the variables with missing values or impute the missing values.")
											 #stop()
											}
		
			#different number of variables in my.data.test and my.data.train 								
			if(!is.null(my.data.train)){
								if (dim(my.data.test)[2]!=dim(my.data.train)[2]) {stop("The number of variables included in the training and test data is not equal.")
											 #stop()
											}								
									}
											
											
			#the number of classes is different than 2
			#if( length(names(table(y)))!=2)  {cat("The function handles binary classification problems. Samples of two different classes should be included in the training data set.")
			#								return()
			#								}								

			if(!is.null(y.test)){
						#different number of samples in y and my.data.train 								
						if (dim(my.data.test)[1]!=length(y.test)) {stop("The number of samples included in the test data is not equal to the length of the class membership vector.")
														 #stop()
														}								

						
						
						#coding of class membership, must be 0 or 1
						if (all(names(table(y.test))!=c("0","1")))  {stop("The class membership of the test data must be coded with 0 and 1.")
														 #stop()
															}
															
						#coding of class membership, must be 0 or 1
						if (class(y.test)!="numeric" & class(y.test)!="integer")  { stop("The class membership of the test data must be given  as a numeric or integer vector (not a factor or a character vector).")
														 #stop()
															}
												

						#any missing values in the y
						#if( any(is.na(y)))  {cat("There are missing values in the class membership vector.")
						#								return()
						#									}

			}					
										
	return()

	}#end .check.test.data



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

.cart.for.cv.ee.Icv.internal.with.arguments <-
function(my.fold, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, num.ee, num.boost, nfolds.Icv, verbose=FALSE, my.env=NULL){

#internal function that performs cross-validated class prediction using easy ensamble - WITH internal cross-validation for the evaluation of the prediction error within the boosting step



	#################################################################################
	#################### step 1: reduce the data set ################################
	#################################################################################
	
	if (verbose==TRUE) {	#cat("Fold: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Fold: ", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}
							
	####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out

		#########which samples to keep in the cross-validated step
		my.index.keep=c(1:N)[-my.fold]

	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N
	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)
	

	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#my.data.test=my.data.train[-my.index.keep,]
	#23/4/2013: fixed bug: keep a matrix format even if only one sample is included in the test set - as happens in LOO-CV
	my.data.test=matrix(my.data.train[-my.index.keep,], ncol=num.var, byrow=F)
	#class membership of the left out samples
	y.test=y[-my.index.keep]	
	
	num.samples.test=nrow(my.data.test)
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed for EE, re-evaluation of the order of the variables ###################################
	
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again
	#24/4/2013: check if it is needed or the original indexes are needed!!!!


	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################

			########################## fit cart on original data ##########################		
			#results.cart.cv=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		





####################### prediction for left out samples ##########################


#prediction.cv=my.predict.cart(results.cart.cv, my.data.test)
#return(list(prediction.cv, results.cart.cv))






	#################################################################################
	#################### step 2: Perform EE ################################
	#################################################################################
	


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
		
		#10/6/2013: added for debugging
		#set.seed(1)
		
		my.folds.ee=lapply(1:num.ee, function(i) sample(which.maj.class, abs(N1-N0)))
		#### 18/3/2013: added the sorting of the folds, the indexes within each fold are sorted
		my.folds.ee=lapply(my.folds.ee, sort)

		##### 18/3/2013: here it is not necessary to re-order the output because the test set is separate
		#the output score is diff.alpha
		#8/5/2013: added the version for Icv
		my.scores=matrix(unlist(lapply(my.folds.ee, .cart.for.ee.Icv.internal.with.arguments, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, nfolds.Icv)), ncol=num.samples.test, byrow=T)
		#.cart.for.ee.interal.with.arguments=function(my.fold, 								    	N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s){

		
		#derive the global score for the test - left out - samples, combining the 
		my.scores.final=apply(my.scores, 2, sum)

		## predicted class membership for the samples that were left out from cross-validation
		y.pred=ifelse(my.scores.final>0, 1, 0)

		#10/6/2013: added for debugging
		#set.seed(1)
		
		num.ties=sum(my.scores.final==0, na.rm=T)
		if(num.ties>0)  y.pred[which(my.scores.final==0)]=sample(c(0,1), num.ties, replace=T)
	
		
		
		## predicted class membership for the samples that were left out from cross-validation
		#### removed 7/5/2013: y.pred=ifelse(my.scores.final>num.ee/2, 1, 0) the classification rule is based on the sum of the diff.alpha only
		
		#are there any ties?  if there are ties assign the class at random
		#num.ties=sum(my.scores.final==num.ee/2, na.rm=T)
		#if(num.ties>0)  y.pred[which(my.scores.final==num.ee/2)]=sample(c(0,1), num.ties, replace=T)
	
		
		
		############## end added for multiple downsizing ####################



#return(y.pred)
#updated 7/5/2013: returns also the scores to the main function, not only the class prediction
return(list(y.pred=y.pred, scores.final=my.scores.final))
}# end of .cart.for.cv.ee.Icv.internal.with.arguments







.cart.for.ee.Icv.internal.with.arguments <-
#browser()
function(my.fold, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, nfolds.Icv, verbose=FALSE, my.env=NULL){

		####same as .cart.for.ee.internal.with.arguments but uses internal CV
		#my.fold: IDs of the samples to leave out

		#print the EE progress on screen and increment the internal index
		if(verbose==TRUE) {	#cat("Easy ensemble step: ", .internal.index, "\n")
							#.internal.index<<-.internal.index+1
							cat("Easy ensemble step:", my.env$.internal.index, "\n")
							my.env$.internal.index=my.env$.internal.index+1
							}

		
		#########which samples to keep in the cross-validated step
		
		my.index.keep=c(1:N)[-my.fold]
	
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)


	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#here the test 
	#my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			


		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################


	
	
			###################### step 1: evaluate cart############################### 
		
			########################## fit cart on original data ##########################		
			results.cart.ee=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
			########################## fit cart on original data ##########################		

	

			#prediction on the same data used to develop the predictor, the balanced training set, it is! cross-validated
			#folds for internal cross-validation
			
			#10/6/2013: added for debugging
			#set.seed(1)
		
			if(nfolds.Icv<num.samples) {my.folds.Icv=.balanced.folds(y, nfolds.Icv) #
						my.folds.Icv=lapply(my.folds.Icv, sort)} else my.folds.Icv=as.list(1:N)

			
			#prediction.ee=unlist(lapply(my.folds.Icv, .cart.for.cv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
			prediction.ee=unlist(lapply(my.folds.Icv, .cart.for.cv.prediction.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]

			#prediction.ee=my.predict.cart(results.cart.ee, my.data.train) - version without CV




#################################### boosting part ###########################



 #matrix storing the predictions for new data
 prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
#alpha from boosting
	alpha=numeric(num.boost)

#saving the cart fit from the original data
	results.cart.original=results.cart.ee

	#results.cart.boosting=results.cart
 
 #initializing the indexes for the boosting

#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.ee.boosting=prediction.ee
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

		
#cat(m, " boosting \n")	


#cat("Boosting iteration #: ", m, "\n") 



	############ naive error of the classifier
	#indicator variable, indicating if the EE-prediction was correct or not, logical indicator, 1 if wrong prediction, 0 otherwise
	incorrect.class.ee=prediction.ee.boosting!=y

	#calculates weighted error from the EE-estimate
	#weighted.error=sum(my.weights.boosting[incorrect.class.ee])
	weighted.error=sum(my.weights[incorrect.class.ee])
	
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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
	prediction.new[,m]=my.predict.cart(results.cart.ee, my.data.test)

		
	#go on to the next step of boosting
	if(!my.stop.boosting){
				
						#cat("Estimating\n")										
						#calculate alpha, 										
						#alpha, or 2*beta from Izenman, Adaboost algorith from page 513, 
						#alpha is positive if the error is small (<0.5), negative if it is large (>0.5)
						#alpha[m]=log((1-weighted.error)/weighted.error), set to a small constant if the error is large
						
												
						#update the weights
						#my.weights.boosting=ifelse(incorrect.class.cv, my.weights.boosting*(1-weighted.error)/weighted.error, my.weights.boosting)
						my.weights=ifelse(incorrect.class.ee, my.weights*(1-weighted.error)/weighted.error, my.weights)
						#renormalize the weights
						#my.weights.boosting=my.weights.boosting/sum(my.weights.boosting)
						my.weights=my.weights/sum(my.weights)
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
						############### fit the cart model and obtain the (non cross validated) class prediction #####################											
							
							#using name results.cart because it is the default argument of cart.for.cv.internal
							#results.cart=cart.internal()
							results.cart.ee=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
															          

							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.ee.boosting=my.predict.cart(results.cart.ee, my.data.train) -version without CV
							
							#CV-class prediction
							
							#10/6/2013: added for debugging
							### correct to sample at this level? 
							#set.seed(1+m)
							#set.seed(1)
							
							if(nfolds.Icv<num.samples) {my.folds.Icv=.balanced.folds(y, nfolds.Icv) #
								my.folds.Icv=lapply(my.folds.Icv, sort)} else my.folds.Icv=as.list(1:N)


							#prediction.ee.boosting=unlist(lapply(my.folds.Icv, .cart.for.cv.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
							prediction.ee.boosting=unlist(lapply(my.folds.Icv, .cart.for.cv.prediction.internal.with.arguments, N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.Icv))]
			

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
	
	###prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)
	
	#are there any ties?  if there are ties assign the class at random
	###num.ties=sum(diff.alpha==0, na.rm=T)
	###if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	
	
#calculates the sum of alpha for those predicted 1, molti
#	prediction.new.data.boosting=unlist(lapply(c(1:num.samples.test), function(i) sum(prediction.new[i,]*alpha, na.rm=T)))
#prediction.new.data.boosting=ifelse(prediction.new.data.boosting>all.alpha-prediction.new.data.boosting, 1, 0)





########### if the class membership of the test set is provided, calculate the error rate for each number of boosting iterations
########## removed: not used in the output #####################
#if(!is.null(y.test)) {

 #	cumsum.alpha.half=cumsum(alpha)/2
	
	#######boosting error on the test set, one value for each possible value of boosting iterations, for example: the second value would be the average error obtained using 2 boosting iterations
# if(m==1) error.boosting=mean(y.test!=prediction.new[,1]) else{
#	error.boosting=c(#first step 
#			mean(y.test!=prediction.new[,1]),
#			#next steps
#			unlist(lapply(2:num.boost, function(m) {

#								diff.alpha=rowSums(mat.prediction.new.y1.times.alpha[,1:m])-cumsum.alpha.half[m]
#								prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)

								#are there any ties?  if there are ties assign the class at random
#								num.ties=sum(diff.alpha==0, na.rm=T)
#								if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
								
#								mean(prediction.new.data.boosting!=y.test)

 
#							}#end function(m) 
#				)#end lapply
#			)#end unlist
#		)#end error.boosting
#}# end else m==1
#} else {error.boosting<-NA}


########### added for EE, calculates the score, defined as diff.alpha, a score is given for each sample 
#return(list(score=diff.alpha,ff=mat.prediction.new.y1.times.alpha))
return(score=diff.alpha)


#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
##return(list(results.cart=results.cart.original, prediction=prediction.cv, ##prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#prediction.boosting=prediction.new, 
#m: number of performed boosting iterations
##weights.after.boosting=my.weights, ##error.boosting=error.boosting,prob=diff.alpha,pred.new.i=prediction.new, m=m))







}# end .cart.for.ee.Icv.internal.with.arguments

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#################################### added 19/6/2013: try to see if different weighting update might help in EE
	
	

.cart.for.ee.cv.internal.w2.with.arguments <-
function(my.fold, N, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain, my.ranks, y.test, num.samples.test, num.boost, my.data.test, my.ind.s, nfolds){

		####function that performs cross-validation within cart function 
		#my.fold: IDs of the samples to leave out
		
		#### the used data are balanced undersampled training sets, since we are using easy ensemble

		#########which samples to keep in the EE 	step
		my.index.keep=c(1:N)[-my.fold]

	
	
	################### step 0: resize the data#############################################
	
	############ using the resizing script, originally used within the cart function (for nodes>1)
	
	### indexes of the observations that need to be removed, referring to the left out observations
	#which.index.rc=t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N

	###  28/2/2013: fixed bug: transformed which.index.rc from matrix (num.var by number of samples to remove) to vector form; an error would appear if the matrix had EXCACTLY two columns (i.e., if exactly two samples had to be removed), 
	### this happened because the matrix would be interpreted as the row, column indeces to remove, and not as the individual observations to be removed
	which.index.rc=as.numeric(t(my.ranks[-my.index.keep,])+seq(0, num.var-1)*N)


	############## updating the quantities already calculated for the complete data
	
	#remove the observations
	my.tr.s=matrix(my.tr.s[-which.index.rc], ncol=num.var, byrow=F)
	my.y.s=matrix(my.y.s[-which.index.rc], ncol=num.var, byrow=F)

	#uses the complete weights!
	#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var)
	#subsets
	#saving also the original weights, without renormalization
	my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
	my.norm.weights=sum(my.w.s[,1])  
	#renormalized weights
	my.w.s=my.w.s/my.norm.weights
	
	my.y.0=matrix(my.y.0[-which.index.rc], ncol=num.var, byrow=F)

	#recalculate "weighted sample size", weighted relative frequencies
	n.l<-apply(my.w.s,2,cumsum)
	n.r<-1-n.l
	
	my.w.times.y.0=my.w.s*my.y.0
	my.w.times.y.s=my.w.s*my.y.s
		
	
	
	
	############ reduce the other quantites needed; 
	#test set
	#here the test 
	#my.data.test=my.data.train[-my.index.keep,]
	
	
	#training set
	
	my.data.train=my.data.train[my.index.keep,]
	#weights
	my.weights=my.weights[my.index.keep]/my.norm.weights
	#class membership
	y=y[my.index.keep]	
	
	#### saving the number of samples and number of variables
	#number of samples, recalculating
	N=num.samples=nrow(my.data.train)
		#number of variables
		#num.var=ncol(my.data.train)

	
	############empirical prior
	my.prior=c(sum(my.weights[y==1]), sum(my.weights[y==0]))

	############## to include the prior in the calculation of the gini gain
	#weighted relative frequencies, for class 1 and class 2
	N.l.w=sum(my.weights[y==1])
	N.r.w=1-N.l.w

	
	
			
		
		
		
			
	################################### needed???????????????????? ###################################
	#producing and saving the indexes only for ordering the training data - each variable (column) is sorted. (efficient coding)
	###############my.ind.s<-matrix(unlist(lapply(1:num.var , function(i) order(my.data.train[,i]))), ncol=num.var) #not needed again



	## derives the ranks - to backtrasform the sorted data into the original data, needed! needs to be recalculated
	#my.ranks=t(colRanks(my.data.train))
	my.ranks=apply(my.data.train, 2, rank, ties="first")
	################################### needed???????????????????? ###################################



####################### beginning of the cart computation on the CV=training set #############################


		
		#if(length(table(y))==1) { #unnecessary check for EE! data are class-balanced
		#return(list(var=NULL, thr=NULL, #nNode1=NULL, nNode2=NULL, #p1Node1=NULL,p2Node1=NULL, p1Node2=NULL,p2Node2=NULL, 
		#classNode.l=NULL, classNode.r=NULL, 
		#		#QT.max=NULL, 
		#		which.n.l=NULL, which.n.r=NULL, my.stop=NA, depth=depth, parent=parent, parent.direction=parent.direction, 
		#		parent.history=parent.history, giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL,
		#		gini.improvement=NULL, gini.l=NULL, gini.r=NULL, pA=NULL, pA.l=NULL, pA.r=NULL,
		#		QT.rpart.max=NULL
		#		#, which.variables.used=NULL
		#		)#,change=change)
		#)#}



#step 1

#my.stump.mat1.Lara.duplicates.faster.internal<-function(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var,  
#		N, n.l, n.r, my.weights, my.index.keep=NULL, 
#		QT.max.old=NULL,  min.samples=2, min.gain=0.01, which.variables.used=NULL, dep
		
		
		
		
		
		########## perform cross-validated CART con the subset of the data, in this case on a balanced (easy ensemble) subset
		
		

		
		############################################################################
		############ step 1: estimate CART on the complete easy ensemble data
		############################################################################


		results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
		
												  
		


#saving the cart fit from the original data
	results.cart.original=results.cart


################################### step 2: obtain cross-validated class predictions #########################################
	
			#my.folds=.balanced.folds(y, nfolds)
			########## add an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list
			##### my.folds.internal: a list containing the IDs of the samples to be used in the 
			
			#10/6/2013: added for debugging
			#set.seed(1)
			
			if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
			######## convert the selected indexes in the original indexes, problem: .balanced.folds does produce indexes that are in the range: 1 to num.samples (reduced by ee), must be re-transfromed in the original indexes to be used together with the old indeces.
#			my.folds.internal.converted=lapply(my.folds.internal, function(x) my.index.keep[x])

			
			#### add to the samples to remove for CV, those that really need to be cross-validated are reported first 
#			my.folds.internal.all=lapply(my.folds.internal.converted, function(fold) c( fold, my.fold )) 
			

			################ obtain cross-validated prediction using nfolds-CV
			##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
			#prediction.cv=unlist(lapply(my.folds.internal, cart.for.cv.internal))[order(unlist(my.folds.internal))]
			###### updated 4/3/2013: the function cart.for.cv.internal was upadted to allow to arguments: the samples to leave out (all those that were not used in the easy-ensemble fold (my.fold[[i]]) and those that are left out because of the internal CV (my.folds.interal[[i]]), the evaluation is then made only for the samples left out for CV
			
#			prediction.cv=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]])[1:length(my.folds.internal[[i]])]))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(my.folds.internal, cart.for.cv.internal))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]], my.folds.internal[[i]] )))[order(unlist(my.folds.internal))]
			#prediction.cv=unlist(lapply(my.folds, cart.for.cv.internal))

			#prediction.cv=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
			prediction.cv=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]

			
			
			if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
			#prediction.cv.boosting2=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
			prediction.cv.boosting2=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]


	if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
			#prediction.cv.boosting3=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
			prediction.cv.boosting3=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
			

#################################### step 3: boosting part ######################################################################



 
	#matrix storing the predictions for new data
	prediction.new=matrix(NA, ncol=num.boost, nrow=num.samples.test)
	#alpha from boosting
	alpha=numeric(num.boost)
	#results.cart.boosting=results.cart
 
	#initializing the indexes for the boosting

	#cat("Boosting interation #: ", 1, "\n") 



 m=1
 my.stop.boosting=FALSE
 prediction.cv.boosting=prediction.cv
 
 #saving the original weights
 my.weights.original=my.weights
 
 ##################my.weights=my.weights[]
 
 
 #initialize the value of the weights, set equal to the initial weights at the beginning of boosting
 #my.weights.boosting=my.weights
 

################## step one of boosting, saving the results

#saving the prediction for new samples for the first step	
#prediction.new[,1]=my.predict.cart(results.cart.original, my.data.test)

 
 
############################# boosting ##################################
while(m<=num.boost & my.stop.boosting==FALSE){
#for(m in 1:num.boost){

	


#cat("Boosting interation #: ", m, "\n") 



	############ naive error of the classifier
	#indicator variable, indicating if the CV-prediction was correct or not, logical indicator, 1 if wrong prediction, 0 otherwise
	#incorrect.class.cv=prediction.cv.boosting!=y
	#incorrect.class.cv=prediction.cv.boosting!=y & prediction.cv.boosting2!=y
	incorrect.class.cv=prediction.cv.boosting!=y | prediction.cv.boosting2!=y | prediction.cv.boosting3!=y

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
	

		

		
		
		
	#saving the prediction for new samples for the first step	
	#prediction.new[,m]=my.predict.cart(results.cart, my.data.test)
	#results.cart.cv contains the cart on the original, balanced data of this step
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
						#my.weights=ifelse(incorrect.class.cv, my.weights*((1-weighted.error)/weighted.error)^(1/2), my.weights)
						my.weights=ifelse(incorrect.class.cv, my.weights*((1-weighted.error)/weighted.error), my.weights)
						#renormalize the weights
						#my.weights.boosting=my.weights.boosting/sum(my.weights.boosting)
						my.weights=my.weights/sum(my.weights)
																	
											
						#### create a vector of weights of the same length of the original ones, 
						###### changed: 20/2/2013: to fix a bug present 
						
						my.weights.original.length=numeric(num.samples+length(my.fold))
						my.weights.original.length[my.index.keep]=my.weights
						

					
							############### update of the matrices used in cart#####################
																	
																	
							#new weights, in matrix format
							#my.w.s<-matrix(my.weights.boosting[my.ind.s], ncol=num.var)
							#my.w.s<-matrix(my.weights[my.ind.s], ncol=num.var), not correct because it uses my.ind.s, which was created for the complete data set
							
							my.w.s=matrix(my.weights.original.length[my.ind.s], ncol=num.var)
							my.w.s=my.w.s=matrix(my.w.s[-which.index.rc], ncol=num.var, byrow=F)
								#my.norm.weights=sum(my.w.s[,1])  
								#renormalized weights
								#my.w.s=my.w.s/my.norm.weights
							
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
							
							#using name results.cart because it is the default argument of cart.for.cv.internal, 
							#also here we need to restrict the attention to the subset of data
							results.cart=.cart.internal.with.arguments(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var, N, n.l, n.r, my.weights,  min.samples, min.gain, y, my.data.train,  my.prior, N.l.w, N.r.w, max.depth, my.y.0, my.ranks)
																	  

							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							### removed, does not do CV if(nfolds<num.samples) my.folds=.balanced.folds(y, nfolds) else my.folds=as.list(1:num.samples)

							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
							#prediction.cv.boosting=unlist(lapply(my.folds, cart.for.cv.internal))[order(unlist(my.folds))]   
						

							#prediction.cv.boosting=my.predict.cart(results.cart, my.data.train)

							
							#my.folds=.balanced.folds(y, nfolds), 
							########## added an option for loo-cv, in this case it does not use balaced folds to produce the folds, but makes the complete list; new folds are derived for each boosting step
							#10/6/2013: added for debugging
							#set.seed(1+m)
							#set.seed(1)
							
							if(nfolds<num.samples) my.folds.boosting=.balanced.folds(y, nfolds) else my.folds.boosting=as.list(1:num.samples)
							my.folds.boosting=lapply(my.folds.boosting, sort)	
								
							######## convert the selected indexes in the original indexes, problem: .balanced.folds does produce indexes that are in the range: 1 to num.samples (reduced by ee), must be re-transfromed in the original indexes to be used together with the old indeces.
	#						my.folds.internal.converted=lapply(my.folds.boosting, function(x) my.index.keep[x])
							
							#### add to the samples to remove for CV, those that really need to be cross-validated are reported first 
	#						my.folds.internal.all=lapply(my.folds.internal.converted, function(fold) c( fold, my.fold )) 
	
							
							################ obtain cross-validated prediction using nfolds-CV
							##############  needs reordering at the end of the CV, to obtain the same order of the prediction as in the original data set, 
							####prediction.cv.boosting=unlist(lapply(my.folds.boosting, cart.for.cv.internal))[order(unlist(my.folds.boosting))]   
							
							###### updated 4/3/2013: the function cart.for.cv.internal was upadted to allow to arguments: the samples to leave out (all those that were not used in the easy-ensemble fold (my.fold[[i]]) and those that are left out because of the internal CV (my.folds.interal[[i]]), the evaluation is then made only for the samples left out for CV
	#						prediction.cv.boosting=unlist(lapply(1:nfolds, function(i) cart.for.cv.internal(my.folds.internal.all[[i]])[1:length(my.folds.boosting[[i]])]))[order(unlist(my.folds.boosting))]


					#prediction.cv.boosting=unlist(lapply(my.folds.boosting, .cart.for.cv.internal.with.arguments, N=N, my.ranks=my.ranks, num.var=num.var, my.tr.s=my.tr.s, my.y.s=my.y.s, my.w.s=my.w.s, my.y.0=my.y.0,my.data.train=my.data.train,  my.weights=my.weights, y=y, min.samples=min.samples, max.depth=max.depth, min.gain))[order(unlist(my.folds.boosting))] 

					prediction.cv.boosting=unlist(lapply(my.folds.boosting, .cart.for.cv.prediction.internal.with.arguments, N=N, my.ranks=my.ranks, num.var=num.var, my.tr.s=my.tr.s, my.y.s=my.y.s, my.w.s=my.w.s, my.y.0=my.y.0,my.data.train=my.data.train,  my.weights=my.weights, y=y, min.samples=min.samples, max.depth=max.depth, min.gain))[order(unlist(my.folds.boosting))] 

	

			if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
			#	prediction.cv.boosting2=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
		prediction.cv.boosting2=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]

	
	if(nfolds<num.samples) my.folds.internal=.balanced.folds(y, nfolds) else my.folds.internal=as.list(1:num.samples)
			my.folds.internal=lapply(my.folds.internal, sort)
			
		#	prediction.cv.boosting3=unlist(lapply(my.folds.internal, .cart.for.cv.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]
		prediction.cv.boosting3=unlist(lapply(my.folds.internal, .cart.for.cv.prediction.internal.with.arguments,  N, my.ranks, num.var, my.tr.s, my.y.s, my.w.s, my.y.0, my.data.train,  my.weights, y, min.samples, max.depth, min.gain))[order(unlist(my.folds.internal))]

	
							
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
	
	#10/6/2013: added for debugging
	#set.seed(1)
	
	num.ties=sum(diff.alpha==0, na.rm=T)
	if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
	
	
#calculates the sum of alpha for those predicted 1, molti
#	prediction.new.data.boosting=unlist(lapply(c(1:num.samples.test), function(i) sum(prediction.new[i,]*alpha, na.rm=T)))
#prediction.new.data.boosting=ifelse(prediction.new.data.boosting>all.alpha-prediction.new.data.boosting, 1, 0)





########### if the class membership of the test set is provided, calculate the error rate for each number of boosting iterations
############ removed, it is not used in the output
#if(!is.null(y.test)) {

# 	cumsum.alpha.half=cumsum(alpha)/2
	
#	if(m==1) error.boosting=mean(y.test!=prediction.new[,1]) else{
	
	#######boosting error on the test set, one value for each possible value of boosting iterations, for example: the second value would be the average error obtained using 2 boosting iterations
#	error.boosting=c(#first step 
#			mean(y.test!=prediction.new[,1]),
			#next steps
#			unlist(lapply(2:num.boost, function(mm) {
#
#								diff.alpha=rowSums(mat.prediction.new.y1.times.alpha[,1:mm])-cumsum.alpha.half[mm]
#								prediction.new.data.boosting=ifelse(diff.alpha>0, 1, 0)
#
#								#are there any ties?  if there are ties assign the class at random
#								num.ties=sum(diff.alpha==0, na.rm=T)
#								if(num.ties>0)  {prediction.new.data.boosting[diff.alpha==0]=sample(c(0,1), num.ties, replace=T)}
#								
#								mean(prediction.new.data.boosting!=y.test)
#
#
#							}#end function(m) 
#				)#end lapply
#			)#end unlist
#		)#end error.boosting
#}#end else m==1
#} else {error.boosting<-NA}


########### added for EE, calculates the score, defined as diff.alpha, a score is given for each sample 
return(score=diff.alpha)


#return(list(results.cart=results.cart, prediction=prediction.cv, my.folds=my.folds))
##return(list(results.cart=results.cart.original, prediction=prediction.cv, ##prediction.new.boosting=prediction.new.data.boosting, alpha=alpha, 
#prediction.boosting=prediction.new, 
#m: number of performed boosting iterations
##weights.after.boosting=my.weights, ##error.boosting=error.boosting,prob=diff.alpha,pred.new.i=prediction.new, m=m))







}

	
	
	
	
	
	
	
	
	
	
	
	
	