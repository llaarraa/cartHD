######### function that fits cart	
############# Roxygen prelude ########################
#' CART for binary outcomes
#' 
#' Fit of a classification and regression tree (CART) for binary outcomes
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @param max.depth maximum number of levels of the tree
#' @param check.data logical indicator, if set to TRUE the function checks the formal correctness of the objects passed as arguments to the function. Set to FALSE to save computational time - for example when running simulations.
#' @return results.cart an object of class CART (contains the fitted CART model)
#' @note \code{cart} doesn't work if there are missing values 
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @seealso \code{\link{print.CART}}, \code{\link{plot.CART}}, \code{\link{predict.CART}}
#' @export
#' @examples 
#' set.seed(1)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000) #simulate data from 10 samples with 1000 variables
#' y=rbinom(10, size=1, prob=0.5) #simulate class membership of the 10 samples (0/1 vector)
#' my.cart.fit=cart(my.data.train, y) #fit the CART model
#' print(my.cart.fit) # print the fitted model



cart <-
function(my.data.train, y, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01,  max.depth=5, check.data=TRUE){

############ function that derives cart


####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train)
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
####max.depth: maximum number of levels of the tree


############## uses internal functions:
#.my.stump.forCART.oneStep.stepOne
#.my.stump.forCART.oneStep.stepNotOne


####################### step 00: check data #############################
			
			#run checks on data only if the functions are not used for simulations		
			if(check.data){
				.check.training.data(my.data.train, y, my.weights, my.prior, max.depth)
				}

###################### step 0: initialization/evaluation of objects needed later ########################

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


				######################## rescaled weights ########################
				############ not needed at the moment, useful for adjusting the results for non-empirical priors
				#### defining rescaled weights, the sum of the rescaled weights 
				#my.weights.prob=rep(0, N)
				#my.weights.prob[y==1]=my.weights[y==1]*my.prior[1]/sum(my.weights[y==1])
				#my.weights.prob[y==0]=my.weights[y==0]*my.prior[2]/sum(my.weights[y==0])
				#my.weights.prob=my.weights.prob/sum(my.weights.prob)
				 
				##### other quantities 
				#my.w.p.s<-matrix(my.weights.prob[my.ind.s], ncol=num.var)
				#my.w.p.times.y.0=my.w.p.s*my.y.0
				#my.w.p.times.y.s=my.w.p.s*my.y.s

				#n.l.p<-apply(my.w.p.s,2,cumsum)
				#n.r.p<-1-n.l.p

				############## to include the prior in the calculation of the gini gain
				#N1.w.p=sum(my.weights.prob[y==1])
				#N2.w.p=1-N1.w.p

				######################## end rescaled weights ########################

				##########################################################



###############################step 1: fit of cart ##################################

#my.stump.mat1.Lara.duplicates.faster.internal<-function(my.tr.s, my.y.s, my.w.s, my.w.times.y.0, my.w.times.y.s,  num.var,  
#		N, n.l, n.r, my.weights, my.index.keep=NULL, 
#		QT.max.old=NULL,  min.samples=2, min.gain=0.01, which.variables.used=NULL, depth=1, parent=0)


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
							#if(is.null(results.cart[[k]]$var))  return(results.cart) #return the list with the tree results

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
																											parent.history=c(results.cart[[which.k.parent]]$parent.history, k), giniP=NULL, p1Node.l.new=NULL, p1Node.r.new=NULL,
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

			
		names(results.cart)=unlist(lapply(1:length(results.cart), function(k) .convert.node(k, results.cart[[k]]$depth )))

		#for the output: add the vector with the class membership of the training set and define the class of the object (class: cart)
		Results.cart=list(tree=results.cart, y=y)
			
			
			#define the class CART for the outout object	
			class(Results.cart)="CART"							
									
									
#return the list with the tree results
return(Results.cart)

}
