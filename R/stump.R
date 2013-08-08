#' Stumps (CART with a single node) for binary outcomes
#' 
#' Function that fits stumps for binary outcomes. Stumps are CART objects with a single split.
#' @param my.data.train training set data; samples by rows, variables by columns
#' @param y vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
#' @param my.weights weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
#' @param my.prior vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
#' @param min.samples minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
#' @param min.gain minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)
#' @return Results.cart CART object containing the fitted stump (CART with one node)
#' @note \code{stump} doesn't work if there are missing values 
#' @references L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone  (1984) "Classification and Regression Trees." Chapman and Hall/CRC
#' @seealso \code{cart}, \code{print.CART}, \code{predict.CART}, \code{plot.CART}
#' @export
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000) #simulated data for 10 samples and 1000 variables
#' y=rbinom(10, size=1, prob=0.5) #simulated class membership 
#' my.stump.fit=stump(my.data.train, y) #stump fit
#' print(my.stump.fit) #display of the result


stump <-
function(my.data.train, y, my.weights=rep(1/nrow(my.data.train), nrow(my.data.train)), my.prior=NULL, min.samples=3, min.gain=0.01){

######### stump: function that fits stumps, a cart with only one node

####my.data.train: training set data; samples by rows, variables by columns
####y: vector with class membership: must be 1 or 0 and its length must be equal to the number of samples (number of rows of my.data.train
####my.weights: weights for each observation, vector with length equal to the number of observations; the default is equal weights to all the observations
####my.prior: vector containing the prior probability for class 1 and class 0; if NULL it is set equal to the empirical prior (weighted frequencies of the classes)
####min.samples: minimum number of samples in a node; if the number of samples is lower or equal the algorithm stops splitting
####min.gain: minimum improvement in the Gini index (calculated for FUTURE observations, the calculation includes the prior information, similarly as in rpart)


############## uses internal functions:
#.my.stump.forCART.oneStep.stepOne


############################ step 0: initialize the variables/save objects used for further calculations ###########################


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
						#8/4/2013: updated to use functin rank, colRanks only possible treatment of 
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



#################################	step 1: evaluate stump ################################

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
					
				
			
		#for the output: add the vector with the class membership of the training set and define the class of the object (class: cart)
		Results.cart=list(tree=results.cart, y=y)
			
			
			#define the class CART for the outout object	
			class(Results.cart)="CART"							

			
				
			#return the list with the tree results
			return(Results.cart)

			}
