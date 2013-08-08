#' CART prediction of class membership for new samples, produces also a classification score for each sample
#' 
#' The function predicts the class membreship of new samples and returns the class membership and a numerical score expressing the probability of the sample of being in Class 1. 	
#' @param results.cart CART object fitted with \code{cart} or \code{stump} functions 
#' @param new.data data set containing the new observations: rows: samples, cols= variables, must have the same variables as my.data.train
#' @return prediction 
#' \item{my.prediction}{predicted class membership obtained for the new sample}
#' \item{my.score}{estimated probability of being in Class 1 for new sample}
#' @note \code{my.predict.cart.scpres} doesn't work if there are missing values 
#' @references Breiman, L., J. H. Friedman, and R. A. Olshen. "Stone., CJ (1984) Classification and Regression Trees." Chapman and Hall/CRC
#' @export
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000)
#' y=rbinom(10, size=1, prob=0.5)
#' my.cart.fit=cart(my.data.train, y)
#' my.data.test=matrix(rnorm(10*1000), ncol=1000)
#' y.test=rbinom(10, size=1, prob=0.5)
#' my.predict.cart=my.predict.cart.with.score(my.cart.fit$tree, my.data.test)


my.predict.cart.with.score <-
function(results.cart, new.data){
#function that predicts the class membership for new data, given a CART object	

			####### results.cart: list containing the cart output, obtained on training data
			######## new data: data set containing the new observations: rows: samples, cols= variables, must have the same variables as my.data.train

			######## modified: 13/11: the output is NA if the tree had only one split, and no variable was selected, can happen if no variable is "good enough" to make a first split

					######### check if the data are a vector instead of a matrix, transforms the vector in a matrix, if necessary
					if(is.null(dim(new.data))) new.data=t(as.matrix(new.data))


					############## supporting function, prediction for a single sample
					my.predict=function(x, results.cart){
							######### relies on the assumption that the stop nodes were correctly assigned by the cart funcion
							k=1
							my.stop.node=is.element(k, which.stop.nodes)

							#continue until reaching a stop node
							while(!my.stop.node){

							#select the next node, left or right depending on the value of the variables, left if less than threshold, right otherwise
							k=ifelse(x[results.cart[[k]]$var]<results.cart[[k]]$thr, results.cart[[k]]$k.son.l, results.cart[[k]]$k.son.r)
							my.stop.node=is.element(k, which.stop.nodes)

							}

							#obtain the class prediction, defined as the assigned class membership assigned to the left or right node
							which.direction=ifelse(x[results.cart[[k]]$var]<results.cart[[k]]$thr, 1, 0)
						
							#my.prediction=ifelse(x[results.cart[[k]]$var]<results.cart[[k]]$thr, results.cart[[k]]$classNode.l, results.cart[[k]]$classNode.r)
							my.prediction=ifelse(which.direction==1, results.cart[[k]]$classNode.l, results.cart[[k]]$classNode.r)
							#my.score=ifelse(x[results.cart[[k]]$var]<results.cart[[k]]$thr, results.cart[[k]]$p1Node.l.new, results.cart[[k]]$p1Node.r.new)
							
							#my.p=ifelse(results.cart[[k]]$var==1 & results.cart[[k]]$thr==99999, 	ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), NA)
							#my.score=ifelse(my.prediction==1, ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), 1-ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ))

							#score: probability of being in class 1
							which.parent=ifelse(k>1, results.cart[[k]]$parent, 1)
							#my.score=ifelse(results.cart[[k]]$parent.direction==0, round(results.cart[[which.parent]]$p1Node.r.new,2), round(results.cart[[which.parent]]$p1Node.l.new,2) )
							#my.score=ifelse(which.direction==0, round(results.cart[[which.parent]]$p1Node.r.new,2), round(results.cart[[which.parent]]$p1Node.l.new,2) )
							#6/5/2013: fixed bug, fixed the problem that appeared when the evaluated node was of no further split (var=1, thr=99999)
							my.score=ifelse(which.direction==0, results.cart[[k]]$p1Node.r.new, results.cart[[k]]$p1Node.l.new )

							
						return(c(my.prediction, my.score))	
					} #end of my.predict, internal function


					
					#obtain the indexes of the stop nodes  
					which.stop.nodes=which(unlist(lapply(results.cart, function(x) x$my.stop==TRUE))==TRUE)
					#which.not.stop.nodes=which(unlist(lapply(results.cart, function(x) x$my.stop==FALSE))==TRUE)


					################# check if there was at least one step that was performed, can be a problem only if the tree stopped at the first iteration
					if(length(which.stop.nodes)<=1 & is.null(results.cart[[1]]$var)) return(rep(NA, nrow(new.data)))
					

					#unlist(apply(new.data, 1, my.predict))
					
					#evaluate for each new sample - rows of the new.data matrix - the class membership and the scores (estimated probability of belonging to class 1)
					as.numeric(apply(new.data, 1, my.predict, results.cart))




}
