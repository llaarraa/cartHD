#' Prediction method for CART objects
#' 
#' The function predicts the class membership of new samples using a fitted CART model (obtained with \code{cart} or \code{stump} function). It also computes some accuracy measures if the class membership of new data is provided.
#' @param results.cart CART object fitted with \code{cart} or \code{stump} functions 
#' @param new.data data set containing the new observations: rows: samples, cols= variables, must have the same variables as my.data.train
#' @param y.new class membership of the test set (must be a numeric vector were the classes are coded with 0 and 1, using the same coding used for the training set class membership). 
#' @return Prediction obtained using the CART fit and new data
#' \item{prediction.test}{Predicted class membership for the new samples}
#' \item{class.test}{True class membership of the new samples}
#' \item{score.test.class1}{Estimated probability of belonging to class 1 for each sample}
#' \item{accuracy.all}{Overall predictive accuracy}
#' \item{accuracy.class1}{Predictive accuracy for class 1}
#' \item{accuracy.class0}{Predictive accuracy for class 0}
#' @method predict CART
#' @note \code{predict.CART} doesn't work if there are missing values 
#' @references Breiman, L., J. H. Friedman, and R. A. Olshen. "Stone., CJ (1984) Classification and Regression Trees." Chapman and Hall/CRC
#' @seealso \code{\link{cart}}, \code{\link{stump}} 
#' @note The accuracy measures are computed only if the class membership of the new samples is passed as a parameter (y.new). The scores calculated for each new sample are the estimated probabilities of belonging to class 1.
#' @export
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(10*1000), ncol=1000)
#' y=rbinom(10, size=1, prob=0.5)
#' my.cart.fit=cart(my.data.train, y)
#' my.data.test=matrix(rnorm(10*1000), ncol=1000)
#' y.test=rbinom(10, size=1, prob=0.5)
#' my.predictions=predict(my.cart.fit, my.data.test, y.test)
#' @export 


predict.CART <-
function(results.cart, new.data, y.new=NULL){
#function that predicts the class membership for new data, given a CART object	


			####### results.cart: list containing the cart output, obtained on training data, class CART
			#use the tree
			results.cart=results.cart$tree

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

							#calculate the score: probability of being in class 1
							which.parent=ifelse(k>1, results.cart[[k]]$parent, 1)
							#my.score=ifelse(results.cart[[k]]$parent.direction==0, round(results.cart[[which.parent]]$p1Node.r.new,2), round(results.cart[[which.parent]]$p1Node.l.new,2) )
							#my.score=ifelse(which.direction==0, round(results.cart[[which.parent]]$p1Node.r.new,2), round(results.cart[[which.parent]]$p1Node.l.new,2) )
							
							
							#score: probability of being in class 1, 
							#6/5/2013: fixed bug, added the score (equal for both splits) also for the "no further splits" nodes (var=1, thr=99999)
							my.score=ifelse(which.direction==0, results.cart[[k]]$p1Node.r.new, results.cart[[k]]$p1Node.l.new )
							

							#my.p=ifelse(my.cart[[k]]$my.stop==T & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, 	ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), NA)
							#it evaluates the score differently, depending on whether the node is a "no further split node" (var=1, thr=99999) or a node where a further split occurs
							#my.score=ifelse(my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, 	ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), ifelse(which.direction==0, round(results.cart[[which.parent]]$p1Node.r.new,2), round(results.cart[[which.parent]]$p1Node.l.new,2) ))

							
							
							
						return(c(my.prediction, my.score))	
					} #end of my.predict, internal function


					
					#obtain the indexes of the stop nodes  
					which.stop.nodes=which(unlist(lapply(results.cart, function(x) x$my.stop==TRUE))==TRUE)
					#which.not.stop.nodes=which(unlist(lapply(results.cart, function(x) x$my.stop==FALSE))==TRUE)


					################# check if there was at least one step that was performed, can be a problem only if the tree stopped at the first iteration
					if(length(which.stop.nodes)<=1 & is.null(results.cart[[1]]$var)) return(rep(NA, nrow(new.data)))
					

					#unlist(apply(new.data, 1, my.predict))
					
					#evaluate for each new sample - rows of the new.data matrix - the class membership (first row of the matrix output) and the associated score of being in Class 1 (second row of the matrix output)
					pred.results=apply(new.data, 1, my.predict, results.cart)
					#saving the dimnames
					dimnames(pred.results)[[2]]=dimnames(new.data)[[1]]
													

	return(list(prediction.test=pred.results[1,], class.test=ifelse(!is.null(y.new), y.new, "Test set class membership not provided"), score.test.class1=pred.results[2,], accuracy.all=ifelse(!is.null(y.new), mean(y.new==pred.results[1,]), "Test set class membership not provided"), accuracy.class1=ifelse(!is.null(y.new), mean(y.new[y.new==1]==pred.results[1,y.new==1]), "Test set class membership not provided"), accuracy.class0=ifelse(!is.null(y.new), mean(y.new[y.new==0]==pred.results[1,y.new==0]), "Test set class membership not provided") ))

}
