#' Print method for CART object
#'
#' The function prints a CART object. 
#' @param my.cart CART object of class CART obtained with \code{cart} function
#' @return prints the CART fit
#' @method print CART
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(100*1000), ncol=1000)
#' y=rbinom(100, size=1, prob=0.5)
#' my.cart.fit=cart(my.data.train, y)
#' print(my.cart.fit)
#' @export

print.CART=function(my.cart, ...){
#print.CART=function(my.cart, reorder.levels=TRUE, ...){

			
			my.y=my.cart$y
			my.cart=my.cart$tree
			
			
			#number of nodes included in the cart tree, some nodes can be empty
			num.nodes=length(my.cart)
			#maximum depth of the tree, at most it is equal to the set max.depth
			max.depth=ceiling(log2(num.nodes+.1))
			#number of nodes in the final level
			num.nodes.max.depth=2^(max.depth-1)
			#depth of each node of the tree
			my.depth=unlist(lapply(my.cart, function(x) x$depth))

		
			kk=1:num.nodes

			#mat.layout = matrix(rep(kk, 2^(max.depth-my.depth)), ncol=num.nodes.max.depth, byrow=T)
			#layout(mat.layout)
			my.output=NULL 
			for(k in 1:num.nodes) {
					#par(mar=c(1,1,1,1))
					
					#check if the node was used, only the nodes that were effectively evaluated are inspected
					if(!is.null(my.cart[[k]]$var)){
					#sample size within the node 
					n=ifelse(k==1, length(my.y), 
					#parent.direction: 1 for left, 0 for right
					#it is necessary to calculate it in this way beacuse the final nodes that require no further split do not have explicitly the number of samples included
					ifelse(my.cart[[k]]$parent.direction==0, length(my.cart[[my.cart[[k]]$parent]]$which.n.r), length(my.cart[[my.cart[[k]]$parent]]$which.n.l) )  )

					#how many samples from class 1
					n1=ifelse(k==1, sum(my.y==1), 
					#parent.direction: 1 for left, 0 for right
					ifelse(my.cart[[k]]$parent.direction==0, sum(my.y[my.cart[[my.cart[[k]]$parent]]$which.n.r]==1), sum(my.y[my.cart[[my.cart[[k]]$parent]]$which.n.l]==1) )  )

					#gini
					gini=round(ifelse(k==1, my.cart[[k]]$giniP, 
					#parent.direction: 1 for left, 0 for right
					ifelse(my.cart[[k]]$parent.direction==0, my.cart[[my.cart[[k]]$parent]]$gini.r,  my.cart[[my.cart[[k]]$parent]]$gini.l)),2)

					#my.class=ifelse(my.cart[[k]]$my.stop==T, my.cart[[k]]$classNode.0, NA)
					#my.p=ifelse(my.cart[[k]]$my.stop==T, 
					#ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), NA)

					
					
				#	my.thr=ifelse(my.cart[[k]]$my.stop==T, NA, round(my.cart[[k]]$thr,2))
				#	my.var=ifelse(my.cart[[k]]$my.stop==T, NA, my.cart[[k]]$var)
				
				#	my.thr=ifelse(my.cart[[k]]$my.stop==TRUE | ( my.cart[[k]]$my.stop==TRUE & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999), NA, round(my.cart[[k]]$thr,2))
				#	my.var=ifelse(my.cart[[k]]$my.stop==TRUE | ( my.cart[[k]]$my.stop==TRUE & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999), NA, my.cart[[k]]$var)
				
					my.thr=ifelse( my.cart[[k]]$my.stop==TRUE & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, NA, round(my.cart[[k]]$thr,2))
					my.var=ifelse( my.cart[[k]]$my.stop==TRUE & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, NA, my.cart[[k]]$var)

					my.class=ifelse(my.cart[[k]]$my.stop==T & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, my.cart[[k]]$classNode.0, NA)
					my.p=ifelse(my.cart[[k]]$my.stop==T & my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999, 	ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) ), NA)

					
					
				
				#my.stop=ifelse(my.cart[[k]]$my.stop, "**", "")
					#a full stop is given if my.stop=T and (var!=1 & thr!=99999) - this coding was used to indicate the nodes for which a complete stop is required, without a further split. 
					
					my.stop=ifelse(my.cart[[k]]$my.stop &  (my.cart[[k]]$var==1 & my.cart[[k]]$thr==99999), "**", NA)
					
					#my.output=rbind.data.frame(my.output, data.frame(Level=my.cart[[k]]$depth, Node=.convert.node(k, my.cart[[k]]$depth),  n=n, n1=n1, n0=n-n1, Gini=gini, STOP=my.stop, "Assigned Class"=my.class, "P(Class1)"=my.p, Var=my.var, Threshold=my.thr,  Path=paste(.convert.node.depth.unknown(my.cart[[k]]$parent.history), collapse="-")))
					my.output=rbind.data.frame(my.output, data.frame(Level=my.cart[[k]]$depth, Node=.convert.node.depth.unknown.code2(k),  NodeID=k, n=n, n1=n1, n0=n-n1, Gini=gini, STOP=my.stop, "Class"=my.class, "Prob"=my.p, Var=my.var, Threshold=my.thr  #, Path=paste(.convert.node.depth.unknown(my.cart[[k]]$parent.history), 	collapse="-"))
					))
					
					
					#add an addtional final node if the current node is a stop node but requires a further split
					if(my.cart[[k]]$my.stop==TRUE & (my.cart[[k]]$var!=1 | my.cart[[k]]$thr!=99999)) {

							my.stop="**"
							my.thr=NA
							my.var=NA
							
							#left node 
							
								#sample size within the node 
								#parent.direction: 1 for left, 0 for right
								n=length(my.cart[[k]]$which.n.l) 

								#how many samples from class 1
								n1=sum(my.y[my.cart[[k]]$which.n.l]==1)   

								#gini
								gini=round(my.cart[[k]]$gini.l,2)

								my.class=my.cart[[k]]$classNode.l
								my.p=round(my.cart[[k]]$p1Node.l.new,2) 
						
							#derive the index for the left son node
							k.node.l=2*(k-2^(my.cart[[k]]$depth-1))+2^(my.cart[[k]]$depth)
							node.l=.convert.node.depth.unknown.code2(k.node.l)

								
								my.output=rbind.data.frame(my.output, data.frame(Level=my.cart[[k]]$depth+1, Node=node.l, NodeID=k.node.l, n=n, n1=n1, n0=n-n1, Gini=gini, STOP=my.stop, "Class"=my.class, "Prob"=my.p, Var=my.var, Threshold=my.thr  
								#Path=paste(.convert.node.depth.unknown(c(my.cart[[k]]$parent.history, node.l)), 	collapse="-")
								))
							
							#right node 
							
								#sample size within the node 
								#parent.direction: 1 for left, 0 for right
								n=length(my.cart[[k]]$which.n.r) 

								#how many samples from class 1
								n1=sum(my.y[my.cart[[k]]$which.n.r]==1)   

								#gini
								gini=round(my.cart[[k]]$gini.r,2)

								my.class=my.cart[[k]]$classNode.r
								my.p=round(my.cart[[k]]$p1Node.r.new,2) 

							#node.r=.convert.node.depth.unknown.code2(2*(k-2^(my.cart[[k]]$depth-1))+2^(my.cart[[k]]$depth)+1)
							node.r=.convert.node.depth.unknown.code2(k.node.l+1)
							
								my.output=rbind.data.frame(my.output, data.frame(Level=my.cart[[k]]$depth+1, Node=node.r, NodeID=k.node.l+1, n=n, n1=n1, n0=n-n1, Gini=gini, STOP=my.stop, "Class"=my.class, "Prob"=my.p, Var=my.var, Threshold=my.thr  
								#Path=paste(.convert.node.depth.unknown(c(my.cart[[k]]$parent.history, node.l+1)), collapse="-")
								))
							
							
							
							
					}# end if() a further split is required before the final stop
					
					
					
					
					

			}#end if(!is.null(my.cart[[k]]$var))
			
		}#end for k
		
		
		
	#calculate the probability of the assigned class
	my.output$Prob=ifelse(my.output$Class==1, my.output$Prob, 1-my.output$Prob)
	
	#remove the missing values from the output		
	my.output[is.na(my.output)]	= ""
	
	#re-order the output by the level of the node in the tree
	#if(reorder.levels) my.output=my.output[order(my.output$NodeID),]
	
	########## add an empty line at the beginning of each level
	
	#newrow=rep("", dim(my.output)[2])
	#level.prev=my.output$Level[1]
	#not.end=TRUE
	#i=2
	#while(not.end) {
	#level.next=my.output$Level[i]
	#	if(level.prev!=level.next) {my.output=rbind.data.frame(my.output[1:(i-1),], newrow, my.output[-c(1:(i-1)),])
	#								i=i+2}
	#	lev.prev=level.next
	#	if(i==dim(my.output)[2]) not.end=FALSE
	#}#end for
	
	#print(my.output[my.output$Level==1,], row.names=FALSE, right=FALSE)
	#for(i in unique(my.output$Level[-1])) {
	#	cat("\n")
	#print(my.output[my.output$Level==i,], row.names=FALSE, right=FALSE, col.names=FALSE)
	#}
	
	#		cat("\n\n")
print(my.output, row.names=FALSE, right=FALSE) 
}# end draw.cart


#print(my.cart.fit)
