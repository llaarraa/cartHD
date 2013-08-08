#' Plot method for CART objects
#'
#' The function plots the fitted CART object. The plot can be displayed in the internal graphical device or save in a PDF file.
#' @param my.cart CART object of class CART obtained with \code{cart} function
#' @param file.out PDF file where the plot will be saved, if NULL the plot will be displayed in the internal graphical device. 
#' @return plots the CART fit
#' @method plot CART
#' @note Large trees are difficult to visualize. Display of large trees in the graphical device might produce errors or unclear graphs. Try saving the output in a PDF file, which will be automatically resized according to the depth of the tree. Specify the name of the file using the file.out parameter, the name of the file must include the pdf subfix (for example file.out="CartResult.pdf"); the plot is saved in the working directory, unless the (valid) full path is specified.
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(100*1000), ncol=1000)
#' y=rbinom(100, size=1, prob=0.5)
#' my.cart.fit=cart(my.data.train, y)
#' plot(my.cart.fit, file.out="cartPlot.pdf")
#' @export 

plot.CART=function(my.cart, file.out=NULL){

			################## obtain the output of the CART object ############
			#my.output=print.CART(my.cart, reorder.levels=FALSE) #### 22/4/2013 removed the option on reordering of the levels as it was unnecessary.
			my.output=print.CART(my.cart)
			
									
			
			############## set the plotting region ######################
			
			#largest index of the nodes
			max.k.node=max(my.output$NodeID)
			
			#number of nodes used
			num.nodes.used=dim(my.output)[1]		
		
			#maximum depth of the displayed tree, at most it is equal to the chosen max.depth + 1
			max.depth=max(my.output$Level)
			
			#maximum number of POSSIBLE nodes in the deepest level
			num.nodes.max.depth=2^(max.depth-1)

			#indexes of all the possible nodes to display - add the extra line beyond max.depth		
			kk.all= 1:c(2^max.depth-1)
			
			#level of each possible node (evaluates the depth for the indexes from 1 to num.nodes.max.depth)
			my.depth = rep(1:max.depth, 2^(c(1:max.depth)-1))
	

			
			my.indexes.layout=rep(kk.all, 2^(max.depth-my.depth))
			#my.indexes.layout[!is.element(my.indexes.layout, my.output$NodeID)]=0
			
			My.indexes.layout=numeric(length(my.indexes.layout))
			
			#substitute the original indexes of the nodes with the row number in which the node appears in the output matrix		
			for(i in 1:num.nodes.used) {My.indexes.layout[my.indexes.layout==my.output$NodeID[i]] <-  i}
	
	
			#convert into a matrix that is used to define the 
			mat.layout = matrix(My.indexes.layout, ncol=num.nodes.max.depth, byrow=T)
			
			
			
				

			################### plot ##########################
			

			######## open the file and use the correct size
			
	
			if(!is.null(file.out))  pdf(file.out, width=3*num.nodes.max.depth, height=max.depth*6) 
			#pdf(file.out, width=30, height=40) 
			layout(mat.layout)
			par(cex=1)
			
			par(mar=c(0,0,0,0), oma=c(0,0,0,0))

			
			for(k in 1:num.nodes.used) {
					
							#frame()
							par(cex=1-my.output$Level[k]/10+.1)
							plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
		
		
						
							my.labels=paste(#my.var, my.thr,
									"n=", my.output$n[k], " \n(n1=", my.output$n1[k], ", n0=", my.output$n[k]-my.output$n1[k], ")\n", "Gini=", my.output$Gini[k],
									ifelse(my.output$STOP[k]=="**", paste("\n Class=", my.output$Class[k], "; Prob=", my.output$Prob[k]), "") 
									, sep="" )#end paste labels								
								
								xh=1.5*strwidth(my.labels, units = "figure", cex = 2, font = NULL, vfont = NULL)
								yh=1.5*strheight(my.labels, units = "figure", cex = 2, font = NULL, vfont = NULL)
								
								#rect(0.5 - xh/2, 0.5 - yh/2, 0.5 + xh/2, 0.5 + yh/2, col = "white")
								#rect(0.5 - xh/2 - .1 , 1 - max(0, yh+.2), 0.5 + xh/2 + .1, 1 , col = "white")
								col.box=ifelse(my.output$STOP[k]!="**", "white", "red")
								if( my.output$STOP[k]=="**" & my.output$Class[k]==0)  col.box = "orange"
								
								rect(0.5 - xh/2 , 1 - yh, 0.5 + xh/2 , 1 , col = col.box)
								
								
								# information about the node, on y-scale it is alligned at the top - almost 1-0.5 , centered on the x-axes (adj=c(0.5, 1-0.05))								
								text(0.5, .95, label=my.labels, cex=2, col=1, adj=c(0.5,1-.05))
								#end text
								
								#draw arrows if it is not a final node
								if(my.output$STOP[k]!="**") {
										#segments(0.5, 0.5, 0.25, 0)	
										#segments(0.5, 0.5, .75, 0)	
										#text(0.25, 0.25, labels=paste("Var ", my.output$Var[k], "<", (my.output$Threshold[k])), cex=par("cex")*3)
										#text(0.75, 0.25, labels=paste("Var ", my.output$Var[k], ">=", (my.output$Threshold[k])), cex=par("cex")*3)
									
										segments(0.5, 1-yh-.1, 0.25,  0.05)	
										segments(0.5, 1-yh-.1, 0.75, 0.05)	
										#text(0.25, 1-yh-.2, labels=paste("Var ", my.output$Var[k], "<", (my.output$Threshold[k])), cex=par("cex")*3)
										#text(0.75, 1-yh-.3, labels=paste("Var ", my.output$Var[k], ">=", (my.output$Threshold[k])), cex=par("cex")*3)
										
										#right alignment on the left side, as the levels go deeper the text is more centered (further away from 1)
										text(.25, 0, labels=paste("Var ", my.output$Var[k], "<", (my.output$Threshold[k])), cex=par("cex")*3, adj=1-my.output$Level[k]/7+1/7)
										#left alignment on the left side, as the levels go deeper the text is more centered (further away from 0)
										text(.75, 0, labels=paste("Var ", my.output$Var[k], ">=", (my.output$Threshold[k])), cex=par("cex")*3, adj=0+my.output$Level[k]/7-1/7)
										
									#abline(h=0, lty=2)
	
								}#end if(my.cart[[k]]$my.stop!=TRUE)

					
				}#end for k
	
	
	
	
	########### close the file
	if(!is.null(file.out)) { tmp=dev.off() 
		cat("\n\n The results were plotted in the ", file.out, " file")
		}
	
	
}# end plot.CART


