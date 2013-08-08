#' Displays graphically the CART fit
#' @param my.cart CART object of class cart obtained with \code{cart} function
#' @return plots the CART fit
#' @export
#' @examples
#' set.seed(1)
#' my.data.train=matrix(rnorm(100*1000), ncol=1000)
#' y=rbinom(100, size=1, prob=0.5)
#' my.cart.fit=cart(my.data.train, y)
#' pdf("cartPlot.pdf", width=20, height=30)
#' draw.cart(my.cart.fit)
#' tmp=dev.off()

draw.cart=function(my.cart){

		my.y=my.cart$y
		my.cart=my.cart$tree
		
			
			num.nodes=length(my.cart)
			max.depth=ceiling(log2(num.nodes+.1))
			num.nodes.max.depth=2^(max.depth-1)
			my.depth=unlist(lapply(my.cart, function(x) x$depth))

		
			kk=1:num.nodes

			mat.layout = matrix(rep(kk, 2^(max.depth-my.depth)), ncol=num.nodes.max.depth, byrow=T)
			layout(mat.layout)
			 
			for(k in 1:num.nodes) {
			#par(mar=c(1,1,1,1))
			frame()
			if(!is.null(my.cart[[k]]$var)){
			#sample size within the node 
			n=ifelse(k==1, length(my.y), 
			#parent.direction: 1 for left, 0 for right
			ifelse(my.cart[[k]]$parent.direction==0, length(my.cart[[my.cart[[k]]$parent]]$which.n.r), length(my.cart[[my.cart[[k]]$parent]]$which.n.l) )  )

			#how many samples from class 1
			n1=ifelse(k==1, sum(my.y==1), 
			#parent.direction: 1 for left, 0 for right
			ifelse(my.cart[[k]]$parent.direction==0, sum(my.y[my.cart[[my.cart[[k]]$parent]]$which.n.r]==1), sum(my.y[my.cart[[my.cart[[k]]$parent]]$which.n.l]==1) )  )

			#gini
			gini=round(ifelse(k==1, my.cart[[k]]$giniP, 
			#parent.direction: 1 for left, 0 for right
			ifelse(my.cart[[k]]$parent.direction==0, my.cart[[my.cart[[k]]$parent]]$gini.r,  my.cart[[my.cart[[k]]$parent]]$gini.l)),2)

			my.class=ifelse(my.cart[[k]]$my.stop==T, paste("Class=", my.cart[[k]]$classNode.0, " "), "")
			my.p=ifelse(my.cart[[k]]$my.stop==T, 
			paste("; P(Class=1)=",  ifelse(my.cart[[k]]$parent.direction==0, round(my.cart[[my.cart[[k]]$parent]]$p1Node.r.new,2), round(my.cart[[my.cart[[k]]$parent]]$p1Node.l.new,2) )),  
			"")

			my.thr=ifelse(my.cart[[k]]$my.stop==T, "", paste("Thr=", round(my.cart[[k]]$thr,2)))
			my.var=ifelse(my.cart[[k]]$my.stop==T, "", paste("Var=", my.cart[[k]]$var))


			mtext(text=paste(#my.var, my.thr,
			"n=", n, "(n1=", n1, ", n0=", n-n1, ")\n", "Gini=", gini,
			my.class, my.p
			))
			
		#draw arrows if it is not a final node
		if(my.cart[[k]]$my.stop!=TRUE) 
			{segments(0.5, 0.5, 0.25, 0)	
			segments(0.5, 0.5, .75, 0)	
			text(0.25, 0.25, labels=paste(my.var, "<", round(my.cart[[k]]$thr,2)), cex=1.5)
			text(0.75, 0.25, labels=paste(my.var, ">", round(my.cart[[k]]$thr,2)), cex=1.5)

	}#end if


}
}
 
}# end draw.cart

