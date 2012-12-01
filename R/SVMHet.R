#This file contains the functions used for analyzing causal heterogeneity using support vector machines under 
#a double LASSO constraint.

#The functions consist of 
	#scale.func: demeans and standardizes variables.
	#scale.func.2: gives variables a standard deviation of one
	#FindIt: a wrapper function for SVM.func
	#maketwoway<-rescales covariates and makes all twoway interactions.
	#SVM.func: the workhorse function.  Given values of y, control, treatment, and lambdas, 
		#it returns betas, the loss (GCV), fitted values, and size of the margin.
	#find.lambda: Conducts the alternating line search.  A large value of lambda.c is selected (10), and alternating
		#line searches of decreasing radius are conducted.
	#lars: An adaptation of the LARS algorithm of Efron, et al. 2004.  Due to the large numbers being considered in the fitting,
		#the code is adapted to not crash.  The original code is commented out.  The remainder of the code is untouched.
	
##################################################################
##################################################################
##################################################################

scale.func<-function(x) {
	x2<-x-mean(x)
		x2/sd(x2)}

scale.func.2<-function(x) {x2<-x
		x2*sd(x2)}

#####################################################

FindIt<-function(y,X.c, treat ,type="single", rescale.c=FALSE, search.lambdas=TRUE,lambdas=NULL,wts=1,scale.c=1,scale.t=1,n.highlow=10){
	
	if(type%in%c("single","multiple")==F){
		print("Type must be either single or multiple")
		break
	}


	if(sum(wts<0)>0) {
		print("Only non-negative weights allowed")
		break
	}
	
	#Give wts mean 1.
	wts<-wts/mean(wts)

	y<-y*wts^.5
	X.c<-X.c*wts^.5
	treat<-treat*wts^.5
	
	#Remove all elements with no variance and add an intercept.
	X.c<-X.c[,apply(X.c,2,sd)>0]
	X.c<-unique(X.c,MARGIN=2)
	if(sd(X.c[,1])>0){
		X.c<-cbind(1,X.c)
		colnames(X.c)[1]<-"Intercept"
	}
	
	
	if(rescale.c==T) {
		scale.c<-c(1,apply(X.c[,-1],2,sd))
		X.c[,-1]<-apply(X.c[,-1],2,scale.func)

	}
	
	if(is.vector(treat)) {
		X.t<-cbind(treat*1e4,(treat>0)*X.c[,-1])
		X.t[treat!=0,-1]<-apply(X.t[treat!=0,-1],2, FUN=function(x) x-mean(x))
		colnames(X.t)<-c("treat",paste("treat",colnames(X.c)[-1],sep=":"))
	} else {
		X.t<-treat
	}
	
	
	if(search.lambdas==T) lambdas<-search.lambda(y,X.c,X.t)
	A<-SVM.func(y,X.c[,-1],X.t,lambdas[1],lambdas[2])


	#Find ten highest and lowest	
	if(type=="single"){
		#Generate outcome treatment and no treatment.
		preds.treat<-cbind(X.c,cbind(1e4,X.c[,-1]))%*%A$beta
		preds.treat<-sign(preds.treat)*pmin(abs(preds.treat),1)
		preds.control<-cbind(X.c,0*cbind(1e4,X.c[,-1]))%*%A$beta	
		preds.control<-sign(preds.control)*pmin(abs(preds.control),1)
		preds.diff<-preds.treat-preds.control

		#Gather highest and lowest treated values.
		highlow.ind<-c(sort(preds.diff,decreasing=T,ind=T)$ix[1:n.highlow],rev(sort(preds.diff,decreasing=F,ind=T)$ix[1:10] ))
		highlow.diff<-c(sort(preds.diff,decreasing=T)[1:n.highlow],rev(sort(preds.diff,decreasing=F)[1:10] ))
		scale.out<-c(1,scale.c,1e-4,scale.c)
		ATE<-mean(preds.treat-preds.control)/2
	} else {
		highlow.ind<-highlow.diff<-NULL
		scale.out<-c(1,scale.c,rep(1,dim(X.t)[2]))

		#Make estimated effect table for every unique treatment.
		unique.t1<-unique.t1.orig<-unique(X.t>0,M=1)
		unique.t1<-unique.t1[(rowSums(unique.t1)>0),]
		names.all<-colnames(unique.t1)
		which.treat<-(1:length(names.all))
		
		X.treat.make<-matrix(0,nrow=length(which.treat),ncol=dim(X.t)[2])
		
		for(i in 1: length(which.treat)) {
			X.treat.make[i,]<-(X.t>0)[which((X.t>0)[,which.treat[i] ]==1)[1],]
		}
		
		preds.treat<-(X.treat.make)%*%A$beta[-c(1:(dim(X.c)[2]))]
		preds.control<-X.c%*%A$beta[1:(dim(X.c)[2])]
		
		
		make.post.prob<-function(x){
			preds<-preds.control
			probs.control<-pmin(pmax(preds,-1),1) 
			
			
			preds<-preds.control+x
			probs.treat<-pmin(pmax(preds,-1),1)
			(probs.treat-probs.control)/2
			
		}
		
		treat.out<-sapply(preds.treat,FUN=function(x) mean(make.post.prob(x)))
		
		out1<-data.frame(names.all,treat.out*100)
		ATE<-out1[sort(out1[,2],ind=T,decreasing=TRUE)$ix,]
		names(ATE)<-c("Treatment","Effect")
		
		
	}
	
	

	
	names.out<-c(colnames(X.c),colnames(X.t))

		output<-(list("coefs"=A$beta,"coefs.orig"=A$beta/scale.out,"fit"=A$fit,"highlow"=cbind(highlow.ind,highlow.diff),
					  "scale"=scale.out,"names.out"=names.out,"y"=y,"X.c"=X.c,"X.t"=X.t,"GCV"=exp(A$loss),"ATE"=ATE,"lambdas"=A$lambdas,"n.highlow"=n.highlow))
		class(output)<-c("FindIt","list")
		output
}

####################################################
print.FindIt<-summary.FindIt<-function(object,...){
  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)
  
	not.zero<-object$coefs.orig!=0
	coef.table<-cbind(object$name,object$coefs.orig/2)[not.zero,]
	coef.print<-cbind(object$name,signif(object$coefs.orig/2,3))[not.zero,]
	rownames(coef.print)<-rownames(coef.table)<-object$name[not.zero]	
  	coef.print[coef.print=="0"]<-"0.000"	
  	colnames(coef.print)<-c("Coefficient","Estimate")

	null.gcv<-var(object$y)/(length(object$y)-1)*length(object$y)
	model.gcv<-(object$GCV)*length(object$y)


	 cat("\nCoefficients:\n")

  print(noquote(coef.print[,-1]))  
  
    cat("\n---------")
  cat("\nModel Fit Statistics: \n")
  cat(c("GACV:\n"))
  	cat(c("  Null: ",round(null.gcv,3)))
  	cat(c("  Model: ",round(model.gcv,3),"\n"))
    cat(c("Percent Misclassified:\n"))
  	cat(c("   Null: ",round(min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2)))
  	cat(c("   Model: ",round(mean(sign(object$y)!=sign(object$fit)),2),"\n"))
  	cat(c("   Percent Improvement, vs. NULL: ",round(100-100*mean(sign(object$y)!=sign(object$fit))/min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2),"% \n"))
   	cat(c("Percent Outside Margin: ",round(mean((object$y^2- object$y*object$fit)<=0 )*100,3),"%, n =",sum((1- object$y*object$fit)<=0 ), "\n")) 	
  out<-list("coefficients"=data.frame(coef.table[,-1]),"GCV"=c(null.gcv,model.gcv),"misclass"=c(min(c(mean(sign(object$y)==-1),mean(object$y==1))),mean(sign(object$y)!=sign(object$fit))  ))
  invisible(out)
  
  
  
}




highlow<-function(object,X0=object$X.c,n.highlow=object$n.highlow, ...){
  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)
  
	X0<-X0[,apply(X0,2,sd)>0]
  	names.ord<-c("1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th",paste(11:500,"th",sep=""))[1:n.highlow]
  	which.notinter<-unique(c(grep(":",colnames(X0),fixed=T),grep(".2",colnames(X0),fixed=T)))
	data.highlow<-cbind(X0[object$highlow[,1],],object$highlow[,2]/2)
	if(length(which.notinter)>0) data.highlow<-cbind(X0[object$highlow[,1],],object$highlow[,2]/2)[,-which.notinter]
	colnames(data.highlow)[dim(data.highlow)[2]]<-"Estimate"
	data.high<-data.highlow[1:n.highlow,]
	data.high.print<-apply(data.high,2,FUN=function(x) if(mean(round(x)==x)==1  ) x else round(x,2)  )
	data.low<-data.highlow[(n.highlow+1):(2*n.highlow),]
	data.low.print<-apply(data.low,2,FUN=function(x) if(mean(round(x)==x)==1  ) x else round(x,2)  )

	rownames(data.high.print)<-rownames(data.low.print)<-names.ord

 	 cat("\nHighest Estimated Effect:\n")
 	 print(noquote(data.high.print))

 	 cat("\nLowest Estimated Effect:\n")
 	 print(noquote(data.low.print))
 	 
 	 out<-list("highlowtable"=data.highlow)
 	 invisible(out)
  
}



#####################################################
##makeallway

makeallway<-function(X,wts=1){
	X<-data.frame(X)
	X2<-NULL
	for(i in 1:dim(X)[2]) {
		X2[[i]]<-cbind(sapply(sort(unique(X[,i])),FUN=function(j) 1*(X[,i]==j)  )    )
		colnames(X2[[i]])<-c(paste(names(X[i]), sort(unique(X[,i])),sep="_"))
		}

	X2<-data.frame(X2)
	X.next<-as.matrix(X2)

	for(i in 1:dim(X)[2])	{
		X.next<-unique(maketwoway(X.next,center=FALSE)$X,MAR=2)
#X.next<-X.next[,apply(X.next,MAR=2,FUN=sd)>0]
		}
	X.next
}



#####################################################
##maketwoway


maketwoway<-function(X, wts=1,center=T){

	X.orig<-X<-X[,apply(X,2,sd)>0]
	rescale.maineff<-c(1,apply(X.orig,2,sd))
	if(center==T) X<-apply(X,2,scale.func)
	X.orig<-cbind(1,X.orig)
	X<-cbind(1,X)
	colnames(X)[1]<-"int"
	X.maineff<-X

	if(center==T) X.maineff[,-1]<-apply(X.maineff[,-1],MARGIN=2, FUN=function(x) 
	 (x-mean(wts^.5*x))/sd(x)
	)
	colnames.X<-X.intereff<-NULL
	rescale.intereff<-NULL
	valid.inter<-NULL


	#Generate interaction effects.
	for(i in 1:(dim(X.maineff)[2]-1)) {
	X.intereff<-cbind(X.maineff[,i]*X.maineff[,-c(1:i)],X.intereff)
		rescale.intereff<-c(rescale.maineff[i]*rescale.maineff[-c(1:i)],rescale.intereff)
	colnames.X<-c(paste(colnames(X.maineff)[i],colnames(X.maineff[,-c(1:i)]),sep=":"),colnames.X)
	valid.inter<-c(apply( as.matrix(X.orig[,i]*X.orig[,-c(1:i)]),2,FUN=function(x) length(unique(x))) ,valid.inter)
#if(center==FALSE) {
			keeps<-colSums(X.intereff^2)!=0
			X.intereff<-X.intereff[,keeps]
			colnames.X<-colnames.X[keeps]
			valid.inter<-valid.inter[keeps]
			rescale.intereff<-rescale.intereff[keeps]
#		}
			
			
		
	}
	
	#Gets lost in the loop, add manually.
	colnames.X[1]<-paste(rev(colnames(X.maineff))[2],rev(colnames(X.maineff))[1],sep=":")

	#Remove dummy-dummy interactions
	X.intereff<-X.intereff[,valid.inter>1]
	colnames.X<-colnames.X[valid.inter>1]


	X.sq.orig<-X.sq<-as.matrix(X.maineff^2)
		rescale.sq<-rescale.maineff^2
	colnames(X.sq)<-colnames(X.maineff)
	X.sq<-as.matrix(X.sq[,apply(X.sq,MARGIN=2,FUN=function(x) length(unique(x))>2)])
	rescale.sq<-rescale.sq[apply(X.sq.orig,MARGIN=2,FUN=function(x) length(unique(x))>2)]
	if(dim(X.sq)[2]==1) colnames(X.sq)[1]<-paste(names(rescale.sq),sep="")
	if(length(rescale.sq)==0) X.sq<-rescale.sq<-NULL
	X.c<-cbind(X.sq,X.intereff)
		rescale.c<-c(rescale.sq,rescale.intereff)
	if(length(colnames(X.sq))>0) colnames.X<-c(paste(colnames(X.sq),".2",sep=""),colnames.X)
	X.c<-X.c[,(dim(X.c)[2]):1]
	rescale.c<-rescale.c[(dim(X.c)[2]):1]
	colnames.X<-sub("int:","",colnames.X)
	colnames.X<-sub(":int","",colnames.X)
	colnames(X.c)<-rev(colnames.X)
	invisible(list("X"=X.c, "scale.X"=rescale.c ))
	
	
}

#####################################################
SVM.func<-function(y,X.c,X.t,lambda.c,lambda.t,wts=1){

	n<-length(y)
	X<-X2<-cbind(X.c,X.t)

	#Originally written to handle a string of lambda.c and lambda.t
		#loss.run creates a matrix of losses, across values of lambda.
	loss.run<-matrix(0,nrow=length(lambda.c),ncol=length(lambda.t))
	min.params<-c(lambda.c[1],lambda.t[1])
	for(i in 1:length(lambda.c)){
	for(j in 1:length(lambda.t)){

	#Scale the X matrix by lambdas.
		X3<-X2
		X3[,1:dim(X.c)[2]]<-1/exp(lambda.c[i])*X2[,1:dim(X.c)[2]]
		X3[,-c(1:dim(X.c)[2])]<-1/exp(lambda.t[j])*X2[,-c(1:dim(X.c)[2])]

	#Declare original values.
		X.new<-X3
		which.use<-rep(T,n)
		beta.new<-beta.curr<-rep(0,dim(X2)[2]+1)

	#The loop to fit the LASSO.  Center on the margin, fit the lasso, update coefficients,
		#then drop all outside the margin (y*yhat>1)
		for(beta.run in 1:100){
			X.new.2<-apply(X.new[which.use,],MARGIN=2,FUN=function(x) x-mean(x))
			X.new.2<-apply(X.new.2,MARGIN=2,FUN=function(x) if(length(unique(x))==1) rnorm(sum(which.use)) else x)
			glmnet.1<-glmnet(X.new.2,y[which.use]-mean(y[which.use]),family="gaussian",lambda=c(5,4,3,2,seq(1.5,1,-.1)),standardize=F)
			#This was the old L1 optimizer.
			#lasso1<-lars(X.new.2,(y[which.use]-mean(y[which.use])),max.steps=dim(X)[2]+2,normalize=F,type="lasso",eps=0)
			#beta.new[-1]<-as.vector(predict(lasso1,s=1,type="coef",mode="lambda")$coef)
			beta.new[-1]<-as.vector(glmnet.1$beta[,10])
			#if(log(sum(which.use))/2*sum(beta.new[-1]!=0)>.9*sum(which.use)) 
			#	beta.new[-1]<-as.vector(predict(lasso1,s=min(c(floor(.9*sum(which.use)*2/log(n)),
			#	dim(X.new)[2]*.8)),type="coef",mode="step")$coef)
			#beta.curr<-.5*beta.new+.5*beta.curr
			beta.new[1]<-mean(y[which.use])-mean(X.new[which.use,]%*%beta.new[-1])
			run.diff<-(mean((beta.new[-1]-beta.curr[-1])^2)/sum(beta.new[-1]^2+1e-10))
			if(run.diff<1e-6) break
			which.use<-(y*cbind(1,X.new)%*%beta.new)<=y^2
		}
		
	#Find fitted values.
		fits<-pmin(abs(cbind(1,X.new)%*%beta.new),y^2)*sign((cbind(1,X.new)%*%beta.new))
		fits2<-cbind(1,X.new)%*%beta.new
		fits2[sign(fits2)==sign(y)]<-pmin(abs(fits2[sign(fits2)==sign(y)]),y[sign(fits2)==sign(y)]^2)*sign(fits2[sign(fits2)==sign(y)])
	

	#Calculate degrees of freedom
		edf<-  1+sum(beta.new[-1]!=0)#*log(n)/2#sum(diag(hat.mat))#+1-sum(beta.new[-c(1:(dim(X.c)[2]+1))]!=0)
	#Left in for estimates of standard deviation.  	
		stdev<-NA
	
	#GCV statistic
		loss.run[i,j]<-sum((y^2-y*fits2)^2)/(sum(which.use)-edf)^2*(mean(which.use))^2
		
	#Gather minimum loss function.
		if(i*j>1) if(loss.run[i,j]==min(loss.run[loss.run!=0])) min.params<-c(lambda.c[i],lambda.t[j])

	}
	}
	
	#Scale betas back.
		beta.new[c(1:dim(X.c)[2]+1)]<-beta.new[c(1:dim(X.c)[2]+1)]/exp(lambda.c[i])
		beta.new[-c(1:(dim(X.c)[2]+1))]<-beta.new[-c(1:(dim(X.c)[2]+1))]/exp(lambda.t[i])
		beta.new<-as.vector(beta.new)
		X<-as.matrix(X)

	#Calculate intercept
		beta.new[1]<-mean(y[which.use])-mean((X%*%beta.new[-1])[which.use])

	output<-list("lambdas"=min.params,"beta"=beta.new,"fits"=fits,"loss"=log(loss.run),"marg"=mean(which.use),"edf"=edf,
		"sd"=stdev)
	invisible(output)

	}
	
	
#####################################################

#####################################################


	search.lambda<-function(y=y,X.c=X.c,X.t=X.t){

	lambda.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,lambda[1],lambda[2])$loss)
	lambda.c.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,lambda,lambda.t)$loss)
	lambda.t.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,lambda.c,lambda)$loss)
	lambda.c.old<-lambda.t.old<-999

	lambda.c.seek<-seq(-15,10,1)
	lambda.t.seek<-seq(-15,10,1)
	lambda.t<-25
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))

	#lambda.c<- -15
	#lambda.t<-5
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t

	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-5,lambda.c+5,1)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-5,lambda.t+5,1)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}



	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-2.5,lambda.c+2.5,.5)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-2.5,lambda.t+2.5,.5)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}

	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-1.5,lambda.c+1.5,.25)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-1.5,lambda.t+1.5,.25)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.5,lambda.c+.5,.1)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.5,lambda.t+.5,.1)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}



	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.25,lambda.c+.25,.05)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.25,lambda.t+.25,.05)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.15,lambda.c+.15,.025)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.15,lambda.t+.15,.025)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.1,lambda.c+.1,.01)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.1,lambda.t+.1,.01)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}

	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.05,lambda.c+.05,.005)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.05,lambda.t+.05,.005)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.025,lambda.c+.025,.0025)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.025,lambda.t+.025,.0025)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.01,lambda.c+.01,.001)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.01,lambda.t+.01,.001)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}


	for(lambda.loop in 1:5){
	lambda.c.seek<-seq(lambda.c-.0025,lambda.c+.0025,.0005)
	lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.0025,lambda.t+.0025,.0005)
	lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	print(c(lambda.c,lambda.t))
	print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
	lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
	}




	
	optim.lambda<-c(lambda.c,lambda.t)


	invisible(optim.lambda)
	}
	
#####################################################
