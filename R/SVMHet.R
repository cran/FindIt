## This file contains the functions used for analyzing causal heterogeneity
## using support vector machines under a double LASSO constraint.

## The functions consist of 
## scale.func: demeans and standardizes variables.
## scale.func.2: gives variables a standard deviation of one
## FindIt: a wrapper function for SVM.func
## maketwoway<-rescales covariates and makes all twoway interactions.
## SVM.func: the workhorse function.  Given values of y, control, treatment, and lambdas, 
## it returns betas, the loss (GCV), fitted values, and size of the margin.
## find.lambda: Conducts the alternating line search.
## A large value of lambda.c is selected (10), and alternating
## line searches of decreasing radius are conducted.
## lars: An adaptation of the LARS algorithm of Efron, et al. 2004.
## Due to the large numbers being considered in the fitting,
## the code is adapted to not crash.  The original code is commented out.
## The remainder of the code is untouched.

##################################################################
##################################################################
##################################################################

scale.func<-function(x) {
    x2<-x-mean(x)
    x2/sd(x2)}

scale.func.2<-function(x) {x2<-x
                           x2*sd(x2)}

#####################################################

##change here. 
FindIt <- function(model.treat, model.main,model.int,data=NULL,
                   search.lambdas=TRUE,lambdas=NULL,
                   wts=1,scale.c=1,scale.int=1,
                   fit.glmnet=TRUE,
                   make.twoway=TRUE,make.allway=TRUE,
                   reference.main=NULL,
                   type="binary",
                   treat.type="multiple",
                   threshold=0.999999,
                   make.reference=TRUE,
		   nway=2){

    ## if(make.allway==FALSE & is.null(reference.main)){
    ##     warning("Need to specify the corresponding reference of treatment matrix")
    ##     break
    ## }
    unique <- FALSE

    if(sum(wts<0)>0) {
        print("Only non-negative weights allowed")
        break
    }

    ## Checke whether the model has the main effect covariates or not. 
    if(!missing(model.main)|!missing(model.int)){
        main <- TRUE
    }    
    if(missing(model.int) & !missing(model.main)){
        model.int <- model.main
        main <- TRUE
    }    
    if(!missing(model.int) & missing(model.main)){
        model.main <- model.int
        main <- TRUE
    }    
    if(missing(model.main) & missing(model.int)){
        main <- FALSE
    }

    ## Extract the data
    treat.frame     <- model.frame(model.treat,data, na.action=NULL)
    if(main){
        main.frame  <- model.frame(model.main, data, na.action=NULL)
        int.frame   <- model.frame(model.int,  data, na.action=NULL)
    }

    ## NA : List-wise deletion.
    data.which <- c()
    if(main){
        for(i in 1:nrow(treat.frame)){
            if(any(is.na(main.frame[i,])) | any(is.na(int.frame[i,])) |
               any(is.na(treat.frame[i,])) )
                {data.which[i] <- 0}
            else{data.which[i] <- 1}
        }
    }
    if(main==FALSE){
        for(i in 1:nrow(treat.frame)){
            if(any(is.na(treat.frame[i,])))
                {data.which[i] <- 0}
            else{data.which[i] <- 1}
        }
    }
    treat.frame <- treat.frame[data.which==1,]
    treat.orig <- treat.frame[,-1]
    treat <- treat.frame[,-1]
    treat <- as.matrix(treat)
    y     <- y.orig      <- treat.frame[,1]
    if(type=="binary"){
        y.orig <- y
        y <-(2*y-1)        
    }
    
    if(main){
        X.c   <- X.c.orig    <- main.frame[data.which==1,]
        X.int <- X.int.orig  <- int.frame[data.which==1,]
        X.c   <- X.c.orig <- as.data.frame(X.c)
        X.int <- X.int.orig <- as.data.frame(X.int)

        if(ncol(X.c)>1){
            if(any(apply(X.c,2,FUN=function(x) length(unique(x)) == 1))){
                warning("There is a column with no variation
                         in the main effect matrix.")
            }
        }
        if(ncol(X.int)>1){
            if(any(apply(X.int,2,FUN=function(x) length(unique(x)) == 1))){
                warning("There is a column with no variation
                         in the interaction matrix.")
            }   
        }
        if(treat.type=="single" & main==FALSE){
            warning("Need to specify the covariates to interact with treatment.")
        }
    }

    ## Make multiple-treatment matrix
    if(treat.type == "multiple"){
        if(make.allway==TRUE){
            Allway <- makeallway(treat,threshold,make.reference=make.reference,nway=nway)
            X.t    <- Allway$FinalData
            reference.main <- Allway$reference
        }
        if(make.allway==FALSE){
            X.t    <- as.matrix(treat)
            ## reference.main <- reference.main
            reference.main <- "No reference"
        }
    }

    ## print("X.t")
    ## print(head(X.t))

    ## Maketwoway    
    if(main){
        X.c.old <- X.c
        X.int.old <- X.int
        if(make.twoway==TRUE){
            if(ncol(X.c.old)>1){
                XC         <- maketwoway(X.c,center=TRUE,wts=wts)
                X.c        <- XC$X
                scale.c    <- XC$scale.X
                XC.deleted <- XC$D
                frame.meanC <- XC$frame.mean
            }
            if(ncol(X.c.old)==1){
                form <- ~ .
                X.c.old <- as.data.frame(X.c.old)
                X.c1 <- model.matrix(form,data=X.c.old)
                X.c1 <- as.data.frame(X.c1)
                X.c <- X.c1[,-1]
                if(is.vector(X.c)){
                    scale.c <- sd(X.c)
                    frame.meanC <- mean(X.c)
                    X.c <- (X.c-mean(wts^.5*X.c))/sd(X.c)    
                    X.c <- as.data.frame(X.c)
                    colnames(X.c) <- colnames(X.c1)[-1]                   
                }else{
                    scale.c <- apply(X.c,2,sd)
                    frame.meanC <- apply(X.c,2,function(x) mean(wts^.5*x))
                    X.c <- apply(X.c,2,FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                }
            }
            if(ncol(X.int.old)>1){
                XInt       <- maketwoway(X.int,center=TRUE,wts=wts)
                X.int      <- XInt$X
                scale.int  <- XInt$scale.X
                XInt.deleted <- XInt$D
                frame.meanInt <- XInt$frame.mean
            }
            if(ncol(X.int.old)==1){
                form <- ~ .
                X.int.old<- as.data.frame(X.int.old)
                X.int1 <- model.matrix(form,data=X.int.old)
                X.int1 <- as.data.frame(X.int1)
                X.int <- X.int1[,-1]
                if(is.vector(X.int)){
                    scale.int <- sd(X.int)
                    frame.meanInt <- mean(X.int)
                    X.int <- (X.int-mean(wts^.5*X.int))/sd(X.int) 
                    X.int <- as.data.frame(X.int)
                    colnames(X.int) <- colnames(X.int1)[-1]
                }
                else{
                    scale.int <- apply(X.int,2,sd)
                    frame.meanInt <- apply(X.int,2,function(x) mean(wts^.5*x))
                    X.int <- apply(X.int,2,FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                }
            }
        }
        
        ## Remove all elements with no variance and add an intercept.
        if(is.matrix(X.c)){
            scale.c <- scale.c[apply(X.c,2,sd)>0]
            X.c<-X.c.1 <- X.c[,apply(X.c,2,sd)>0]
            ## X.c<-unique(X.c,MARGIN=2)
            scale.c <- scale.c[colnames(X.c.1)==colnames(X.c)]
            frame.meanC <- frame.meanC[colnames(X.c.1)==colnames(X.c)]
        }
        if(is.vector(X.c)){
            X.c<-cbind(1,X.c)
            colnames(X.c)[1]<-"Intercept"
        }
        else{
            if(sd(X.c[,1])>0){
                X.c<-cbind(1,X.c)
                colnames(X.c)[1]<-"Intercept"
            }
        }
        name.c <- colnames(X.c)
        ## This X.c is used in the SVM.
        
        ## Remove all elements with no variance and add an intercept.
        if(is.matrix(X.int)){
            scale.int <- scale.int[apply(X.int,2,sd)>0]
            frame.meanInt <- frame.meanInt[apply(X.int,2,sd)>0]
            X.int<- X.int.1<- X.int[,apply(X.int,2,sd)>0]
            scale.int <- scale.int[colnames(X.int.1)==colnames(X.int)]
            frame.meanInt <- frame.meanInt[colnames(X.int.1)==colnames(X.int)]
        }
        if(is.vector(X.int)){
            X.int<-cbind(1,X.int)
            colnames(X.int)[1]<-"Intercept"
        }
        else{
            if(sd(X.int[,1])>0){
                X.int<-cbind(1,X.int)
                colnames(X.int)[1]<-"Intercept"
            }
        }
        name.int <- colnames(X.int)
        ## This X.int is used in the SVM

        ## Make X.t
        ## I changed X.c to X.int        
        if(treat.type == "single") {
            X.t<-cbind(treat*1e4,(treat>0)*X.int[,-1])
            X.t[treat!=0,-1]<-apply(X.t[treat!=0,-1],2, FUN=function(x) x-mean(x))
            ## This X.t is used.            
            colnames(X.t)<-c("treat",paste("treat",colnames(X.int)[-1],sep=":"))
        }
    }

###################################
    ## Create the FindIt.fit function.
################################### 

    FindIt.fit<-function(y,X.c,X.t,X.int,treat,treat.type,                     
                         search.lambdas,lambdas,
                         wts,scale.c,scale.int,
                         reference.main,
                         fit.glmnet,
                         type,
                         X.c.orig,X.int.orig,treat.orig,
                         threshold,
                         unique,
                         main){

        ## Give wts mean 1.
	wts<-wts/mean(wts)
        
	y <- y*wts^.5
        X.t <- X.t*wts^.5
        name.t <- colnames(X.t)
        if(main){
            X.c <- X.c*wts^.5
            ## X.int <- X.int*wts^.5
            name.c <- colnames(X.c)            
            name.int <- colnames(X.int)           
        }
        ## treat<-treat*wts^.5

        

        X.t <- as.matrix(X.t)
        if(main){
            X.c <- as.matrix(X.c)
            X.int <- as.matrix(X.int)
        }
	
	if(search.lambdas==TRUE){
            if(main){
                lambdas <-search.lambda(y,X.c,X.t,
                                        fit.glmnet=fit.glmnet,
                                        main=main,
                                        type=type)
            }
            if(main==FALSE){
                lambdas <- search.lambda.nocov(y=y,X.t=X.t,
                                               fit.glmnet=fit.glmnet,
                                               main=main,
                                               type=type)
            }
        }
        

        if(main){
            A<-SVM.func(y,X.c[,-1],X.t,
                        lambdas[1],lambdas[2],
                        fit.glmnet=fit.glmnet,
                        type=type,
                        main=main)
        }
        if(main==FALSE){
            A<-SVM.func(y=y,X.c[,-1],X.t=X.t,
                        lambda.c=1,lambda.t=lambdas,
                        fit.glmnet=fit.glmnet,
                        type=type,
                        main=main)
        }
        ## When we input X.c doesnt have intercept. 
        
	if(treat.type=="single"){
            ## Keep the names of coefficients.
            reference.coef <- as.data.frame(A$beta[-c(1:(dim(X.c)[2]))])
            rownames(reference.coef) <- c("treat",colnames(X.int[,-1]))
            if(type=="binary"){
                reference.coef <- reference.coef/2
            }
            ## Generate outcome treatment and no treatment.
            ## I changed the code.
            preds.treat<-cbind(X.c,cbind(1e4,X.int[,-1]))%*%A$beta
            preds.control<-cbind(X.c,0*cbind(1e4,X.int[,-1]))%*%A$beta
            ## preds <- cbind(X.c,X.t) %*% A$beta
            ## diff  <- mean(y.inter)- mean(preds)
            
            ## cbind(1e4,X.int[,-1]) is X.t when treat =1.
            ## preds.treat<-cbind(X.c,cbind(1e4,X.c[,-1]))%*%A$beta
            ## in here, X.c includes the intercept.
            ## 1e4 means treated group         
            ## preds.control<-cbind(X.c,0*cbind(1e4,X.c[,-1]))%*%A$beta

            if(type=="binary"){
                preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                ATE <- mean(preds.treat-preds.control)/2
                preds.diff<- (preds.treat-preds.control)/2
            }
            if(type=="continuous"){
                preds.diff<-preds.treat-preds.control
                ATE<-mean(preds.treat-preds.control)
            }

            scale.out<-c(1,scale.c,1e-4,scale.int)
            Treatment.version <- NULL
            
	}
        if(treat.type == "multiple"){
            ## make the reference matrix.
            if(main){
                ## preds <- cbind(X.c,X.t) %*% A$beta
                ## diff <- mean(y.inter) - mean(preds)
                scale.out<-c(1,scale.c,rep(1,dim(X.t)[2]))
                coef.orig.ref <- A$beta[-c(1:ncol(X.c))]/scale.out[-c(1:ncol(X.c))]
            }
            if(main==FALSE){
                ## preds <- cbind(1,X.t) %*% A$beta
                ## diff <- mean(y.inter) - mean(preds)
                scale.out <- c(1,rep(1,ncol(X.t)))
                ## I need one more for intercept.
                coef.orig.ref <- A$beta[-1]/scale.out[-1]
            }
            
            if(type=="binary"){
                reference.coef <- as.data.frame(signif(coef.orig.ref/2,6))
            }
            if(type=="continuous"){
                reference.coef <- as.data.frame(signif(coef.orig.ref,6))
            }
            
            rownames(reference.coef) <- colnames(X.t)

###############################
            ## create the Treatment.version
###############################
            coefs <- A$beta
            
            if(type=="binary"){
                if(main){
                    preds.treat    <- cbind(X.c,X.t) %*% coefs
                    preds.control  <- X.c %*% coefs[c(1:ncol(X.c))]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                    preds.diff <- (preds.treat-preds.control)/2
                }
                if(main==FALSE){
                    preds.treat <- X.t %*% coefs[-1]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.diff <- preds.treat/2
                }
                ATE <- mean(preds.diff)
            }
            if(type=="continuous"){
                if(main){
                    preds.diff    <- X.t %*% coefs[-c(1:(dim(X.c)[2]))]
                }
                if(main==FALSE){
                    preds.diff <- X.t %*% coefs[-1]
                }
                ATE <- mean(preds.diff)
            }
            
            if(unique==FALSE){
                ##preds.control <- X.c %*% coefs[1:ncol(X.c)]
                ##preds <- cbind(X.c,X.t)%*%coefs
                ##preds.diff <- preds - preds.control
                preds.diff <- preds.diff
                Treatment.version <- cbind(preds.diff,treat.orig)
                colnames(Treatment.version) <- c("Treatment.effect",
                                                 colnames(treat.orig))
                
                Treatment.version <- as.data.frame(Treatment.version)
                
            }
            if(unique==TRUE){
                X.t.u <- X.t
                preds.diff.u <- as.data.frame(preds.diff)
                rownames(preds.diff.u) <- rownames(X.t.u) <- seq(1:nrow(X.t.u))
                treat.unique  <- unique(treat.orig,MARGIN=1)
                X.t.unique    <- unique(X.t.u, MARGIN=1)
                X.t.unique2   <- unique(X.t, MARGIN=1)
                preds.diff   <- preds.diff.u[rownames(preds.diff.u) %in%
                                             rownames(X.t.unique),]
                
                preds.diff <- as.data.frame(preds.diff)
                ##name back
                rownames(X.t.unique) <- rownames(X.t.unique2)
                rownames(preds.diff) <- rownames(X.t.unique2)

                Treatment.version <- cbind(preds.diff,treat.unique)
                colnames(Treatment.version) <- c("Treatment.effect",
                                                 colnames(treat.orig))
                
                Treatment.version <- as.data.frame(Treatment.version)
                Treatment.version<-
                    Treatment.version[order(Treatment.version$Treatment.effect,
                                            decreasing=TRUE),]
            }            
	}
        
        if(main){
            names.out <- c(colnames(X.c),colnames(X.t))
        }
        else{
            names.out <- c("Intercept",colnames(X.t))
        }

        if(main){
            output<-(list("coefs"=A$beta,"coefs.orig"=A$beta/scale.out,
                          "fit"=A$fit,
                          "scale.out"=scale.out,
                          "scale.c" = scale.c,
                          "scale.int"=scale.int,
                          "frame.meanC" = frame.meanC,
                          "frame.meanInt" = frame.meanInt,
                          "names.out"=names.out,
                          "y"=y,"X.c"=X.c,"X.t"=X.t,"X.int"=X.int,
                          "GCV"=exp(A$loss),
                          "name.c" =name.c,
                          "name.int" = name.int,
                          "name.t" =name.t,
                          "ATE" = ATE,
                          "Treatment.version" = Treatment.version,
                          "lambdas"=lambdas,
                          "reference.main"=reference.main,
                          "reference.coef"=reference.coef)
                     )
        }
        if(main==FALSE){
            output<-(list("coefs"=coefs,"coefs.orig"=A$beta/scale.out,
                          "fit"=A$fit,
                          "scale.out"=scale.out,
                          "names.out"=names.out,
                          "y"=y,"X.t"=X.t,
                          "GCV"=exp(A$loss),
                          "name.t" =name.t,
                          "ATE" = ATE,
                          "Treatment.version" = Treatment.version,
                          "lambdas"=lambdas,
                          "reference.main"=reference.main,
                          "reference.coef"=reference.coef)
                     )
        }
	class(output)<-c("FindIt","list")
	output
    }

###########################################
    ## The FindIt.fit function is finished
##########################################
    
    if(main==FALSE){
        X.c <- X.int <- scale.c <- scale.int <- frame.meanC <- frame.meanInt <- NULL
        X.c.orig <- X.int.orig <- NULL
    }

    
    Fit.output <- FindIt.fit(y=y,X.c=X.c,X.t=X.t,X.int=X.int,treat=treat,
                             treat.type=treat.type,
                             search.lambdas = search.lambdas,
                             lambdas = lambdas,
                             wts = wts,
                             scale.c = scale.c,
                             scale.int = scale.int,
                             reference.main = reference.main,
                             fit.glmnet = fit.glmnet,
                             type = type,
                             X.c.orig = X.c.orig,
                             X.int.orig = X.int.orig,
                             treat.orig = treat.orig,
                             threshold = threshold,
                             unique = unique,
                             main=main)

#########################################################
#### Make the reference matrix when treat.type is multiple.
#########################################################
    if(make.reference==TRUE & make.allway==TRUE){
        if(treat.type=="multiple"){
            reference.main <- Fit.output$reference.main
            reference.coef <- Fit.output$reference.coef
            coef   <- c()
            
            for(i in 1:nrow(reference.main)){
                if(rownames(reference.main)[i] %in% rownames(reference.coef)){
                    coef[i] <-
                        reference.coef[rownames(reference.coef)==
                                       rownames(reference.main)[i],1]
                }
                else{coef[i] <- NA
                 }
            }
            
            reference.Final <- as.data.frame(cbind(reference.main,coef))
            colnames(reference.Final) <- c(colnames(reference.main),
                                           "Coefficients"
                                           )
        }
    }else{
        reference.Final <- "No reference"
    }
    if(treat.type=="single"){
        reference.Final <- "No reference"
    }

    if(main){
        output<-(list("coefs"=Fit.output$coefs,"coefs.orig"=Fit.output$coefs.orig,
                      "fit"=Fit.output$fit,                     
                      "scale.out"=Fit.output$scale.out,
                      "scale.c" = Fit.output$scale.c,
                      "scale.int" = Fit.output$scale.int,
                      "frame.meanC" = frame.meanC,
                      "frame.meanInt" = frame.meanInt,
                      "names.out"=Fit.output$names.out,
                      "y"=Fit.output$y,"X.c"=Fit.output$X.c,
                      "X.t"=Fit.output$X.t, "X.int"=Fit.output$X.int,
                      "treat" = Fit.output$treat,
                      "GCV"=Fit.output$GCV,
                      "y.orig"=y.orig, "X.c.orig" = X.c.orig,
                      "X.int.orig" = X.int.orig, "treat.orig" = treat.orig,
                      "name.c" = Fit.output$name.c,
                      "name.int" = Fit.output$name.int,
                      "name.t" = Fit.output$name.t,
                      "treat.type" = treat.type,
                      "model.main"=model.main,
                      "model.int" = model.int,
                      "model.treat"=model.treat,
                      "ATE"=Fit.output$ATE,
                      "data"=data,
                      "Treatment.version" = Fit.output$Treatment.version,
                      "lambdas"=Fit.output$lambdas,
                      "threshold" = threshold,
                      "make.twoway" = make.twoway,
                      "make.allway" = make.allway,
                      "model.main"=model.main,
                      "model.int"=model.int, "model.treat"=model.treat,
                      "reference"= reference.Final,"type" = type,
                      "main" =main,
		      "nway" =nway
                      ))
    }
    if(main==FALSE){
        output<-(list("coefs"=Fit.output$coefs,"coefs.orig"=Fit.output$coefs.orig,
                      "fit"=Fit.output$fit,
                      "scale.out"=Fit.output$scale.out,
                      "names.out"=Fit.output$names.out,
                      "y"=Fit.output$y,                    
                      "X.t"=Fit.output$X.t,                     
                      "treat" = Fit.output$treat,
                      "GCV"=Fit.output$GCV,
                      "y.orig"=y.orig,                     
                      "treat.orig" = treat.orig,
                      "name.t" = Fit.output$name.t,
                      "treat.type" = treat.type,
                      "model.treat"=model.treat,
                      "ATE"=Fit.output$ATE,
                      "data"=data,
                      "Treatment.version" = Fit.output$Treatment.version,
                      "lambdas"=Fit.output$lambdas,                      
                      "threshold" = threshold,
                      "make.twoway" = make.twoway,
                      "make.allway" = make.allway,
                      "reference"= reference.Final,
                      "type" = type,
                      "main" =main,
		      "nway" =nway
                      ))
    }
    
    class(output)<-c("FindIt","list")
    output

}

predict.FindIt<-
    function(object,newdata,sort=TRUE,decreasing =TRUE,
             wts=1,unique=FALSE,...){
        treat.type   <- object$treat.type
        type <- object$type
        threshold <- object$threshold
        make.twoway <- object$make.twoway
        make.allway  <- object$make.allway
        main <- object$main
	nway <- object$nway
        model.treat <- object$model.treat
        if(main){
            model.main <- object$model.main
            model.int  <- object$model.int 
        }

        if(missing(newdata)){
            different <- 0
        }
        if(!missing(newdata)){
            different <- 1
            newdata <- newdata
        }
        
        if(different==0){
            coefs2       <- object$coefs.orig
            y            <- object$y
            treat        <- object$treat
            y.orig     <- object$y.orig
            treat.orig   <- object$treat.orig
            X.t          <- object$X.t
            if(main){
                X.c          <- object$X.c        
                X.int        <- object$X.int
                X.int.orig <- object$X.int.orig
                scale.c      <- object$scale.c
                scale.int    <- object$scale.int
                X.c.orig   <- object$X.c.orig
            }                      
        }
        if(different==1){
            terms.treat     <- terms(model.treat)
            treat.frame     <- model.frame(terms.treat,data=newdata, na.action=NULL) 
            if(main){
                terms.main  <- terms(model.main)
                terms.int   <- terms(model.int)
                main.frame      <- model.frame(terms.main,data=newdata,na.action=NULL)
                int.frame       <- model.frame(terms.int,data=newdata,na.action=NULL)
            }
            ## NA manipulation: List-wise deletion
            data.which <- c()
            if(main){
                for(i in 1:nrow(treat.frame)){
                    if(any(is.na(main.frame[i,])) | any(is.na(int.frame[i,])) |
                       any(is.na(treat.frame[i,])) )
                        {data.which[i] <- 0}
                    else{data.which[i] <- 1}
                }
            }
            if(main==FALSE){
                for(i in 1:nrow(treat.frame)){
                    if(any(is.na(treat.frame[i,])))
                        {data.which[i] <- 0}
                    else{data.which[i] <- 1}
                }
            }

            treat.frame <- treat.frame[data.which==1,]
            treat <- treat.orig <- treat.frame[,-1]
            
            y     <- y.orig      <- treat.frame[,1]
            if(type=="binary"){
                y <- 2*y-1
            }
            y <- y*wts^.5
            ## treat<-treat*wts^.5
            
            if(main==TRUE){
                X.c   <- main.frame[data.which==1,]
                X.int <- int.frame[data.which==1,]
                X.int.orig <- X.int
                
                if(make.twoway==TRUE){
                    frame.meanPreC <- object$frame.meanC
                    XC      <- maketwoway(X.c,wts=wts,center=TRUE,frame.meanPre=frame.meanPreC,predict=TRUE)
                    X.c     <- XC$X
                    scale.c <- XC$scale.X
                    frame.meanPreInt <- object$frame.meanInt
                    XInt    <- maketwoway(X.int,wts=wts,center=TRUE,frame.meanPre=frame.meanPreInt,predict=TRUE)
                    X.int   <- XInt$X
                    scale.int <- XInt$scale.X
                    X.c <- as.matrix(X.c)
                    X.int <- as.matrix(X.int)

                    ## print("check after fix")
                    ## print(dim(X.c))
                    ## print(colnames(X.c))
                    ## print(dim(X.int))
                    ## print(colnames(X.c))
                }
                
                
                ## Remove all elements with no variance and add an intercept.
                ## X.c<-X.c[,apply(X.c,2,sd)>0]
                ## in the prediction model, it is allowed.
                X.c.m <- X.c
                ## X.c<-unique(X.c,MARGIN=2)
                scale.c <- scale.c[colnames(X.c.m)==colnames(X.c)]
                if(all(X.c[,1]==1)==FALSE){
                    X.c<-cbind(1,X.c)
                    colnames(X.c)[1]<-"Intercept"
                }
                X.c <- X.c*wts^.5
                ## X.int <- X.int*wts^.5  

                if(treat.type=="single"){
                    ## Remove all elements with no variance and add an intercept.
                    ## X.int<-X.int[,apply(X.int,2,sd)>0]
                    ## in the prediction, it is OK.
                    X.int.m <- X.int
                    ## X.int<-unique(X.int,MARGIN=2)
                    scale.int <- scale.int[colnames(X.int.m)==colnames(X.int)] 
                    if(all(X.int[,1]==1)==FALSE){
                        X.int<-cbind(1,X.int)
                        colnames(X.int)[1]<-"Intercept"
                    }       
                    X.t<-cbind(treat*1e4,(treat>0)*X.int[,-1])
                    X.t[treat!=0,-1]<-
                        apply(X.t[treat!=0,-1],2,
                              FUN=function(x) x-mean(x))
                    colnames(X.t)<-c("treat",paste("treat",
                                                   colnames(X.int)[-1],
                                                   sep=":"))
                }
                
            }
            
            if(treat.type=="multiple"){
                if(make.allway==TRUE){
                    treat <- as.matrix(treat)
                    Allway <- makeallway(treat,threshold,nway=nway)
                    X.t    <- Allway$FinalData
                    X.t    <- as.matrix(X.t)
                    reference.main <- Allway$reference
                }
                if(make.allway==FALSE){
                    X.t    <- as.matrix(treat)
                    ## reference.main <- reference.main
                    reference.main <- "No Reference"
                }
            }
            X.t <- X.t*wts^.5
        }

        
        ## Adjust the coefficients.
        ## When we have smaller variation in the new data.
        if(different==1){
            if(treat.type=="single"){
                coefs.orig         <- object$coefs.orig
                names(coefs.orig)  <- c(object$name.c,"treat",object$name.int[-1])      
                coefs.orig.c   <- coefs.orig[1:ncol(object$X.c)]
                coefs.orig.int <- coefs.orig[-c(1:ncol(object$X.c))]                  
                ## X.c include intercept.
                ## I use the same formula, so the order is the same.
                coefs.c2    <- coefs.orig.c[is.element(object$name.c,
                                                       colnames(X.c))]
                coefs.int2  <- coefs.orig.int[is.element(object$name.int,
                                                         colnames(X.int))]
                coefs2      <- c(coefs.c2,coefs.int2)
            }
            if(treat.type=="multiple"){
                coefs.orig <- object$coefs.orig
                if(main==FALSE){
                    names(coefs.orig)  <- c("Intercept", object$name.t)
                    coefs.orig.t <- coefs.orig[-1]
                    coefs.t1   <- coefs.orig.t[is.element(object$name.t,
                                                          colnames(X.t))]
                    coefs2   <- c(coefs.orig[1], coefs.t1)
                    names(coefs2) <- c("Intercept", names(coefs.t1))
                }                
                if(main){
                    names(coefs.orig)  <- c(object$name.c,object$name.t)   
                    coefs.orig.c   <- coefs.orig[1:ncol(object$X.c)]
                    coefs.orig.t   <- coefs.orig[-c(1:ncol(object$X.c))]              
                    ## X.c include intercept.
                    ## I use the same formula, so the order is the same.
                    coefs.c2    <- coefs.orig.c[is.element(object$name.c,
                                                           colnames(X.c))]
                    coefs.t2    <- coefs.orig.t[is.element(object$name.t,
                                                           colnames(X.t))]
                    coefs2      <- c(coefs.c2,coefs.t2)
                }
                X.t <- X.t[,is.element(colnames(X.t), names(coefs.t2))]
            }
        }

        if(main){
            if(any(is.element(colnames(X.c),object$name.c)==FALSE)){
                warning("The new data set has more variation
                         than the original model.")
            }
            
            if(any(is.element(colnames(X.int), object$name.int)==FALSE)){
                warning(
                    "The new data set has more variation
                     than the original model."
                    )
            }            
        }

        if(any(is.element(colnames(X.t), object$name.t)==FALSE)){
            warning(
                "The new data set has more variation than the original model."
                )
        }

        
        if(treat.type=="single"){        
            scale.out<-c(1,scale.c,1e-4,scale.int)
            coefs <- coefs2 * scale.out            
            ## Generate outcome treatment and no treatment.
            preds.treat<-cbind(X.c,cbind(1e4,X.int[,-1]))%*%coefs
            ## preds.treat<-cbind(X.c,cbind(1e4,X.c[,-1]))%*%coefs
            preds.control<-cbind(X.c,0*cbind(1e4,X.int[,-1]))%*%coefs
            ## preds.control<-cbind(X.c,0*cbind(1e4,X.c[,-1]))%*%coefs
            preds <- cbind(X.c,X.t)%*%coefs
            
            if(type=="binary"){
                preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                preds         <- sign(preds)*pmin(abs(preds),1)
                ATE <- mean(preds.treat-preds.control)/2
                preds.diff<- (preds.treat-preds.control)/2
            }
            if(type=="continuous"){
                preds     <- preds
                preds.diff<-preds.treat-preds.control
                ATE<-mean(preds.treat-preds.control)
            }

            pred.data <- cbind(preds.diff,y.orig,treat.orig,X.int.orig)
            pred.data <- as.data.frame(pred.data)
            colnames(pred.data) <- c("Treatment.effect",
                                     "outcome",
                                     "treatment",
                                     colnames(X.int.orig))
            
            if(sort==TRUE){
                pred.data.out <- pred.data[order(pred.data$Treatment.effect,
                                                 decreasing=decreasing),]
            }
            else{pred.data.out <- pred.data}
            
        }
        if(treat.type=="multiple"){            
            if(main){
                scale.out<-c(1,scale.c,rep(1,dim(X.t)[2]))
            }
            if(main==FALSE){
                scale.out <- c(1,rep(1,ncol(X.t)))
            }

            coefs <- coefs2 * scale.out

            if(type=="binary"){
                if(main){
                    X.t1 <- as.matrix(X.t)
                    preds.treat    <- cbind(X.c,X.t1) %*% coefs
                    preds.control  <- X.c %*% coefs[c(1:ncol(X.c))]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.control <- sign(preds.control)*pmin(abs(preds.control),1)
                    preds.diff <- (preds.treat-preds.control)/2
                }
                if(main==FALSE){
                    preds.treat <- X.t %*% coefs[-1]
                    preds.treat   <- sign(preds.treat)*pmin(abs(preds.treat),1)
                    preds.diff <- preds.treat/2
                }
                ATE <- mean(preds.diff)
            }
            if(type=="continuous"){
                if(main){
                    preds.diff    <- X.t %*% coefs[-c(1:(dim(X.c)[2]))]
                }
                if(main==FALSE){
                    preds.diff <- X.t %*% coefs[-1]
                }
                ATE <- mean(preds.diff)
            }

            ## Make the SATE for Multiple treatment version.            
            if(unique==FALSE){
                ##preds.control <- X.c %*% coefs[1:ncol(X.c)]
                ##preds <- cbind(X.c,X.t)%*%coefs
                ##preds.diff <- preds - preds.control
                preds.diff <- preds.diff
                pred.data <- cbind(preds.diff,treat.orig)
                ##print(head(pred.data))
                ##print(head(treat.orig))
                colnames(pred.data) <- c("Treatment.effect",                          
                                         colnames(treat.orig)
                                         )
                pred.data <- as.data.frame(pred.data)
                
                if(sort==TRUE){
                    pred.data.out <- pred.data[order(pred.data$Treatment.effect,
                                                     decreasing=decreasing),]
                }
                else{pred.data.out <- pred.data}
                
            }
            if(unique==TRUE){
                X.t <- as.matrix(X.t)
                X.t.u <- X.t
                preds.diff.u <- as.data.frame(preds.diff)
                rownames(preds.diff.u) <- rownames(X.t.u) <- seq(1:nrow(X.t.u))
                treat.unique  <- unique(treat.orig,MARGIN=1)
                X.t.unique    <- unique(X.t.u, MARGIN=1)
                X.t.unique2   <- unique(X.t, MARGIN=1)
                preds.diff   <- preds.diff.u[rownames(preds.diff.u) %in%
                                             rownames(X.t.unique),]
                
                preds.diff <- as.data.frame(preds.diff)
                ##name back
                rownames(X.t.unique) <- rownames(X.t.unique2)
                rownames(preds.diff) <- rownames(X.t.unique2)

                pred.data <- cbind(preds.diff,treat.unique)
                colnames(pred.data) <- c("Treatment.effect",
                                         colnames(treat.orig)
                                         )
                pred.data <- as.data.frame(pred.data)
                pred.data.out <-
                    pred.data[order(pred.data$Treatment.effect,
                                    decreasing=decreasing
                                    ),]
            }
            
        }
        pred.data.out <- as.data.frame(pred.data.out)
        out <- list("treat.type"=treat.type, "ATE"=ATE, "data" =pred.data.out,
                    "coefs"=coefs,"orig.coef"=coefs2)
        class(out) <- "PredictFindIt"
        invisible(out)
    }


plot.PredictFindIt <- function(x,labels="index",...){
    object <- x
    treat.type <- object$treat.type
    if(treat.type=="single"){
        pred.data <- object$data
        ATE <- object$ATE
        pred.data.out.p <- pred.data[order(pred.data$Treatment.effect,
                                           decreasing=FALSE),]   
        xp <- seq(1:nrow(pred.data.out.p))
        zero.1 <- min(which(min(abs(pred.data.out.p$Treatment.effect))==
                            abs(pred.data.out.p$Treatment.effect)))
        zero <- xp[zero.1]
        low <- min(pred.data.out.p$Treatment.effect)
        plot(xp,pred.data.out.p$Treatment.effect,type="l",col="red",
             main ="Causal Moderation: Heterogeneous Treatment Effect ",
             ylab="Treatment effect",xlab="index of observation")
        text(zero,low, labels=as.character(zero))
        abline(h=0,lty="dotdash")
        abline(h=ATE, col="blue")
        abline(v=zero, col="grey")                                    
        if(labels=="index"){
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=rownames(pred.data.out.p)
                              ),silent=TRUE)
        }
        else{
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=pred.data.out.p[,labels]
                              ),silent=TRUE)
        }
    }
    if(treat.type=="multiple"){
        pred.data <- object$data
        ## pred.data <- data[,3:ncol(data)]
        ATE <- object$ATE
        pred.data.out.p <- pred.data[order(pred.data$Treatment.effect,
                                           decreasing=FALSE),]   
        xp <- seq(1:nrow(pred.data.out.p))
        zero.1 <- min(which(min(abs(pred.data.out.p$Treatment.effect))==
                            abs(pred.data.out.p$Treatment.effect)))
        zero <- xp[zero.1]
        low <- min(pred.data.out.p$Treatment.effect)
        plot(xp,pred.data.out.p$Treatment.effect,col="red",
             main ="Causal Interaction: Heterogeneous Treatment Effect",
             ylab="Treatment effect",xlab="index of observation")
        text(zero,low, labels=as.character(zero))
        abline(h=0,lty="dotdash")
        abline(h=ATE, col="blue")
        abline(v=zero, col="grey")
        
        if(labels=="index"){
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=rownames(pred.data.out.p)
                              ),silent=TRUE)
        }
        else{
            p <- try(identify(xp,pred.data.out.p$Treatment.effect,
                              labels=pred.data.out.p[,labels]
                              ),silent=TRUE)
        }
    }
    
}

####################################################
summary.FindIt<-function(object,...){
    ## x <- summary.glm(object, dispersion = dispersion,
    ## correlation = correlation, symbolic.cor = symbolic.cor, ...)
    treat.type <- object$treat.type
    type <- object$type
    main <- object$main
    not.zero<-object$coefs.orig!=0
    
    if(type=="binary"){
        coef.table<-cbind(object$names.out,object$coefs.orig/2)[not.zero,]
        coef.print<-cbind(object$names.out,signif(object$coefs.orig/2,3))[not.zero,]
    }
    if(type=="continuous"){
        coef.table<-cbind(object$names.out,object$coefs.orig)[not.zero,]
        coef.print<-cbind(object$names.out,signif(object$coefs.orig,3))[not.zero,]
    }
    
    rownames(coef.print)<-rownames(coef.table)<-object$names.out[not.zero]	
    coef.print[coef.print=="0"]<-"0.000"	
    colnames(coef.print)<-c("Coefficient","Estimate")

    null.gcv<-var(object$y)/(length(object$y)-1)*length(object$y)
    model.gcv<-(object$GCV)*length(object$y)

    if(main){
        model.main  <- object$model.main
    }
    if(main & treat.type=="single"){
        model.int   <- object$model.int
    }
    
    model.treat <- object$model.treat
    
    cat("\nCall:\n")
    cat(" Treatment Model: ")
    print( model.treat)
    
    if(main){
        cat(" Main Model : ")
        print(model.main)
    }
    
    if(main & treat.type=="single"){
        cat(" Interaction Covariates: ")
        print(model.int)
    }
    

    cat(" Treatment type: ")
    print(treat.type)

    cat(" Outcome type: ")
    print(type)
    
    
    cat("\nATE:\n")
    print(object$ATE)

    
    cat("\nCoefficients:\n")

    print(noquote(coef.print[,-1]))  
    
    cat("\n---------")
    
    cat("\nModel Fit Statistics:\n")
    cat(c("GCV:\n"))
    cat(c("  Null: ",round(null.gcv,3)))
    cat(c("  Model: ",round(model.gcv,3),"\n"))
    cat(c("Percent Misclassified:\n"))
    cat(c("   Null: ",round(min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2)))
    cat(c("   Model: ",round(mean(sign(object$y)!=sign(object$fit)),2),"\n"))
    cat(c("   Percent Improvement, vs. NULL: ",round(100-100*mean(sign(object$y)!=sign(object$fit))/min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2),"% \n"))
    cat(c("Percent Outside Margin:\n ",round(mean((object$y^2- object$y*object$fit)<=0 )*100,3),"%, n =",sum((1- object$y*object$fit)<=0 ), "\n")) 	
    out<-list("coefficients"=noquote(coef.print[,-1]),
              "GCV"=c(null.gcv,model.gcv),
              "misclass"=c(min(c(mean(sign(object$y)==-1),mean(object$y==1))),mean(sign(object$y)!=sign(object$fit))))
    invisible(out)  
}




## highlow<-function(object,X0=object$X.c,n.highlow=object$n.highlow, ...){
##     ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)

##     X0<-X0[,apply(X0,2,sd)>0]
##     names.ord<-c("1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th",paste(11:500,"th",sep=""))[1:n.highlow]
##     which.notinter<-unique(c(grep(":",colnames(X0),fixed=T),grep(".2",colnames(X0),fixed=T)))
##     data.highlow<-cbind(X0[object$highlow[,1],],object$highlow[,2]/2)
##     if(length(which.notinter)>0) data.highlow<-cbind(X0[object$highlow[,1],],object$highlow[,2]/2)[,-which.notinter]
##     colnames(data.highlow)[dim(data.highlow)[2]]<-"Estimate"
##     data.high<-data.highlow[1:n.highlow,]
##     data.high.print<-apply(data.high,2,FUN=function(x) if(mean(round(x)==x)==1  ) x else round(x,2)  )
##     data.low<-data.highlow[(n.highlow+1):(2*n.highlow),]
##     data.low.print<-apply(data.low,2,FUN=function(x) if(mean(round(x)==x)==1  ) x else round(x,2)  )

##     rownames(data.high.print)<-rownames(data.low.print)<-names.ord

##     cat("\nHighest Estimated Effect:\n")
##     print(noquote(data.high.print))

##     cat("\nLowest Estimated Effect:\n")
##     print(noquote(data.low.print))

##     out<-list("highlowtable"=data.highlow)
##     invisible(out)

## }



#####################################################
##makeallway

###New version
makeallway<-function(X,threshold=0.999999,deletion=TRUE,
                     make.reference=TRUE,
                     sparse.use=FALSE,
		     nway){
    if(!is.vector(X) & !is.matrix(X)){
        warning("X should be a vector or matrix.")
    }
    X2<-NULL
    ## Make indicator variables for all columns.
    if(is.matrix(X)){
        matrix <- 1
        X <- data.frame(X)
        for(i in 1:dim(X)[2]) {
            X2[[i]]<-cbind(sapply(sort(unique(X[,i])),FUN=function(j) 1*(X[,i]==j)))
            colnames(X2[[i]])<-c(paste(names(X[i]), sort(unique(X[,i])),sep="_"))
        }
        n <- ncol(X)
    }

    ## if(is.vector(X)){
    ##     matrix <- 0
    ##     n <- length(unique(X))
    ##     X2 <- matrix(NA, ncol=n, nrow=length(X))
    ##     name <- c()
    ##     for(z in 1:n){
    ##         for(i in 1:length(X)){
    ##             X2[i,z]<- as.numeric(X[i]==sort(unique(X))[z])
    ##         }
    ##         name[z]<-c(paste(names(X),sort(unique(X))[z],sep="_"))              
    ##     }
    ##     X2 <- as.data.frame(X2)
    ##     colnames(X2) <- name
    ## }

    ## print("makeallway")
    ## print(head(X2[[1]]))

    
    ## NEW data sets
    ##One-way 
    one.way.data <- NULL
    one.way.data <- as.data.frame(X2)
    if(sparse.use==TRUE){
        one.way.data <- as.matrix(one.way.data)
        one.way.data <- Matrix(one.way.data,sparse=TRUE)
    }
    ## print(head(one.way.data))
    
    
    if(matrix == 1){
        ##Two-way
        if(ncol(X)>=2){
            two.way.data <- NULL
            two.way.name.w <- NULL
            for(j in 1:choose(n,2)){
                formula  <- ~ X2[[combn(n,2)[1,j]]]:X2[[combn(n,2)[2,j]]]
                two.way.data[[j]] <- model.matrix(formula)[,-1]
                if(sparse.use==TRUE){
                    two.way.data[[j]] <- Matrix(two.way.data[[j]],sparse=TRUE)
                }
                two.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,2)[1,j]]]),
                                                   colnames(X2[[combn(n,2)[2,j]]]))
                name <- c()
                for(i in 1:nrow(two.way.name.w[[j]])){
                    name[i] <- paste(two.way.name.w[[j]][i,1],
                                     two.way.name.w[[j]][i,2],sep=":")
                }
                colnames(two.way.data[[j]]) <- name
            }
            if(sparse.use==TRUE){
                two.way.data <- do.call(cBind,two.way.data)
            }else{
                two.way.data <- as.data.frame(two.way.data)
            }
            ## two.way.data <- as.matrix(two.way.data)
            ## two.way.data <- Matrix(two.way.data,sparse=TRUE)
        }
        
        ## Three-way
         if(nway>=3){
             three.way.data <- NULL
             three.way.name.w <- NULL
             for(j in 1:choose(n,3)){
                 formula  <- ~ X2[[combn(n,3)[1,j]]]:X2[[combn(n,3)[2,j]]]:
                     X2[[combn(n,3)[3,j]]]
                 three.way.data[[j]] <- model.matrix(formula)[,-1]
                 three.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,3)[1,j]]]),
                                                      colnames(X2[[combn(n,3)[2,j]]]),
                                                      colnames(X2[[combn(n,3)[3,j]]]))
                 name <- c()
                 for(i in 1:nrow(three.way.name.w[[j]])){
                     name[i] <- paste(three.way.name.w[[j]][i,1],
                                      three.way.name.w[[j]][i,2],
                                      three.way.name.w[[j]][i,3],
                                      sep=":")
                 }
                 colnames(three.way.data[[j]]) <- name
             }
             three.way.data <- as.data.frame(three.way.data)
             ## three.way.data <- as.matrix(three.way.data)
             ## three.way.data <- Matrix(three.way.data,sparse=TRUE)
         }
        
        
        ## Four-way
        if(nway>=4){
            four.way.data <- NULL
            four.way.name.w <- NULL
            for(j in 1:choose(n,4)){
                formula  <- ~ X2[[combn(n,4)[1,j]]]:X2[[combn(n,4)[2,j]]]:
                    X2[[combn(n,4)[3,j]]]:X2[[combn(n,4)[4,j]]]
                four.way.data[[j]] <- model.matrix(formula)[,-1]
                four.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,4)[1,j]]]),
                                                    colnames(X2[[combn(n,4)[2,j]]]),
                                                    colnames(X2[[combn(n,4)[3,j]]]),
                                                    colnames(X2[[combn(n,4)[4,j]]])
                                                    )
                name <- c()
                for(i in 1:nrow(four.way.name.w[[j]])){
                    name[i] <- paste(four.way.name.w[[j]][i,1],
                                     four.way.name.w[[j]][i,2],
                                     four.way.name.w[[j]][i,3],
                                     four.way.name.w[[j]][i,4],
                                     sep=":")
                }
                colnames(four.way.data[[j]]) <- name
            }
            four.way.data <- as.data.frame(four.way.data)
            ## four.way.data <- as.matrix(four.way.data)
            ## four.way.data <- Matrix(four.way.data,sparse=TRUE)
        }
        
        
        ## ## Five-way
        ## if(ncol(X)>=5){
        ##     five.way.data <- NULL
        ##     five.way.name.w <- NULL
        ##     for(j in 1:choose(n,5)){
        ##         formula  <- ~ X2[[combn(n,5)[1,j]]]:X2[[combn(n,5)[2,j]]]:
        ##             X2[[combn(n,5)[3,j]]]:X2[[combn(n,5)[4,j]]]:X2[[combn(n,5)[5,j]]]
        
        ##         five.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         five.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,5)[1,j]]]),
        ##                                             colnames(X2[[combn(n,5)[2,j]]]),
        ##                                             colnames(X2[[combn(n,5)[3,j]]]),
        ##                                             colnames(X2[[combn(n,5)[4,j]]]),
        ##                                             colnames(X2[[combn(n,5)[5,j]]])
        ##                                             )
        ##         name <- c()
        ##         for(i in 1:nrow(five.way.name.w[[j]])){
        ##             name[i] <- paste(five.way.name.w[[j]][i,1],
        ##                              five.way.name.w[[j]][i,2],
        ##                              five.way.name.w[[j]][i,3],
        ##                              five.way.name.w[[j]][i,4],
        ##                              five.way.name.w[[j]][i,5],
        ##                              sep=":")
        ##         }
        ##         colnames(five.way.data[[j]]) <- name
        ##     }
        ##     five.way.data <- as.data.frame(five.way.data)
        ##     five.way.data <- as.matrix(five.way.data)
        ##     five.way.data <- Matrix(five.way.data,sparse=TRUE)
        ## }
        

        ## ## Six-way
        ## if(ncol(X)>=6){
        ##     six.way.data <- NULL
        ##     six.way.name.w <- NULL
        ##     for(j in 1:choose(n,6)){
        ##         formula  <- ~ X2[[combn(n,6)[1,j]]]:X2[[combn(n,6)[2,j]]]:
        ##             X2[[combn(n,6)[3,j]]]:X2[[combn(n,6)[4,j]]]:
        ##                 X2[[combn(n,6)[5,j]]]:X2[[combn(n,6)[6,j]]]
        ##         six.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         six.way.name.w[[j]] <- expand.grid( colnames(X2[[combn(n,6)[1,j]]]),
        ##                                            colnames(X2[[combn(n,6)[2,j]]]),
        ##                                            colnames(X2[[combn(n,6)[3,j]]]),
        ##                                            colnames(X2[[combn(n,6)[4,j]]]),
        ##                                            colnames(X2[[combn(n,6)[5,j]]]),
        ##                                            colnames(X2[[combn(n,6)[6,j]]])
        ##                                            )
        ##         name <- c()
        ##         for(i in 1:nrow(six.way.name.w[[j]])){
        ##             name[i] <- paste(six.way.name.w[[j]][i,1],
        ##                              six.way.name.w[[j]][i,2],
        ##                              six.way.name.w[[j]][i,3],
        ##                              six.way.name.w[[j]][i,4],
        ##                              six.way.name.w[[j]][i,5],
        ##                              six.way.name.w[[j]][i,6],
        ##                              sep=":")
        ##         }
        ##         colnames(six.way.data[[j]]) <- name
        ##     }
        ##     six.way.data <- as.data.frame(six.way.data)
        ##     six.way.data <- as.matrix(six.way.data)
        ##     six.way.data <- Matrix(six.way.data,sparse=TRUE)
        ## }
        

        ## ## Seven-way
        ## if(ncol(X)>=7){
        ##     seven.way.data <- NULL
        ##     seven.way.name.w <- NULL
        ##     for(j in 1:choose(n,7)){
        ##         formula  <- ~ X2[[combn(n,7)[1,j]]]:X2[[combn(n,7)[2,j]]]:
        ##             X2[[combn(n,7)[3,j]]]:X2[[combn(n,7)[4,j]]]:
        ##                 X2[[combn(n,7)[5,j]]]:X2[[combn(n,7)[6,j]]]:
        ##                     X2[[combn(n,7)[7,j]]]
        ##         seven.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         seven.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,7)[1,j]]]),
        ##                                              colnames(X2[[combn(n,7)[2,j]]]),
        ##                                              colnames(X2[[combn(n,7)[3,j]]]),
        ##                                              colnames(X2[[combn(n,7)[4,j]]]),
        ##                                              colnames(X2[[combn(n,7)[5,j]]]),
        ##                                              colnames(X2[[combn(n,7)[6,j]]]),
        ##                                              colnames(X2[[combn(n,7)[7,j]]])
        ##                                              )
        ##         name <- c()
        ##         for(i in 1:nrow(seven.way.name.w[[j]])){
        ##             name[i] <- paste(seven.way.name.w[[j]][i,1],
        ##                              seven.way.name.w[[j]][i,2],
        ##                              seven.way.name.w[[j]][i,3],
        ##                              seven.way.name.w[[j]][i,4],
        ##                              seven.way.name.w[[j]][i,5],
        ##                              seven.way.name.w[[j]][i,6],
        ##                              seven.way.name.w[[j]][i,7],
        ##                              sep=":")
        ##         }
        ##         colnames(seven.way.data[[j]]) <- name
        ##     }
        ##     seven.way.data <- as.data.frame(seven.way.data)
        ##     seven.way.data <- as.matrix(seven.way.data)
        ##     seven.way.data <- Matrix(seven.way.data,sparse=TRUE)
        ## }
        

        ## ## Eight-way
        ## if(ncol(X)>=8){
        ##     eight.way.data <- NULL
        ##     eight.way.name.w <- NULL
        ##     for(j in 1:choose(n,8)){
        ##         formula  <- ~ X2[[combn(n,8)[1,j]]]:X2[[combn(n,8)[2,j]]]:
        ##             X2[[combn(n,8)[3,j]]]:X2[[combn(n,8)[4,j]]]:
        ##                 X2[[combn(n,8)[5,j]]]:X2[[combn(n,8)[6,j]]]:
        ##                     X2[[combn(n,8)[7,j]]]:X2[[combn(n,8)[8,j]]]
        ##         eight.way.data[[j]] <- model.matrix(formula)[,-1]
        ##         eight.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,8)[1,j]]]),
        ##                                              colnames(X2[[combn(n,8)[2,j]]]),
        ##                                              colnames(X2[[combn(n,8)[3,j]]]),
        ##                                              colnames(X2[[combn(n,8)[4,j]]]),
        ##                                              colnames(X2[[combn(n,8)[5,j]]]),
        ##                                              colnames(X2[[combn(n,8)[6,j]]]),
        ##                                              colnames(X2[[combn(n,8)[7,j]]]),
        ##                                              colnames(X2[[combn(n,8)[8,j]]])
        ##                                              )
        ##         name <- c()
        ##         for(i in 1:nrow(eight.way.name.w[[j]])){
        ##             name[i] <- paste(eight.way.name.w[[j]][i,1],
        ##                              eight.way.name.w[[j]][i,2],
        ##                              eight.way.name.w[[j]][i,3],
        ##                              eight.way.name.w[[j]][i,4],
        ##                              eight.way.name.w[[j]][i,5],
        ##                              eight.way.name.w[[j]][i,6],
        ##                              eight.way.name.w[[j]][i,7],
        ##                              eight.way.name.w[[j]][i,8],
        ##                              sep=":")
        ##         }
        ##         colnames(eight.way.data[[j]]) <- name
        ##     }
        ##     eight.way.data <- as.data.frame(eight.way.data)
        ##     eight.way.data <- as.matrix(eight.way.data)
        ##     eight.way.data <- Matrix(eight.way.data,sparse=TRUE)
        ## }
    }
    
    ## if(is.vector(X)){FinalData <- as.data.frame(one.way.data)}
    if(ncol(X)==1){FinalData <- one.way.data}
    if(matrix==1){
        if(sparse.use==TRUE){
            if(ncol(X)>=2){FinalData <- cBind(one.way.data,two.way.data)}
        }else{
            if(ncol(X)>=2){FinalData <- cbind(one.way.data,two.way.data)}
        }
        
        if(nway==3){FinalData <- cbind(one.way.data,two.way.data,three.way.data)}
        if(nway==4){FinalData <- cbind(one.way.data,two.way.data,three.way.data,
                   four.way.data)}
        ## if(ncol(X)==5){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data)}
        ## if(ncol(X)==6){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data)}
        ## if(ncol(X)==7){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data,
        ##            seven.way.data)}
        ## if(ncol(X)==8){FinalData <- cBind(one.way.data,two.way.data,three.way.data,
        ##            four.way.data,five.way.data,six.way.data,
        ##            seven.way.data,eight.way.data)}
    }
    
    ## FinalData <- as.data.frame(FinalData)
    ## print(head(FinalData))
    
    if(deletion==TRUE){
        ## Discard the columns with no variation or the same variation
        ## with other columns,
        ## while keeping the column names
        No.variation <- colnames(FinalData[,apply(FinalData,2,sd)==0])
        T2 <- FinalData[,apply(FinalData,2,sd)>0]
        T2.sparse <- T2
        
        T2 <- as.matrix(T2)
        ## T2 <- t(T2)
        coldrop      <- sapply(1:ncol(T2),
                               FUN=function(x)
                               which(cor(T2[,x] ,T2)^2 > threshold))
        
        col.names.w  <- sapply(coldrop, FUN=function(x) colnames(T2)[x])
        col.names    <- list(Perfect.Correlated.Columns=col.names.w,
                             no.variation.columns=No.variation)
        
        ## it could be three columns are the same. 
        coldrop.cols <- sapply(coldrop,FUN=function(x) x[1])
        Keep.cols <- unique(coldrop.cols)
        ## if there are more than two, the first one is gotten.
        ## so the second or the third one is the discarded one.
        ## if there are only one, then it is the same.
        
        ## colnames(T2)[coldrop[[Keep.cols]][-1]]
        ## This gives the corresponding dropped coefficients.
        T3 <- T2.sparse[,Keep.cols]
        
        if(make.reference==TRUE){
            Corresponding <- list()
            ## Keep Discarded Data.
            Discarded.cols <- seq(1:ncol(T2))[-Keep.cols]
            
            for(i in Keep.cols){
                if(length(coldrop[[i]][-1])==0){
                    Corresponding[[i]] <- "No Match"
                }
                if(length(coldrop[[i]][-1])!=0){
                    Corresponding[[i]] <- colnames(T2)[coldrop[[i]][-1]]
                }
            }
            for(i in Discarded.cols){
                Corresponding[[i]] <- "Discarded"
            }
            
            C <- max(sapply(Corresponding, FUN=function(x) length(x)))
            
            CorrespondingM <- matrix(NA,nrow=ncol(T2), ncol=C)
            Variation <- colnames(T2)
            rownames(CorrespondingM) <- colnames(T2)
            for(j in 1:C){
                for(i in 1:ncol(T2)){
                    if(length(Corresponding[[i]])>=j){
                        CorrespondingM[i,j] <- Corresponding[[i]][j]
                    }
                }
            }
            
            Reference <- matrix(NA,nrow=ncol(FinalData), ncol=C)
            rownames(Reference) <- colnames(FinalData)
            
            for(i in 1:nrow(Reference)){
                if(rownames(Reference)[i] %in% No.variation){
                    Reference[i,1] <- "No Variation"
                }
                if(rownames(Reference)[i] %in% Variation){
                    Reference[i,]  <- CorrespondingM[rownames(CorrespondingM)==
                                                     rownames(Reference)[i]]
                }
            }
            
            names <- c("Matched Variables", rep("Other matched variables",C-1))
            colnames(Reference) <- names
            Reference <- as.data.frame(Reference)
        }
        else{
            Reference <- NULL}
    }
    if(deletion==FALSE){
        T3 <- FinalData
        Reference <- NULL
    }
    
    return(list(FinalData=T3, reference=Reference))
}



#####################################################
##maketwoway

maketwoway<-function(X,wts=1,center=TRUE,deletion=TRUE,threshold=0.99999,frame.meanPre,predict=FALSE){
    ## We have the basic code.
    ## I need to fix the column number and scaling.
    X <- as.data.frame(X)
    c <- c(lapply(X,class))
    num <- is.element(c, c("integer","numeric"))
    X.num <- X[,num==TRUE]
    X.categ <- X[,num==FALSE]
    formula <- ~ .
    ##print(head(X.categ))
    ##print(head(X.num))
    X.num <- as.data.frame(X.num)
    X.categ <- as.data.frame(X.categ)

    ## The deletion for constant categorical variables.
    ## Deletion tweak is with predict function.
    if(sum(num)<ncol(X)){
        frame.categ <- model.frame(formula,data=X.categ)
        ## The first deletion for categorical.
        if(ncol(frame.categ)==1){
            constant <- length(unique(frame.categ)) < 2
            if(constant==TRUE){
                delete1 <- colnames(frame.categ)
                Cat <- FALSE
                if(deletion==FALSE){
                    frame.categ.c <- frame.categ
                }
            }
            else{
                delete1 <- NULL
                Cat <- TRUE
                if(deletion==FALSE){
                    frame.categ.c <- frame.categ
                }
                
            }
        }
        if(ncol(frame.categ)>=2){
            constant <- apply(frame.categ,2,
                              FUN=function(x)
                              length(unique(x)) < 2)
            if(sum(constant)>0){
                delete1 <- colnames(frame.categ)[constant]
            }
            else{delete1 <- NULL}
            ## deletion finish for categorical.
            if(sum(constant)< ncol(frame.categ)){
                if(deletion==TRUE){
                    frame.categ <- frame.categ[,constant==FALSE]
                    Cat <- TRUE
                }
                if(deletion==FALSE){
                    frame.categ.n <- frame.categ[,constant==FALSE]
                    frame.categ.c <- frame.categ[,constant==TRUE]
                }
            }
            if(sum(constant) == ncol(frame.categ)){
                Cat <- FALSE
                frame.categ.c <- frame.categ
            }
        }

        Cat <- TRUE
        
        if(Cat==TRUE){
            if(deletion==TRUE){
                mat.categ   <- model.matrix(formula,
                                            data=frame.categ)[,-1]
            }
            if(deletion==FALSE){
                mat.categ.c <- frame.categ.c
                mat.categ.n <- model.matrix(formula,
                                            data=frame.categ.n)[,-1]
                mat.categ <- cbind(mat.categ.c,mat.categ.n)
            }
            
            ##print("mat.categ")
            ##print(head(mat.categ))
            
            mat.categ.s <- matrix(NA,ncol=ncol(mat.categ),nrow=nrow(mat.categ))
            colnames(mat.categ.s) <- colnames(mat.categ)

            if(center==TRUE){
                if(predict==FALSE){
                    mat.categ.s <-
                        apply(mat.categ,
                              MARGIN=2,
                              FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                    mean.categ<-
                        apply(mat.categ,
                              MARGIN=2,
                              FUN=function(x) mean(wts^.5*x))
                }else{                    
                    mean.categ <- frame.meanPre[1:ncol(mat.categ)]
                    mat.categ.s1 <- matrix(NA,ncol=ncol(mean.categ),nrow=nrow(mean.categ))
                    for(j in 1:ncol(mat.categ)){
                        mat.categ.s1[,j] <- mat.categ[,j] - mean.categ[j]
                    }
                    mat.categ.s <-
                        apply(mat.categ.s1,
                              MARGIN=2,
                              FUN=function(x) x/sd(x))
                }
                scale.categ<-
                    apply(mat.categ,
                          MARGIN=2,
                          FUN=function(x) sd(x))                       
                frame.categ.s <- as.data.frame(mat.categ.s)
                colnames(frame.categ.s) <- colnames(mat.categ)
            }
            else{
                mat.categ.s <- mat.categ
                frame.categ.s <- frame.categ
            }
        }
    }else{Cat <- FALSE}
    
    ##print(head(frame.categ.s))
    
    if(sum(num)>0){
        Num <- TRUE
        frame.num <- model.frame(formula,data=X.num)
        mat.num   <- model.matrix(formula,frame.num)[,-1]
        
        mat.num.s <- matrix(NA, ncol=ncol(mat.num),nrow=nrow(mat.num))
        colnames(mat.num.s) <- colnames(mat.num)
        
        ## Scaling
        ## Binary and Categorical Variable are also scaled.
        if(center==TRUE){
            if(predict==FALSE){
                mat.num.s <-
                    apply(mat.num,
                          MARGIN=2,
                          FUN=function(x) (x-mean(wts^.5*x))/sd(x))
                mean.num <-
                    apply(mat.num,
                          MARGIN=2,
                          FUN=function(x) mean(wts^.5*x))
            }else{
                if(Cat){
                    mean.num   <- frame.meanPre[(ncol(mat.categ)+1):length(frame.meanPre)]
                    mat.num.s1 <- matrix(NA,ncol=ncol(mat.num),nrow=nrow(mat.num))
                    for(j in 1:ncol(mat.num)){
                        mat.num.s1[,j] <- mat.num[,j] - mean.num[j]
                    }
                    colnames(mat.num.s1) <- colnames(mat.num)
                }else{                    
                    mean.num   <- frame.meanPre
                    mat.num.s1 <- matrix(NA,ncol=ncol(mat.num),nrow=nrow(mat.num))
                    for(j in 1:ncol(mat.num)){
                        mat.num.s1[,j] <- mat.num[,j] - mean.num[j]
                    }
                    colnames(mat.num.s1) <- colnames(mat.num)
                }
                    
                mat.num.s <-
                    apply(mat.num.s1,
                          MARGIN=2,
                          FUN=function(x) x/sd(x))                
            }
            scale.num<-
                apply(mat.num,
                      MARGIN=2,
                      FUN=function(x) sd(x))                       
            frame.num.s <- as.data.frame(mat.num.s)
            colnames(frame.num.s) <- colnames(frame.num)
        }
        else{
            mat.num.s <- mat.num
            frame.num.s <- frame.num
        }
    }
    if(sum(num)==0){
        Num <- FALSE
    }
    

    if(center==TRUE){
        if(Cat & Num){
            frame.name<-merge(frame.categ,frame.num,sort=FALSE,
                              by=c("row.names",
                                  intersect(names(frame.categ),
                                            names(frame.num))))[,-1]
            frame <-merge(frame.categ.s,frame.num.s,sort=FALSE,
                          by=c("row.names",
                              intersect(names(frame.categ.s),
                                        names(frame.num.s))))[,-1]
            frame.scale <- c(scale.categ,scale.num)
            frame.mean <- c(mean.categ,mean.num)
        }
        if(Cat==TRUE & Num==FALSE){
            frame.name <- frame.categ
            frame <- frame.categ.s
            frame.scale <- scale.categ
            frame.mean  <- mean.categ
        }
        if(Cat==FALSE & Num==TRUE){
            frame.name <- frame.num
            frame <- frame.num.s
            frame.scale <- scale.num
            frame.mean <- mean.num
        }
        frame.s <- frame
    }
    else{
        frame   <- model.frame(formula, data=X)
    }

    ##print(head(frame.s))
    ## interaction term
    formula2 <- ~ .*.
    ## int.matrix.s <- model.matrix(formula2, frame.s)[,-1]
    int.matrix <- model.matrix(formula2, frame)[,-1]
    int.matrix.name <- model.matrix(formula2, frame.name)[,-1]
    


    frame.scale.mat <- matrix(frame.scale,
                              ncol=length(frame.scale),
                              nrow=3,
                              byrow=TRUE)
    frame.scale.mat <- as.data.frame(frame.scale.mat)
    if(center==TRUE){
        scale.int.mat   <- model.matrix(formula2, frame.scale.mat)[,-1]
        scale.int <- scale.int.mat[1,]
    }else{
        scale.int <- rep(1,ncol(int.matrix))
    }
    names(scale.int) <- colnames(int.matrix)

    ## Reduce the column number for INT matrix.
    ## Categorical variables are treated properly.
    int.in <- is.element(colnames(int.matrix),colnames(int.matrix.name))
    int.matrix <- int.matrix[,int.in]
    scale.int  <- scale.int[int.in]
    ## no variation columns. Try to remove the no-variation in the original scale.
    no.variation.int <- apply(int.matrix.name,2, FUN=function(x) sd(x)==0)
    delete.novar.int <- colnames(int.matrix)[no.variation.int]
    int.matrix <- int.matrix[,no.variation.int==FALSE]
    scale.int <- scale.int[no.variation.int==FALSE]
    
    ## Squared Matrix
    if(Num==TRUE){
        nonbin <- apply(mat.num.s,2,FUN=function(x) length(unique(x)) > 2)
        ##print(nonbin)
        if(sum(nonbin)>0){
            mat.sq   <- matrix(NA, ncol=sum(nonbin),nrow=nrow(mat.num.s))
            mat.sq <- mat.num.s[,nonbin]^2
            scale.sq <- scale.num[nonbin]^2
            mat.sq <- as.data.frame(mat.sq)
            if(sum(nonbin)>1){
                colnames(mat.sq) <- c(paste(colnames(mat.num.s)[nonbin],
                                            ".2",sep=""))
                names(scale.sq) <- colnames(mat.sq)
            }
            if(sum(nonbin)==1){
                colnames(mat.sq) <- c(paste(colnames(mat.num.s)[nonbin==TRUE],
                                            ".2",sep=""))
                names(scale.sq) <- colnames(mat.sq)
            }
            ##print(dim(mat.sq))
            ##print(dim(int.matrix))
            Xtwo <- cbind(int.matrix, mat.sq)
            scale.out <- c(scale.int, scale.sq)
        }else{
            Xtwo <- int.matrix
            scale.out <- scale.int 
        }
    }else{
        Xtwo <- int.matrix
        scale.out <- scale.int
    }

    ## remove no variation columns.(deletion 2) for Squared Matrix.
    no.variation <- apply(Xtwo,2, FUN=function(x) sd(x)==0)
    delete.novar <- colnames(Xtwo)[no.variation]
    Xtwo <- Xtwo[,no.variation==FALSE]
    scale.out <- scale.out[no.variation==FALSE]

    ## remove correlation is almost 1. (deletion 3).
    T <- Xtwo
    multcor <-  sapply(1:ncol(T),
                       FUN=function(x)
                       which(cor(T[,x] ,T)^2 > threshold))
    multcor.keep <- unique(sapply(multcor,
                                  FUN=function(x) x[1]))
    delete.cor   <- colnames(Xtwo)[-multcor.keep]
    Xtwo <- Xtwo[,multcor.keep]
    scale.out <- scale.out[multcor.keep]

    ## Deleted <- c(delete1,delete.novar, delete.cor)
    Deleted <- c(delete.novar, delete.cor)
    out <- list(X=Xtwo,scale.X=scale.out,delete=Deleted,frame.mean=frame.mean)
    invisible(out)
}







#####################################################
SVM.func<-function(y,X.c,X.t,lambda.c,lambda.t,wts=1,fit.glmnet,type,main){

    ### From the construction, this X.c doesnt include the intercept. 

    n<-length(y)
    if(main){
        X<-X2<-cbind(X.c,X.t)
    }
    else{
        X <- X2 <- X.t 
    }

    ##Originally written to handle a string of lambda.c and lambda.t
    ##loss.run creates a matrix of losses, across values of lambda.
    ## loss.run<-matrix(0,nrow=length(lambda.c),ncol=length(lambda.t))
    ## min.params<-c(lambda.c[1],lambda.t[1])
    ## for(i in 1:length(lambda.c)){
    ## for(j in 1:length(lambda.t)){

    ##Scale the X matrix by lambdas.
    X3<-X2
    if(main){
        if(is.matrix(X.c)){
            X3[,1:dim(X.c)[2]]<-1/exp(lambda.c)*X2[,1:dim(X.c)[2]]
            X3[,-c(1:dim(X.c)[2])]<-1/exp(lambda.t)*X2[,-c(1:dim(X.c)[2])]
        }
        if(is.vector(X.c)){
            X3[,1]<-1/exp(lambda.c)*X2[,1]
            X3[,-1]<-1/exp(lambda.t)*X2[,-1]
        }
    }

    if(main==FALSE){
        X3 <- 1/exp(lambda.t)*X2
    }
    ##X3.c <- 1/exp(lambda.c[i])*X.c
    ##X3.t <- 1/exp(lambda.t[j])*X.t
    ##X3   <- cbind(X.c,X.t)
    

    ##Declare original values.
    X.new<-X3
    which.use<-rep(T,n)
    beta.new<-beta.curr<-rep(0,dim(X2)[2]+1)
    ## this +1 is for the intercept.

    ## The loop to fit the LASSO.
    ## Center on the margin, fit the lasso, update coefficients,
    ##then drop all outside the margin (y*yhat>1)
    if(type=="binary"){
        for(beta.run in 1:100){
            X.new.2<-apply(X.new[which.use,],MARGIN=2,FUN=function(x) x-mean(x))
            X.new.2<-apply(X.new.2,MARGIN=2,
                           FUN=function(x)
                           if(length(unique(x))==1) rnorm(sum(which.use))
                           else x)
            if(fit.glmnet==TRUE){
                glmnet.1<-glmnet(X.new.2,y[which.use]-mean(y[which.use]),
                                 family="gaussian",
                                 lambda=c(5,4,3,2,seq(1.5,1,-.1)),
                                 standardize=F)
                beta.new[-1]<-as.vector(glmnet.1$beta[,10])
            }else{
                ##This was the old L1 optimizer.
                lasso1<-lars(X.new.2,(y[which.use]-mean(y[which.use])),
                             max.steps=dim(X)[2]+2,
                             normalize=F,
                             type="lasso",
                             eps=0)
                beta.new[-1]<-as.vector(predict(lasso1,s=1,type="coef",
                                                mode="lambda")$coef)
                if(log(sum(which.use))/2*sum(beta.new[-1]!=0)>.9*sum(which.use)){
                    beta.new[-1]<-
                        as.vector(predict(lasso1,
                                          s=min(c(floor(.9*sum(which.use)*2/
                                              log(n)),dim(X.new)[2]*.8)),
                                          type="coef",mode="step")$coef)
                }
            }
            beta.curr<-.5*beta.new+.5*beta.curr
            
            beta.new[1]<-mean(y[which.use])-mean(X.new[which.use,]%*%beta.new[-1])
            ## This is for the intercept. This corrects the demean tricks.
            run.diff<-(mean((beta.new[-1]-beta.curr[-1])^2)/sum(beta.new[-1]^2+1e-10))
            if(run.diff<1e-6) break
            which.use<-(y*cbind(1,X.new)%*%beta.new)<=y^2
        }
        
        ##Find fitted values.
        fits<-pmin(abs(cbind(1,X.new)%*%beta.new),y^2)*sign((cbind(1,X.new)%*%beta.new))
        fits2<-cbind(1,X.new)%*%beta.new
        fits2[sign(fits2)==sign(y)]<-pmin(abs(fits2[sign(fits2)==sign(y)]),
                      y[sign(fits2)==sign(y)]^2)*sign(fits2[sign(fits2)==sign(y)])
        
    }
    if(type=="continuous"){
        X.new.2<-apply(X.new,MARGIN=2,FUN=function(x) x-mean(x))
        X.new.2<-apply(X.new.2,MARGIN=2,
                       FUN=function(x) if(length(unique(x))==1) rnorm(sum(which.use)) else x)
        if(fit.glmnet==TRUE){
            glmnet.1<-glmnet(X.new.2,y-mean(y),
                             family="gaussian",
                             lambda=c(5,4,3,2,seq(1.5,1,-.1)),
                             standardize=F)
            beta.new[-1]<-as.vector(glmnet.1$beta[,10])
        }else{
            ##This was the old L1 optimizer.
            lasso1<-lars(X.new.2,(y-mean(y)),
                         max.steps=dim(X)[2]+2,normalize=F,type="lasso",eps=0)
            beta.new[-1]<-as.vector(predict(lasso1,s=1,type="coef",mode="lambda")$coef)
            if(log(length(y))/2*sum(beta.new[-1]!=0)>.9*length(y)) 
                beta.new[-1]<-
                    as.vector(predict(
                        lasso1,
                        s=min(c(floor(.9*length(y)*2/log(n)),dim(X.new)[2]*.8)),
                        type="coef",mode="step")$coef)
        }
        beta.curr<-.5*beta.new+.5*beta.curr
        
        beta.new[1]<-mean(y)-mean(X.new%*%beta.new[-1])
        ## This is for the intercept
        ##run.diff<-(mean((beta.new[-1]-beta.curr[-1])^2)/sum(beta.new[-1]^2+1e-10))
        ##if(run.diff<1e-6) break
        ##which.use<-(y*cbind(1,X.new)%*%beta.new)<=y^2
        ##}
        ## We dont need loop. Because this is not SVM.
        
        ##Find fitted values.
        fits2 <- cbind(1,X.new)%*%beta.new
        fits  <- fits2
        ##fits<-pmin(abs(cbind(1,X.new)%*%beta.new),y^2)*sign((cbind(1,X.new)%*%beta.new))
        ##fits2<-cbind(1,X.new)%*%beta.new
        ##fits2[sign(fits2)==sign(y)]<-pmin(abs(fits2[sign(fits2)==sign(y)]),
        ##              y[sign(fits2)==sign(y)]^2)*sign(fits2[sign(fits2)==sign(y)])
    }
    
    ##Calculate degrees of freedom
    edf<-  1+sum(beta.new[-1]!=0)
    ##*log(n)/2#sum(diag(hat.mat))#+1-sum(beta.new[-c(1:(dim(X.c)[2]+1))]!=0)
    ##Left in for estimates of standard deviation.  	
    stdev<-NA
    
    ##GCV statistic
    if(type=="binary"){
        loss.run<-sum((y^2-y*fits2)^2)/(sum(which.use)-edf)^2*(mean(which.use))^2*n
    }
    if(type=="continuous"){
        loss.run<-sum((y-fits2)^2)/(n*((1-edf/n)^2))
    }
    
    ## Gather minimum loss function.
    ## if(i*j>1) if(loss.run[i,j]==min(loss.run[loss.run!=0]))
    ## min.params<-c(lambda.c[i],lambda.t[j])

    ##     }
    ## }
    
    ##Scale betas back.
    if(main){
        if(is.matrix(X.c)){
            beta.new[c(1:dim(X.c)[2]+1)]<-beta.new[c(1:dim(X.c)[2]+1)]/exp(lambda.c)
            ##for the main covariates and intercepts
            beta.new[-c(1:(dim(X.c)[2]+1))]<-beta.new[-c(1:(dim(X.c)[2]+1))]/exp(lambda.t)
            ##for the interaction covariates
        }
        if(is.vector(X.c)){
            beta.new[c(1:2)]<-beta.new[c(1:2)]/exp(lambda.c)
            ##for the main covariates and intercepts
            beta.new[-c(1:2)]<-beta.new[-c(1:2)]/exp(lambda.t)
            ##for the interaction covariates
        }
    }
    if(main==FALSE){
        beta.new <- beta.new/exp(lambda.t)
    }
    
    beta.new<-as.vector(beta.new)
    X<-as.matrix(X)

    ##Calculate intercept
    if(type=="binary"){
        beta.new[1]<-mean(y[which.use])-mean((X%*%beta.new[-1])[which.use])
    }
    if(type=="continuous"){
        beta.new[1]<-mean(y)-mean((X%*%beta.new[-1]))
    }

    output<-list(##"lambdas"=min.params,
                 "beta"=beta.new,"fits"=fits,
                 "loss"=log(loss.run),##"marg"=mean(which.use),
                 "edf"=edf,"sd"=stdev)
    invisible(output)
}


#####################################################

#####################################################


search.lambda<-function(y=y,X.c=X.c,X.t=X.t,fit.glmnet,main,type){
    lambda.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                      lambda[1],lambda[2],
                                                      fit.glmnet=fit.glmnet,
                                                      main=main,
                                                      type=type)$loss)
    lambda.c.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                        lambda,lambda.t,
                                                        fit.glmnet=fit.glmnet,
                                                        main=main,
                                                        type=type)$loss)
    lambda.t.find<-function(lambda) as.numeric(SVM.func(y,X.c,X.t,
                                                        lambda.c,lambda,
                                                        fit.glmnet=fit.glmnet,
                                                        main=main,
                                                        type=type)$loss)
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

    ##lambda.c<- -15
    ##lambda.t<-5
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
search.lambda.nocov<-function(y=y,X.t=X.t,fit.glmnet,main,type){
    lambda.t.find<-function(lambda) as.numeric(SVM.func(y,X.c=NULL,X.t,
                                                        lambda.c=1,lambda,
                                                        fit.glmnet=fit.glmnet,
                                                        main=FALSE,
                                                        type=type)$loss)
    ## This lambda.c will not be used. 
    ## lambda.c.old<-lambda.t.old<-999
    lambda.t.old<-999

##############################
    ## The Initial Try
############################## 
    ## lambda.c.seek<-seq(-15,10,1)
    lambda.t.seek<-seq(-15,10,1)
    ## lambda.t<-25
    ## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
    ## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
    lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
    lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
    ## print(c(lambda.c,lambda.t))
    ## print(range(lambda.c.seek))
    print(range(lambda.t.seek))

    ##lambda.c<- -15
    ##lambda.t<-5
    ## lambda.c.old<-lambda.c
    lambda.t.old<-lambda.t

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-5,lambda.c+5,1)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-5,lambda.t+5,1)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-2.5,lambda.c+2.5,.5)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-2.5,lambda.t+2.5,.5)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
        ## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-1.5,lambda.c+1.5,.25)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-1.5,lambda.t+1.5,.25)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.5,lambda.c+.5,.1)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.5,lambda.t+.5,.1)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }



    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.25,lambda.c+.25,.05)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.25,lambda.t+.25,.05)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.15,lambda.c+.15,.025)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.15,lambda.t+.15,.025)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.1,lambda.c+.1,.01)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.1,lambda.t+.1,.01)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }

    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.05,lambda.c+.05,.005)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.05,lambda.t+.05,.005)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.025,lambda.c+.025,.0025)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.025,lambda.t+.025,.0025)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.01,lambda.c+.01,.001)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.01,lambda.t+.01,.001)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }


    for(lambda.loop in 1:5){
	## lambda.c.seek<-seq(lambda.c-.0025,lambda.c+.0025,.0005)
	## lambda.c.seek<-unique(pmax(lambda.c.seek,-15))
	lambda.t.seek<-seq(lambda.t-.0025,lambda.t+.0025,.0005)
	## lambda.c.out<-sapply(lambda.c.seek,lambda.c.find)
	## lambda.c<-lambda.c.seek[min(which(lambda.c.out==min(lambda.c.out)))]
	## lambda.c<-max(-15,lambda.c)
	lambda.t.out<-sapply(lambda.t.seek,lambda.t.find)
	lambda.t<-lambda.t.seek[min(which(lambda.t.out==min(lambda.t.out)))]
	## print(c(lambda.c,lambda.t))
	## print(range(lambda.c.seek))
	print(range(lambda.t.seek))
	## if(lambda.c.old==lambda.c&lambda.t.old==lambda.t) break
        if(lambda.t.old==lambda.t) break
	## lambda.c.old<-lambda.c
	lambda.t.old<-lambda.t
    }
    
    ## optim.lambda<-c(lambda.c,lambda.t)
    optim.lambda<- lambda.t


    invisible(optim.lambda)
    
}


## Effect<- function(object,unif.data,column,Restriction,sort=TRUE,base=1,dist="unique",
##                   compare=FALSE,order=2,measure="TIE"){

##     if(compare==TRUE){
##         if(missing(order)){
##             warning("Need to specify the order")
##         }
##     }
    
##     ## if(all(class(object)!="FindIt")){
##     ##     warning("Effect function can take only class-FindIt.")
##     ## }
    


##     coefs <- object$coefs.orig
##     names(coefs) <- c("Intercept",object$name.t)
##     outcome.type <- object$type

    
##     if(dist=="unique"){
##         data <- object$Treatment.version
##         data <- unique(data,MARGIN=1)
##     }
##     if(dist=="sample"){
##         data <- object$Treatment.version
##     }

##     if(dist=="uniform"){
##         data <- unif.data
##         print("Using Data from Uniform Distribution or Researcher-Specify Data.")
##     }else{
##         print("Using the same Data as the one used to fit the model.")
##     }

##     if(!missing(column)){
##         ind <- c()
##         for(i in 1:length(column)){
##             ind[i] <- which(colnames(data)==column[i])
##         }
##         column <- column[order(ind,decreasing=FALSE)]
##     }

    
    
##     ##Marginal Effect
##     Marginal <- function(data,base,Restriction){
##     Effect.list <- list()
##     Name.list <- list()
##     range.mar <- c()
##     for(i in 1:(ncol(data)-1)){
##         columnM <- colnames(data)[i+1]
##         A <- tapply(data[,1],data[,columnM],mean,simplify=FALSE)
##         Effect1 <- unlist(A)
##         if(base=="min"){
##             base <- which(Effect1==min(Effect1))
##         }
##         TEffect <- Effect1 - Effect1[base]

##         if(!missing(Restriction)){
##             for(res in 1:length(Restriction)){
##                 restriction <- Restriction[[res]]
##                 if(is.element(columnM,restriction$col.restricted)){
##                     ## base.ind<- which(is.element(names(Effect1),restriction$fac.restricted))
##                     ## base <-min(seq(1:length(Effect1))[-base.ind])
                    
##                     ## if(base=="min"){
##                     ##     base <- which(Effect1==min(Effect1))
##                     ## }
                    
##                     TEffect <- Effect1 - Effect1[base]            
##                     ind.AC <- which(is.element(data[,restriction$col.restricting],
##                                                restriction$fac.restricting))
##                     dataAC <- data[ind.AC,]
##                     AC <- tapply(dataAC[,1],dataAC[,columnM],mean,simplify=FALSE)
##                     Effect1AC <- unlist(AC)
##                     TEffectAC <- Effect1AC - Effect1AC[base]
##                     for(z in 1:length(Effect1)){
##                         if(is.element(names(Effect1)[z],restriction$fac.restricted)){
##                             ## print(names(Effect1)[z])
##                             TEffect[z] <- TEffectAC[z]
##                         }
##                     }            
##                 }
##             }
##         }

##         range.mar[i] <- max(TEffect) - min(TEffect)
##         Effect.list[[i]] <- TEffect
##         Name.list[[i]] <- paste(colnames(data)[(i+1)],names(TEffect),sep="_")
##     }
##     names(range.mar) <- colnames(data)[2:ncol(data)]
##     name <- unlist(Name.list)
##     Treatment.Effect1 <- unlist(Effect.list)
##     Treatment.Effect1 <- as.data.frame(Treatment.Effect1)
##     colnames(Treatment.Effect1) <- "Marginal.Effect"
##     rownames(Treatment.Effect1) <- name
##     ## print(rownames(Treatment.Effect1))
##     ## print(names(Treatment.Effect1))
##     output <- list("Treatment.Effect1"=Treatment.Effect1,"range"=range.mar)
##     return(output)
## }

##     if(compare==FALSE & missing(column)){
##         Treatment.Effect1 <- Marginal(data=data,
##                                       base=base,
##                                       Restriction=Restriction)$Treatment.Effect1
##         Treatment.Effect <- Treatment.Effect1
##     }
##     if(compare==TRUE & order==1){
##         Marginal.range <- Marginal(data=data,
##                                       base=base,
##                                       Restriction=Restriction)$range
##         print("Range of Marginal Effects")
##         return(Marginal.range)
##     }
    
    
    
##     comb.prem <- function(data,column,Restriction,order){
##         if(!missing(Restriction)){
##             ##Restriction
##             ## Col.Restricted
##             Col.Restricted <- c()    
##             for(i in 1:length(Res)){
##                 Col.Restricted <- c(Col.Restricted,Res[[i]]$col.restricted)
##             }
##             ## Col.Restricting
##             Col.Restricting <- c()    
##             for(i in 1:length(Res)){
##                 Col.Restricting <- c(Col.Restricting,Res[[i]]$col.restricting)
##             }

##             ##Pattern
##             ## Now we only support until order 2.
##             Res.type <- "Normal"
##             ## One 
##             if(sum(is.element(column,Col.Restricted))==1 &
##                sum(is.element(column,Col.Restricting))==0){Res.type <- "Restricted1"}
##             if(sum(is.element(column,Col.Restricted))==0 &
##                sum(is.element(column,Col.Restricting))==1){Res.type <- "Restrict1"}
##             ## two
##             if(sum(is.element(column,Col.Restricted))==2){Res.type <- "Restricted2"}
##             if(sum(is.element(column,Col.Restricting))==2){Res.type <- "Restrict2"}
##             ## Mix
##             if(sum(is.element(column,Col.Restricted))==1 &
##                sum(is.element(column,Col.Restricting))==1){
##                 col.restricted.simple <- column[is.element(column,Col.Restricted)==TRUE]
##                 col.restricting.simple <- column[is.element(column,Col.Restricting)==TRUE]
##                 ## Simple Case
##                 for(resC in 1:length(Restriction)){
##                     restriction.check <- Restriction[[resC]]
##                     if((restriction.check$col.restricted==col.restricted.simple)&
##                        (restriction.check$col.restricting==col.restricting.simple)){
##                         Res.type <- "Restriction.mix.simple"
##                     }else{
##                         Res.type <- "Restriction.mix.hard"
##                     }
##                 }
                
##             }
            
##             ## Make the data for Full common Support
##             if(Res.type=="Restrict1"| Res.type=="Restriction.mix.simple"){
##                 column.restrict1 <- column[is.element(column,Col.Restricting)==TRUE]
##                 for(z in 1:length(Restriction)){
##                     restrictionR1 <- Restriction[[z]]
##                     if(restrictionR1$col.restricting==column.restrict1){
##                         use.restriction <- restrictionR1
##                         break
##                     }
##                 }
##                 ind.res1 <- which(is.element(data[,use.restriction$col.restricted],
##                                              use.restriction$fac.restricted))
##                 data.use <- data[-ind.res1,]
##             }
##             if(Res.type=="Restrict2"){
##                 ind.use <- c()
##                 column.restrict2 <- column
##                 for(z in 1:length(Restriction)){
##                     restrictionR2 <- Restriction[[z]]
##                     if(is.element(restrictionR2$col.restricting,column.restrict2)){
##                         use.restriction <- restrictionR2
##                         ind.res2 <- which(is.element(data[,use.restriction$col.restricted],
##                                                      use.restriction$fac.restricted))
##                         ind.use <- union(ind.use,ind.res2)
##                     }
##                 }       
##                 data.use <- data[-ind.use,]
##             }
##             if(Res.type=="Restricted1"){
##                 column.restricted1 <- column[is.element(column,Col.Restricted)==TRUE]
##                 for(z in 1:length(Restriction)){
##                     restrictionRt1 <- Restriction[[z]]
##                     if(restrictionRt1$col.restricted==column.restricted1){
##                         use.restriction <- restrictionRt1
##                         break
##                     }
##                 }
##                 ind.rest1 <- which(is.element(data[,use.restriction$col.restricting],
##                                               use.restriction$fac.restricting))
##                 data.use <- data[ind.rest1,]
##             }
##             if(Res.type=="Restricted2"){
##                 ind.use <- seq(1:nrow(data))
##                 column.restricted2 <- column
##                 for(z in 1:length(Restriction)){
##                     restrictionRt2 <- Restriction[[z]]
##                     if(is.element(restrictionRt2$col.restricted,column.restricted2)){
##                         use.restriction <- restrictionRt2
##                         ind.rest2 <- which(is.element(data[,use.restriction$col.restricting],
##                                                       use.restriction$fac.restricting))
##                         ind.use <- intersect(ind.use,ind.rest2)
##                     }
##                 }    
##                 data.use <- data[ind.use,]
##             }
##             if(Res.type=="Normal"){
##                 data.use <- data
##             }
##             if(Res.type=="Restriction.mix.hard"){
##                 column.restrictedM <- column[is.element(column,Col.Restricted)==TRUE]
##                 column.restrictM <- column[is.element(column,Col.Restricting)==TRUE]
##                 for(z in 1:length(Restriction)){
##                     restrictionRtM <- Restriction[[z]]
##                     if(restrictionRtM$col.restricted==column.restrictedM){
##                         use.restricted <- restrictionRtM               
##                     }
##                     if(restrictionRtM$col.restricting==column.restrictM){
##                         use.restricting <- restrictionRtM               
##                     }
##                 }
##                 ## Keep ind.restrictedM
##                 ind.restrictedM <- which(is.element(data[,use.restricted$col.restricting],
##                                                     use.restricted$fac.restricting))
##                 data.use1 <- data[ind.restrictedM,]
##                 ## Remove ind.restrictM
##                 ind.restrictM <- which(is.element(data.use1[,use.restricting$col.restricted],
##                                                   use.restricting$fac.restricted))
##                 data.use <- data.use1[-ind.restrictM,]                
##             }

##             data <- data.use
##         }

##         ## Column name and Restriction Type
##         ## print("Column name")
##         ## print(column)
##         ## print("Restriction Type")
##         ## print(Res.type)
        
##         ##Treatment Effect
##         A <- tapply(data[,1],data[,column],mean,simplify=FALSE)
##         A2 <- tapply(data[,1],data[,column],mean,simplify=TRUE)
##         A[is.na(A2)] <- NA
##         Effect1 <- unlist(A)
##         if(base=="min"){
##             base <- which(Effect1==min(Effect1))
##         }        
        
##         ## Combination Effect
##         Comb.Effect <- Effect1 - Effect1[base]
##         Treatment.Effect <- cbind(Comb.Effect,
##                                   expand.grid(dimnames(A)))
##         Treatment.Effect <- na.omit(Treatment.Effect)
##         Comb.Effect <- Treatment.Effect[,1]
##         Combination.name <- Treatment.Effect[,-1]
##         Treatment.Effect.print <- Treatment.Effect
        
##         for(j in 2:ncol(Treatment.Effect)){
##             Treatment.Effect[,j] <- paste(colnames(Treatment.Effect)[j],
##                                           Treatment.Effect[,j],
##                                           sep="_")
##         }

##         ## Sum of Marginals
##         Treatment.Effect1 <- Marginal(data=data,
##                                       base=base,
##                                       Restriction=Restriction)$Treatment.Effect1
##         Main.comb <- c()
##         for(i in 1:nrow(Treatment.Effect)){
##             main <- c()
##             for(j in 1:(ncol(Treatment.Effect)-1)){
##                 main[j] <- Treatment.Effect1[rownames(Treatment.Effect1)==
##                                              Treatment.Effect[i,(j+1)],1]
##             }
##             Main.comb[i] <- sum(main)
##         }
##         Sum.Mar.Effect <- Main.comb - Main.comb[base]

##         if(order==2){
##             TIE <- round(Comb.Effect - Sum.Mar.Effect, digits=8)

##         }
##         if(order==3){
##             ##Two combination Effect List
##             data.column.three <- data[,column]
##             Two.way.comb.list <- list()
##             for(j in 1:3){
##                 column.com.three <- c()
##                 Two.way.comb <- list()
##                 column.com.three <-
##                     colnames(data.column.three)[c(combn(length(data.column.three),2)[,j])]
##                 A.Three <- tapply(data[,1],data[,column.com.three],mean,simplify=FALSE)
##                 A2.Three <- tapply(data[,1],data[,column.com.three],mean,simplify=TRUE)
##                 A.Three[is.na(A2.Three)] <- NA
##                 EffectTwo <- unlist(A.Three)
##                 Two.way.comb <- EffectTwo - EffectTwo[base]
##                 Two.way.comb2 <- cbind(Two.way.comb,
##                                           expand.grid(dimnames(A.Three)))
##                 for(k in 2:3){
##                     Two.way.comb2[,k] <- paste(colnames(Two.way.comb2)[k],
##                                                Two.way.comb2[,k],sep="_")
##                 }
##                 colnames(Two.way.comb2) <- c("Two.way.comb","first","second")
##                 Two.way.comb.list[[j]]  <- na.omit(Two.way.comb2)
##                 ## print(j)
##                 ## print("th")
##                 ## print(Two.way.comb.list[[j]])
##             }
##             Two.way.comb.Final <- do.call(rbind,Two.way.comb.list)


##             Two.Way.sum <- c()
##             for(i in 1:nrow(Treatment.Effect)){
##                 two.sum <- c()
##                 for(j in 1:nrow(Two.way.comb.Final)){
##                     ind.two <- sum(is.element(Two.way.comb.Final[j,2:3],
##                                          Treatment.Effect[i,2:4]))
##                     if(ind.two==2){
##                         two.sum[j] <- Two.way.comb.Final[j,1]
##                     }else{
##                         two.sum[j] <- 0
##                     }
##                 }
##                 Two.Way.sum[i] <- sum(two.sum)
##             }
##             Sum.Two.Effect <- Two.Way.sum - Two.Way.sum[base]
##             TIE <- round(Comb.Effect - Sum.Two.Effect + Sum.Mar.Effect, digits=8)
##         }
       
##         ## print("The range of TIE")
##         range.change <- max(TIE) - min(TIE)
##         ## Order diff
##         Order.data <- cbind(Sum.Mar.Effect,Comb.Effect)
##         Order.data <- Order.data[order(Order.data[,1],decreasing=TRUE),]
##         Order.data <- as.data.frame(Order.data)
##         Order.data$order.main <- seq(1:nrow(Order.data))
##         Order.comb <- Order.data[order(Order.data[,2],decreasing=TRUE),"order.main"]
##         Order.diff <- Order.comb - seq(1:nrow(Order.data))
##         max.order.diff <- max(abs(Order.diff)) 
##         ## print(range.change)
##         return(list("A"=A,
##                     "Treatment.Effect"=Treatment.Effect,
##                     "Comb.Effect"=Comb.Effect,
##                     "Combination.name"=Combination.name,
##                     "Sum.Mar.Effect"=Sum.Mar.Effect,
##                     "TIE"=TIE,
##                     "range.change"=range.change,
##                     "max.order.diff"=max.order.diff))
##     }

##     if(compare==FALSE){
##         if(!missing(column)){
##             X <- comb.prem(data=data,column=column,Restriction=Restriction,order=order)
##             A <- X$A
##             Treatment.Effect <- X$Treatment.Effect
##             Comb.Effect <- X$Comb.Effect
##             Combination.name <- X$Combination.name
##             Sum.Mar.Effect <- X$Sum.Mar.Effect
##             TIE <- X$TIE
##             range.change <- X$range.change     
            
##         }
##     }else{        
##         data.column <- data[,-1]
##         column.com <- list()
##         Range.com <- list()
##         for(j in 1:choose(length(data.column),order)){
##             column.com[[j]] <-
##                 colnames(data.column)[c(combn(length(data.column),order)[,j])]
##             if(measure=="TIE"){
##                 Range.com[[j]] <- comb.prem(data=data,
##                                             column=column.com[[j]],
##                                             Restriction=Restriction,
##                                             order=order)$range.change
##             }
##             if(measure=="Order.Diff"){
##                 Range.com[[j]] <- comb.prem(data=data,
##                                             column=column.com[[j]],
##                                             Restriction=Restriction,
##                                             order=order)$max.order.diff
##             }
##             ## print(Range.com[[j]])
##             ## print(column.com[[j]])
##         }
##         column.com <- as.data.frame(matrix(unlist(column.com),ncol=order,byrow=TRUE))
##         Compare <- as.data.frame(cbind(column.com,unlist(Range.com)))
##         Compare <- Compare[order(Compare[,(order+1)],decreasing=TRUE),]
##         if(measure=="TIE"){
##             colnames(Compare) <- c(colnames(Compare)[1:order],"range.TIE")
##         }
##         if(measure=="Order.Diff"){
##             colnames(Compare) <- c(colnames(Compare)[1:order],"max.Order.diff")
##         }
##         return(Compare)
##     }


##     ##Put the coefficients together
##     if(!missing(column)){
##         if(length(column)==2){
##             Treatment.coefs.name <- Treatment.Effect[,2]
##             if(ncol(Treatment.Effect)>=3){
##                 for(i in 1:nrow(Treatment.Effect)){
##                     for(j in 3:ncol(Treatment.Effect)){
##                         Treatment.coefs.name[i] <-
##                             paste(Treatment.coefs.name[i],
##                                   Treatment.Effect[i,j],
##                                   sep=".")
##                     }
##                 }
##             }
##             Treatment.coefs.name <- gsub(" ", ".", Treatment.coefs.name)
##         }
##     }

##     if(missing(column)){
##         Treatment.coefs.name <- rownames(Treatment.Effect1)
##         Treatment.coefs.name <- gsub(" ", ".", Treatment.coefs.name)
##     }

##     ##print("This is point")
##     ##print(Treatment.coefs.name)
##     ## print(Treatment.coefs.name2)
##     ##print(names(coefs))
##     ##print(Treatment.coefs.name[1])
##     ##print(dim(Treatment.Effect))
##     ##print(is.element(names(coefs), Treatment.coefs.name))
##     ##print(coefs)
##     ##print(coefs[names(coefs)==Treatment.coefs.name[1]])

##     if(missing(column)){
##         Coefficients <- c()
##         for(i in 1:length(Treatment.coefs.name)){
##             if(Treatment.coefs.name[i] %in% names(coefs)){
##                 Coefficients[i] <- coefs[names(coefs)==Treatment.coefs.name[i]]                
##             }else{
##                 Coefficients[i] <- NA
##             }
##         }
##         if(outcome.type=="binary"){
##             Coefficients <- Coefficients/2
##         }
##     }
##     if(!missing(column)){
##         if(length(column)==2){
##             Coefficients <- c()
##             for(i in 1:length(Treatment.coefs.name)){
##                 if(Treatment.coefs.name[i] %in% names(coefs)){
##                     Coefficients[i] <- coefs[names(coefs)==Treatment.coefs.name[i]]
##                 }else{
##                     Coefficients[i] <- NA
##                 }
##             }
        
##             if(outcome.type=="binary"){
##                 Coefficients <- Coefficients/2
##             }
##         }
##     }



##     ## print(Treatment.coefs)

##     if(!missing(column)){
##         if(length(column)==2){
##             Treatment.Effect.print <- cbind(Comb.Effect,
##                                             TIE,
##                                             Coefficients,
##                                             Combination.name)
##             Main.matrix <- cbind(Sum.Mar.Effect,
##                                  TIE,
##                                  Coefficients,
##                                  Combination.name)
##             Inequality.Cor <- cor(Sum.Mar.Effect,TIE)
##             Relative.change <- cbind(TIE,Combination.name)
##             Relative.change <- Relative.change[order(Relative.change[,1],decreasing=TRUE),]
##         }
##         if(length(column)>=3){
##             Treatment.Effect.print <- cbind(Comb.Effect,
##                                             TIE,
##                                             Combination.name)
##             Main.matrix <- cbind(Sum.Mar.Effect,
##                                  TIE,
##                                  Combination.name)
##             Relative.change <- cbind(TIE,Combination.name)
##             Relative.change <- Relative.change[order(Relative.change[,1],decreasing=TRUE),]
##         }
##     }else{
##         Treatment.Effect.print <- cbind(Treatment.Effect1,
##                                         Coefficients
##                                         )
##         Man.matrix <- NULL
##     }

##     Treatment.Effect.print <- as.data.frame(Treatment.Effect.print)
##     if(!missing(column)){
##         Main.matrix <- as.data.frame(Main.matrix)
##     }

##     if(sort==TRUE){
##         Treatment.Effect.print1 <- Treatment.Effect.print
##         Treatment.Effect.print <-
##             Treatment.Effect.print[order(Treatment.Effect.print[,1],
##                                          decreasing=TRUE),]
        
##         if(!missing(column)){
##             ## Order.data <- cbind(Main.matrix$Sum,Treatment.Effect.print1[,1])
##             ## Order.data <- Order.data[order(Order.data[,1],decreasing=TRUE),]
##             ## Order.data <- as.data.frame(Order.data)
##             ## Order.data$order.main <- seq(1:nrow(Order.data))
##             ## Order.comb <- Order.data[order(Order.data[,2],decreasing=TRUE),"order.main"]
##             ## Order.diff <- Order.comb - seq(1:nrow(Order.data)) 
            
##             Main.matrix <-
##                 Main.matrix[order(Main.matrix$Sum.Mar.Effect,
##                                   decreasing=TRUE),]
            
##             Treatment.Effect.print <- cbind(
##                 Treatment.Effect.print[,1:2],
##                 ## Order.diff,
##                 Treatment.Effect.print[,3:ncol(Treatment.Effect.print)])
##         }
##     }

##     if(!missing(column)){
##         if(range.change>0){
##             print("Interaction matters")
##         }else{
##             print("No Interaction Effect")
##         }
##     }

##     ## if(!missing(column)){
##     ##     if(all(Treatment.Effect.print[,"Coefficients"]==Main.matrix[,"Coefficients"])){
##     ##         print("No Significant Interaction Effect")
##     ##     }else{
##     ##         print("Interaction Effect changes the order")
##     ##     }
##     ## }
##     if(!missing(column)){
##         if(sort){
##             return(list("Range Treatment Interactions"=range.change,
##                         "Treatment Interactions"=Relative.change,
##                         ## "Inequality Correlation"=Inequality.Cor,
##                         ## "Absolute Order change"=max(abs(Order.diff)),
##                         "Combination Effect"=Treatment.Effect.print,
##                         "Sum of Marginal Effect"=Main.matrix
##                         ))
##         }else{
##             return(list("Range Treatment Interactions"=range.change,
##                         "Treatment Interactions"=Relative.change,
##                         ## "Inequality Correlation"=Inequality.Cor,
##                         "Combination Effect"=Treatment.Effect.print,
##                         "Sum of Marginal Effect"=Main.matrix
##                         ))
##         }
##     }else{
##         return(list("Treatment.Effect"=Treatment.Effect.print))
##     }
## }

## Restriction should look like this.
## restriction1 <- list("col.restricted"="FeatReason",
##                     "fac.restricted"="Escape",
##                     "col.restricting"="FeatCountry",
##                     "fac.restricting"=c("China","Sudan","Somalia","iraq"))

## restriction2 <- list("col.restricted"="FeatJob",
##                     "fac.restricted"=c(
##                         "Doctor","Research scientist",
##                         "Computer programmer","Financial analyst"),
##                     "col.restricting"="FeatEd",
##                     "fac.restricting"=c("Twoyears.college","College","Graduate"))

## Res <- list(restriction1,restriction2)
