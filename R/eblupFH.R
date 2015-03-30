eblupFH <-
function(formula,vardir,method="REML",MAXITER=100,PRECISION=0.0001,data) {

   result <- list(eblup=NA, 
                  fit=list(method=method, convergence=TRUE, iterations=0, estcoef=NA, 
                  refvar=NA, goodness=NA)
                 ) 

   if (method!="REML" & method!="ML" & method!="FH")
       stop(" method=\"",method, "\" must be \"REML\", \"ML\" or \"FH\".")

   namevar     <- deparse(substitute(vardir))
   if (!missing(data))
   {
      formuladata <- model.frame(formula,na.action = na.omit,data)
      X           <- model.matrix(formula,data)        
      vardir      <- data[,namevar]
   } else
   {
      formuladata <- model.frame(formula,na.action = na.omit)
      X <- model.matrix(formula)        
   }
   y <- formuladata[,1]            

   if (attr(attributes(formuladata)$terms,"response")==1)
      textformula <- paste(formula[2],formula[1],formula[3])
   else
      textformula <- paste(formula[1],formula[2])

   if (length(na.action(formuladata))>0)
      stop("Argument formula=",textformula," contains NA values.")
   if (any(is.na(vardir)))
      stop("Argument vardir=",namevar," contains NA values.")

  m<-length(y) # Sample size or number of areas
  p<-dim(X)[2] # Num. of X columns of num. of auxiliary variables (including intercept)
  Xt<-t(X)
  
  # Fisher-scoring algorithm for ML estimator of variance A starts
  
  if (method=="ML") {    

    # Initial value of variance A is fixed to the median of sampling variances vardir
    Aest.ML<-0
    Aest.ML[1]<-median(vardir)
    
    k<-0
    diff<-PRECISION+1
    while ((diff>PRECISION)&(k<MAXITER))
    {
      k<-k+1
      Vi<-1/(Aest.ML[k]+vardir)
      XtVi<-t(Vi*X)
      Q<-solve(XtVi%*%X)
      P<-diag(Vi)-t(XtVi)%*%Q%*%XtVi
      Py<-P%*%y
      # Score function obtained from restricted log-likelihood
      s<-(-0.5)*sum(Vi)+0.5*(t(Py)%*%Py) 
      # Fisher information obtained from restricted log-likelihood 
      F<-0.5*sum(Vi^2)                  
      # Updating equation
      Aest.ML[k+1]<-Aest.ML[k]+s/F
      # Relative difference of estimators in 2 iterations for stopping condition
      diff<-abs((Aest.ML[k+1]-Aest.ML[k])/Aest.ML[k])
    } # End of while

    # Final estimator of variance A
    A.ML<-max(Aest.ML[k+1],0)
    # print(Aest.ML)
    
    # Indicator of convergence
    result$fit$iterations  <- k  
    #if(k<MAXITER) {conv<-TRUE} else {conv<-FALSE}
    if(k>=MAXITER && diff>=PRECISION) 
    {
       result$fit$convergence <- FALSE
       return(result)
    }

  
    # Computation of the coefficients' estimator beta  
  
    Vi<-1/(A.ML+vardir)
    XtVi<-t(Vi*X)
    Q<-solve(XtVi%*%X)
    beta.ML<-Q%*%XtVi%*%y
  
    # Significance of the regression coefficients
  
    varA<-1/F
  
    std.errorbeta<-sqrt(diag(Q))
    tvalue<-beta.ML/std.errorbeta
    pvalue<-2*pnorm(abs(tvalue),lower.tail=FALSE)
  
    # Goodness of fit measures: loglikelihood, AIC, BIC
    
    Xbeta.ML<-X%*%beta.ML
    resid<-y-Xbeta.ML
    
    loglike<-(-0.5)*(sum(log(2*pi*(A.ML+vardir))+(resid^2)/(A.ML+vardir)))
    AIC<-(-2)*loglike+2*(p+1)
    BIC<-(-2)*loglike+(p+1)*log(m)
  
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
    
    # Computation of the empirical best (EB) predictor
    coef     <- data.frame(beta=beta.ML,std.error=std.errorbeta,tvalue,pvalue)
    variance <- A.ML
    EBLUP    <- Xbeta.ML+A.ML*Vi*resid

    # Fisher-scoring algorithm for REML estimator of variance A starts
 
  } else if (method=="REML") {  

    # Initial value of variance A is fixed to the median of sampling variances vardir
    Aest.REML<-0
    Aest.REML[1]<-median(vardir)
    
    k<-0
    diff<-PRECISION+1
    while ((diff>PRECISION)&(k<MAXITER))
    {
      k<-k+1
      Vi<-1/(Aest.REML[k]+vardir)
      XtVi<-t(Vi*X)
      Q<-solve(XtVi%*%X)
      P<-diag(Vi)-t(XtVi)%*%Q%*%XtVi
      Py<-P%*%y
      # Score function obtained from restricted log-likelihood
      s<-(-0.5)*sum(diag(P))+0.5*(t(Py)%*%Py) 
      # Fisher information obtained from restricted log-likelihood 
      F<-0.5*sum(diag(P%*%P))                  
      # Updating equation
      Aest.REML[k+1]<-Aest.REML[k]+s/F
      # Relative difference of estimators in 2 iterations for stopping condition
      diff<-abs((Aest.REML[k+1]-Aest.REML[k])/Aest.REML[k])
    } # End of while

    # Final estimator of variance A

    A.REML<-max(Aest.REML[k+1],0)
    # print(Aest.REML)
    
    # Indicator of convergence
    result$fit$iterations  <- k  
    #if(k<MAXITER) {conv<-TRUE} else {conv<-FALSE}
    if(k>=MAXITER && diff>=PRECISION) 
    {
       result$fit$convergence <- FALSE
       return(result)
    }
  
    # Computation of the coefficients' estimator beta
    Vi<-1/(A.REML+vardir)
    XtVi<-t(Vi*X)
    Q<-solve(XtVi%*%X)
    beta.REML<-Q%*%XtVi%*%y
  
    # Significance of the regression coefficients
    varA<-1/F
  
    std.errorbeta<-sqrt(diag(Q))
    tvalue<-beta.REML/std.errorbeta
    pvalue<-2*pnorm(abs(tvalue),lower.tail=FALSE)
  
    # Goodness of fit measures: loglikelihood, AIC, BIC
    
    Xbeta.REML<-X%*%beta.REML
    resid<-y-Xbeta.REML
    
    loglike<-(-0.5)*(sum(log(2*pi*(A.REML+vardir))+(resid^2)/(A.REML+vardir)))
    AIC<-(-2)*loglike+2*(p+1)
    BIC<-(-2)*loglike+(p+1)*log(m)
  
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
    
    # Computation of the empirical best (EB) predictor
  
    coef     <- data.frame(beta=beta.REML,std.error=std.errorbeta,tvalue,pvalue)
    variance <- A.REML
    EBLUP    <- Xbeta.REML+A.REML*Vi*resid
 
 
    # Fisher-scoring algorithm for REML estimator of variance A starts
 
  } else   #FH
  {
  
    # Initial value of variance A is fixed to the median of sampling variances vardir
    Aest.FH<-NULL
    Aest.FH[1]<-median(vardir)
  
    k<-0
    diff<-PRECISION+1
    while ((diff>PRECISION)&(k<MAXITER)){


      k<-k+1
      Vi<-1/(Aest.FH[k]+vardir)
      XtVi<-t(Vi*X)
      Q<-solve(XtVi%*%X)
      betaaux<-Q%*%XtVi%*%y
      resaux<-y-X%*%betaaux   
      # Left-hand side of equation for FH estimator
      s<-sum((resaux^2)*Vi)-(m-p)
      # Expectation of negative derivative of s 
      F<-sum(Vi)
      # Updating equation
      Aest.FH[k+1]<-Aest.FH[k]+s/F
      # Relative difference of estimators in 2 iterations for stopping condition
      diff<-abs((Aest.FH[k+1]-Aest.FH[k])/Aest.FH[k])
  
    } # End of while
 
    A.FH<-max(Aest.FH[k+1],0)
    #print(Aest.FH)
  
    # Indicator of convergence
    result$fit$iterations  <- k  
    #if(k<MAXITER) {conv<-TRUE} else {conv<-FALSE}
    if(k>=MAXITER && diff>=PRECISION) 
    {
       result$fit$convergence <- FALSE
       return(result)
    }
    
    # Computation of the coefficients' estimator beta
  
    Vi<-1/(A.FH+vardir)
    XtVi<-t(Vi*X)
    Q<-solve(XtVi%*%X)
    beta.FH<-Q%*%XtVi%*%y
  
    # Significance of the regression coefficients
  
    varA<-1/F
    varbeta<-diag(Q)
    std.errorbeta<-sqrt(varbeta)
    zvalue<-beta.FH/std.errorbeta
    pvalue<-2*pnorm(abs(zvalue),lower.tail=FALSE)
  
    # Goodness of fit measures: loglikelihood, AIC, BIC
    
    Xbeta.FH<-X%*%beta.FH
    resid<-y-Xbeta.FH
  
    loglike<-(-0.5)*(sum(log(2*pi*(A.FH+vardir))+(resid^2)/(A.FH+vardir)))
    AIC<-(-2)*loglike+2*(p+1)
    BIC<-(-2)*loglike+(p+1)*log(m)
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
  
    # Computation of the empirical best (EB) predictor
    coef     <- data.frame(beta=beta.FH,std.error=std.errorbeta,tvalue=zvalue,pvalue)
    variance <- A.FH
    EBLUP    <- Xbeta.FH+A.FH*Vi*resid
  }  

  result$fit$estcoef  <- coef
  result$fit$refvar    <- variance
  result$fit$goodness <- goodness

  result$eblup <- EBLUP

  return(result)

}

