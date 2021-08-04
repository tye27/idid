#' Simulate dataset
#'
#' @param case Simulation scenarios used in Tables 1-2 in Ye et al., (2020)
#' @param n Sample size
#'
#' @return A data frame
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2020). Instrumented Difference-in-Differences.
#'
#' @import stats truncnorm
#' @export
#' @examples
#' df<-data_gen("case1",1e5)
data_gen<-function(case,n){
  if(case=="case1"){
    z<-rbinom(n,1,0.5)
    u0<-rtruncnorm(n,a=-1,b=1,mean = 0,sd = 1)
    u1<-1+rtruncnorm(n,a=-1,b=1,mean = 0,sd = 1)
    x<-rnorm(n)
    pi_d0<-(z+1)/8*u0+0.5
    pi_d1<-(z+1)/8*u1+0.5
    d0<-rbinom(n,1,pi_d0)
    d1<-rbinom(n,1,pi_d1)
    y0<-(1+x)*d0+2+2*u0+z+x+rnorm(n) # observed outcome at t=0
    y1<-(1+x)*d1+2+2*u1+z+x+rnorm(n) # observed outcome at t=1
  }else if (case=="case2"){
    x<-rnorm(n)
    z<-rbinom(n,1,exp(x/2)/(1+exp(x/2)))
    u0<-rnorm(n,0,1)
    u1<-rnorm(n,1,1)
    pi_d0<-(-1)*z*u0-0.5+1.5*u0
    pi_d1<-(-1)*z*u1-0.5+1.5*u1
    d0<-rbinom(n,1,exp(pi_d0)/(1+exp(pi_d0)))
    d1<-rbinom(n,1,exp(pi_d1)/(1+exp(pi_d1)))
    y0<-(1+x)*d0+2+2*u0+z+x+rnorm(n) # observed outcome at t=0
    y1<-(1+x)*d1+2+2*u1+z+x+rnorm(n) # observed outcome at t=1
  }
  df_full<-data.frame(z,x,d0,d1,y0,y1)
  t<-rbinom(n,1,.5)
  df_full$t<-t
  df0<-df_full[df_full$t==0,c(1,2,3,5,7)]
  df1<-df_full[df_full$t==1,c(1,2,4,6,7)]
  names(df0)<-c("z","x","d","y","t")
  names(df1)<-c("z","x","d","y","t")
  df<-rbind(df0,df1)
  return(df)
}


#' One-sample Instrumented DID Wald estimator
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#'
#' @return
#' \describe{
#' \item{beta.est}{One-sample instrumented DID Wald estimate.}
#' \item{beta.se}{Standard error of \code{beta.est}.}
#' \item{f.statistic}{F statistic from the first-stage regression. It is preferably larger than 10.}
#' }
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2020). Instrumented Difference-in-Differences.
#'
#' @import stats
#' @export
#' @examples
#' df<-data_gen("case1",1e5)
#' idid_wald1(df$y,df$d,df$t,df$z)
idid_wald1<-function(yv,dv,tv,zv){
  n<-length(yv)
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)
  delta_D<-(mean(df$d[df$t==1 & df$z==1],na.rm=TRUE)- mean(df$d[df$t==0 & df$z==1],na.rm=TRUE) )-
    (mean(df$d[df$t==1 & df$z==0],na.rm=TRUE)- mean(df$d[df$t==0 & df$z==0],na.rm=TRUE))
  delta_Y<-(mean(df$y[df$t==1 & df$z==1],na.rm=TRUE)- mean(df$y[df$t==0 & df$z==1],na.rm=TRUE) )-
    (mean(df$y[df$t==1 & df$z==0],na.rm=TRUE)- mean(df$y[df$t==0 & df$z==0],na.rm=TRUE))
  beta.est<-delta_Y/delta_D
  tmp<-numeric(4)
  i<-1
  for(t in 0:1){
    for(z in 0:1){
      tmp[i]<-var(df$y[df$t==t& df$z==z]-beta.est*df$d[df$t==t& df$z==z],na.rm = TRUE)/
        (length(which(df$t==t & df$z==z))/n)
      i<-i+1
    }
  }
  beta.se<-sqrt(sum(tmp)/delta_D^2)/sqrt(n)
  df$ind<-1*(df$z==1 & df$t==1)
  fit1<-lm(ind~t+z,data=df)
  fit2<-lm(d~fit1$residuals,data=df)
  f.statistic<-anova(fit2)$F[1]
  list(beta.est=beta.est,beta.se=beta.se,f.statistic=f.statistic)
}


#' Two-sample Summary-Data Instrumented DID Wald estimator
#'
#' @param mu.y A 4-dimensional vector of outcome means, for every (t,z) combination
#' @param mu.d A 4-dimensional vector of exposure rates, for every (t,z) combination
#' @param se.y A 4-dimensional vector of the standard error of \code{mu.y}, for every (t,z) combination
#' @param se.d A 4-dimensional vector of the standard error of \code{mu.d}, for every (t,z) combination
#' @param tv A 4-dimensional vector for time
#' @param zv A 4-dimensional vector for the level of IV for DID
#'
#' @return
#' \describe{
#' \item{beta.est}{Two-sample summary-data instrumented DID Wald estimate.}
#' \item{beta.se}{Standard error of \code{beta.est}.}
#' \item{f.statistic}{F statistic from the first-stage regression. It is preferably larger than 10.}
#' }
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2020). Instrumented Difference-in-Differences.
#'
#' @import stats
#' @export
#'
#' @examples
#'dfa<-data_gen(case="case1",n=1e5) # outcome dataset
#'dfb<-data_gen(case="case1",n=1e5) # exposure dataset
#'res<-data.frame(mu.y=0,mu.d=0, se.y=0, se.d=0, t=c(0,1,0,1),z=c(0,0,1,1))
#'for(i in 1:4){
#' t<-res$t[i]
#' z<-res$z[i]
#' res$mu.y[i]<-mean(dfa$y[dfa$t==t & dfa$z==z])
#' res$mu.d[i]<-mean(dfb$d[dfb$t==t & dfb$z==z])
#' res$se.y[i]<-sqrt(var(dfa$y[dfa$t==t & dfa$z==z])/length(dfa$y[dfa$t==t & dfa$z==z]))
#' res$se.d[i]<-sqrt(var(dfb$d[dfb$t==t & dfb$z==z])/length(dfb$d[dfb$t==t & dfb$z==z]))
#'}
#'idid_wald2(mu.y=res$mu.y,mu.d=res$mu.d,se.y=res$se.y,se.d=res$se.d,tv=res$t,zv=res$z)
#'
idid_wald2<-function(mu.y,mu.d,se.y,se.d,tv,zv){
  delta_Y<-mu.y[tv==1 & zv==1]-mu.y[tv==0 & zv==1]-mu.y[tv==1 & zv==0]+mu.y[tv==0 & zv==0]
  delta_D<-mu.d[tv==1 & zv==1]-mu.d[tv==0 & zv==1]-mu.d[tv==1 & zv==0]+mu.d[tv==0 & zv==0]
  beta.est<-delta_Y/delta_D
  se_delta_D<-sqrt(sum(se.d^2))
  f.statistic<-(delta_D/se_delta_D)^2
  wald.var<-(sum(se.y^2)+beta.est^2*sum(se.d^2))/delta_D^2
  return(list(beta.est=beta.est, beta.se=sqrt(wald.var), f.statistic=f.statistic))
}


# the efficient influence function for the iDID semiparametric estimator
IF<-function(psi, yv,dv,tv,zv, V, delta_Y, delta_D, pi, y_fit, d_fit){
  res<-as.vector(delta_Y/delta_D - (V %*% psi) +
                   (2*zv-1)*(2*tv-1)/pi/delta_D*( yv-y_fit- delta_Y/delta_D*(dv-d_fit) ))* V
  res_mean<-apply(res,2,mean)
  return(res_mean)
}

#' Find the root based on the efficient influence function
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#'
#' @details This function can be modified if one would like to use alternative working models.
#' @return A vector of roots
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2020). Instrumented Difference-in-Differences.
#'
#' @import rootSolve
#' @export
#'
#' @examples
#'df<-data_gen("case2",1e5)
# constant working model
#'V<-matrix(1,ncol=1,nrow=dim(df)[1])
#'idid_semi_root(df$y,df$d,df$t,df$z,df$x,V)
# linear working model
#'V<-cbind(matrix(1,ncol=1,nrow=dim(df)[1]),df$x)
#'idid_semi_root(df$y,df$d,df$t,df$z,df$x,V)
idid_semi_root<-function(yv,dv,tv,zv,X,V){
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)
  fit_y<-lm(y~as.factor(t)*as.factor(z)*X,data=df)
  fit_d<-lm(d~as.factor(t)*as.factor(z)*X,data=df)
  mu_y_fit<-matrix(nrow=dim(df)[1], ncol=4)
  mu_d_fit<-matrix(nrow=dim(df)[1], ncol=4)
  df_tmp<-df
  b<-1
  for(t in 0:1){
    for(z in 0:1){
      df_tmp$t<-as.factor(t)
      df_tmp$z<-as.factor(z)
      mu_y_fit[,b]<-predict(fit_y,df_tmp)
      mu_d_fit[,b]<-predict(fit_d,df_tmp)
      b<-b+1
    }
  }
  delta_Y<-(mu_y_fit[,1]-mu_y_fit[,2]-mu_y_fit[,3]+mu_y_fit[,4])
  delta_D<-(mu_d_fit[,1]-mu_d_fit[,2]-mu_d_fit[,3]+mu_d_fit[,4])
  pi_z<-glm(z~X+t,data=df,family = "binomial")$fitted.values
  pi_t<-glm(t~X,data=df,family = "binomial")$fitted.values
  pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
    ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
    (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
    ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)
  y_fit<-fit_y$fitted.values
  d_fit<-fit_d$fitted.values
  res<-multiroot(IF,start=rep(0,dim(V)[2]),  yv=yv,dv=dv,tv=tv,zv=zv, V=V, delta_Y=delta_Y,
                 delta_D=delta_D, pi=pi, y_fit=y_fit, d_fit=d_fit)$root
  return(res)
}


idid_semi_bootstrap_sd<-function(yv,dv,tv,zv,X,V, n.boot, func_root, alpha=.05){
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)
  npar<-dim(V)[2]
  boot_res<-matrix(nrow=n.boot,ncol=npar)
  for(i in 1:n.boot){
    id_boot<-sample.int(n=dim(df)[1],replace = TRUE)
    df_boot<-df[id_boot,]
    V_boot<-V[id_boot,,drop=FALSE]
    X_boot<-X[id_boot,,drop=FALSE]
    boot_res[i,]<-func_root(df_boot$y,df_boot$d,df_boot$t,df_boot$z,X_boot,V_boot)
  }
  boot_se<-numeric(npar)
  boot_se<-apply(boot_res,2,function(x){(quantile(x,probs = 1-alpha/2)-quantile(x,probs = alpha/2))/2/qnorm(alpha/2,lower.tail = FALSE)
  })
  boot_se
}


#' Instrumented DID semiparametric estimator
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#' @param func_root A function that calculates the root, e.g., idid_semi_root
#' @param n.boot Number of bootstrap interations
#'
#' @return
#' \describe{
#' \item{psi.est}{Semiparametric estimators.}
#' \item{psi.se}{Standard error of \code{psi.est}.}
#' \item{f.statistic}{F statistic from the first-stage regression. It is preferably larger than 10.}
#' }
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2020). Instrumented Difference-in-Differences.
#'
#' @import rootSolve
#' @export
#'
#' @examples
#'df<-data_gen("case2",1e5)
# constant working model
#'V<-matrix(1,ncol=1,nrow=dim(df)[1])
#'idid_semi(df$y,df$d,df$t,df$z,matrix(df$x,ncol=1),V,idid_semi_root,n.boot = 10)
#'# n.boot is set to a small number for illustrative purpose
# linear working model
#'V<-cbind(matrix(1,ncol=1,nrow=dim(df)[1]),df$x)
#'idid_semi(df$y,df$d,df$t,df$z,matrix(df$x,ncol=1),V,idid_semi_root,n.boot = 10)
#'# n.boot is set to a small number for illustrative purpose
idid_semi<-function(yv,dv,tv,zv,X,V,func_root,n.boot=100){
  psi.est<-func_root(yv,dv,tv,zv,X,V)
  psi.se<-idid_semi_bootstrap_sd(yv,dv,tv,zv,X,V,n.boot=n.boot,func_root = func_root)
  ind<-1*(zv==1 & tv==1)
  fit1<-lm(ind~tv+zv+X)
  fit2<-lm(dv~fit1$residuals)
  f.statistic<-anova(fit2)$F[1]
  return(list(psi.est=psi.est, psi.se=psi.se,f.statistic=f.statistic))
}





