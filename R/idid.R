#' Simulate dataset (used in Table 1 of Ye et al. (2021))
#'
#' @param n Sample size
#'
#' @return A data frame
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @import stats
#' @export
#' @examples
#' df<-data_gen(1e5)
data_gen<-function(n){
  x1<-rnorm(n)
  x2<-rnorm(n)
  z<-rbinom(n,1,exp((1*(x1>0)+1*(x2>0))/2)/(1+exp((1*(x1>0)+1*(x2>0))/2)))
  u0<-rnorm(n,-1,1)
  u1<-rnorm(n,1,1)
  pi_d0<-(-1)*z*u0-0.5+1.5*u0
  pi_d1<-(-1)*z*u1-0.5+1.5*u1
  d0<-rbinom(n,1,exp(pi_d0)/(1+exp(pi_d0)))
  d1<-rbinom(n,1,exp(pi_d1)/(1+exp(pi_d1)))
  y0<-(1+x1+x2)*d0+2+2*u0+z+(1+x1+x2)+rnorm(n) # observed outcome at t=0
  y1<-(1+x1+x2)*d1+2+2*u1+z+(1+x1+x2)+rnorm(n) # observed outcome at t=1
  df_full<-data.frame(z,x1,x2,d0,d1,y0,y1)
  t<-rbinom(n,1,.5)
  df_full$t<-t
  df0<-df_full[df_full$t==0,c(1,2,3,4,6,8)]
  df1<-df_full[df_full$t==1,c(1,2,3,5,7,8)]
  names(df0)<-c("z","x1","x2","d","y","t")
  names(df1)<-c("z","x1","x2","d","y","t")
  df<-rbind(df0,df1)
  return(df)
}

expit<-function(x){
  return(1/(1+exp(-x)))
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
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @import stats
#' @export
#' @examples
#' df<-data_gen(1e5)
#' with(df[df$x1<=0 & df$x2<=0,],idid_wald1(y,d,t,z))
#' with(df[df$x1<=0 & df$x2>0,],idid_wald1(y,d,t,z))
#'
#'
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
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @import stats
#' @export
#'
#' @examples
#'dfa<-data_gen(1e5) # outcome dataset
#'dfb<-data_gen(1e5) # exposure dataset
#'dfa<-dfa[dfa$x1<=0&dfa$x2<=0,]
#'dfb<-dfb[dfb$x1<=0&dfb$x2<=0,]
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
IF.n<-function(psi,yv,dv,tv,zv,V,delta, delta_D, pi, y_fit, d_fit){
  res<-as.vector(delta - (V %*% psi) +
                   (2*zv-1)*(2*tv-1)/pi/delta_D*( yv-y_fit- delta*(dv-d_fit) ))* V
  return(res)
}

IF<-function(psi,yv,dv,tv,zv,V,delta, delta_D, pi, y_fit, d_fit){
  return(apply(IF.n(psi,yv,dv,tv,zv,V,delta, delta_D, pi, y_fit, d_fit),2,mean))
}

#' Regression-based semiparametric estimator
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#'
#' @details Here, all the imposed nuisance parameters are linear in X. beta(v;psi) is also linear in v. This function can be modified if one would like to use alternative nuisance parameter specifications.
#' @return A list
#' \describe{
#' \item{est}{Estimated psi.}
#' \item{se}{Standard error of \code{est}.}
#' }
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @export
#'
#' @examples
#'df<-data_gen(1e5)
#'X<-as.matrix(cbind(matrix(data=1,nrow = dim(df)[1],ncol=1),df[,grepl("x",names(df))]))
# constant working model
#'V<-matrix(data=1,nrow = dim(df)[1],ncol=1)
#'idid_reg(df$y,df$d,df$t,df$z,X,V)
# linear working model
#'V<-as.matrix(cbind(V,df[,grepl("x1",names(df))]))
#'idid_reg(df$y,df$d,df$t,df$z,X,V)
idid_reg<-function(yv,dv,tv,zv,X,V){
  V<-as.matrix(V)
  x.reg<-X
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)

  g<-function(all.param,df){
    psi<-all.param[1:dim(V)[2]]
    param<-all.param[-(1:dim(V)[2])]
    bY<-param[1:dim(x.reg)[2]]
    mYZ<-param[(1+dim(x.reg)[2]):(2*dim(x.reg)[2])]
    mYT<-param[(1+2*dim(x.reg)[2]):(3*dim(x.reg)[2])]
    alpha<-param[(1+3*dim(x.reg)[2]):(4*dim(x.reg)[2])]
    bD<-param[(1+4*dim(x.reg)[2]):(5*dim(x.reg)[2])]
    mDZ<-param[(1+5*dim(x.reg)[2]):(6*dim(x.reg)[2])]
    mDT<-param[(1+6*dim(x.reg)[2]):(7*dim(x.reg)[2])]

    x.reg.all<-cbind(x.reg,x.reg*df$z,x.reg*df$t)*(1-df$t*df$z)
    g.d<-as.vector(df$d-x.reg.all%*%c(bD,mDZ,mDT))*x.reg.all
    g.y<-as.vector(df$y-x.reg.all%*%c(bY,mYZ,mYT))*x.reg.all
    d_fit<-cbind(x.reg,x.reg*df$z,x.reg*df$t)%*% c(bD,mDZ,mDT)
    y_fit<-cbind(x.reg,x.reg*df$z,x.reg*df$t)%*% c(bY,mYZ,mYT)

    g.alpha<-as.vector(df$y-y_fit-x.reg %*% alpha *(df$d-d_fit))*x.reg

    g<-as.vector(x.reg%*% alpha- (V %*% psi))*V

    g<-cbind(g, g.d, g.y, g.alpha)
    return(g)
  }

  g.mean<-function(all.param,df){
    return(apply(g(all.param,df),2,mean))
  }

  n.all.param<-dim(V)[2]+dim(x.reg)[2]*7
  est<-rootSolve::multiroot(g.mean,start=rep(0.1,n.all.param),df=df)$root
  score.deriv<-numDeriv::jacobian(fun=g.mean,x=est,df=df)

  meat<-g(est,df)
  var<-diag(solve(score.deriv)%*%  t(meat) %*% meat %*%
              t(solve(score.deriv))/dim(df)[1]/dim(df)[1])

  return(list(est=est[1:dim(V)[2]],se=sqrt(var)[1:dim(V)[2]]))
}



#' IPW semiparametric estimator
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#' @param fit_z An glm object of Z using \code{X}
#' @param fit_t An glm object of T using \code{X}
#'
#' @details Here, delta_D(x) is linear in x and beta(v;psi) is linear in v. This function can be modified if one would like to use alternative nuisance parameter specifications.
#' @return A list
#' \describe{
#' \item{est}{Estimated psi.}
#' \item{se}{Standard error of \code{est}.}
#' }
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.

#' @export
#'
#' @examples
#'df<-data_gen(1e5)
#'fit_z<-glm(z~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'fit_t<-glm(t~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'X<-as.matrix(cbind(matrix(data=1,nrow = dim(df)[1],ncol=1),df[,grepl("x",names(df))]))
# constant working model
#'V<-matrix(data=1,nrow = dim(df)[1],ncol=1)
#'idid_ipw(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
# linear working model
#'V<-as.matrix(cbind(V,df[,grepl("x1",names(df))]))
#'idid_ipw(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
idid_ipw<-function(yv,dv,tv,zv,X,V,fit_z,fit_t){
  V<-as.matrix(V)
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)
  x.theta<-X
  pi_z<-fit_z$fitted.values
  pi_t<-fit_t$fitted.values

  theta_score.n<-function(theta, df, pi){
    res<-as.vector((2*df$z-1)*(2*df$t-1)*df$d/pi-x.theta %*% theta)* x.theta
    return(res)
  }

  ipw_score.n<-function(psi,df,V,delta_D,pi){
    res<-as.vector((2*df$z-1)*(2*df$t-1)/pi/delta_D*df$y- V %*% psi)* V
    return(res)
  }

  g<-function(all.param,df){
    n.psi<-dim(V)[2]
    psi<-all.param[1:n.psi]
    param<-all.param[-(1:n.psi)]
    x.model<-model.matrix(fit_z)
    gamma_z<-param[1:dim(x.model)[2]]
    gamma_t<-param[(1+dim(x.model)[2]):(2*dim(x.model)[2])]
    theta<-param[(2*dim(x.model)[2]+1):length(param)]
    pi_z<-as.vector(expit(x.model%*%gamma_z))
    pi_t<-as.vector(expit(x.model%*%gamma_t))
    pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
      ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
      (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
      ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)
    gz<-as.vector(df$z-pi_z)*x.model # n*dim(x.model)[2] matrix
    gt<-as.vector(df$t-pi_t)*x.model # n*dim(x.model)[2] matrix
    gtheta<-theta_score.n(theta,df,pi)
    delta_D<-x.theta %*% theta
    g.psi<-ipw_score.n(psi,df,V,delta_D,pi)
    g.all<-cbind(g.psi,gz,gt,gtheta)
    return(g.all)
  }

  g.mean<-function(all.param,df){
    return(apply(g(all.param,df),2,mean))
  }

  n.all.param<-dim(V)[2]+length(fit_z$coefficients)+length(fit_t$coefficients)+dim(x.theta)[2]
  est<-rootSolve::multiroot(g.mean,start=rep(0.1,n.all.param),df=df)$root
  score.deriv<-numDeriv::jacobian(fun=g.mean,x=est,df=df)

  meat<-g(est,df)
  se<-sqrt(diag(solve(score.deriv)%*%  t(meat) %*% meat %*%
                  t(solve(score.deriv))/dim(df)[1]/dim(df)[1])[1:dim(V)[2]])

  return(list(est=est[1:dim(V)[2]],se=se))
}


#' G computation semiparametric estimator
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#' @param fit_z An glm object of Z using \code{X}
#' @param fit_t An glm object of T using \code{X}
#'
#' @details Here, delta(x) is linear in x and beta(v;psi) is linear in v. This function can be modified if one would like to use alternative nuisance parameter specifications.
#' @return A list
#' \describe{
#' \item{est}{Estimated psi.}
#' \item{se}{Standard error of \code{est}.}
#' }
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @export
#'
#' @examples
#'df<-data_gen(1e5)
#'fit_z<-glm(z~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'fit_t<-glm(t~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'X<-as.matrix(cbind(matrix(data=1,nrow = dim(df)[1],ncol=1),df[,grepl("x",names(df))]))
# constant working model
#'V<-matrix(data=1,nrow = dim(df)[1],ncol=1)
#'idid_g_comp(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
# linear working model
#'V<-as.matrix(cbind(V,df[,grepl("x1",names(df))]))
#'idid_g_comp(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
idid_g_comp<-function(yv,dv,tv,zv,X,V,fit_z,fit_t){
    V<-as.matrix(V)
    df<-data.frame(y=yv,d=dv,t=tv,z=zv)
    x.alpha<-X
    pi_z<-fit_z$fitted.values
    pi_t<-fit_t$fitted.values

  alpha_score.n<-function(alpha, df, pi){
    res<-as.vector((2*df$z-1)*(2*df$t-1)/pi*
                     (df$y- x.alpha %*% alpha * df$d))* x.alpha
    return(res)
  }

  g_comp_score.n<-function(psi,df,V,delta){
    res<-as.vector( delta- (V %*% psi))*V
    return(res)
  }

  g<-function(all.param,df){
    psi<-all.param[1:dim(V)[2]]
    param<-all.param[-(1:dim(V)[2])]
    x.model<-model.matrix(fit_z)
    gamma_z<-param[1:dim(x.model)[2]]
    gamma_t<-param[(1+dim(x.model)[2]):(2*dim(x.model)[2])]
    alpha<-param[(2*dim(x.model)[2]+1):length(param)]
    pi_z<-as.vector(expit(x.model%*%gamma_z))
    pi_t<-as.vector(expit(x.model%*%gamma_t))
    pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
      ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
      (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
      ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)
    gz<-as.vector(df$z-pi_z)*x.model # n*dim(x.model)[2] matrix
    gt<-as.vector(df$t-pi_t)*x.model # n*dim(x.model)[2] matrix
    galpha<-alpha_score.n(alpha,df,pi)
    delta<-x.alpha%*% alpha
    g.psi<-g_comp_score.n(psi,df,V,delta)
    g.all<-cbind(g.psi,gz,gt,galpha)
    return(g.all)
  }

  g.mean<-function(all.param,df){
    return(apply(g(all.param,df),2,mean))
  }

  n.all.param<-dim(V)[2]+length(fit_z$coefficients)+length(fit_t$coefficients)+dim(x.alpha)[2]
  est<-rootSolve::multiroot(g.mean,start=rep(0.1,n.all.param),df=df)$root
  score.deriv<-numDeriv::jacobian(fun=g.mean,x=est,df=df)

  meat<-g(est,df)
  var<-diag(solve(score.deriv)%*%  t(meat) %*% meat %*%
              t(solve(score.deriv))/dim(df)[1]/dim(df)[1])

  return(list(est=est[1:dim(V)[2]],se=sqrt(var)[1:dim(V)[2]]))
}


#' Multiply robust semiparametric estimator
#'
#' @param yv A vector of outcome for each subject, missing is coded as NA
#' @param dv A vector of exposure for each subject, missing is coded as NA, length(dv)==length(yv)
#' @param tv A vector of time for each subject, missing is coded as NA, length(tv)==length(yv)
#' @param zv A vector of IV for DID for each subject, missing is coded as NA, length(zv)==length(yv)
#' @param X A matrix of observed covariates for each subject, missing is coded as NA, dim(X)[1]==length(yv)
#' @param V A matrix of effect modifiers for each subject, missing is coded as NA, dim(V)[1]==length(yv)
#' @param fit_z An glm object of Z using \code{X}
#' @param fit_t An glm object of T using \code{X}
#'
#' @details Here, all nuisance parameters except for \code{fit_z} adn \code{fit_t} are linear in x and beta(v;psi) is linear in v. This function can be modified if one would like to use alternative nuisance parameter specifications.
#' @return A list
#' \describe{
#' \item{est}{Estimated psi.}
#' \item{se}{Standard error of \code{est}.}
#' }
#'
#' @references Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
#'
#' @export
#'
#' @examples
#'df<-data_gen(1e5)
#'fit_z<-glm(z~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'fit_t<-glm(t~I(x1>0)+I(x2>0),data=df,family = "binomial") # one can use other models
#'X<-as.matrix(cbind(matrix(data=1,nrow = dim(df)[1],ncol=1),df[,grepl("x",names(df))]))
# constant working model
#'V<-matrix(data=1,nrow = dim(df)[1],ncol=1)
#'idid_multiR(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
# linear working model
#'V<-as.matrix(cbind(V,df[,grepl("x1",names(df))]))
#'idid_multiR(df$y,df$d,df$t,df$z,X,V,fit_z,fit_t)
idid_multiR<-function(yv,dv,tv,zv,X,V,fit_z,fit_t){
  V<-as.matrix(V)
  df<-data.frame(y=yv,d=dv,t=tv,z=zv)
  x.alpha<-X
  x.reg<-X
  x.theta<-X
  pi_z<-fit_z$fitted.values
  pi_t<-fit_t$fitted.values

  reg_score<-function(param,df){
    bY<-param[1:dim(x.reg)[2]]
    mYZ<-param[(1+dim(x.reg)[2]):(2*dim(x.reg)[2])]
    mYT<-param[(1+2*dim(x.reg)[2]):(3*dim(x.reg)[2])]
    bD<-param[(1+3*dim(x.reg)[2]):(4*dim(x.reg)[2])]
    mDZ<-param[(1+4*dim(x.reg)[2]):(5*dim(x.reg)[2])]
    mDT<-param[(1+5*dim(x.reg)[2]):(6*dim(x.reg)[2])]

    x.reg.all<-cbind(x.reg,x.reg*df$z,x.reg*df$t)*(1-df$t*df$z)
    g.d<-as.vector(df$d-x.reg.all%*%c(bD,mDZ,mDT))*x.reg.all
    g.y<-as.vector(df$y-x.reg.all%*%c(bY,mYZ,mYT))*x.reg.all
    g.reg<-cbind(g.d, g.y)
    return(g.reg)
  }

  reg_score.mean<-function(param,df){
    return(apply(reg_score(param,df),2,mean))
  }

  param.reg<-rootSolve::multiroot(reg_score.mean,start=rep(0.1,(6*dim(x.reg)[2])),df=df)$root

  pi.score<-function(param,df){
    z.model<-model.matrix(fit_z)
    gamma_z<-param[1:dim(z.model)[2]]
    gamma_t<-param[(1+dim(z.model)[2]):(2*dim(z.model)[2])]
    pi_z<-as.vector(expit(z.model%*%gamma_z))
    pi_t<-as.vector(expit(z.model%*%gamma_t))
    pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
      ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
      (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
      ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)
    gz<-as.vector(df$z-pi_z)*z.model # n*dim(z.model)[2] matrix
    gt<-as.vector(df$t-pi_t)*z.model # n*dim(z.model)[2] matrix
    return(cbind(gz,gt))
  }

  pi.score.mean<-function(param,df){
    return(apply(pi.score(param,df),2,mean))
  }

  z.model<-model.matrix(fit_z)
  gamma.hat<-rootSolve::multiroot(pi.score.mean,start=rep(0.1,2*dim(z.model)[2]),df=df)$root
  gamma_z.hat<-gamma.hat[1:dim(z.model)[2]]
  gamma_t.hat<-gamma.hat[(1+dim(z.model)[2]):(2*dim(z.model)[2])]

  pi_z<-as.vector(expit(z.model%*%gamma_z.hat))
  pi_t<-as.vector(expit(z.model%*%gamma_t.hat))
  pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
    ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
    (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
    ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)

  x.reg.all<-cbind(x.reg,x.reg*df$z,x.reg*df$t)
  y_fit<-x.reg.all%*% param.reg[1:(3*dim(x.reg)[2])]
  d_fit<-x.reg.all%*% param.reg[(1+3*dim(x.reg)[2]):(6*dim(x.reg)[2])]

  alpha_score.dr.mean<-function(alpha, df, pi, y_fit, d_fit){
    res<-as.vector((2*df$z-1)*(2*df$t-1)/pi*(df$y-y_fit-x.alpha %*% alpha *(df$d-d_fit)))* x.alpha
    return(apply(res,2,mean))
  }

  theta_score.dr.mean<-function(theta, df, pi,d_fit){
    res<-as.vector((2*df$z-1)*(2*df$t-1)/pi*(df$d-d_fit-x.theta %*% theta *df$z*df$t))* x.theta
    return(apply(res,2,mean))
  }

  alpha.hat<-rootSolve::multiroot(alpha_score.dr.mean,start=rep(0.1,dim(x.alpha)[2]),df=df,pi=pi,y_fit=y_fit,d_fit=d_fit)$root
  theta.hat<-rootSolve::multiroot(theta_score.dr.mean,start=rep(0.1,dim(x.theta)[2]),df=df,pi=pi,d_fit=d_fit)$root
  delta=x.alpha %*% alpha.hat
  delta_D=x.theta %*% theta.hat

  IF.n<-function(psi,df,V,delta, delta_D, pi, y_fit, d_fit){
    res<-as.vector(delta - (V %*% psi) +
                     (2*df$z-1)*(2*df$t-1)/pi/delta_D*( df$y-y_fit- delta*(df$d-d_fit) ))* V
    return(res)
  }

  IF<-function(psi,df,V,delta, delta_D, pi, y_fit, d_fit){
    return(apply(IF.n(psi,df,V,delta, delta_D, pi, y_fit, d_fit),2,mean))
  }

  psi.est<-rootSolve::multiroot(IF,start=rep(0,dim(V)[2]),df=df,V=V,delta=delta,delta_D=delta_D, pi=pi,y_fit=y_fit,d_fit=d_fit)$root

  est<-c(psi.est,param.reg,gamma_z.hat,gamma_t.hat,theta.hat,alpha.hat)

  g<-function(all.param,df){
    psi<-all.param[1:dim(V)[2]]
    param<-all.param[-(1:dim(V)[2])]
    param.reg<-param[1:(6*dim(x.reg)[2])]
    z.model<-model.matrix(fit_z)
    gamma_z<-param[(1+6*dim(x.reg)[2]):(6*dim(x.reg)[2]+dim(z.model)[2])]
    gamma_t<-param[(1+6*dim(x.reg)[2]+dim(z.model)[2]):(6*dim(x.reg)[2]+2*dim(z.model)[2])]
    theta<-param[(1+6*dim(x.reg)[2]+2*dim(z.model)[2]):(6*dim(x.reg)[2]+2*dim(z.model)[2]+dim(x.theta)[2])]
    alpha<-param[(1+6*dim(x.reg)[2]+2*dim(z.model)[2]+dim(x.theta)[2]):
                   (6*dim(x.reg)[2]+2*dim(z.model)[2]+dim(x.theta)[2]+dim(x.alpha)[2])]
    g.reg<-reg_score(param.reg,df) # matrix with n rows

    pi_z<-as.vector(expit(z.model%*%gamma_z))
    pi_t<-as.vector(expit(z.model%*%gamma_t))
    pi<-(pi_z*pi_t)*(df$z==1 & df$t==1)+
      ((1-pi_z)*pi_t)*(df$z==0 & df$t==1)+
      (pi_z*(1-pi_t))*(df$z==1& df$t==0)+
      ((1-pi_z)*(1-pi_t))*(df$z==0& df$t==0)
    gz<-as.vector(df$z-pi_z)*z.model # n*dim(z.model)[2] matrix
    gt<-as.vector(df$t-pi_t)*z.model # n*dim(z.model)[2] matrix

    x.reg.all<-cbind(x.reg,x.reg*df$z,x.reg*df$t)
    y_fit<-x.reg.all%*% param.reg[1:(3*dim(x.reg)[2])]
    d_fit<-x.reg.all%*% param.reg[(1+3*dim(x.reg)[2]):(6*dim(x.reg)[2])]

    gtheta<-as.vector((2*df$z-1)*(2*df$t-1)/pi*(df$d-d_fit-x.theta %*% theta *df$z*df$t))* x.theta
    galpha<-as.vector((2*df$z-1)*(2*df$t-1)/pi*(df$y-y_fit-x.alpha %*% alpha *(df$d-d_fit)))* x.alpha

    delta<-x.alpha%*% alpha
    delta_D<-x.theta%*%theta
    g.psi<-IF.n(psi,df,V,delta,delta_D,pi,y_fit,d_fit)

    g.all<-cbind(g.psi,gz,gt,gtheta,galpha,g.reg)
    return(g.all)

  }

  g.mean<-function(all.param,df){
    return(apply(g(all.param,df),2,mean))
  }

  n.all.param<-dim(V)[2]+6*dim(x.reg)[2]+2*dim(z.model)[2]+dim(x.theta)[2]+dim(x.alpha)[2]
  #est<-multiroot(g.mean,start=rep(0.1,n.all.param),df=df)$root # this also works but is slower

  score.deriv<-numDeriv::jacobian(fun=g.mean,x=est,df=df)

  meat<-g(est,df)
  var<-diag(solve(score.deriv)%*%  t(meat) %*% meat %*%
              t(solve(score.deriv))/dim(df)[1]/dim(df)[1])

  return(list(est=est[1:dim(V)[2]],se=sqrt(var)[1:dim(V)[2]]))
}



