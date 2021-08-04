
#' Generate simulation datasets.
#'
#' @param I Number of the broad case matched sets.
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @return A data frame.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#' data_gen(I=1e3, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15)
#'
data_gen<-function(I, J, pi, bT, bC, etaT, etaC){
  dat<-data.frame(id=1:(I*J), matchvec=rep(1:I,each=J), broad.case=rep(c(1,rep(0,J-1)), I))
  pZ1<-bT*pi/(bT*pi+bC*(1-pi))
  pZ0<-(1-bT)*pi/((1-bT)*pi+(1-bC)*(1-pi))
  dat$z<-rbinom(I*J,size = 1,prob = pZ1*dat$broad.case + pZ0*(1-dat$broad.case))
  dat$narrow.case<-rbinom(I*J,size = 1,prob =
                            etaT*dat$broad.case*dat$z+ etaC*dat$broad.case*(1-dat$z) + 0*(1-dat$broad.case))
  return(dat)
}


#' Design sensitivity formulas for the broad case and narrow case tests when data is generated as described in Section 4.2
#'
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#'
#' @return A matrix
#' \describe{
#' \item{broad.case}{The design sensitivity for the broad case test.}
#' \item{narrow.case}{The design sensitivity for the narrow case test.}
#' }
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @export
#' @examples
#' design_sensitivity(J=6,pi=1/3,bT=0.3,bC=0.1,etaT=0.3,etaC=0.15)
#'
#'
design_sensitivity<-function(J, pi, bT, bC, etaT, etaC){
  res<-matrix(nrow=1,ncol=2)
  colnames(res)<-c("broad.case","narrow.case")
  res[1,1]<-bT/(1-bT)/bC*(1-bC)
  res[1,2]<-bT/(1-bT)/bC*(1-bC)*etaT/etaC
  return(res)
}


#' Sensitivity Analysis for Mantel-Haenszel Broad Case Test
#'
#' @param broad.case A vector of broad case status with no missing data.
#' @param matchvec Matched set indicator, 1, 2, ..., sum(broad.case) with length(matchvec)==length(broad.case). Matched set indicators should be either integers or a factor.
#' @param z Treatment indicator, z=1 for treated, z=0 for control with length(z)==length(broad.case).
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param two.sided One-sided or two-sided test. Default is FALSE.
#'
#' @details Sensitivity analysis for the Mantel-Haenszel broad case test \eqn{Y_b} of no treatment effect based on the broad case definition for matched case-control studies.
#'
#' @return Upper bound on the p-value for all distributions of treatment assignment consistent with the sensitivity parameter Gamma.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @export
#'
#' @examples
#'
#' dat<-data_gen(I=1e3, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15)
#' broadcase.p.value(dat$broad.case,dat$matchvec,dat$z,Gamma=3)
#'
broadcase.p.value<-function(broad.case,matchvec,z,Gamma,two.sided=FALSE){
  df<-data.frame(broad.case=broad.case,z=z,matchvec=matchvec)
  Y_b<-sum(subset(df,broad.case==1)$z)
  m<-numeric(length(unique(df$matchvec)))
  for(i in 1:dim(df)[1]){
    m[i]<-sum(df$z[df$matchvec==i])
  }
  gdf<-data.frame(matchvec=1:length(unique(df$matchvec)),m=m)
  J<-table(df$matchvec)[1]
  p2<-gdf$m*Gamma/(gdf$m*Gamma+J-gdf$m)
  p.value2<-pnorm((Y_b-sum(p2))/sqrt(sum(p2*(1-p2))),lower.tail = FALSE)
  p1<-gdf$m/(gdf$m+(J-gdf$m)*Gamma)
  p.value1<-pnorm((Y_b-sum(p1))/sqrt(sum(p1*(1-p1))))
  if(two.sided==TRUE){
    p.value<-min(2*(min(p.value1,p.value2)),1)
  }else{
    p.value<-p.value2
  }
  return(p.value)
}


#' Sensitivity Analysis for Mantel-Haenszel Narrow Case Test
#'
#' @param narrow.case A vector of narrow case status with no missing data.
#' @param matchvec Matched set indicator, 1, 2, ..., sum(narrow.case) with length(matchvec)==length(narrow.case). Matched set indicators should be either integers or a factor.
#' @param z Treatment indicator, z=1 for treated, z=0 for control with length(z)==length(narrow.case).
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param Theta The sensitivity parameter \eqn{\Theta} at which the test is conducted, where \eqn{\Theta\ge 1}. \eqn{1-\Theta} can be interpreted as the proportion of narrow cases when treated changing case definition when untreated. Setting \eqn{\Theta = 1} means that the treatment on average does not change case definition among always-cases.
#' @param two.sided One-sided or two-sided test. Default is FALSE.
#'
#' @details Sensitivity analysis for the Mantel-Haenszel broad case test \eqn{Y_b} of no treatment effect based on the broad case definition for matched case-control studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @return Upper bound on the p-value for all distributions of treatment assignment consistent with the sensitivity parameter Gamma.
#'
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @export
#'
#' @examples
#'
#' dat<-data_gen(I=1e3, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15)
#' narrow.case.matched.index<-unique(subset(dat,narrow.case==1)$matchvec)
#' dat.narrow<-subset(dat,matchvec %in% narrow.case.matched.index)
#' narrowcase.p.value(dat.narrow$narrow.case,dat.narrow$matchvec,dat.narrow$z,Gamma=3,Theta=1.5)
#'
narrowcase.p.value<-function(narrow.case,matchvec,z,Gamma,Theta,two.sided=FALSE){
  df<-data.frame(narrow.case=narrow.case,z=z,matchvec=matchvec)
  Y_n<-sum(subset(df,narrow.case==1)$z)
  m<-numeric(length(unique(df$matchvec)))
  for(i in 1:dim(df)[1]){
    m[i]<-sum(df$z[df$matchvec==i])
  }
  gdf<-data.frame(matchvec=1:length(unique(df$matchvec)),m=m)
  J<-table(df$matchvec)[1]
  p2<-gdf$m*Gamma*Theta/(gdf$m*Gamma*Theta+J-gdf$m)
  p.value2<-pnorm((Y_n-sum(p2))/sqrt(sum(p2*(1-p2))),lower.tail = FALSE)
  p1<-gdf$m/(gdf$m+(J-gdf$m)*Gamma)
  p.value1<-pnorm((Y_n-sum(p1))/sqrt(sum(p1*(1-p1))))
  if(two.sided==TRUE){
    p.value<-min(2*(min(p.value1,p.value2)),1)
  }else{
    p.value<-p.value2
  }
  return(p.value)
}

#' Simulated Power of Sensitivity Analysis
#'
#' @param I Number of the broad case matched sets.
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @param n_sim Number of simulation repetition.
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param Theta The sensitivity parameter \eqn{\Theta} at which the test is conducted, where \eqn{\Theta\ge 1}. \eqn{1-\Theta} can be interpreted as the proportion of narrow cases when treated changing case definition when untreated. Setting \eqn{\Theta = 1} means that the treatment on average does not change case definition among always-cases.
#' @param alpha Significance level, usually 0.05.
#' @param two.sided One-sided or two-sided test. Default is FALSE.
#' @return A matrix with simulation results.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#'
#' sim_power(I=1e3, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3,
#' etaC=0.15, n_sim=10, Gamma=3, Theta=1.5, alpha=0.05)
#'
sim_power<-function(I, J, pi, bT, bC, etaT, etaC, n_sim, Gamma, Theta, alpha, two.sided=FALSE){
  pval_res<-matrix(nrow=n_sim,ncol=3)
  colnames(pval_res)<-c("broad.case","narrow.case","bonferroni")
  for(i in 1:n_sim){
    dat<-data_gen(I, J, pi, bT, bC, etaT, etaC)
    pval_res[i,1]<-broadcase.p.value(dat$broad.case,dat$matchvec,dat$z,Gamma)
    narrow.case.matched.index<-unique(subset(dat,narrow.case==1)$matchvec)
    dat.narrow<-subset(dat,matchvec %in% narrow.case.matched.index)
    pval_res[i,2]<-narrowcase.p.value(dat.narrow$narrow.case,dat.narrow$matchvec,dat.narrow$z,Gamma=3,Theta=1.5)
    pval_res[i,3]<-min(pval_res[i,1],pval_res[i,2])*2
  }
  res<-apply(pval_res,2,function(x) length(which(x<alpha))/n_sim)
  return(res)
}

#' Calculate Power of Sensitivity Analysis (by formula)
#'
#' @param I Number of the broad case matched sets.
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param Theta The sensitivity parameter \eqn{\Theta} at which the test is conducted, where \eqn{\Theta\ge 1}. \eqn{1-\Theta} can be interpreted as the proportion of narrow cases when treated changing case definition when untreated. Setting \eqn{\Theta = 1} means that the treatment on average does not change case definition among always-cases.
#' @param alpha Significance level, usually 0.05.
#' @return A matrix with calculated power results.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#'
#' cal_power(I=1e3, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15, Gamma=3, Theta=1.5, alpha=0.05)
#'
cal_power<-function(I, J, pi, bT, bC, etaT, etaC, Gamma, Theta, alpha){
  res<-matrix(nrow=1,ncol=2)
  colnames(res)<-c("broad.case","narrow.case")
  #broad.case
  m_prob<-matrix(nrow=J+1,ncol=2)
  m_prob[,1]<-0:J
  m_prob[1,2]<-(1-pi)^J*bC*(1-bC)^(J-1)/(bT*pi+bC*(1-pi))/((1-bT)*pi+(1-bC)*(1-pi))^(J-1)
  for(t in 1:J){
    m_prob[t+1,2]<-pi^t*(1-pi)^(J-t)*(bT*choose(J-1,t-1)*(1-bT)^(t-1)*(1-bC)^(J-t)+
                                        bC*choose(J-1,t)*(1-bT)^t*(1-bC)^(J-1-t))/(bT*pi+bC*(1-pi))/
      ((1-bT)*pi+(1-bC)*(1-pi))^(J-1)
  }
  m_exp<-sum(m_prob[,1]*Gamma/(m_prob[,1]*Gamma+ J-m_prob[,1])*m_prob[,2])
  m_var<-sum(m_prob[,1]*Gamma*(J-m_prob[,1])/(m_prob[,1]*Gamma+ J-m_prob[,1])^2*m_prob[,2])
  tmp<-(m_exp-bT*pi/(bT*pi+bC*(1-pi))+qnorm(1-alpha)/sqrt(I)*sqrt(m_var))/sqrt(bT*pi*bC*(1-pi)/(bT*pi+bC*(1-pi))^2/I)
  res[1,1]<-1-pnorm(tmp)
  #narrow.case
  I.narrow<-I*(etaT*bT*pi/(bT*pi+bC*(1-pi)) + etaC*bC*(1-pi)/(bT*pi+bC*(1-pi)))
  m_prob<-matrix(nrow=J+1,ncol=2)
  m_prob[,1]<-0:J
  m_prob[1,2]<-(1-pi)^J*bC*etaC*(1-bC)^(J-1)/(bT*etaT*pi+bC*etaC*(1-pi))/((1-bT)*pi+(1-bC)*(1-pi))^(J-1)
  for(t in 1:J){
    m_prob[t+1,2]<-pi^t*(1-pi)^(J-t)*(bT*etaT*choose(J-1,t-1)*(1-bT)^(t-1)*(1-bC)^(J-t)+
                                        bC*etaC*choose(J-1,t)*(1-bT)^t*(1-bC)^(J-1-t))/(bT*etaT*pi+bC*etaC*(1-pi))/
      ((1-bT)*pi+(1-bC)*(1-pi))^(J-1)
  }
  m_exp<-sum(m_prob[,1]*Gamma*Theta/(m_prob[,1]*Gamma*Theta+ J-m_prob[,1])*m_prob[,2])
  m_var<-sum(m_prob[,1]*Gamma*Theta*(J-m_prob[,1])/(m_prob[,1]*Gamma*Theta+ J-m_prob[,1])^2*m_prob[,2])
  tmp<-(m_exp-bT*etaT*pi/(bT*etaT*pi+bC*etaC*(1-pi))+qnorm(1-alpha)/sqrt(I.narrow)*sqrt(m_var))/
    sqrt(bT*etaT*pi*bC*etaC*(1-pi)/(bT*etaT*pi+bC*etaC*(1-pi))^2/I.narrow)
  res[1,2]<-1-pnorm(tmp)
  return(res)
}

#' Calculate the Number of Broad Case Matched Sets to Achieve the Target Power Using the Broad Case Test
#'
#' @param pw Target power.
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param alpha Significance level, usually 0.05.
#' @param upper Upper bound on \eqn{I}. If the calculated power is below pw at the upper bound, output Inf.
#' @return The number of broad case matched sets needed.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#'
#' cal_size_b(pw=0.8, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15, Gamma=3, alpha=0.05)
#'
cal_size_b<-function(pw, J, pi, bT, bC, etaT, etaC, Gamma, alpha, upper=1e7){
  if(cal_power(upper, J, pi, bT, bC, etaT, etaC, Gamma, Theta=1, alpha)[1,1]<pw){
    res<-Inf
  }else{
    res<-round(uniroot(function(I){cal_power(I,J,pi,bT,bC,etaT,etaC,Gamma,Theta=1,alpha)[1,1]-pw} ,c(1,upper))$root)
  }
  return(res)
}

#' Calculate the Expected Number of Narrow Case Matched Sets to Achieve the Target Power Using the Narrow Case Test
#'
#' @param pw Target power.
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @param Gamma The sensitivity parameter \eqn{\Gamma} at which the test is conducted, where \eqn{\Gamma\ge 1}. Setting \eqn{\Gamma = 1} is equivalent to assuming exchangeability given the matched sets, and it performs a within-set randomization test.
#' @param Theta The sensitivity parameter \eqn{\Theta} at which the test is conducted, where \eqn{\Theta\ge 1}. \eqn{1-\Theta} can be interpreted as the proportion of narrow cases when treated changing case definition when untreated. Setting \eqn{\Theta = 1} means that the treatment on average does not change case definition among always-cases.
#' @param alpha Significance level, usually 0.05.
#' @param upper Upper bound on \eqn{I}. If the calculated power is below pw at the upper bound, output Inf.
#' @return The expected number of narrow case matched sets needed.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#'
#' cal_size_n(pw=0.8, J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15, Gamma=3, Theta=1.3, alpha=0.05)
#'
cal_size_n<-function(pw, J, pi, bT, bC, etaT, etaC, Gamma, Theta, alpha, upper){
  if(cal_power(1e7, J, pi, bT, bC, etaT, etaC, Gamma, Theta, alpha)[1,2]<pw){
    res<-Inf
  }else{
    res<-round(uniroot(function(I){cal_power(I,J,pi,bT,bC,etaT,etaC,Gamma,Theta,alpha)[1,2]-pw} ,c(1,1e7))$root)
  }
  res<-res*(etaT*bT*pi/(bT*pi+bC*(1-pi))+etaC*bC*(1-pi)/(bT*pi+bC*(1-pi)))
  return(res)
}


#' Plot Power Against the Number of Broad Case Matched Sets
#'
#' @param J Number of subjects in each matched sets.
#' @param pi Probability of receiving treatment.
#' @param bT Probability of being a broad case when treated.
#' @param bC Probability of being a broad case when not treated.
#' @param etaT Probability of a broad case being a narrow case when treated.
#' @param etaC Probability of a broad case being a narrow case when not treated.
#' @param alpha Significance level, usually 0.05.
#' @param n_sim Number of simulation repetition.
#' @param two.sided One-sided or two-sided test. Default is FALSE.
#' @param method Whether to only include the formula power or both the formula and simulation power.
#' @return A power plot saved at the working directory.
#' @references Ting Ye and Dylan S. Small (2021). Combining Broad and Narrow Case Definitions in Matched Case-Control Studies.
#' @importFrom reshape2 melt
#' @importFrom  ggplot2 ggplot aes geom_line facet_grid ylab xlab theme_bw scale_linetype_manual scale_color_manual ggsave theme label_parsed
#' @importFrom stats pnorm qnorm rbinom uniroot
#' @importFrom utils write.csv
#' @export
#' @examples
#'
#' plot_power_byI(J=6, pi=1/3, bT=0.3, bC=0.1, etaT=0.3, etaC=0.15, alpha=0.05, method="formula")
#'
plot_power_byI<-function(J, pi, bT, bC, etaT, etaC, alpha, n_sim=NULL,
                     two.sided=FALSE,method=c("formula","both")){
  ds.broad<-design_sensitivity(J,pi,bT,bC,etaT,etaC)[1,1]
  ds.narrow<-design_sensitivity(J,pi,bT,bC,etaT,etaC)[1,2]
  Theta_seq<-c(1,1.5,2)
  Gamma_seq<-c(3,3.5,4)
  B_seq<-round(seq(10,1000,length.out = 20))
  if(method=="formula"){
    res.broad<-matrix(nrow=length(Gamma_seq)*length(Theta_seq)*length(B_seq),ncol=5)
    colnames(res.broad)<-c("Method","Theta","Gamma","I","value")
    res.narrow<-matrix(nrow=length(Gamma_seq)*length(Theta_seq)*length(B_seq),ncol=5)
    colnames(res.narrow)<-c("Method","Theta","Gamma","I","value")
    res.broad<-as.data.frame(res.broad)
    res.narrow<-as.data.frame(res.narrow)
    res.broad[,1]<-"Broad case (formula)"
    res.narrow[,1]<-"Narrow case (formula)"
    b<-1
    for(i in 1:length(Theta_seq)){
      for(j in 1:length(Gamma_seq)){
        for(k in 1:length(B_seq)){
          res.broad[b,2:4]<-c(Theta_seq[i],Gamma_seq[j],B_seq[k])
          res.narrow[b,2:4]<-c(Theta_seq[i],Gamma_seq[j],B_seq[k])
          res.broad[b,5]<-cal_power(B_seq[k],J, pi, bT, bC, etaT, etaC, Gamma_seq[j], Theta_seq[i], alpha)[1,1]
          res.narrow[b,5]<-cal_power(B_seq[k],J, pi, bT, bC, etaT, etaC, Gamma_seq[j], Theta_seq[i], alpha)[1,2]
          b<-b+1
        }
      }
    }
    gdf<-rbind(res.broad,res.narrow)
    gdf$Theta<-as.factor(gdf$Theta)
    levels(gdf$Theta)<-c(expression(paste(Theta, " = 1")), expression(paste(Theta, " = 1.5")),
                         expression(paste(Theta, " = 2")))
    gdf$Gamma<-as.factor(gdf$Gamma)
    levels(gdf$Gamma)<-c(expression(paste(Gamma, " = 3")),expression(paste(Gamma, " = 3.5")),
                         expression(paste(Gamma, " = 4")))
    ggplot(gdf,aes(x=I,y=value,group=Method))+ geom_line(aes(linetype=Method,color=Method))+
      facet_grid(Gamma ~ Theta,labeller = label_parsed)+
      ylab("Power")+
      xlab("Number of matched sets")+
      theme_bw()+
      theme(legend.position="top")+
      scale_linetype_manual(values=c("solid","dashed"))+
      scale_color_manual(values=c("black",'black'))
    ggsave(paste0("byB_Method",method,"J",J,"pi",pi,"bT",bT,"bC",bC,"etaT",
                  etaT,"etaC",etaC,"n_sim",n_sim,".pdf"),width=6,height=5)
  }else if (method=="both"){
    res<-matrix(nrow=length(Gamma_seq)*length(Theta_seq)*length(B_seq),ncol=8)
    colnames(res)<-c("Theta","Gamma","I","Broad case (formula)","Narrow case (formula)",
                     "Broad case (simulation)","Narrow case (simulation)","Bonferroni (simulation)")
    b<-1
    for(i in 1:length(Theta_seq)){
      for(j in 1:length(Gamma_seq)){
        for(k in 1:length(B_seq)){
          res[b,1:3]<-c(Theta_seq[i],Gamma_seq[j],B_seq[k])
          res[b,4]<-cal_power(B_seq[k],J, pi, bT, bC, etaT, etaC, Gamma_seq[j], Theta_seq[i], alpha)[1,1]
          res[b,5]<-cal_power(B_seq[k],J, pi,
                              bT, bC, etaT, etaC, Gamma_seq[j], Theta_seq[i], alpha)[1,2]
          res[b,c(6,8)]<-sim_power(round(B_seq[k]),J, pi, bT, bC,etaT, etaC, n_sim,Gamma_seq[j],Theta_seq[i], alpha)[c(1,3)]
          res[b,7]<-sim_power(round(B_seq[k]),J, pi,
                              bT, bC,etaT, etaC, n_sim,Gamma_seq[j],Theta_seq[i], alpha)[2]
          b<-b+1
        }
      }
    }
    write.csv(res,paste0("by_B_","J",J,"bT",bT,"bC",bC,"etaT",
                         etaT,"etaC",etaC,"n_sim",n_sim,".csv"))
    res<-data.frame(res)
    res$id<-1:dim(res)[1]
    gdf<- melt(res,id.vars=c("id","Gamma","Theta","I"))
    gdf$Theta<-as.factor(gdf$Theta)
    levels(gdf$Theta)<-c(expression(paste(Theta, " = 1")), expression(paste(Theta, " = 1.5")),expression(paste(Theta, " = 2")))
    gdf$Gamma<-as.factor(gdf$Gamma)
    levels(gdf$Gamma)<-c(expression(paste(Gamma, " = 3")),expression(paste(Gamma, " = 3.5")), expression(paste(Gamma, " = 4")))
    names(gdf)[names(gdf)=="variable"]<-"Method"
    levels(gdf$Method)<-c("Broad case (formula)", "Narrow case (formula)",
                          "Broad case (simulation)","Narrow case (simulation)","Combined (simulation)")
    ggplot(gdf,aes(x=I,y=value,group=Method))+ geom_line(aes(linetype=Method,color=Method))+
      facet_grid(Gamma ~ Theta,labeller = label_parsed)+
      ylab("Power")+
      xlab("Number of broad case matched sets")+
      theme_bw()+
      scale_linetype_manual(values=c("solid","dashed","solid","dashed","dotted"))+
      scale_color_manual(values=c("black",'black',"#999999","#999999","#999999"))
    ggsave(paste0("byB_Method",method,"J",J,"bT",bT,"bC",bC,"etaT",
                  etaT,"etaC",etaC,"n_sim",n_sim,".pdf"),width=8,height=5)

    ggplot(subset(gdf,Method%in%c("Broad case (formula)", "Narrow case (formula)",
                                  "Combined (simulation)")),aes(x=I,y=value,group=Method))+ geom_line(aes(linetype=Method,color=Method))+
      facet_grid(Gamma ~ Theta,labeller = label_parsed)+
      ylab("Power")+
      xlab("Number of broad case matched sets")+
      theme_bw()+
      theme(legend.position="top")+
      scale_linetype_manual(values=c("solid","dashed","dotted"))+
      scale_color_manual(values=c("black",'black',"#999999"))
    ggsave(paste0("byB_Method","formula","J",J,"bT",bT,"bC",bC,"etaT",
                  etaT,"etaC",etaC,"n_sim",n_sim,".pdf"),width=6,height=5)

      }
}

