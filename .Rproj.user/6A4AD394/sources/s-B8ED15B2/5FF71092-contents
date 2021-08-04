#' Calculate the DID parameters used to construct the brackets.
#'
#' @param df A long-form data frame with three columns (g, t, Y), where column g is the group indicator, column t is the time indicator, column Y is the observed outcome. Specifically, g takes values trt, a, b, where a and b represent the two control groups; t takes values 1, 2, ..., where the treated group adopts the treatment between t=1 and t=2.
#' @param l Estimate the treatment effect for the treated group l time periods after adopting the treatment.
#' @export
#' @examples
#'
#' df<-data_gen(1000,"case1")
#' did_list(df,l=2)
#'
did_list<-function(df,l){
  tau<-list()
  for(tmp in 0:l){
    tau_a<-(mean(df$Y[df$g=="trt" & df$t==2+tmp])-mean(df$Y[df$g=="trt" & df$t==1+tmp]))-
      (mean(df$Y[df$g=="a" & df$t==2+tmp])-mean(df$Y[df$g=="a" & df$t==1+tmp]))
    tau_b<-(mean(df$Y[df$g=="trt" & df$t==2+tmp])-mean(df$Y[df$g=="trt" & df$t==1+tmp]))-
      (mean(df$Y[df$g=="b" & df$t==2+tmp])-mean(df$Y[df$g=="b" & df$t==1+tmp]))
    tau[[tmp+1]]<-c(tau_a, tau_b)
  }
  names(tau)<-paste0("t=",2:(2+l))
  return(tau)
}

union_CI<-function(df,l,alpha=0.05,n_boot=300){
  did<-did_list(df,l)
  did_all_combn<-expand.grid(did) # number of all combinations= 2^(l+1)
  did_vec<-rowSums(did_all_combn)
  did_vec_boot<-matrix(nrow=n_boot, ncol=length(did_vec))
  for(i in 1:n_boot){
    df_boot<-df[sample.int(dim(df)[1],dim(df)[1],replace = TRUE),] # Efron bootstrap
    did_boot<-did_list(df_boot,l)
    did_vec_boot[i,]<-rowSums(expand.grid(did_boot))
  }
  did_vec_CI<-matrix(nrow=length(did_vec),ncol=2)
  did_vec_CI[,1]<-apply(did_vec_boot, 2, function(x) quantile(x, probs = alpha/2))
  did_vec_CI[,2]<-apply(did_vec_boot, 2, function(x) quantile(x, probs = 1-alpha/2))
  ci<-c(min(did_vec_CI[,1]), max(did_vec_CI[,2]))
  return(ci)
}

#' Bootstrap confidence interval for DID brackets with a union bounds form.
#' @param df A long-form data frame with three columns (g, t, Y), where column g is the group indicator, column t is the time indicator, column Y is the observed outcome.
#' Specifically, g takes values trt, a, b, where a and b represent the two control groups; t takes values 1, 2, ..., where the treated group adopts the treatment between t=1 and t=2.
#' @param l Estimate the treatment effect for the treated group l time periods after adopting the treatment.
#' @param alpha 1-alpha level confidence interval. Default is 0.05.
#' @param n_boot Number of bootstrap iterations. Default is 300.
#'
#' @details Construct bootstrap confidence interval for the average treatment effect for the treated at time 2+l (i.e., l times periods after adopting the treatment). See Ye et al., (2020) Section 4 for details.
#'
#' @return A matrix
#' \describe{
#' \item{CI lower (identified set)}{Lower end of the 1-alpha bootstrap confidence interval for the identified set.}
#' \item{CI upper (identified set)}{upper end of the 1-alpha bootstrap confidence interval for the identified set.}
#' \item{CI lower (ATT)}{lower end of the 1-alpha bootstrap confidence interval for ATT.}
#' \item{CI upper (ATT)}{upper end of the 1-alpha bootstrap confidence interval for ATT.}
#' }
#'
#' @references Ting Ye, Luke Keele, Raiden Hasegawa, Dylan S. Small (2020). A Negative Correlation Strategy for Bracketing in Difference-in-Differences with Application to the Effect of Voter Identification Laws on Voter Turnout.
#'
#' @import stats
#' @export
#'
#' @examples
#' df<-data_gen(1000,"case1")
#' bootstrap_CI(df,l=0)
#'
bootstrap_CI<-function(df,l,alpha=0.05,n_boot=300){
  did<-did_list(df,l)
  did_all_combn<-expand.grid(did) # number of all combinations= 2^(l+1)
  did_vec<-rowSums(did_all_combn)
  m<-length(did_vec)
  boot_res_min<-numeric(n_boot)
  boot_res_max<-numeric(n_boot)
  for(i in 1:n_boot){
    df_boot<-df[sample.int(dim(df)[1],dim(df)[1],replace = TRUE),] # Efron bootstrap
    did_boot<-did_list(df_boot,l)
    did_vec_boot<-rowSums(expand.grid(did_boot))
    boot_res_min[i]<-min(did_vec_boot)
    boot_res_max[i]<-max(did_vec_boot)
  }
  res<-matrix(NA, nrow = 1, ncol=6)
  boot_res_min<-boot_res_min-min(did_vec)
  boot_res_max<-boot_res_max-max(did_vec)

  res[1,1]<-min(did_vec)-sort(boot_res_min)[ceiling(n_boot*1/2)]
  res[1,2]<-max(did_vec)-sort(boot_res_max)[floor(n_boot*1/2)]
  res[1,3]<-min(did_vec)-sort(boot_res_min)[ceiling(n_boot*(1-alpha/2))]
  res[1,4]<-max(did_vec)-sort(boot_res_max)[floor(n_boot*alpha/2)]
  # confidence interval for ATT_t
  w<-res[1,2]-res[1,1]
  w_plus<-max(0,w)
  rho<-max(sort(boot_res_max)[floor(n_boot*3/4)]-(sort(boot_res_max)[floor(n_boot*1/4)]),
           sort(boot_res_min)[ceiling(n_boot*3/4)]-sort(boot_res_min)[ceiling(n_boot*1/4)])
  rho<-1/(rho*log(dim(df)[1]))
  pn<-1-pnorm(rho*w_plus)*alpha
  res[1,5]<-min(did_vec)-sort(boot_res_min)[ceiling(n_boot*(pn))]
  res[1,6]<-max(did_vec)-sort(boot_res_max)[floor(n_boot*(1-pn))]
  colnames(res)<-c("half_median_L","half_median_U","CI lower (identified set)","CI upper (identified set)",
                   "CI lower (ATT)","CI upper (ATT)")
  return(res[1,3:6])
}

#' Generate simulation datasets.
#'
#' @param N Sample size.
#' @param case Take values "case1", "case2", "S_t", used to indicate three simulation scenarios in Ye et al., (2020).
#' @param att The true treatment effects for the treated group respectively at t=1,2,3,4. The first element needs to be 0 since this corresponds to the pre-treatment period. Default is c(0,2,3,1).
#' @references Ting Ye, Luke Keele, Raiden Hasegawa, Dylan S. Small (2020). A Negative Correlation Strategy for Bracketing in Difference-in-Differences with Application to the Effect of Voter Identification Laws on Voter Turnout.
#' @export
#' @examples
#' data_gen(1000,"case1")
#'
data_gen<-function(N,case=c("case1","case2","S_t"), att=c(0,2,3,1)){
  if(case %in% c("case1","case2")){
    if(case=="case1"){ # case 1: parallel trend #
      mu_a<-c(10,11,9,8)
      mu_b<-c(4,5,3,2)
      mu_trt<-c(3,4,2,1)
    }else if (case=="case2"){ # case 2: partial parallel #
      mu_a<-c(10,11,10,11)
      mu_b<-c(4,6,2,3)
      mu_trt<-c(3,4,0,1)
    }
    mu_trt<-att+mu_trt
    n_a<-N*0.2
    n_b<-N*0.5
    n_trt<-N*0.3
    df<-data.frame(g=c(rep("trt",n_trt*4), rep("b",n_b*4), rep("a",n_a*4)),
                   t=c(rep(1:4,each=n_trt), rep(1:4, each=n_b), rep(1:4,each=n_a)))
    df$Y<-0
    for(t in 1: 4){
      df$Y[df$t==t& df$g=="trt"]<-rnorm(n_trt,mu_trt[t],1)
      df$Y[df$t==t& df$g=="a"]<-rnorm(n_a,mu_a[t],1)
      df$Y[df$t==t& df$g=="b"]<-rnorm(n_b,mu_b[t],1)
    }
  }else if (case=="S_t"){
    mu_a<-c(3,5,1,1) # 2, -4, 0
    mu_b<-c(10,12,8,11) # 2, -4, 3
    mu_trt<-c(4,6,2,3)# 2, -4, 1
    n_a<-N*0.2
    n_b<-N*0.5
    n_trt<-N*0.3
    df<-data.frame(g=c(rep("trt",n_trt*4), rep("b",n_b*4), rep("a",n_a*4)),
                   t=c(rep(1:4,each=n_trt), rep(1:4, each=n_b), rep(1:4,each=n_a)))
    df$Y<-0
    df$S<-0
    for(t in 1: 4){
      # time varying covariates
      df$S[df$t==t & df$g=="trt"]<-rnorm(n_trt, mu_trt[t], 0.5)
      df$S[df$t==t & df$g=="a"]<-rnorm(n_a, mu_a[t],0.5)
      df$S[df$t==t & df$g=="b"]<-rnorm(n_b, mu_b[t],0.5)
    }
    for(t in 1:4){
      df$Y[df$t==t]<-2-t+0.2*t*df$S[df$t==t]+rnorm(N)
      df$Y[df$g=="trt" & df$t==t]<-df$Y[df$g=="trt" & df$t==t]+att[t]
    }
  }
  return(df)
}

#' Reproduce simulation results in Section 6 (Ye et al., 2020).
#'
#' @param full_df A data frame with all the simulated datasets. Use data(df_case1), data(df_case2), data(df_S_t) to call the simulated datasets across 1000 simulation runs.
#' @param n_rep The simulation results in Section 6 are based on n_rep=1000, but n_rep can be set to a smaller number to show the results based on the first n_rep simulated datasets.
#' @param case Take values "case1", "case2", "S_t", used to indicate three simulation scenarios. Should agree with the specified dataset.
#'
#' @return A matrix with simulation results.
#' @references Ting Ye, Luke Keele, Raiden Hasegawa, Dylan S. Small (2020). A Negative Correlation Strategy for Bracketing in Difference-in-Differences with Application to the Effect of Voter Identification Laws on Voter Turnout.
#'
#' @import  stats
#' @import sandwich
#' @export
#' @examples
#'
#' data(df_case1)
#' simu(df_case1,"case1",n_rep=5) # Case I in Table 1 can be obtained by setting n_rep=1000
#' data(df_case2)
#' simu(df_case2,"case2",n_rep=5) # Case II in Table 1 can be obtained by setting n_rep=1000
#' data(df_S_t)
#' simu(df_S_t,"S_t",n_rep=5) # Table 2 can be obtained by setting n_rep=1000
#'
simu<-function(full_df,case,n_rep=1000){
  res_union0<-matrix(nrow=n_rep,ncol=2)
  res_union1<-matrix(nrow=n_rep,ncol=2)
  res_union2<-matrix(nrow=n_rep,ncol=2)
  res_boot0<-matrix(nrow=n_rep,ncol=4)
  res_boot1<-matrix(nrow=n_rep,ncol=4)
  res_boot2<-matrix(nrow=n_rep,ncol=4)
  did0<-matrix(nrow=n_rep, ncol=2)
  did1<-matrix(nrow=n_rep, ncol=2)
  did2<-matrix(nrow=n_rep, ncol=2)
  for(k in 1:n_rep){
    df<-full_df[full_df$simu==k,]
    res_boot0[k,]<-bootstrap_CI(df = df, l = 0)
    res_boot1[k,]<-bootstrap_CI(df = df, l = 1)
    res_boot2[k,]<-bootstrap_CI(df = df, l = 2)
    if(case %in% c("case1","case2")){
      res_union0[k,]<-union_CI(df = df,l = 0)
      res_union1[k,]<-union_CI(df = df,l = 1)
      res_union2[k,]<-union_CI(df = df,l = 2)
    }else if(case=="S_t"){
      df$trt<-0
      df$trt<-1*(df$t==2 & df$g=="trt")+ 2*(df$t==3 & df$g=="trt")+3*(df$t==4 & df$g=="trt")
      df$trt<-as.factor(df$trt)
      fit<-lm(Y~S+g+as.factor(t)+trt,data=df)
      v_tmp<-diag(vcovHC(fit,type="HC3"))
      did0[k,]<-c(fit$coefficients[8],sqrt(v_tmp[8]))
      did1[k,]<-c(fit$coefficients[9],sqrt(v_tmp[9]))
      did2[k,]<-c(fit$coefficients[10],sqrt(v_tmp[10]))
    }
  }
  res<-table_summary(res_boot0,res_boot1,res_boot2,res_union0,res_union1,res_union2,did0,did1,did2,case)
  return(res)
}

table_summary<-function(res_boot0,res_boot1,res_boot2,res_union0,res_union1,res_union2,did0,did1,did2,case){
  z_alpha<-qnorm(0.025,lower.tail = FALSE)
  if(case %in% c("case1","case2")){
    res<-matrix(NA, nrow=3, 6)
    rownames(res)<-c("t=2","t=3","t=4")
    colnames(res)<-c("CI_len (set)","CP (set)","CI_len (ATT)","CP (ATT)","CI_len (naive)","CP (naive)")
    res_mean_boot0<-apply(res_boot0, 2, mean)
    res[1,1]<-round(res_mean_boot0[2]-res_mean_boot0[1],3)
    res[1,2]<-mean(1*(2>res_boot0[,1] & 2<res_boot0[,2]))
    res[1,3]<-round(res_mean_boot0[4]-res_mean_boot0[3],3)
    res[1,4]<-mean(1*(2>res_boot0[,3] & 2<res_boot0[,4]))
    res_mean_boot1<-apply(res_boot1, 2, mean)
    res[2,1]<-round(res_mean_boot1[2]-res_mean_boot1[1],3)
    res[2,2]<-mean(1*(3>res_boot1[,1] & 3<res_boot1[,2]))
    res[2,3]<-round(res_mean_boot1[4]-res_mean_boot1[3],3)
    res[2,4]<-mean(1*(3>res_boot1[,3] & 3<res_boot1[,4]))
    res_mean_boot2<-apply(res_boot2, 2, mean)
    res[3,1]<-round(res_mean_boot2[2]-res_mean_boot2[1],3)
    res[3,2]<-mean(1*(1>res_boot2[,1] & 1<res_boot2[,2]))
    res[3,3]<-round(res_mean_boot2[4]-res_mean_boot2[3],3)
    res[3,4]<-mean(1*(1>res_boot2[,3] & 1<res_boot2[,4]))
    #union
    res_mean_union0<-apply(res_union0, 2, mean)
    res[1,5]<-round(res_mean_union0[2]-res_mean_union0[1],3)
    res[1,6]<-mean(1*(2>res_union0[,1] & 2<res_union0[,2]))
    res_mean_union1<-apply(res_union1, 2, mean)
    res[2,5]<-round(res_mean_union1[2]-res_mean_union1[1],3)
    res[2,6]<-mean(1*(3>res_union1[,1] & 3<res_union1[,2]))
    res_mean_union2<-apply(res_union2, 2, mean)
    res[3,5]<-round(res_mean_union2[2]-res_mean_union2[1],3)
    res[3,6]<-mean(1*(1>res_union2[,1] & 1<res_union2[,2]))
  } else if(case=="S_t"){
    res<-matrix(NA, nrow=3, ncol=12)
    rownames(res)<-c("t=2","t=3","t=4")
    colnames(res)<-c("CI_lower (set)","CI_upper (set)", "CP (set)",
                     "CI_lower (ATT)","CI_upper (ATT)", "CP (ATT)",
                     "mean (fixed effect)","SD (fixed effect)",
                     "SE (fixed effect)", "CI_lower (fixed effect)","CI_upper (fixed effect)",
                     "CP (fixed effect")
    res[1,1:2]<-round(apply(res_boot0[,1:2], 2, mean),3)
    res[1,3]<-mean(1*(2>res_boot0[,1] & 2<res_boot0[,2]))
    res[1,4:5]<-round(apply(res_boot0[,3:4], 2, mean),3)
    res[1,6]<-mean(1*(2>res_boot0[,3] & 2<res_boot0[,4]))
    res[2,1:2]<-round(apply(res_boot1[,1:2], 2, mean),3)
    res[2,3]<-mean(1*(3>res_boot1[,1] & 3<res_boot1[,2]))
    res[2,4:5]<-round(apply(res_boot1[,3:4], 2, mean),3)
    res[2,6]<-mean(1*(3>res_boot1[,3] & 3<res_boot1[,4]))
    res[3,1:2]<-round(apply(res_boot2[,1:2], 2, mean),3)
    res[3,3]<-mean(1*(1>res_boot2[,1] & 1<res_boot2[,2]))
    res[3,4:5]<-round(apply(res_boot2[,3:4], 2, mean),3)
    res[3,6]<-mean(1*(1>res_boot2[,3] & 1<res_boot2[,4]))
    # linear fixed effect
    res[1,7]<-round(mean(did0[,1]),3)
    res[1,8]<-round(sd(did0[,1]),3)
    res[1,9]<-round(mean(did0[,2]),3)
    res[1,10:11]<-round(c(mean(did0[,1])-z_alpha*mean(did0[,2]),
                          mean(did0[,1])+z_alpha*mean(did0[,2]) ),3)
    cover_ind<-1*(2>(did0[,1]-z_alpha*did0[,2]) & 2<(did0[,1]+z_alpha*did0[,2]) )
    res[1,12]<-round(mean(cover_ind),3)
    res[2,7]<-round(mean(did1[,1]),3)
    res[2,8]<-round(sd(did1[,1]),3)
    res[2,9]<-round(mean(did1[,2]),3)
    res[2,10:11]<-round(c(mean(did1[,1])-z_alpha*mean(did1[,2]),
                          mean(did1[,1])+z_alpha*mean(did1[,2]) ),3)
    cover_ind<-1*(3>(did1[,1]-z_alpha*did1[,2]) & 3<(did1[,1]+z_alpha*did1[,2]) )
    res[2,12]<-round(mean(cover_ind),3)
    res[3,7]<-round(mean(did2[,1]),3)
    res[3,8]<-round(sd(did2[,1]),3)
    res[3,9]<-round(mean(did2[,2]),3)
    res[3,10:11]<-round(c(mean(did2[,1])-z_alpha*mean(did2[,2]),
                          mean(did2[,1])+z_alpha*mean(did2[,2]) ),3)
    cover_ind<-1*(1>(did2[,1]-z_alpha*did2[,2]) & 1<(did2[,1]+z_alpha*did2[,2]) )
    res[3,12]<-round(mean(cover_ind),3)
  }
  return(res)
}


