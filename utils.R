# para_python <- input_priors_csv('/Users/kaihangzhao/Documents/GitHub/targeting-sample-size/simulated_priors.csv')
# s_c <- para_python$s_c
# w_c <- para_python$w_c
# mu_c <- para_python$mu_c
# sigma_c <- para_python$sigma_c
# sanity_check(para_python)
# 
# DELTA <- 1.5 # Expected Performance Requirement (see Section 3)
# GAMMA <- 0.3 # Probability Requirement (see Section 3.5)
# 
# ALPHA <- 0.05 # (1-ALPHA) indicates confidence in a statistical test (see Section 4)
# BETA <- 0.3 # (1-BETA) indicates power in a statistical test (see Section 4)
# N_MAX <- 100000 # Maximum experiment size (Algorithms 2-5)
# B <- 10000 # Simulation iterations (Algorithms 2 and 4)
# N_PARAL <- 1000 # Values of Ntr for parallel evaluation (memory constraint: B x N_PARAL x C)

input_priors_csv<- function(input_priors){
 
  C <- dim(input_priors)[1]
  s_c <- input_priors$s_c
  w_c <- input_priors$w_c
  mu_c <- input_priors$mu_c
  sigma_c <- input_priors$sigma_c
  stopifnot((1-sum(w_c))<1e-3, sum(w_c>=0) == C)
  print(paste("Number of Segments = ", C))
  print(paste("Average Treatment Effect = ", round(sum(w_c*mu_c), 2)))
  output <- list(w_c = w_c, mu_c = mu_c, sigma_c = sigma_c, s_c = s_c)
  return(output)
}

sanity_check <- function(data){
  s_c <- na.omit(data$s_c)
  w_c <- na.omit(data$w_c)
  mu_c <- na.omit(data$mu_c)
  sigma_c <- na.omit(data$sigma_c)
  if (length(s_c)!= length(w_c) | length(s_c)!=length(mu_c) | length(s_c)!= length(sigma_c)){
    return('Length is the not the same for all parameters after removing null values, please check your dataset uploaded.')
  }
  if (abs(sum(w_c)-1)>0.01){
    return('The sum of w_c is not equal to 1, please check your dataset uploaded.')
  }
  if (sum(w_c<0)>0) {
    return('There is at least one negative value in w_c, please check your dataset uploaded.')
  }
  if (sum(s_c<0)>0) {
    return('There is at least one negative value in s_c, please check your dataset uploaded.')
  }
  if (sum(sigma_c<0)>0) {
    return('There is at least one negative value in sigma_c, please check your dataset uploaded.')
  }
  
}

get_expected_perf_limit<- function(mu_c, sigma_c, s_c, w_c){
  pdf_c <- dnorm(mu_c/sigma_c)
  cdf_c <-pnorm(mu_c/sigma_c)  
  mu_limit <- sum(w_c*(sigma_c*pdf_c+mu_c*cdf_c))
  return(mu_limit)
}

get_expected_perf<- function(mu_c,sigma_c,s_c,w_c,n_tr) {
  if (n_tr==0){
    mu_delta <- sum(w_c*mu_c*(1*(mu_c>0)))
    var_delta <- 0
  }
  else {
    k_c <- sqrt((1+4*s_c**2/w_c/n_tr/(sigma_c**2)))
    # dnorm is pdf
    pdf_c <- dnorm(mu_c*k_c/sigma_c)
    # pnorm is cdf
    cdf_c <-pnorm(mu_c*k_c/sigma_c)  
    mu_delta <-sum(w_c*(sigma_c/k_c*pdf_c+mu_c*cdf_c))
    var_delta<-sum((w_c**2*(mu_c**2*cdf_c*(1-cdf_c)+sigma_c**2/k_c**2*(cdf_c-pdf_c**2)+mu_c*sigma_c/k_c*pdf_c*(1-2*cdf_c))))
  }
  output <- list ('mu_delta'= mu_delta, 'var_delta'= var_delta)
  return(output)
}

solve_d_expected_improvement_search<- function(delta,mu_c,sigma_c,s_c,w_c,N_MIN,N_MAX){
  
  if ((N_MAX%/%1) == (N_MIN%/%1)){
    return(N_MIN)
  }
  N_MID <- (N_MIN+N_MAX)/2
  mid_val <- get_expected_perf(mu_c,sigma_c,s_c,w_c,N_MID)$mu_delta-delta
  print(mid_val)
  if (mid_val<=0){
    return(solve_d_expected_improvement_search(delta,mu_c,sigma_c,s_c,w_c,N_MID,N_MAX))
  }
  else{
    return(solve_d_expected_improvement_search(delta,mu_c,sigma_c,s_c,w_c,N_MIN,N_MID))
  }
}

solve_d_expected_improvement<- function(delta, mu_c, sigma_c, s_c, w_c){
  # Algorithm 1
  output_0 <- get_expected_perf(mu_c,sigma_c,s_c,w_c,0)
  mu_delta_0 <- output_0$mu_delta
  
  if(mu_delta_0 >= delta){return(0)}
  
  mu_delta_max <- get_expected_perf_limit(mu_c, sigma_c, s_c, w_c)
  if(mu_delta_max<=delta){return(-1)}
  
  step_ <- 1e7
  ntr_curr <- step_
  while(get_expected_perf(mu_c, sigma_c, s_c, w_c, ntr_curr)$mu_delta<=delta){
    ntr_curr <- ntr_curr + step_
  }
  
  n_opt <- solve_d_expected_improvement_search(delta,mu_c,sigma_c,s_c,w_c, ntr_curr-step_, ntr_curr)
  return (ceiling(n_opt))
  
}

chunker<- function(se, size){
  output<- c()
  for (pos in seq(0,length(se), size)){
    #output<- se[pos:(pos+size-1)]
    output<- c(output, list(se[pos:(pos+size-1)]))
  }
  return(output)
}

stats_g_prob_d_expected_improvement<- function(delta, gamma, mu_c, sigma_c, s_c, w_c, ntr_set, B){
  if(B==0){
    k_c <- sqrt(1+4*s_c**2/outer(w_c,(ntr_set+1e-30))/(sigma_c**2))
    pdf_c <- dnorm(mu_c*k_c/sigma_c)
    cdf_c <- pnorm(mu_c*k_c/sigma_c)    
    mu_delta<- apply(w_c*(sigma_c/k_c*pdf_c+mu_c*cdf_c), 2, sum)
    var_delta<-apply(w_c**2*(mu_c**2*cdf_c*(1-cdf_c)+sigma_c**2/k_c**2*(cdf_c-pdf_c**2)+mu_c*sigma_c/k_c*pdf_c*(1-2*cdf_c)),2,sum)
    prob <- pnorm((mu_delta-delta)/sqrt(var_delta))
    quant <- mu_delta - qnorm(1-gamma)*sqrt(var_delta)
  }
  else{
    C <- length(mu_c)
    N <- length(ntr_set)
    mu_c_post_var <-sigma_c**4/(sigma_c**2+4*s_c**2/outer(w_c, (ntr_set+1e-30)))
    mu_c_post_sample <- array(rnorm(C*N*B, mean = mu_c, sd = sqrt(mu_c_post_var)), dim = c(C, N, B))
    
    
    v_post<- colSums(w_c*mu_c_post_sample*(1*(mu_c_post_sample>=0)))
    prob <- apply(1*(v_post>=delta), 1, mean)
    quant <- apply(v_post, 1, quantile , probs = gamma)
  }
  output<- list('prob'=prob, 'quant'=quant)
  return(output)
}

solve_g_prob_d_expected_improvement<- function(delta, gamma, mu_c, sigma_c, s_c, w_c, n_max, parallel=1000000, B=0){
  # B>0: Algorithm 4 (simulation); requires Parallel & B to avoid memory overflow
  # B=0: Algorithm 5 (analytical approximation); Parallel & B are irrelevants
  
  ntr_full_set <- seq(0, n_max, 1)
  ntr_opt <- -1
  for (ntr_group in chunker(ntr_full_set, parallel)){
    output<- stats_g_prob_d_expected_improvement(delta, gamma, mu_c, sigma_c, s_c, w_c, ntr_group, B)
    quant <- output$quant
    cond_check <- quant >= delta
    if (sum(1*cond_check, na.rm = TRUE)>0){
      ntr_opt<- ntr_group[cond_check][1]
      break
    }
  }
  return(ntr_opt)
}

get_prob_ab_cert_sim<- function(n_tr,n_ce,alpha,mu_c,sigma_c,s_c,w_c,B=10000){
  C <- length(mu_c)
  mu_c_post_var<- sigma_c**4/(sigma_c**2+4*s_c**2/w_c/(n_tr+1e-30))
  mu_c_post_mean<- mu_c
  mu_c_post_sample<-array(rnorm(mean=mu_c_post_mean, sd=sqrt(mu_c_post_var),n=C*B), dim=c(C,B))
  v_post <-apply(w_c*mu_c_post_sample*(1*(mu_c_post_sample>=0)),2, sum)
  var_c_post <- 1/(1/sigma_c**2+w_c*n_tr/4/s_c**2)
  w_U <- w_c**2*(var_c_post+4*s_c**2/w_c/n_ce)
  U <- apply(w_U*(1*(mu_c_post_sample>0)), 2, sum)
  Y <- v_post/sqrt(U+1e-10)
  prob <- mean(pnorm(Y-qnorm(1-alpha)))
  return(prob)
}

get_prob_ab_cert_aprx <- function(n_tr, n_ce, alpha, mu_c, sigma_c, s_c, w_c, taylor=0){
  # Estimate LHS Equation 41
  k_c <- sqrt(1+4*s_c**2/w_c/(n_tr+1e-30)/(sigma_c**2))
  # dnorm is pdf
  pdf_c <- dnorm(mu_c*k_c/sigma_c)
  # pnorm is cdf
  cdf_c <-pnorm(mu_c*k_c/sigma_c)  
  var_c_post <- 1/(1/sigma_c**2+w_c*n_tr/4/s_c**2)
  mu_D <- sum(w_c*(sigma_c/k_c*pdf_c+mu_c*cdf_c))
  
  
  p_c <- cdf_c
  mu_U <- sum(w_c**2*(var_c_post+4*s_c**2/w_c/n_ce)*p_c)
  
  prob<- pnorm((mu_D/sqrt(mu_U) - qnorm(1-alpha)))
  if (taylor){
    # Higher-order Taylor series (Appendix D2) 
    var_D <- sum(w_c**2*(mu_c**2*cdf_c*(1-cdf_c)+sigma_c**2/k_c**2*(cdf_c-pdf_c**2)+mu_c*sigma_c/k_c*pdf_c*(1-2*cdf_c)))
    var_U <- sum(p_c*(1-p_c)*w_c**4*(var_c_post+4*s_c**2/w_c/n_ce)**2)
    var_DU <- var_D/mu_D**2+var_U/(4*mu_U**2)-1/mu_D/mu_U*sum(w_c**3*(1-p_c)*(sigma_c/k_c*pdf_c+mu_c*cdf_c)*(var_c_post+4*s_c**2/w_c/n_ce))
    prob<- pnorm((mu_D-qnorm(1-alpha)*sqrt(mu_U))/sqrt(mu_U+mu_D**2*var_DU))
  }
  return(prob)}

get_prob_ab_cert_aprx_array<- function(ntr, nce_set, alpha, mu_c, sigma_c, s_c, w_c, taylor=0){
  ntr <- ntr + 1e-30
  k_c <- sqrt(1+4*s_c**2/w_c/ntr/(sigma_c**2))
  # dnorm is pdf
  pdf_c <- dnorm(mu_c*k_c/sigma_c)
  # pnorm is cdf
  cdf_c <-pnorm(mu_c*k_c/sigma_c) 
  var_c_post <- 1/(1/sigma_c**2+w_c*ntr/4/s_c**2)
  
  mu_D <- sum(w_c*(sigma_c/k_c*pdf_c+mu_c*cdf_c))
  p_c <- cdf_c
  mu_U <-apply(w_c**2*(var_c_post+4*s_c**2/outer(w_c,nce_set))*p_c, 2, sum)
  prob <- pnorm((mu_D/sqrt(mu_U)-qnorm(1-alpha)))
  
  if(taylor){
    var_D <- sum(w_c**2*(mu_c**2*cdf_c*(1-cdf_c)+sigma_c**2/k_c**2*(cdf_c-pdf_c**2)+mu_c*sigma_c/k_c*pdf_c*(1-2*cdf_c)))
    var_U <-apply(p_c*(1-p_c)*w_c**4*(var_c_post+4*s_c**2/outer(w_c,nce_set))**2,2,sum)
    var_DU <-var_D/mu_D**2+var_U/(4*mu_U**2)-1/mu_D/mu_U*apply(w_c**3*(1-p_c)*(sigma_c/k_c*pdf_c+mu_c*cdf_c)*(var_c_post+4*s_c**2/outer(w_c,nce_set)),2,sum)
    prob <- pnorm((mu_D-qnorm(1-alpha)*sqrt(mu_U))/sqrt(mu_U+mu_D**2*var_DU))
  }
  return(prob)
}


solve_nce_aprx <- function(ntr_set, alpha, beta,mu_c, sigma_c, s_c, w_c){
  Z_ab <- qnorm(1-alpha)+qnorm(1-beta)
  sigma_c_p <- sqrt(1/(1/sigma_c**2+outer(w_c, ntr_set)/4/s_c**2))
  k_c <- sqrt(1+4*s_c**2/outer(w_c, (ntr_set+1e-30))/(sigma_c**2))
  # dnorm is pdf
  pdf_c <- dnorm(mu_c*k_c/sigma_c)
  # pnorm is cdf
  cdf_c <-pnorm(mu_c*k_c/sigma_c) 
  mu_delta <- apply(w_c*(sigma_c/k_c*pdf_c+mu_c*cdf_c), 2, sum)
  p_c <- cdf_c
  nce <- 4*Z_ab**2*apply(p_c*w_c*s_c**2, 2, sum)/(mu_delta**2-Z_ab**2*apply(w_c**2*sigma_c_p**2*p_c,2, sum))    
  return (ceiling(nce))
}

solve_nce_sim_inner<- function(ntr,nce_min,nce_max,alpha,beta,mu_c,sigma_c,s_c,w_c,B){
  
  if ((nce_max%/%1) == (nce_min%/%1)){
    return (nce_min)
  }
  
  nce_mid <- (nce_min+nce_max)/2
  mid_val <- get_prob_ab_cert_sim(ntr,nce_mid,alpha,mu_c,sigma_c,s_c,w_c,B)-(1-beta)
  if (mid_val<=0){
    return(solve_nce_sim_inner(ntr,nce_mid,nce_max,alpha,beta,mu_c,sigma_c,s_c,w_c,B))
  }
  else{return(solve_nce_sim_inner(ntr,nce_min,nce_mid,alpha,beta,mu_c,sigma_c,s_c,w_c,B))
  }
}

solve_nce_sim<- function(ntr,nce_max,alpha,beta,mu_c,sigma_c,s_c,w_c,B){
  nce_opt <- -1
  max_prob <- get_prob_ab_cert_sim(ntr,nce_max,alpha,mu_c,sigma_c,s_c,w_c,B)
  if (max_prob>1-beta){
    nce_opt <- solve_nce_sim_inner(ntr,0,nce_max,alpha,beta,mu_c,sigma_c,s_c,w_c,B)
  }
  return (ceiling(nce_opt))
}


solve_nce_taylor<- function(ntr,nce_max,alpha,beta,mu_c,sigma_c,s_c,w_c){
  nce_full_set <- seq(0, nce_max, 1)+1e-30
  nce_opt <- -1
  for (nce_group in chunker(nce_full_set, 150000)){
    prob_group <- get_prob_ab_cert_aprx_array(ntr, nce_group, alpha, mu_c, sigma_c, s_c, w_c, taylor=1)
    cond_check <- prob_group>= (1-beta)
    if (sum(1*cond_check, na.rm = TRUE)>0){
      nce_opt<- nce_group[cond_check][1]
      break
    }
  }
  return(nce_opt)
}


# Key solver for (a,b)-certification
solve_ab_certification<- function(alpha, beta, mu_c, sigma_c, s_c, w_c, n_max, B=0, taylor=0){
  # Key solver for (a,b)-certification
  # B > 0: simulated solution
  # B==0: analytical approximation; taylor = 1 uses a more-precise but slower approximation (Appendix D2)
  
  ntr_full_set <- seq(0, n_max, 1)
  ntr_opt <- -1
  nce_opt <- -1
  ntotal_opt <- n_max+1
  
  if (B>0){
    #Algorithm 2 
    for (ntr in ntr_full_set){
      if (ntr>=ntotal_opt){next}
      nce <- solve_nce_sim(ntr, ntotal_opt-ntr, alpha, beta, mu_c, sigma_c, s_c, w_c, B)
      if (nce>0){
        ntr_opt <- ntr
        nce_opt <- nce
        ntotal_opt <- ntr+nce
      }
    }
  }
  
  if (B==0 & taylor ==1){
    # Using higher-order Taylor series (Appendix D2)
    for (ntr in ntr_full_set){
      if (ntr>=ntotal_opt){next}
      nce <- solve_nce_taylor(ntr, ntotal_opt-ntr, alpha, beta, mu_c, sigma_c, s_c, w_c)
      if (nce>0){
        ntr_opt <- ntr
        nce_opt <- nce
        ntotal_opt <- ntr+nce
      }
    }
  }
  
  if (B==0 & taylor ==0) {
    # Algorithm 3
    for (ntr_group in chunker(ntr_full_set, 100000)){
      nce_group <- solve_nce_aprx(
        ntr_group, alpha, beta, mu_c, sigma_c, s_c, w_c)
      nce_group[nce_group<0] <- n_max+1
      ntotal_group <- ntr_group+nce_group
      ntotal_group_min <- min(ntotal_group)
      if (ntotal_group_min<ntotal_opt){
        idx <- which(ntotal_group==ntotal_group_min)[1]
        ntr_opt <- ntr_group[idx]
        nce_opt <- nce_group[idx]
        ntotal_opt <- ntotal_group_min
      }
      if (ntr_group[length(ntr_group)]>=ntotal_opt){
        break
      }
    }
  }
  
  
  output<-list('ntr_opt'=ntr_opt, 'nce_opt'= round(nce_opt))
  return (output)
}

