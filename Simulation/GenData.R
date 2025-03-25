
# generate data 

# ------------- Assumptions 1–2 and Condition 1 hold  -----------------

gen_dat1 <- function(n, num_trial){
  
  G <- rep(1:num_trial, each = n)     # indicator of trials
  A <- rep(NA, n * num_trial)
  Y0 <- rep(NA, n * num_trial)
  Y1 <- rep(NA, n * num_trial)
  
  prob.y0 <- seq(0.5, 0.8, len = num_trial)
  
  for(g in 1:10){
    A[G == g] <- rbinom(n, size = 1, prob = 0.5)
    Y0[G == g] <- rbinom(n, size = 1, prob = prob.y0[g])
    
    prob.y1 <- exp(Y0[G == g] - 0.5) / (1+ exp(Y0[G == g] - 0.5)) 
    Y1[G == g] <- rbinom(n, size = 1, prob = prob.y1)
  }
  
  Y <- A*Y1 + (1-A)*Y0
  dat <- data.frame(G=G, A=A, Y=Y)
}


gen_dat2 <- function(n, num_trial){
  
  G <- rep(1:num_trial, each = n)     # indicator of trials
  A <- rep(NA, n * num_trial)
  Y0 <- rep(NA, n * num_trial)
  Y1 <- rep(NA, n * num_trial)
  
  prob.y0 <- seq(0.5, 0.8, len = num_trial)
  
  for(g in 1:10){
    A[G == g] <- rbinom(n, size = 1, prob = 0.5)
    Y0[G == g] <- rbinom(n, size = 1, prob = prob.y0[g])
    
    prob.y1 <- exp(Y0[G == g] + 0.5) / (1+ exp(Y0[G == g] + 0.5)) 
    Y1[G == g] <- rbinom(n, size = 1, prob = prob.y1)
  }
  
  Y <- A*Y1 + (1-A)*Y0
  dat <- data.frame(G=G, A=A, Y=Y)
}



# --------- Assumptions 4–5 and Condition 3 hold   --------------------
gen_dat3 <- function(n, num_trial){
  
  G <- rep(1:num_trial, each = n)     # indicator of trials
  A <- rep(NA, n * num_trial)
  S0 <- rep(NA, n * num_trial)
  S1 <- rep(NA, n * num_trial)
  Y0 <- rep(NA, n * num_trial)
  Y1 <- rep(NA, n * num_trial)
  
  prob.y0 <- seq(0.5, 0.8, len = num_trial)
  prob.s0 <- seq(0.5, 0.8, len = num_trial)
  
  for(g in 1:10){
    A[G == g] <- rbinom(n, size = 1, prob = 0.5)
    S0[G == g] <- rbinom(n, size = 1, prob = prob.s0[g])
    Y0[G == g] <- rbinom(n, size = 1, prob = prob.y0[g])
    
    prob.s1 <- exp( 0.5*(Y0[G == g]+S0[G == g]) - 0.5) / (1+ exp( 0.5*(Y0[G == g]+S0[G == g]) - 0.5))
    prob.y1 <- exp( 0.5*(Y0[G == g]+S0[G == g]) + 0.5) / (1+ exp( 0.5*(Y0[G == g]+S0[G == g]) + 0.5))
    S1[G == g] <- rbinom(n, size = 1, prob = prob.s1)
    Y1[G == g] <- rbinom(n, size = 1, prob = prob.y1)
  }
  
  S <- A*S1 + (1-A)*S0
  Y <- A*Y1 + (1-A)*Y0
  dat <- data.frame(G=G, A=A, S=S, Y=Y)
}


# ----------- Assumptions 4, 6–7 and Condition 5 hold   ----------------
gen_dat4 <- function(n, num_trial){

  G <- rep(1:num_trial, each = n)     # indicator of trials
  A <- rep(NA, n * num_trial)
  S0 <- rep(NA, n * num_trial)
  S1 <- rep(NA, n * num_trial)
  Y0 <- rep(NA, n * num_trial)
  Y1 <- rep(NA, n * num_trial)
  
  prob.s0 <- seq(0.3, 0.6, len = num_trial)
  prob.s1 <- seq(0.5, 0.8, len = num_trial)
  
  for(g in 1:10){
    A[G == g] <- rbinom(n, size = 1, prob = 0.5)
    S0[G == g] <- rbinom(n, size = 1, prob = prob.s0[g])
    S1[G == g] <- rbinom(n, size = 1, prob = prob.s1[g])
    S0 <- ifelse(S1 == 1, S0, 0)       # ensure monotonicity
    
    prob.y0 <- exp( 0.5*(S0[G == g]+S1[G == g]) - 0.5) / (1+ exp( 0.5*(S0[G == g]+S1[G == g]) - 0.5))
    prob.y1 <- exp( 0.5*(S0[G == g]+S1[G == g]) + 0.5) / (1+ exp( 0.5*(S0[G == g]+S1[G == g]) + 0.5))
    Y0[G == g] <- rbinom(n, size = 1, prob = prob.y0)
    Y1[G == g] <- rbinom(n, size = 1, prob = prob.y1)
  }
  
  S <- A*S1 + (1-A)*S0
  Y <- A*Y1 + (1-A)*Y0
  dat <- data.frame(G=G, A=A, S=S, Y=Y)
}

