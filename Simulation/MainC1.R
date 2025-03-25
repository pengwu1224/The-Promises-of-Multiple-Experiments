# Encoding: UTF-8

rm(list = ls( )) 
library(dplyr)

# ------------------------ Setup   -----------------------------
# function used to generate data
source('GenData.R')

true.value <- exp(c(0,1)-0.5) / (1+ exp(c(0,1)-0.5)) # true value of \theta


LOOP <- 1000                      # number of simulations
SampleSize <- c(100, 200, 500)   # sample size for each trial 
NumTrial <- 10                    # number of trials 

# used to save results 
Result <- data.frame(matrix(nrow = 2, ncol = 14))
Result[, 1] <- rep('case2', each = 2)  # rep(c(200, 500, 1000), each = 2)
Result[, 2] <- c('pi_10','pi_11')      # rep(c('pi_10','pi_11'), 3)
names(Result) <- c('Case','pi_ab',
                   c( paste0(c('Bias','SSE','ESE', 'CP95'), '.100'), 
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.200'),
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.500')))

# -------------------- Simulation   ---------------------------------


set.seed(3174)


for(s in 1:3){
  n <- SampleSize[s]
  
  PARA <- matrix(nrow = LOOP, ncol = 2)
  ESE <- matrix(nrow = LOOP, ncol = 2)
  COUNT95 <- matrix(nrow = LOOP, ncol = 2)
  
  cat('========  n =', n,'  ========= \n')
  for(loop in 1:LOOP){
    
    # step 1. generate the data 
    dat <- gen_dat1(n, num_trial = NumTrial)
    
    # step 2. estimate the parameter \theta 
    y1.prob.1.est <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
    y0.prob.1.est <- rep(NA, 10)   # P(Y0 = 1 | G=g) = P(Y = 1 | G=g, A=0)
    for(g in 1:10){
      dat.tmp <- filter(dat, G == g)
      A.g <- dat.tmp$A
      Y.g <- dat.tmp$Y
      y1.prob.1.est[g] <- sum( Y.g * A.g) / sum(A.g)
      y0.prob.1.est[g] <- sum( Y.g * (1-A.g)) / sum(1-A.g)
    }
    y0.prob.0.est <- 1- y0.prob.1.est    # P(Y0 = 0 | G=g) 
    
    mod <- lm(y1.prob.1.est ~ y0.prob.0.est + y0.prob.1.est - 1)
    para <- mod$coefficients
    
    # step 3. asymptotic variance estimation (bootstrap)
    B <- 100
    est.tmp <- matrix(nrow = B, ncol  =2)
    for(b in 1:B){
      ind <- sapply(1:NumTrial, function(j){sample(1:n, replace = T)})
      ind <- c(ind)
      ind <- ind + rep( (1:NumTrial)*n - n, each = n) 
      dat.b <- dat[ind, ]
      
      y1.prob.1.b <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
      y0.prob.1.b <- rep(NA, 10)   # P(Y0 = 1 | G=g) = P(Y = 1 | G=g, A=0)
      for(g in 1:10){
        dat.tmp <- filter(dat.b, G == g)
        A.g <- dat.tmp$A
        Y.g <- dat.tmp$Y
        y1.prob.1.b[g] <- sum( Y.g * A.g) / sum(A.g)
        y0.prob.1.b[g] <- sum( Y.g * (1-A.g)) / sum(1-A.g)
      }
      y0.prob.0.b <- 1- y0.prob.1.b    # P(Y0 = 0 | G=g) 
      
      mod <- lm(y1.prob.1.b ~ y0.prob.0.b + y0.prob.1.b - 1)
      est.tmp[b, ] <- mod$coefficients
    }
    ese.boot <- apply(est.tmp, 2, sd) 
    
    # step 5. save the results
    PARA[loop,] <- para 
    ESE[loop,] <-  ese.boot
    COUNT95[loop,] <- abs(para - true.value) / ese.boot  <= 1.96
    
    cat(loop,'\r')
  }
  
  bias <- apply(PARA, 2, mean, na.rm = T) - true.value
  sse <- apply(PARA, 2, sd, na.rm = T)
  ese <- sqrt( apply(ESE*ESE, 2, mean, na.rm = T) )
  cp95 <- apply(COUNT95, 2, mean, na.rm = T)
  
  Result[1:2, 3 + (s-1)*4] <- bias  
  Result[1:2, 4 + (s-1)*4] <- sse
  Result[1:2, 5 + (s-1)*4] <- ese
  Result[1:2, 6 + (s-1)*4] <- cp95
  
  print(Result)
}
write.csv(Result, file = 'result/CaseC1_result.csv')

