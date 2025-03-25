# Encoding: UTF-8

rm(list = ls( )) 
library(dplyr)

# ------------------------ Setup   -----------------------------
# function used to generate data
source('GenData.R')

true.value.s1 <- exp(c(0,0.5,0.5,1)-0.5) / (1+ exp(c(0,0.5,0.5,1)-0.5)) # true value of \theta
true.value.y1 <- exp(c(0,0.5,0.5,1)+0.5) / (1+ exp(c(0,0.5,0.5,1)+0.5)) # true value of \theta
true.value <- c(true.value.s1, true.value.y1)

LOOP <- 1000                      # number of simulations
SampleSize <- c(100, 200, 500)   # sample size for each trial 
NumTrial <- 10                    # number of trials 

# used to save results 
Result <- data.frame(matrix(nrow = 8, ncol = 14))
Result[, 1] <- rep('case4', each = 8)  # rep(c(200, 500, 1000), each = 2)
Result[, 2] <- c('pi_s_00','pi_s_01', 'pi_s_10', 'pi_s_11',
                 'pi_y_00','pi_y_01', 'pi_y_10', 'pi_s_11')     
names(Result) <- c('Case','pi_ab',
                   c( paste0(c('Bias','SSE','ESE', 'CP95'), '.100'), 
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.200'),
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.500')))

# -------------------- Simulation   ---------------------------------


set.seed(1434)

for(s in 1:3){
  n <- SampleSize[s]
  
  PARA <- matrix(nrow = LOOP, ncol = 8)
  ESE <- matrix(nrow = LOOP, ncol = 8)
  COUNT95 <- matrix(nrow = LOOP, ncol = 8)
  
  cat('========  n =', n,'  ========= \n')
  for(loop in 1:LOOP){
    
    # step 1. generate the data 
    dat <- gen_dat3(n, num_trial = NumTrial)
    
    # step 2. estimate the parameter \theta 
    y1.prob.1 <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
    s1.prob.1 <- rep(NA, 10)   # P(S1 = 1 | G=g) = P(S = 1 | G=g, A=1)
    
    s0.y0.prob.00 <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
    s0.y0.prob.01 <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
    s0.y0.prob.10 <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
    s0.y0.prob.11 <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)
    
    for(g in 1:10){
      dat.tmp <- filter(dat, G == g)
      A.g <- dat.tmp$A
      S.g <- dat.tmp$S         # surrogate
      Y.g <- dat.tmp$Y         # outcome 
      
      y1.prob.1[g] <- sum( Y.g * A.g) / sum(A.g)
      s1.prob.1[g] <- sum( S.g * A.g) / sum(A.g)
      
      s0.y0.prob.00[g] <- sum( (1-S.g)*(1-Y.g)*(1-A.g)) / sum((1-A.g))
      s0.y0.prob.01[g] <- sum( (1-S.g)*   Y.g *(1-A.g)) / sum((1-A.g))
      s0.y0.prob.10[g] <- sum(    S.g *(1-Y.g)*(1-A.g)) / sum((1-A.g))
      s0.y0.prob.11[g] <- sum(    S.g*    Y.g *(1-A.g)) / sum((1-A.g))
    }
    
    mod1 <- lm(s1.prob.1 ~ s0.y0.prob.00 + s0.y0.prob.01 + 
                 s0.y0.prob.10+ s0.y0.prob.11 -1)
    para1 <- mod1$coefficients
    mod2 <- lm(y1.prob.1 ~ s0.y0.prob.00 + s0.y0.prob.01 + 
                 s0.y0.prob.10+ s0.y0.prob.11 -1)
    para2 <- mod2$coefficients
    
    para <- c(para1, para2)
    
    
    # step 3. asymptotic variance estimation (bootstrap)
    B <- 100
    est.tmp <- matrix(nrow = B, ncol  = 8)
    for(b in 1:B){
      ind <- sapply(1:NumTrial, function(j){sample(1:n, replace = T)})
      ind <- c(ind)
      ind <- ind + rep( (1:NumTrial)*n - n, each = n) 
      dat.b <- dat[ind, ]
      
      y1.prob.1.b <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
      s1.prob.1.b <- rep(NA, 10)   # P(S1 = 1 | G=g) = P(S = 1 | G=g, A=1)
      
      s0.y0.prob.00.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
      s0.y0.prob.01.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
      s0.y0.prob.10.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
      s0.y0.prob.11.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)
      
      for(g in 1:10){
        dat.tmp <- filter(dat.b, G == g)
        A.g <- dat.tmp$A
        S.g <- dat.tmp$S         # surrogate
        Y.g <- dat.tmp$Y         # outcome 
        
        y1.prob.1.b[g] <- sum( Y.g * A.g) / sum(A.g)
        s1.prob.1.b[g] <- sum( S.g * A.g) / sum(A.g)
        
        s0.y0.prob.00.b[g] <- sum( (1-S.g)*(1-Y.g)*(1-A.g)) / sum((1-A.g))
        s0.y0.prob.01.b[g] <- sum( (1-S.g)*   Y.g *(1-A.g)) / sum((1-A.g))
        s0.y0.prob.10.b[g] <- sum(    S.g *(1-Y.g)*(1-A.g)) / sum((1-A.g))
        s0.y0.prob.11.b[g] <- sum(    S.g*    Y.g *(1-A.g)) / sum((1-A.g))
      }
      
      mod1.b <- lm(s1.prob.1.b ~ s0.y0.prob.00.b + s0.y0.prob.01.b + 
                     s0.y0.prob.10.b + s0.y0.prob.11.b -1)
      para1.b <- mod1.b$coefficients
      mod2.b <- lm(y1.prob.1.b ~ s0.y0.prob.00.b + s0.y0.prob.01.b + 
                     s0.y0.prob.10.b + s0.y0.prob.11.b -1)
      para2.b <- mod2.b$coefficients
      
      est.tmp[b, ] <- c(para1.b, para2.b)
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
  
  Result[, 3 + (s-1)*4] <- bias  
  Result[, 4 + (s-1)*4] <- sse
  Result[, 5 + (s-1)*4] <- ese
  Result[, 6 + (s-1)*4] <- cp95
  
  print(Result)
}
write.csv(Result, file = 'result/CaseC3_result.csv')

