# Encoding: UTF-8

rm(list = ls( )) 
library(dplyr)

# ------------------------ Setup   -----------------------------
# function used to generate data
source('GenData.R')

true.value.y0 <- exp(c(0,0.5,1)-0.5) / (1+ exp(c(0,0.5,1)-0.5)) # true value of \theta
true.value.y1 <- exp(c(0,0.5,1)+0.5) / (1+ exp(c(0,0.5,1)+0.5)) # true value of \theta
true.value <- c(true.value.y0, true.value.y1)

LOOP <- 1000                      # number of simulations
SampleSize <- c(100, 200, 500)   # sample size for each trial 
NumTrial <- 10                    # number of trials 

# used to save results 
Result <- data.frame(matrix(nrow = 6, ncol = 14))
Result[, 1] <- rep('case5', each = 6)  
Result[, 2] <- c('pi_y0_00','pi_y0_01', 'pi_y0_11',
                 'pi_y1_00','pi_y1_01', 'pi_y1_11')     
names(Result) <- c('Case','pi_ab',
                   c( paste0(c('Bias','SSE','ESE', 'CP95'), '.100'), 
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.200'),
                      paste0(c('Bias','SSE','ESE',  'CP95'), '.500')))

# -------------------- Simulation   ---------------------------------

set.seed(7724)

for(s in 1:3){
  n <- SampleSize[s]
  
  PARA <- matrix(nrow = LOOP, ncol = 6)
  ESE <- matrix(nrow = LOOP, ncol = 6)
  COUNT95 <- matrix(nrow = LOOP, ncol = 6)
  
  cat('========  n =', n,'  ========= \n')
  for(loop in 1:LOOP){
    
    # step 1. generate the data 
    dat <- gen_dat4(n, num_trial = NumTrial)
    
    # step 2. estimate the parameter \theta 
    # (a) estimate the principal scores P(S0=a, S1=b|G=g)
    y1.prob.1 <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
    y0.prob.1 <- rep(NA, 10)   # P(S0 = 1 | G=g) = P(Y = 1 | G=g, A=0)
    
    s0.s1.prob.11 <- rep(NA, 10)   # P(S0 = 1, S1 = 1 | G=g)
    s0.s1.prob.01 <- rep(NA, 10)   # P(S0 = 0, S1 = 1 | G=g)
    s0.s1.prob.00 <- rep(NA, 10)   # P(S0 = 0, S1 = 0 | G=g)
    
    for(g in 1:10){
      dat.tmp <- filter(dat, G == g)
      A.g <- dat.tmp$A
      S.g <- dat.tmp$S         # surrogate
      Y.g <- dat.tmp$Y         # outcome 
      
      y1.prob.1[g] <- sum( Y.g * A.g) / sum(A.g)
      y0.prob.1[g] <- sum( Y.g * (1-A.g)) / sum(1-A.g)
      
      s0.s1.prob.11[g] <- sum(S.g * (1-A.g)) / sum(1-A.g)
      s0.s1.prob.01[g] <- sum(S.g * A.g) / sum(A.g) - s0.s1.prob.11[g]
      s0.s1.prob.00[g] <- 1- s0.s1.prob.11[g]  - s0.s1.prob.01[g] 
      # sum( (1-S.g) * A.g) / sum(A.g)
    }
    
    # (b) estimate pi_y0_11 = P(Y0=1|S0=1,S1=1) and pi_y1_00=P(Y1=1|S0=0,S1=0)
    pi_y0_11 <- sum( dat$Y * dat$S * (1-dat$A) ) / sum(dat$S * (1-dat$A)) 
    pi_y1_00 <- sum( dat$Y * (1-dat$S) * dat$A ) / sum((1-dat$S) * dat$A) 
    
    # (c) estimate 'pi_y0_00','pi_y0_01' and 'pi_y1_01', 'pi_y1_11'
    # estimate 'pi_y0_00','pi_y0_01'
    mod1 <- lm(y0.prob.1 - pi_y0_11*s0.s1.prob.11 ~ s0.s1.prob.00 + s0.s1.prob.01 -1)
    para1 <- c(mod1$coefficients, pi_y0_11)
    
    # estimate 'pi_y1_01', 'pi_y1_11'
    mod2 <- lm(y1.prob.1 - pi_y1_00*s0.s1.prob.00 ~ s0.s1.prob.01 + s0.s1.prob.11 -1)
    para2 <- c(pi_y1_00, mod2$coefficients)
    
    para <- c(para1, para2)
    
    
    # step 3. asymptotic variance estimation (bootstrap)
    B <- 100
    est.tmp <- matrix(nrow = B, ncol  = 6)
    for(b in 1:B){
      ind <- sapply(1:NumTrial, function(j){sample(1:n, replace = T)})
      ind <- c(ind)
      ind <- ind + rep( (1:NumTrial)*n - n, each = n) 
      dat.b <- dat[ind, ]
      
      # (a) estimate the principal scores P(S0=a, S1=b|G=g)
      y1.prob.1.b <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
      y0.prob.1.b <- rep(NA, 10)   # P(S0 = 1 | G=g) = P(Y = 1 | G=g, A=0)
      
      s0.s1.prob.11.b <- rep(NA, 10)   # P(S0 = 1, S1 = 1 | G=g)
      s0.s1.prob.01.b <- rep(NA, 10)   # P(S0 = 0, S1 = 1 | G=g)
      s0.s1.prob.00.b <- rep(NA, 10)   # P(S0 = 0, S1 = 0 | G=g)
      
      for(g in 1:10){
        dat.tmp <- filter(dat.b, G == g)
        A.g <- dat.tmp$A
        S.g <- dat.tmp$S         # surrogate
        Y.g <- dat.tmp$Y         # outcome 
        
        y1.prob.1.b[g] <- sum( Y.g * A.g) / sum(A.g)
        y0.prob.1.b[g] <- sum( Y.g * (1-A.g)) / sum(1-A.g)
        
        s0.s1.prob.11.b[g] <- sum(S.g * (1-A.g)) / sum(1-A.g)
        s0.s1.prob.01.b[g] <- sum(S.g * A.g) / sum(A.g) - s0.s1.prob.11.b[g]
        s0.s1.prob.00.b[g] <- 1- s0.s1.prob.11.b[g]  - s0.s1.prob.01.b[g] 
        # sum( (1-S.g) * A.g) / sum(A.g)
      }
      
      # (b) estimate pi_y0_11 = P(Y0=1|S0=1,S1=1) and pi_y1_00=P(Y1=1|S0=0,S1=0)
      pi_y0_11.b <- sum( dat.b$Y * dat.b$S * (1-dat.b$A) ) / sum(dat.b$S * (1-dat.b$A)) 
      pi_y1_00.b <- sum( dat.b$Y * (1-dat.b$S) * dat.b$A ) / sum((1-dat.b$S) * dat.b$A) 
      
      # (c) estimate 'pi_y0_00','pi_y0_01' and 'pi_y1_01', 'pi_y1_11'
      # estimate 'pi_y0_00','pi_y0_01'
      mod1.b <- lm(y0.prob.1.b - pi_y0_11.b*s0.s1.prob.11.b ~ s0.s1.prob.00.b 
                   + s0.s1.prob.01.b -1)
      para1.b <- c(mod1.b$coefficients, pi_y0_11.b)
      
      # estimate 'pi_y1_01', 'pi_y1_11'
      mod2.b <- lm(y1.prob.1.b - pi_y1_00.b*s0.s1.prob.00.b ~ s0.s1.prob.01.b 
                   + s0.s1.prob.11.b -1)
      para2.b <- c(pi_y1_00.b, mod2.b$coefficients)
      
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

write.csv(Result, file = 'result/CaseC4_result.csv')
