

rm(list=ls(all=TRUE))
library(jointseg)
library(lsei)
library(dplyr)
setwd('~/Desktop/MultipleExp/Application')


# ===================== Real Data (ACCTs data)   ============================

dat.trial <- vector(mode = 'list', length = 10)  # 10 trials, saved as a list 

# trial 1
dat.trial[[1]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),           # treatment
  S = c(0, 0, 1, 1, 0, 0, 1, 1),           # surrogate
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),           # outcome
  N = c(122.3, 11.7, 22.7, 218.3, 101.2, 5.8, 21.2, 220.8)) # number of units

# trial 2
dat.trial[[2]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(52.0, 5.0, 4.0, 65.0, 36.0, 4.0, 7.0, 74.0))

# trial 3
dat.trial[[3]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(164.1, 28.0, 14.0, 262.9, 126.0, 11.0, 14.0, 306.0))

# trial 4
dat.trial[[4]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(104.5, 13.9, 26.1, 313.5, 88.8, 10.5, 8.7, 330.0))

# trial 5
dat.trial[[5]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(152.5, 16.5, 18.1, 335.9, 110.6, 18.6, 10.3, 379.5))

# trial 6
dat.trial[[6]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(84.3, 8.0, 9.1, 128.6, 58.2, 15.1, 10.2, 141.6))

# trial 7
dat.trial[[7]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(171.7, 23.8, 27.5, 470.0, 150.9, 19.9, 23.6, 502.5))

# trial 8
dat.trial[[8]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(46.3, 5.0, 9.2, 92.5, 63.0, 8.2, 15.2, 168.6))

# trial 9
dat.trial[[9]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(219.5, 35.6, 26.6, 788.4, 208.6, 31.3, 27.9, 798.1))

# trial 10
dat.trial[[10]] <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 1, 1),
  S = c(0, 0, 1, 1, 0, 0, 1, 1),
  Y = c(0, 1, 0, 1, 0, 1, 0, 1),
  N = c(115.2, 17.2, 18.5, 286.1, 111.7, 15.1, 13.2, 301.0))

dat <- data.frame(matrix(nrow = 80, ncol = 5))
names(dat) <- c('G', 'A', 'S', 'Y', 'N')
dat$G <- rep(1:10, each  = 8)

for(i in 1:10){
  dat[8*(i-1)+1:8, 2:5] = dat.trial[[i]]
}

rm(dat.trial, i)


# =====================  The effect of A on Y  =============================
#library(dplyr)

y1.prob.1 <- rep(NA, 10)   # P(Y1 = 1 | G=g)
y0.prob.1 <- rep(NA, 10)   # P(Y0 = 1 | G=g)
for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  Y.g <- dat.tmp$Y
  y1.prob.1[g] <- sum(N.g[(A.g==1) & (Y.g == 1)]) / sum(N.g[A.g==1])
  y0.prob.1[g] <- sum(N.g[(A.g==0) & (Y.g == 1)]) / sum(N.g[A.g==0])
}
y0.prob.0 <- 1- y0.prob.1    # P(Y0 = 0 | G=g) 


## ATE for each trial 
ATEs <- y1.prob.1 - y0.prob.1
par(mfrow = c(1, 3), mar = c(5, 4.5, 2, 2))
names(ATEs) <- 1:10
barplot(ATEs, xlab = 'Trial Number', 
        ylab =  'ATEs of A on Y', 
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)

plot(y0.prob.0, y1.prob.1,
     xlab = expression(hat(P)(Y^0 == 0 ~ "|" ~  G==g)),
     ylab = expression(hat(P)(Y^1 == 1 ~ "|" ~  G==g)),
     cex.lab = 1.2, cex.axis = 1.2, cex = 1.3)
plot(y0.prob.1, y1.prob.1,
     xlab = expression(hat(P)(Y^0 == 1 ~ "|" ~  G==g)),
     ylab = expression(hat(P)(Y^1 == 1 ~ "|" ~  G==g)),
     cex.lab = 1.2,  cex.axis = 1.2, cex = 1.3)


mod <- lm(y1.prob.1 ~ y0.prob.0 + y0.prob.1 - 1)
summary(mod)              # The coefficients are significant 
est <- mod$coefficients  # estimators of $\pi_{1|0}$ and $\pi_{1|1}$ 

est.all <- c(est, 1 - est)  
# the last two are estimators of $\pi_{0|0}$ and $\pi_{0|1}$ 

# estimate the joint distributions P(Y0=0, Y1=1), P(Y0=1, Y1=1)
#                                  P(Y0=0, Y1=0), P(Y0=1, Y1=0)
prob.y0.y1 <- matrix(nrow = 10, ncol = 4)
for(g in 1:10){
  prob.y0.y1[g, ] <- est.all * c(y0.prob.0[g], y0.prob.1[g], y0.prob.0[g], y0.prob.1[g])
}
prob.y0.y1.vec <- c(prob.y0.y1)

## calculate confidence interval with 500 bootstraps
set.seed(123)
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # total sample size for each trial 

B <- 500
EstAll <- matrix(nrow = B, ncol = 4)
Res <- matrix(nrow = B, ncol = 10)
Est.prob.y0.y1.vec <- matrix(nrow = B, ncol = 40)

for(b in 1:B){
  cat('------- bootstrap ', b, ' ---------- \r')
  
  # step 1: generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # step 2: point estimates based on the boostraped sample 
  y1.prob.1.b <- rep(NA, 10)   # P(Y1 = 1 | G=g)
  y0.prob.1.b <- rep(NA, 10)   # P(Y0 = 1 | G=g)
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    ng <- dat.tmp$N
    Y.g <- dat.tmp$Y
    y1.prob.1.b[g] <- sum(ng[(A.g==1) & (Y.g == 1)]) / sum(ng[A.g==1])
    y0.prob.1.b[g] <- sum(ng[(A.g==0) & (Y.g == 1)]) / sum(ng[A.g==0])
    
  }
  y0.prob.0.b <- 1- y0.prob.1.b    # P(Y0 = 0 | G=g) 
  
  mod.b <- lm(y1.prob.1.b ~ y0.prob.0.b + y0.prob.1.b - 1)
  est.b <- mod.b$coefficients
  est.all.b <- c(est.b, 1- est.b) # state transition probabilities
  EstAll[b,] <- est.all.b     
  Res[b, ] <- y1.prob.1.b - est.b[1]*y0.prob.0.b -  est.b[2]*y0.prob.1.b # residual
  
  # joint distributions 
  prob.y0.y1.b <-  matrix(nrow = 10, ncol = 4)
  for(g in 1:10){
    prob.y0.y1.b[g, ] <- est.all.b * c(y0.prob.0.b[g], y0.prob.1.b[g], 
                                       y0.prob.0.b[g], y0.prob.1.b[g])
  }
  prob.y0.y1.vec.b <- c(prob.y0.y1.b)
  Est.prob.y0.y1.vec[b, ] <- prob.y0.y1.vec.b
}
ese.boot <- apply(EstAll, 2, sd) 
res.ese.boot <- apply(Res, 2, sd) 
prob.y0.y1.vec.est.boot <- apply(Est.prob.y0.y1.vec, 2, sd)
prob.y0.y1.est.boot <- matrix(prob.y0.y1.vec.est.boot, nrow = 10, ncol = 4)

# 95% confidence interval
est.all - 1.96*ese.boot
est.all + 1.96*ese.boot

## Test Assumption 2
stat <- sum( (y1.prob.1 - est[1]*y0.prob.0 -  est[2]*y0.prob.1)^2/res.ese.boot^2 )
1 - pchisq(stat, df = 8)  # p-value     

## Plot the joint distributions
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 1, 2))

plot(1:10, prob.y0.y1[, 3], type = 'b', lty = 5,    # P(Y0=0, Y1=0|G=g)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.45), 
     cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(hat(P)(Y^0 == 0, Y^1 ==  0~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.y0.y1[, 3] - 1.96*prob.y0.y1.est.boot[,3],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.y0.y1[, 3] + 1.96*prob.y0.y1.est.boot[,3],  
      lwd = 1.5, col = 'black', lty = 2)


plot(1:10, prob.y0.y1[, 1], type = 'b', lty = 5,    # P(Y0=0, Y1=1|G=g)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.35), 
     cex.axis = 0.9,  lwd = 1.5)    
title(ylab = expression(hat(P)(Y^0 == 0, Y^1 ==  1~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.y0.y1[, 1] - 1.96*prob.y0.y1.est.boot[,1],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.y0.y1[, 1] + 1.96*prob.y0.y1.est.boot[,1],  
      lwd = 1.5, col = 'black', lty = 2)


plot(1:10, prob.y0.y1[, 4], type = 'b', lty = 5,    # P(Y0=1, Y1=0|G=g)
     xlab = '',
     ylab = '', 
     ylim = c(-0.03, 0.30), 
     cex.axis = 0.9,  lwd = 1.5) 
title(xlab = 'Trial Number',
      ylab = expression(hat(P)(Y^0 == 1, Y^1 ==  0~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.y0.y1[, 4] - 1.96*prob.y0.y1.est.boot[,4],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.y0.y1[, 4] + 1.96*prob.y0.y1.est.boot[,4],  
      lwd = 1.5, col = 'black', lty = 2)

plot(1:10, prob.y0.y1[, 2], type = 'b', lty = 5,    # P(Y0=1, Y1=1)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.80), 
     cex.axis = 0.9,  lwd = 1.5)    
title(xlab = 'Trial Number',
      ylab = expression(hat(P)(Y^0 == 1, Y^1 ==  1~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.y0.y1[, 2] - 1.96*prob.y0.y1.est.boot[,2],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.y0.y1[, 2] + 1.96*prob.y0.y1.est.boot[,2],  
      lwd = 1.5, col = 'black', lty = 2)



# ======================  The effect of A on S  =============================


s1.prob.1 <- rep(NA, 10)   # P(S1 = 1 | G=g)
s0.prob.1 <- rep(NA, 10)   # P(S0 = 1 | G=g)
for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  S.g <- dat.tmp$S         # surrogate
  s1.prob.1[g] <- sum(N.g[(A.g==1) & (S.g == 1)]) / sum(N.g[A.g==1])
  s0.prob.1[g] <- sum(N.g[(A.g==0) & (S.g == 1)]) / sum(N.g[A.g==0])
}
s0.prob.0 <- 1- s0.prob.1    # P(Y0 = 0 | G=g) 

## ATE for each trial 
ATEs <- s1.prob.1 - s0.prob.1

par(mfrow = c(1, 3), mar = c(5, 4.5, 2, 2))
names(ATEs) <- 1:10
barplot(ATEs, xlab = 'Trial Number', 
        ylab =  'ATEs of A on S', 
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)

plot(s0.prob.0, s1.prob.1,
     xlab = expression(hat(P)(S^0 == 0 ~ "|" ~  G==g)),
     ylab = expression(hat(P)(S^1 == 1 ~ "|" ~  G==g)),
     cex.lab = 1.2,  cex.axis = 1.2, cex = 1.3)
plot(s0.prob.1, s1.prob.1,
     xlab = expression(hat(P)(S^0 == 1 ~ "|" ~  G==g)),
     ylab = expression(hat(P)(S^1 == 1 ~ "|" ~  G==g)),
     cex.lab = 1.2,  cex.axis = 1.2, cex = 1.3)

mod <- lm(s1.prob.1 ~ s0.prob.0 + s0.prob.1 - 1)
summary(mod)              # The coefficients are significant 
est <- mod$coefficients          # estimators of $\pi_{1|0}$ and $\pi_{1|1}$ 
est.all <- c(est, 1 - est)  # the last two are estimators of $\pi_{0|0}$ and $\pi_{0|1}$ 

# estimate the joint distributions P(S0=0, S1=1), P(S0=1, S1=1)
#                                  P(S0=0, S1=0), P(S0=1, S1=0)
prob.s0.s1 <- matrix(nrow = 10, ncol = 4)
for(g in 1:10){
  prob.s0.s1[g, ] <- est.all * c(s0.prob.0[g], s0.prob.1[g], s0.prob.0[g], s0.prob.1[g])
}
prob.s0.s1.vec <- c(prob.s0.s1)



## calculate confidence interval with bootstrap
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # sample size for each trial 

B <- 500
EstAll <- matrix(nrow = B, ncol = 4)
Res <- matrix(nrow = B, ncol = 10)
Est.prob.s0.s1.vec <- matrix(nrow = B, ncol = 40)

set.seed(123)
for(b in 1:B){
  cat('------- bootstrap ', b, ' ---------- \r')
  # step 1: generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # step 2: point estimate based on the boostraped sample 
  s1.prob.1.b <- rep(NA, 10)   
  s0.prob.1.b <- rep(NA, 10)   
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    ng <- dat.tmp$N
    S.g <- dat.tmp$S  
    s1.prob.1.b[g] <- sum(ng[(A.g==1) & (S.g == 1)]) / sum(ng[A.g==1])
    s0.prob.1.b[g] <- sum(ng[(A.g==0) & (S.g == 1)]) / sum(ng[A.g==0])
  }
  s0.prob.0.b <- 1- s0.prob.1.b    # P(Y0 = 0 | G=g) 
  
  mod.b <- lm(s1.prob.1.b ~ s0.prob.0.b + s0.prob.1.b - 1)
  est.b <- mod.b$coefficients
  est.all.b <- c(est.b, 1- est.b) # state transition probabilities
  EstAll[b,] <- est.all.b     
  Res[b, ] <- s1.prob.1.b - est.b[1]*s0.prob.0.b -  est.b[2]*s0.prob.1.b # residual
  
  # joint distributions 
  prob.s0.s1.b <-  matrix(nrow = 10, ncol = 4)
  for(g in 1:10){
    prob.s0.s1.b[g, ] <- est.all.b * c(s0.prob.0.b[g], s0.prob.1.b[g], 
                                       s0.prob.0.b[g], s0.prob.1.b[g])
  }
  prob.s0.s1.vec.b <- c(prob.s0.s1.b)
  Est.prob.s0.s1.vec[b, ] <- prob.s0.s1.vec.b
  }
ese.boot <- apply(EstAll, 2, sd) 
res.ese.boot <- apply(Res, 2, sd) 
prob.s0.s1.vec.est.boot <- apply(Est.prob.s0.s1.vec, 2, sd)
prob.s0.s1.est.boot <- matrix(prob.s0.s1.vec.est.boot, nrow = 10, ncol = 4)


# 95% confidence interval
est.all - 1.96*ese.boot
est.all + 1.96*ese.boot

## Test Assumption 2
stat <- sum( (s1.prob.1 - est[1]*s0.prob.0 -  est[2]*s0.prob.1)^2/res.ese.boot^2 )
1 - pchisq(stat, df = 8)  # p-value     


## Plot the joint distributions
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 1, 2))

plot(1:10, prob.s0.s1[, 3], type = 'b', lty = 5,    # P(S0=0, S1=0|G=g)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.45), 
     cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(hat(P)(S^0 == 0, S^1 ==  0~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.s0.s1[, 3] - 1.96*prob.s0.s1.est.boot[,3],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.s0.s1[, 3] + 1.96*prob.s0.s1.est.boot[,3],  
      lwd = 1.5, col = 'black', lty = 2)


plot(1:10, prob.s0.s1[, 1], type = 'b', lty = 5,    # P(S0=0, S1=1)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.35), 
     cex.axis = 0.9,  lwd = 1.5)    
title(ylab = expression(hat(P)(S^0 == 0, S^1 ==  1~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.s0.s1[, 1] - 1.96*prob.s0.s1.est.boot[,1],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.s0.s1[, 1] + 1.96*prob.s0.s1.est.boot[,1],  
      lwd = 1.5, col = 'black', lty = 2)


plot(1:10, prob.s0.s1[, 4], type = 'b', lty = 5,    # P(S0=1, S1=0)
     xlab = '',
     ylab = '', 
     ylim = c(-0.03, 0.30), 
     cex.axis = 0.9,  lwd = 1.5) 
title(xlab = 'Trial Number',
      ylab = expression(hat(P)(S^0 == 1, S^1 ==  0~"|"~G==g)),
      line = 2, cex.lab = 0.75)
lines(1:10, prob.s0.s1[, 4] - 1.96*prob.s0.s1.est.boot[,4],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.s0.s1[, 4] + 1.96*prob.s0.s1.est.boot[,4],  
      lwd = 1.5, col = 'black', lty = 2)

plot(1:10, prob.s0.s1[, 2], type = 'b', lty = 5,    # P(S0=1, S1=1)
     xlab = '',
     ylab = '', 
     ylim = c(0, 0.80), 
     cex.axis = 0.9,  lwd = 1.5)    
title(xlab = 'Trial Number',
      ylab = expression(hat(P)(S^0 == 1, S^1 ==  1~"|"~G==g)), 
      line = 2, cex.lab = 0.75)
lines(1:10, prob.s0.s1[, 2] - 1.96*prob.s0.s1.est.boot[,2],  
      lwd = 1.5, col = 'black', lty = 2)
lines(1:10, prob.s0.s1[, 2] + 1.96*prob.s0.s1.est.boot[,2],  
      lwd = 1.5, col = 'black', lty = 2)




# ===================  Evaluation of principal surrogate  ===================


# ----------  Method 1: Based on Assumptions 5,7, and 8 (S1>=S0)   -----------

# Step 1: estimate the principal scores P(S0=a, S1=b|G=g), P(Y1 = 1 | G=g),P(Y0 = 1 | G=g)
y1.prob.1 <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
y0.prob.1 <- rep(NA, 10)   # P(Y0 = 1 | G=g) = P(Y = 1 | G=g, A=0)

delta.00 <- rep(NA, 10)    # P(S0 = 0, S1 = 0 | G=g)
delta.01 <- rep(NA, 10)    # P(S0 = 0, S1 = 1 | G=g)
delta.10 <- rep(0,  10)    # P(S0 = 1, S1 = 0 | G=g)
delta.11 <- rep(NA, 10)    # P(S0 = 1, S1 = 1 | G=g)

for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  S.g <- dat.tmp$S         # surrogate
  Y.g <- dat.tmp$Y
  
  y1.prob.1[g] <- sum(N.g[(A.g==1) & (Y.g == 1)]) / sum(N.g[A.g==1])
  y0.prob.1[g] <- sum(N.g[(A.g==0) & (Y.g == 1)]) / sum(N.g[A.g==0])
  
  delta.11[g] <- sum(N.g[(A.g==0) & (S.g == 1)]) / sum(N.g[A.g==0])
  delta.01[g] <- sum(N.g[(A.g==1) & (S.g == 1)]) / sum(N.g[A.g==1]) - delta.11[g]
  delta.00[g] <- 1 - delta.11[g] - delta.01[g]
}


# Step 2: estimate pi.y0.11 = P(Y0=1|S0=1,S1=1) and pi.y1.00=P(Y1=1|S0=0,S1=0)
N <- dat$N
A <- dat$A
S <- dat$S
Y <- dat$Y
pi.y1.00 <- sum(N[(A==1) & (S==0) & (Y==1)]) / sum(N[(A==1) & (S==0)])
pi.y0.11 <- sum(N[(A==0) & (S==1) & (Y==1)]) / sum(N[(A==0) & (S==1)])


# Step 3: estimate 'pi.y0.00','pi.y0.01' and 'pi.y1.01', 'pi.y1.11'
# (a) estimate 'pi.y0.00','pi.y0.01'
mod1 <- lm(y0.prob.1 - pi.y0.11*delta.11 ~ delta.00 + delta.01 -1)
para1 <- c(mod1$coefficients, 0, pi.y0.11)   # 'pi.y0.00','pi.y0.01','pi.y0.10','pi.y0.11'

# (b) estimate 'pi.y1.01', 'pi.y1.11'
mod2 <- lm(y1.prob.1 - pi.y1.00*delta.00 ~ delta.01 + delta.11 - 1)
para2 <- c(pi.y1.00, mod2$coefficients[1], 0, mod2$coefficients[2])

psace <- para2 - para1 
names(psace) <- c('psace.00','psace.01', 'psace.10', 'psace.11')


# Step 4: estimate standard error for PSACE_{ab|g} with bootstrap
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # sample size for each trial 

B <- 500
EstPSACE.all <- matrix(nrow = B, ncol = 4)
set.seed(123)
for(b in 1:B){
  cat('------- bootstrap ', b, ' ---------- \r')
  # 1. generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # 2. point estimate based on the bootstrap sample 
  
  # (a) estimate the principal scores P(S0=a, S1=b|G=g), P(Y1 = 1 | G=g),P(Y0 = 1 | G=g)
  y1.prob.1.b <- rep(NA, 10)   # P(Y1 = 1 | G=g) = P(Y = 1 | G=g, A=1) 
  y0.prob.1.b <- rep(NA, 10)   # P(Y0 = 1 | G=g) = P(Y = 1 | G=g, A=0)
  
  delta.00.b <- rep(NA, 10)    # P(S0 = 0, S1 = 0 | G=g)
  delta.01.b <- rep(NA, 10)    # P(S0 = 0, S1 = 1 | G=g)
  delta.10.b <- rep(0,  10)    # P(S0 = 1, S1 = 0 | G=g)
  delta.11.b <- rep(NA, 10)    # P(S0 = 1, S1 = 1 | G=g)
  
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    n.g <- dat.tmp$N
    S.g <- dat.tmp$S         # surrogate
    Y.g <- dat.tmp$Y
    
    y1.prob.1.b[g] <- sum(n.g[(A.g==1) & (Y.g == 1)]) / sum(n.g[A.g==1])
    y0.prob.1.b[g] <- sum(n.g[(A.g==0) & (Y.g == 1)]) / sum(n.g[A.g==0])
    
    delta.11.b[g] <- sum(n.g[(A.g==0) & (S.g == 1)]) / sum(n.g[A.g==0])
    delta.01.b[g] <- sum(n.g[(A.g==1) & (S.g == 1)]) / sum(n.g[A.g==1]) - delta.11.b[g]
    delta.00.b[g] <- 1 - delta.11.b[g] - delta.01.b[g]
  }
  
  
  # (b) estimate pi.y0.11 = P(Y0=1|S0=1,S1=1) and pi.y1.00=P(Y1=1|S0=0,S1=0)
  N <- dat.b$N
  A <- dat.b$A
  S <- dat.b$S
  Y <- dat.b$Y
  pi.y1.00.b <- sum(N[(A==1) & (S==0) & (Y==1)]) / sum(N[(A==1) & (S==0)])
  pi.y0.11.b <- sum(N[(A==0) & (S==1) & (Y==1)]) / sum(N[(A==0) & (S==1)])
  
  # (c) estimate 'pi.y0.00','pi.y0.01' and 'pi.y1.01', 'pi.y1.11'
  ## estimate 'pi.y0.00','pi.y0.01'
  mod1.b <- lm(y0.prob.1.b - pi.y0.11.b*delta.11.b ~ delta.00.b + delta.01.b -1)
  para1.b <- c(mod1.b$coefficients, 0, pi.y0.11.b)   # 'pi.y0.00','pi.y0.01','pi.y0.10','pi.y0.11'
  
  ## estimate 'pi.y1.01', 'pi.y1.11'
  mod2.b <- lm(y1.prob.1.b - pi.y1.00.b*delta.00.b ~ delta.01.b + delta.11.b - 1)
  para2.b <- c(pi.y1.00.b, mod2.b$coefficients[1], 0, mod2.b$coefficients[2])
  
  psace.b <- para2.b - para1.b 
  EstPSACE.all[b, ] <- psace.b
}
psace.ese.boot <- apply(EstPSACE.all, 2, estimateSd) 

# plots 
psace.low <- psace - 1.96 * psace.ese.boot 
psace.up  <- psace + 1.96 * psace.ese.boot 


par(mfrow = c(4, 1), mar = c(3.5, 3.5, 1, 2))
plot(1:10, rep(psace[1],10), type = 'b', lty = 5,    # PSACE_{00|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(PSACE['00|g']), 
      line = 2, cex.lab = 1)
lines(1:10, rep(psace.low[1],10), lwd = 1.5, col = 'black', lty = 2)
lines(1:10, rep(psace.up[1],10), lwd = 1.5, col = 'black', lty = 2)


plot(1:10, rep(psace[2],10), type = 'b', lty = 5,    # PSACE_{01|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.5, 1.1), cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(PSACE['01|g']), 
      line = 2, cex.lab = 1)
lines(1:10, rep(psace.low[2],10), lwd = 1.5, col = 'black', lty = 2)
lines(1:10, rep(psace.up[2],10), lwd = 1.5, col = 'black', lty = 2)

#rep(psace[3],10)
plot(1:10, rep(NA, 10), type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.25, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(PSACE['10|g']), 
      line = 2, cex.lab = 1)
# lines(1:10, rep(psace.low[3],10), lwd = 1.5, col = 'black', lty = 2)
# lines(1:10, rep(psace.up[3],10), lwd = 1.5, col = 'black', lty = 2)


plot(1:10, rep(psace[4],10), type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.15), cex.axis = 0.9,  lwd = 1.5)
title(ylab = expression(PSACE['11|g']),  xlab = 'Trial Number', 
      line = 2, cex.lab = 1)
lines(1:10, rep(psace.low[4],10), lwd = 1.5, col = 'black', lty = 2)
lines(1:10, rep(psace.up[4],10), lwd = 1.5, col = 'black', lty = 2)




# ---- Method 2: Based on Assumptions 5,6, and S1>=S0, Y1>=Y0  -------------

# Step 1: estimate P(S1, Y1 | G=g) and P(S0, Y0 | G=g)  
s1.y1.prob.00 <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
s1.y1.prob.01 <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
s1.y1.prob.10 <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
s1.y1.prob.11 <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)

s0.y0.prob.00 <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
s0.y0.prob.01 <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
s0.y0.prob.10 <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
s0.y0.prob.11 <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)

for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  S.g <- dat.tmp$S         # surrogate
  Y.g <- dat.tmp$Y
  
  s1.y1.prob.00[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.01[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==1])
  s1.y1.prob.10[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.11[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==1])
  
  s0.y0.prob.00[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.01[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==0])
  s0.y0.prob.10[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.11[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==0])
  
}

# Step 2: estimate the invariant parameters pi.cd.ab = P(S1=c,Y1=d|S0=a,Y0=b) 
b <- c(s1.y1.prob.00, s1.y1.prob.01, s1.y1.prob.10, s1.y1.prob.11)
a <- matrix(0, nrow = 40, ncol = 16)
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
a[1:10,1:4] <- tmp
a[11:20,5:8] <- tmp
a[21:30,9:12] <- tmp
a[31:40,13:16] <- tmp

c <- matrix(0, nrow = 11, ncol = 16)
c[1:3, 2:4] <- diag(c(1,1,1))     # constraints: pi.00.01 = 0, pi.00.10 = 0, pi.00.11 = 0 
c[4:5, 7:8] <- diag(c(1,1))       # constraints: pi.01.10 = 0, pi.01.11 = 0
c[6:7, c(10,12)] <- diag(c(1,1))  # constraints: pi.10.01 = 0, pi.10.11 = 0
c[8, 16] <- 1                     # constraints: pi.11.11 = 1
c[9,  c(1,5,9,13)] <- 1
c[10, c(2,6,10,14)] <- 1
c[11, c(3,7,11,15)] <- 1

d <- c(0,0,0,0,0,0,0, 1, 1,1,1)  

e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
f <- c( rep(0, 16), rep(-1, 16) )

coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)

pi.00.ab <- coef[1:4]
pi.01.ab <- coef[5:8]
pi.10.ab <- coef[9:12]
pi.11.ab <- coef[13:16]


# Step 3: estimate joint distributions P(S0,Y0, S1,Y1 |G=g) for g=1,...,10 
# Step 4: estimate PSACE_{ab|g} 
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
PSACE.all <- data.frame(matrix(nrow = 10, ncol = 4))  # 10 trials
names(PSACE.all) <- c('psace.00', 'psace.01', 'psace.10', 'psace.11') 

for(g in 1:10){
  
  # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
  joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
  names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
  joint_dis_g[9:16, 1] <- 1
  joint_dis_g[c(5:8,13:16), 2] <- 1
  joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
  joint_dis_g[seq(2,16,by=2), 4] <- 1
  
  joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab * tmp[g, ]
  joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab * tmp[g, ]
  joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab * tmp[g, ]
  joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab * tmp[g, ]
  
  # -----   estimate PSACE_{ab|g}  ---------  #
  S0 <- joint_dis_g$S0
  Y0 <- joint_dis_g$Y0
  S1 <- joint_dis_g$S1
  Y1 <- joint_dis_g$Y1
  
  prob_tmp <- joint_dis_g$prob 
  
  # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
  psace.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
  
  psace.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
  
  psace.10 <- 0
  
  psace.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
  
  PSACE.all[g,] <- c(psace.00, psace.01, psace.10, psace.11)
}


# step 5: estimate standard error for PSACE_{ab|g} with bootstrap
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # sample size for each trial 

B <- 500
EstPSACE.all <- matrix(nrow = B, ncol = 40)

set.seed(12)
for(loop in 1:B){
  cat('------- bootstrap ', loop, ' ---------- \r')
  # 1. generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # 2. point estimate based on the bootstrap sample 
  s1.y1.prob.00.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
  s1.y1.prob.01.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
  s1.y1.prob.10.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
  s1.y1.prob.11.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)
  
  s0.y0.prob.00.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
  s0.y0.prob.01.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
  s0.y0.prob.10.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
  s0.y0.prob.11.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)
  
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    n.g <- dat.tmp$N
    S.g <- dat.tmp$S         # surrogate
    Y.g <- dat.tmp$Y
    
    s1.y1.prob.00.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.01.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==1])
    s1.y1.prob.10.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.11.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==1])
    
    s0.y0.prob.00.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.01.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==0])
    s0.y0.prob.10.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.11.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==0])
  }
  
  # estimate the invariant parameters pi.cd.ab = P(S1=c,Y1=d|S0=a,Y0=b)
  
  b <- c(s1.y1.prob.00.b, s1.y1.prob.01.b, s1.y1.prob.10.b, s1.y1.prob.11.b)
  a <- matrix(0, nrow = 40, ncol = 16)
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  a[1:10,1:4] <- tmp
  a[11:20,5:8] <- tmp
  a[21:30,9:12] <- tmp
  a[31:40,13:16] <- tmp
  
  c <- matrix(0, nrow = 11, ncol = 16)
  c[1:3, 2:4] <- diag(c(1,1,1))     # constraints: pi.00.01 = 0, pi.00.10 = 0, pi.00.11 = 0 
  c[4:5, 7:8] <- diag(c(1,1))       # constraints: pi.01.10 = 0, pi.01.11 = 0
  c[6:7, c(10,12)] <- diag(c(1,1))  # constraints: pi.10.01 = 0, pi.10.11 = 0
  c[8, 16] <- 1                     # constraints: pi.11.11 = 1
  c[9,  c(1,5,9,13)] <- 1
  c[10, c(2,6,10,14)] <- 1
  c[11, c(3,7,11,15)] <- 1
  
  d <- c(0,0,0,0,0,0,0, 1, 1,1,1)  
  
  e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
  f <- c( rep(0, 16), rep(-1, 16) )
  
  coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)
  
  pi.00.ab.b <- coef[1:4]
  pi.01.ab.b <- coef[5:8]
  pi.10.ab.b <- coef[9:12]
  pi.11.ab.b <- coef[13:16]
  

  # estimate joint distributions P(S0,Y0, S1,Y1 |G=g) and PSACE_{ab|g} 
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  psace.all <- matrix(nrow = 10, ncol = 4)  # 10 trials
  
  for(g in 1:10){
    
    # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
    joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
    names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
    joint_dis_g[9:16, 1] <- 1
    joint_dis_g[c(5:8,13:16), 2] <- 1
    joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
    joint_dis_g[seq(2,16,by=2), 4] <- 1
    
    joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab.b * tmp[g, ]
    joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab.b * tmp[g, ]
    joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab.b * tmp[g, ]
    joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab.b * tmp[g, ]
    
    # -----   estimate PSACE_{ab|g}  ---------  #
    S0 <- joint_dis_g$S0
    Y0 <- joint_dis_g$Y0
    S1 <- joint_dis_g$S1
    Y1 <- joint_dis_g$Y1
    
    prob_tmp <- joint_dis_g$prob 
    
    # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
    psace.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
    
    psace.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
    
    psace.10 <- 0
    
    psace.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
    
    psace.all[g,] <- c(psace.00, psace.01, psace.10, psace.11)
  }
  
  psace.all.vec <- c(psace.all)
  EstPSACE.all[loop, ] <- psace.all.vec
}
psace.vec.ese.boot <- apply(EstPSACE.all, 2, estimateSd) 
psace.ese.boot <- matrix(psace.vec.ese.boot, nrow = 10, ncol = 4)


# plots 
psace <- PSACE.all    # point estimate
psace.low <- PSACE.all - 1.96 * psace.ese.boot 
psace.up  <- PSACE.all + 1.96 * psace.ese.boot 


par(mfrow = c(4, 1), mar = c(3.5, 3.5, 1, 2))
plot(1:10, psace[, 1], type = 'b', lty = 5,    # PSACE_{00|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['00|g']), 
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 1], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 1], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 2], type = 'b', lty = 5,    # PSACE_{01|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.5, 1.1), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['01|g']), 
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 2], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 2], lwd = 1.5, col = 'black', lty = 2)

#psace[, 3]
plot(1:10, rep(NA, 10), type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.25, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['10|g']), 
      line = 2, cex.lab = 1)
# lines(1:10, psace.low[, 3], lwd = 1.5, col = 'black', lty = 2)
# lines(1:10, psace.up[, 3], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 4], type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.15), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['11|g']),  
      xlab = 'Trial Number', 
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 4], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 4], lwd = 1.5, col = 'black', lty = 2)



# ---- Method 3: Based on Assumptions 5,6, and Y1>=Y0  -------------

# Step 1: estimate P(S1, Y1 | G=g) and P(S0, Y0 | G=g)  
s1.y1.prob.00 <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
s1.y1.prob.01 <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
s1.y1.prob.10 <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
s1.y1.prob.11 <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)

s0.y0.prob.00 <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
s0.y0.prob.01 <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
s0.y0.prob.10 <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
s0.y0.prob.11 <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)

for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  S.g <- dat.tmp$S         # surrogate
  Y.g <- dat.tmp$Y
  
  s1.y1.prob.00[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.01[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==1])
  s1.y1.prob.10[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.11[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==1])
  
  s0.y0.prob.00[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.01[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==0])
  s0.y0.prob.10[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.11[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==0])
  
}

# Step 2: estimate the invariant parameters pi.cd.ab = P(S1=c,Y1=d|S0=a,Y0=b)
#library(lsei)
b <- c(s1.y1.prob.00, s1.y1.prob.01, s1.y1.prob.10, s1.y1.prob.11)
a <- matrix(0, nrow = 40, ncol = 16)
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
a[1:10,1:4] <- tmp
a[11:20,5:8] <- tmp
a[21:30,9:12] <- tmp
a[31:40,13:16] <- tmp

c <- matrix(0, nrow = 8, ncol = 16)
c[1:2, c(2,4)]   <- diag(c(1,1))    # constraints: pi.00.01 = 0, pi.00.11 = 0 
c[3:4, c(10,12)] <- diag(c(1,1))    # constraints: pi.10.01 = 0, pi.10.11 = 0

c[5, c(1,5,9,13)] <- 1
c[6, c(2,6,10,14)] <- 1
c[7, c(3,7,11,15)] <- 1
c[8, c(4,8,12,16)] <- 1

d <- c(0,0,0,0, 1,1,1,1)  

e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
f <- c( rep(0, 16), rep(-1, 16) )

coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)

pi.00.ab <- coef[1:4]
pi.01.ab <- coef[5:8]
pi.10.ab <- coef[9:12]
pi.11.ab <- coef[13:16]



# Step 3: estimate joint distributions P(S0,Y0, S1,Y1 |G=g) for g=1,...,10 
# Step 4: estimate PSACE_{ab|g} 
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
PSACE.all <- data.frame(matrix(nrow = 10, ncol = 4))  # 10 trials
names(PSACE.all) <- c('psace.00', 'psace.01', 'psace.10', 'psace.11') 

for(g in 1:10){
  
  # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
  joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
  names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
  joint_dis_g[9:16, 1] <- 1
  joint_dis_g[c(5:8,13:16), 2] <- 1
  joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
  joint_dis_g[seq(2,16,by=2), 4] <- 1
  
  joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab * tmp[g, ]
  joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab * tmp[g, ]
  joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab * tmp[g, ]
  joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab * tmp[g, ]
  
  # -----   estimate PSACE_{ab|g}  ---------  #
  S0 <- joint_dis_g$S0
  Y0 <- joint_dis_g$Y0
  S1 <- joint_dis_g$S1
  Y1 <- joint_dis_g$Y1
  
  prob_tmp <- joint_dis_g$prob 
  
  # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
  psace.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
  
  psace.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
  
  psace.10 <-  sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])-
    sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])
  
  psace.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
  
  PSACE.all[g,] <- c(psace.00, psace.01, psace.10, psace.11)
}


# step 5: estimate standard error for PSACE_{ab|g} with bootstrap
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # sample size for each trial 

B <- 500
EstPSACE.all <- matrix(nrow = B, ncol = 40)

set.seed(12)
for(loop in 1:B){
  cat('------- bootstrap ', loop, ' ---------- \r')
  # 1. generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # 2. point estimate based on the bootstrap sample 
  s1.y1.prob.00.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
  s1.y1.prob.01.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
  s1.y1.prob.10.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
  s1.y1.prob.11.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)
  
  s0.y0.prob.00.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
  s0.y0.prob.01.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
  s0.y0.prob.10.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
  s0.y0.prob.11.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)
  
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    n.g <- dat.tmp$N
    S.g <- dat.tmp$S         # surrogate
    Y.g <- dat.tmp$Y
    
    s1.y1.prob.00.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.01.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==1])
    s1.y1.prob.10.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.11.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==1])
    
    s0.y0.prob.00.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.01.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==0])
    s0.y0.prob.10.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.11.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==0])
  }
  
  # estimate the invariant parameters pi.cd.ab = P(S1=c,Y1=d|S0=a,Y0=b)
  b <- c(s1.y1.prob.00.b, s1.y1.prob.01.b, s1.y1.prob.10.b, s1.y1.prob.11.b)
  a <- matrix(0, nrow = 40, ncol = 16)
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  a[1:10,1:4] <- tmp
  a[11:20,5:8] <- tmp
  a[21:30,9:12] <- tmp
  a[31:40,13:16] <- tmp
  
  c <- matrix(0, nrow = 8, ncol = 16)
  c[1:2, c(2,4)]   <- diag(c(1,1))    # constraints: pi.00.01 = 0, pi.00.11 = 0 
  c[3:4, c(10,12)] <- diag(c(1,1))    # constraints: pi.10.01 = 0, pi.10.11 = 0
  c[5, c(1,5,9,13)] <- 1
  c[6, c(2,6,10,14)] <- 1
  c[7, c(3,7,11,15)] <- 1
  c[8, c(4,8,12,16)] <- 1
  d <- c(0,0,0,0, 1,1,1,1)  
  
  e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
  f <- c( rep(0, 16), rep(-1, 16) )
  
  coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)
  
  pi.00.ab.b <- coef[1:4]
  pi.01.ab.b <- coef[5:8]
  pi.10.ab.b <- coef[9:12]
  pi.11.ab.b <- coef[13:16]
  
  # estimate joint distributions P(S0,Y0, S1,Y1 |G=g) and PSACE_{ab|g} 
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  psace.all <- matrix(nrow = 10, ncol = 4)  # 10 trials
  
  for(g in 1:10){
    
    # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
    joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
    names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
    joint_dis_g[9:16, 1] <- 1
    joint_dis_g[c(5:8,13:16), 2] <- 1
    joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
    joint_dis_g[seq(2,16,by=2), 4] <- 1
    
    joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab.b * tmp[g, ]
    joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab.b * tmp[g, ]
    joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab.b * tmp[g, ]
    joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab.b * tmp[g, ]
    
    # -----   estimate PSACE_{ab|g}  ---------  #
    S0 <- joint_dis_g$S0
    Y0 <- joint_dis_g$Y0
    S1 <- joint_dis_g$S1
    Y1 <- joint_dis_g$Y1
    
    prob_tmp <- joint_dis_g$prob 
    
    # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
    psace.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
    
    psace.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
    
    psace.10 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])-
      sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])
    
    psace.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
    
    psace.all[g,] <- c(psace.00, psace.01, psace.10, psace.11)
  }
  
  psace.all.vec <- c(psace.all)
  EstPSACE.all[loop, ] <- psace.all.vec
}
psace.vec.ese.boot <- apply(EstPSACE.all, 2, estimateSd) 
psace.ese.boot <- matrix(psace.vec.ese.boot, nrow = 10, ncol = 4)


# plots 
psace <- PSACE.all    # point estimate
psace.low <- PSACE.all - 1.96 * psace.ese.boot 
psace.up  <- PSACE.all + 1.96 * psace.ese.boot 


par(mfrow = c(4, 1), mar = c(3.5, 3.5, 1, 2))
plot(1:10, psace[, 1], type = 'b', lty = 5,    # PSACE_{00|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['00|g']), 
  line = 2, cex.lab = 1)
lines(1:10, psace.low[, 1], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 1], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 2], type = 'b', lty = 5,    # PSACE_{01|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.5, 1.1), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['01|g']), 
  line = 2, cex.lab = 1)
lines(1:10, psace.low[, 2], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 2], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 3], type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.25, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['10|g']), 
  line = 2, cex.lab = 1)
lines(1:10, psace.low[, 3], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 3], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 4], type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '', 
     ylim = c(-0.1, 0.15), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['11|g']),  
  xlab = 'Trial Number', 
  line = 2, cex.lab = 1)
lines(1:10, psace.low[, 4], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 4], lwd = 1.5, col = 'black', lty = 2)


# ---------------  Method 4: Based on Assumptions 5,6  --------------------

# step 1: estimate P(S1, Y1 | G=g) and P(S0, Y0 | G=g)  
s1.y1.prob.00 <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
s1.y1.prob.01 <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
s1.y1.prob.10 <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
s1.y1.prob.11 <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)

s0.y0.prob.00 <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
s0.y0.prob.01 <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
s0.y0.prob.10 <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
s0.y0.prob.11 <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)

for(g in 1:10){
  dat.tmp <- filter(dat, G == g)
  A.g <- dat.tmp$A
  N.g <- dat.tmp$N
  S.g <- dat.tmp$S         # surrogate
  Y.g <- dat.tmp$Y
  
  s1.y1.prob.00[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.01[g] <- N.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==1])
  s1.y1.prob.10[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==1])
  s1.y1.prob.11[g] <- N.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==1])
  
  s0.y0.prob.00[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.01[g] <- N.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(N.g[A.g==0])
  s0.y0.prob.10[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(N.g[A.g==0])
  s0.y0.prob.11[g] <- N.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(N.g[A.g==0])
  
}

# step 2: estimate the invariant parameters P(S1,Y1|S0,Y0) by fit four linear models
#library(lsei)
b <- c(s1.y1.prob.00, s1.y1.prob.01, s1.y1.prob.10, s1.y1.prob.11)
a <- matrix(0, nrow = 40, ncol = 16)
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
a[1:10,1:4] <- tmp
a[11:20,5:8] <- tmp
a[21:30,9:12] <- tmp
a[31:40,13:16] <- tmp

c <- matrix(0, nrow = 4, ncol = 16)
c[1, c(1,5,9,13)] <- 1
c[2, c(2,6,10,14)] <- 1
c[3, c(3,7,11,15)] <- 1
c[4, c(4,8,12,16)] <- 1
d <- c(1,1,1,1)  

e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
f <- c( rep(0, 16), rep(-1, 16) )

coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)

pi.00.ab <- coef[1:4]
pi.01.ab <- coef[5:8]
pi.10.ab <- coef[9:12]
pi.11.ab <- coef[13:16]


# step 3: estimate joint distributions P(S0,Y0, S1,Y1 |G=g) for g=1,...,10 
# step 4: estimate PSACE_{ab|g} 
tmp <- cbind(s0.y0.prob.00, s0.y0.prob.01, s0.y0.prob.10, s0.y0.prob.11)
PSACE.all <- data.frame(matrix(nrow = 10, ncol = 4))  # 10 trials
names(PSACE.all) <- c('PSACE.00', 'PSACE.01', 'PSACE.10', 'PSACE.11') 

for(g in 1:10){
  
  # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
  joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
  names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
  joint_dis_g[9:16, 1] <- 1
  joint_dis_g[c(5:8,13:16), 2] <- 1
  joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
  joint_dis_g[seq(2,16,by=2), 4] <- 1

  joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab * tmp[g, ]
  joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab * tmp[g, ]
  joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab * tmp[g, ]
  joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab * tmp[g, ]
  
  # -----   estimate PSACE_{ab|g}  ---------  #
  S0 <- joint_dis_g$S0
  Y0 <- joint_dis_g$Y0
  S1 <- joint_dis_g$S1
  Y1 <- joint_dis_g$Y1
  
  prob_tmp <- joint_dis_g$prob 
  
  # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
  PSACE.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
  
  PSACE.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
  
  PSACE.10 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])-
    sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])
  
  PSACE.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
    sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
  
  PSACE.all[g,] <- c(PSACE.00, PSACE.01, PSACE.10, PSACE.11)
}


# step 5: estimate standard error for PSACE_{ab|g} with bootstrap
tmp <-  group_by(dat, G) %>% summarize(Ng = sum(N))
Ng <- tmp$Ng     # sample size for each trial 

B <- 500
EstPSACE.all <- matrix(nrow = B, ncol = 40)

set.seed(123)
for(loop in 1:B){
  cat('------- bootstrap ', loop, ' ---------- \r')
  # 1. generate the bootstrap sample 
  N.b <- rep(NA, nrow(dat))
  for(g in 1:10){
    ng <- Ng[g]     # total sample size for trial g 
    ng.breaks <- cumsum( dat$N[dat$G == g] )
    ng.breaks <- c(0, ng.breaks[-8], Inf)
    ind.g.boot <- sample(1:ng, replace = T) # indexes of bootstrap sample for trial g
    ng.boot <- unname(table(cut(ind.g.boot, breaks = ng.breaks)))
    N.b[1:8+(g-1)*8] <-  ng.boot 
  }
  dat.b <- dat
  dat.b$N <- N.b
  
  # 2. point estimate based on the bootstrap sample 
  s1.y1.prob.00.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 0 | G=g)
  s1.y1.prob.01.b <- rep(NA, 10)   # P(S1 = 0, Y1 = 1 | G=g)
  s1.y1.prob.10.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 0 | G=g)
  s1.y1.prob.11.b <- rep(NA, 10)   # P(S1 = 1, Y1 = 1 | G=g)
  
  s0.y0.prob.00.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 0 | G=g)
  s0.y0.prob.01.b <- rep(NA, 10)   # P(S0 = 0, Y0 = 1 | G=g)
  s0.y0.prob.10.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 0 | G=g)
  s0.y0.prob.11.b <- rep(NA, 10)   # P(S0 = 1, Y0 = 1 | G=g)
  
  for(g in 1:10){
    dat.tmp <- filter(dat.b, G == g)
    A.g <- dat.tmp$A
    n.g <- dat.tmp$N
    S.g <- dat.tmp$S         # surrogate
    Y.g <- dat.tmp$Y
    
    s1.y1.prob.00.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.01.b[g] <- n.g[(A.g==1) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==1])
    s1.y1.prob.10.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==1])
    s1.y1.prob.11.b[g] <- n.g[(A.g==1) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==1])
    
    s0.y0.prob.00.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.01.b[g] <- n.g[(A.g==0) & (S.g == 0) & (Y.g == 1)] / sum(n.g[A.g==0])
    s0.y0.prob.10.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 0)] / sum(n.g[A.g==0])
    s0.y0.prob.11.b[g] <- n.g[(A.g==0) & (S.g == 1) & (Y.g == 1)] / sum(n.g[A.g==0])
  }
  
  # estimate the invariant parameters P(S1,Y1|S0,Y0)
  b <- c(s1.y1.prob.00.b, s1.y1.prob.01.b, s1.y1.prob.10.b, s1.y1.prob.11.b)
  a <- matrix(0, nrow = 40, ncol = 16)
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  a[1:10,1:4] <- tmp
  a[11:20,5:8] <- tmp
  a[21:30,9:12] <- tmp
  a[31:40,13:16] <- tmp
  
  c <- matrix(0, nrow = 4, ncol = 16)
  c[1, c(1,5,9,13)] <- 1
  c[2, c(2,6,10,14)] <- 1
  c[3, c(3,7,11,15)] <- 1
  c[4, c(4,8,12,16)] <- 1
  d <- c(1,1,1,1)  
  
  e <- rbind( diag(rep(1, 16)), diag(rep(-1, 16)) )  # constraints: 0 <= x <= 1
  f <- c( rep(0, 16), rep(-1, 16) )
  
  coef <- lsei(a = a, b = b, c = c, d = d, e = e, f = f)
  
  pi.00.ab.b <- coef[1:4]
  pi.01.ab.b <- coef[5:8]
  pi.10.ab.b <- coef[9:12]
  pi.11.ab.b <- coef[13:16]
  
  
  # estimate joint distributions P(S0,Y0, S1,Y1 |G=g) for g=1,...,10 
  # estimate PSACE_{ab|g} 
  tmp <- cbind(s0.y0.prob.00.b, s0.y0.prob.01.b, s0.y0.prob.10.b, s0.y0.prob.11.b)
  psace.all <- matrix(nrow = 10, ncol = 4)  # 10 trials
  
  for(g in 1:10){
    
    # -------- estimate joint distributions P(S0,Y0, S1,Y1 |G=g)  -------- #
    joint_dis_g <- data.frame(matrix(0, nrow = 16, ncol = 5))
    names(joint_dis_g) <- c('S0', 'Y0', 'S1', 'Y1', 'prob')
    joint_dis_g[9:16, 1] <- 1
    joint_dis_g[c(5:8,13:16), 2] <- 1
    joint_dis_g[c(3:4,7:8,11:12,15:16), 3] <- 1
    joint_dis_g[seq(2,16,by=2), 4] <- 1

    joint_dis_g[c(1,5,9,13),  5] <- pi.00.ab.b * tmp[g, ]
    joint_dis_g[c(2,6,10,14), 5] <- pi.01.ab.b * tmp[g, ]
    joint_dis_g[c(3,7,11,15), 5] <- pi.10.ab.b * tmp[g, ]
    joint_dis_g[c(4,8,12,16), 5] <- pi.11.ab.b * tmp[g, ]
    
    
    # -----   estimate PSACE_{ab|g}  ---------  #
    S0 <- joint_dis_g$S0
    Y0 <- joint_dis_g$Y0
    S1 <- joint_dis_g$S1
    Y1 <- joint_dis_g$Y1
    
    prob_tmp <- joint_dis_g$prob 
    
    # P(Y1=1|S0, S1) - P(Y0=1|S0, S1) 
    psace.00 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 0])/sum(prob_tmp[S0 == 0 & S1 == 0])
    
    psace.01 <- sum(prob_tmp[Y1==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 0 & S1 == 1])/sum(prob_tmp[S0 == 0 & S1 == 1])
    
    psace.10 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])-
      sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 0])/sum(prob_tmp[S0 == 1 & S1 == 0])
    
    psace.11 <- sum(prob_tmp[Y1==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])-
      sum(prob_tmp[Y0==1 & S0 == 1 & S1 == 1])/sum(prob_tmp[S0 == 1 & S1 == 1])
    
    psace.all[g,] <- c(psace.00, psace.01, psace.10, psace.11)
  }
  
  psace.all.vec <- c(psace.all)
  EstPSACE.all[loop, ] <- psace.all.vec
}
psace.vec.ese.boot <- apply(EstPSACE.all, 2, estimateSd)
psace.ese.boot <- matrix(psace.vec.ese.boot, nrow = 10, ncol = 4)


# plots 
psace <- PSACE.all    # point estimate
psace.low <- PSACE.all - 1.96 * psace.ese.boot 
psace.up  <- PSACE.all + 1.96 * psace.ese.boot 

par(mfrow = c(4, 1), mar = c(3.5, 3.5, 1, 2))
plot(1:10, psace[, 1], type = 'b', lty = 5,    # PSACE_{00|g}
     xlab = '',
     ylab = '',
     ylim = c(-0.25, 0.25), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['00|g']),
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 1], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 1], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 2], type = 'b', lty = 5,    # PSACE_{01|g}
     xlab = '',
     ylab = '',
     ylim = c(-0.5, 1.1), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['01|g']),
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 2], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 2], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 3], type = 'b', lty = 5,    # PSACE_{10|g}
     xlab = '',
     ylab = '',
     ylim = c(-2, 0.8), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['10|g']),
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 3], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 3], lwd = 1.5, col = 'black', lty = 2)


plot(1:10, psace[, 4], type = 'b', lty = 5,    # PSACE_{11|g}
     xlab = '',
     ylab = '',
     ylim = c(-0.1, 0.15), cex.axis = 0.9,  lwd = 1.5)
title(#ylab = expression(PSACE['11|g']),  
      xlab = 'Trial Number',
      line = 2, cex.lab = 1)
lines(1:10, psace.low[, 4], lwd = 1.5, col = 'black', lty = 2)
lines(1:10, psace.up[, 4], lwd = 1.5, col = 'black', lty = 2)

