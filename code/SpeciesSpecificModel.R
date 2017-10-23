rm(list=ls())
setwd("~/GitHub/relicDNA/code")
library("vegan")
library("png")
library("grid")
set.seed(Sys.time())

# Evenness Function:  Simpsons Evenness
simp_even <- function(SAD = " "){
  SAD <- subset(SAD, SAD > 0)
  S <- length(SAD)
  N <- sum(SAD)
  X <- rep(NA, S)
  for (i in 1:S){
    X[i] <- (SAD[i]*(SAD[i] - 1)) / (N * (N - 1))
  }
  D <- sum(X)
  e_d <- (1/D)/S
  return(e_d)
}

# Define OTUs and Rel Abundance in Regional Pool
otus <- paste("OTU", sprintf("%05d", seq(1:5000)), sep = "")
rel <- rep(NA, length(otus))
rel <- rlnorm(length(otus), meanlog = 1, sdlog = 0.5)
dlnorm(i, meanlog = -1, sdlog = 2.5)
rel <- rel[rev(order(rel))]

gamma1 <- rlnorm(n=10000, meanlog = 1, sdlog = 0.98)
gamma1 <- gamma1[rev(order(gamma1))]

# Plot Gamma 
layout(1)
par(mar = c(5, 5, 3, 2.5))
plot(log10(gamma1), col = "forestgreen",
     main = "Regional Species Abundance Distribution",
     xlab = "Species Index", ylab = "Log Sampling Probability",
     las = 1)


relic.pool <- function(n.gamma = 1000, sd.gamma = 2.5, n = 10000, 
                       r = 0.1, m = 0.1, im = 100, im.2 = 0, 
                       d = 0.1, t = 10^3, uniform.decay = FALSE,
                       a.decay = 0.6, b.decay = 0.6){
  
  if (exists("R")){
    rm(R)
  }
  if (exists("N")){
    rm(N)
  }
  
  # Define OTUs and Rel Abundance in Regional Pool
  otus <- paste("OTU", sprintf("%05d", seq(1:n.gamma)), sep = "")
  rel <- rep(NA, length(otus))
  rel <- rlnorm(length(otus), meanlog = -1, sdlog = sd.gamma)
  rel <- rel[rev(order(rel))] / sum(rel)
  
  # Decay Probability
  if(uniform.decay == TRUE){
    sp.decay <- rep(1, length(otus))
  } else { 
    sp.decay <- round(rbeta(length(otus), shape1 = a.decay, shape2 = b.decay), 3)
    # hist(sp.decay)
    #hist(rbeta(length(otus), shape1 = 6, shape2 = 2))
  }

  # Initiate Community Table
  N <- table(otus) * 0
  
  # Initial Population
  comm <- sample(otus, size = n, replace = T, prob = rel)
  ind <- which(otus %in% comm)
  N[ind] <- N[ind] + table(comm)
  
  # Run Simulation
  for (i in 1:t){
    
    N0 <- sum(N)

    # Immigration
    i.comm <- sample(otus, size = im, replace = T, prob = rel)
    i.ind <- which(otus %in% i.comm)
    N[i.ind] <- N[i.ind] + table(i.comm)
    
    # Birth
    b.comm <- sample(otus, size = N0 * r, replace = T, prob = N/sum(N))
    b.ind <- which(otus %in% b.comm)
    N[b.ind] <- N[b.ind] + table(b.comm)
    
    # Death
    if (!exists("R")){
      R <- table(otus) * 0
    }
    d.comm <- sample(otus, size = N0 * m  + length(i.comm), replace = T, prob = N/sum(N))
    d.ind <- which(otus %in% d.comm)
    R[d.ind] <- R[d.ind] + table(d.comm)
    N[d.ind] <- N[d.ind] - table(d.comm)
    while(min(N) < 0){
      bad.ind <- which(N < 0)
      temp <- N[bad.ind]
      l.temp <- abs(sum(temp))
      N[bad.ind] <- N[bad.ind] - temp
      R[bad.ind] <- R[bad.ind] + temp
      d.comm2 <- sample(otus, size = l.temp, replace = T, prob = N/sum(N))
      d.ind2 <- which(otus %in% d.comm2)
      R[d.ind2] <- R[d.ind2] + table(d.comm2)
      N[d.ind2] <- N[d.ind2] - table(d.comm2)
      R[which(R < 0)] <- 0
    }

    # Immigration to Relic
    i.relic <- vector(mode = "character", length = 0)
    if (im.2 > 0){
      i.relic <- sample(otus, size = im.2, replace = T)
      relic.ind <- which(otus %in% i.relic)
      R[relic.ind] <- R[relic.ind] + table(i.relic)
    } 
    
    # Decay
    if (sum(R) > 0){
      d.prob <- (sp.decay * (R/sum(R)))
      d.relic <- sample(otus, size = sum(R) * d + length(i.relic), replace = T, prob = d.prob)
      decay.ind <- which(otus %in% d.relic)
      R[decay.ind] <- R[decay.ind] - table(d.relic)
      while(min(R) < 0){
        bad.ind <- which(R < 0)
        temp <- R[bad.ind]
        l.temp <- abs(sum(temp))
        R[bad.ind] <- R[bad.ind] - temp
        d.relic2 <- sample(otus, size = l.temp, replace = T, prob = (R/sum(R)))
        relic2.ind <- which(otus %in% d.relic2)
        R[relic2.ind] <- R[relic2.ind] - table(d.relic2)
      }
    }
    
    # Upper Cap on Relic Size
    # if (sum(R) > (n * 10^3)){
    #   max.relic <- sample(otus, size = (n * 10^3), replace = T, prob = R/sum(R))
    #   max.ind <- which(otus %in% max.relic)
    #   R <- R * 0
    #   R[max.ind] <- R[max.ind] + table(max.relic)
    # }
 
    # Calculations
    N_intact  <- sum(N)
    S_intact  <- sum(N > 0)
    N_relic   <- sum(R)
    S_relic   <- sum(R > 0)
    Prop_Relic <- round(sum(R) / (sum(R) + sum(N)), 2)
    
    Total <- N + R
    Rich_Ratio_C <- round(sum(Total > 0) / sum(N > 0), 2)
    
    Rich_Ratio_S_raw <- rep(NA, 20)
    for(j in 1:20){
      N_Sample <- sample(otus, size = n * 0.5, replace = T, prob = N/sum(N))
      T_Sample <- sample(otus, size = n * 0.5, replace = T, prob = Total/sum(Total))
      Rich_Ratio_S_raw[j] <- round(length(unique(T_Sample)) / length(unique(N_Sample)), 2)
    }
    Rich_Ratio_S <- round(mean(Rich_Ratio_S_raw), 2)

    # Print Statements
    if (i %in% seq(0, t, 1000)){
      print(paste("N_intact = ", N_intact, sep = ""), quote = F)
      print(paste("S_intact = ", S_intact, sep = ""), quote = F)
      print(paste("N_relic = ", N_relic, sep = ""), quote = F)
      print(paste("S_relic = ", S_relic, sep = ""), quote = F)
      print(paste("Prop = ", Prop_Relic, sep = ""), quote = F)
      print(paste("Rich_Ratio = ", Rich_Ratio_C, sep = ""), quote = F)
      print(paste("Rich_Ratio w/ Sampling = ", Rich_Ratio_S, sep = ""), quote = F)
    }
    
    # Output
    out <- list(Intact = N, Relic = R, 
                N_intact = N_intact, S_intact = S_intact,
                N_relic = N_relic, S_relic = S_relic,
                Prop = Prop_Relic, Rich_Ratio_C = Rich_Ratio_C, 
                Rich_Ratio_S = Rich_Ratio_S)
  }
  return(out)
}

# Test Run
test <- relic.pool(n.gamma = 5000, sd.gamma = 2.5, n = 10000, 
                   r = 0.1, m = 0.1, im = 100, im.2 = 0,
                   d = 0.001, t = 10^3, uniform.decay = FALSE,
                   a.decay = 2, b.decay = 4)


##########
# Shape Options
otus <- paste("OTU", sprintf("%05d", seq(1:1000)), sep = "")
sp.decay <- rbeta(length(otus), shape1 = 0.7, shape2 = 0.7)
hist(sp.decay, xlim = c(0, 1))

d <- c(seq(0.02, 0.08, length = 3), seq(0.1, 0.4, length = 3))


d <- c(0.02, 0.05, 0.075, 0.10, 0.15, 0.25, 0.40)
uniform.test <- matrix(NA, nrow = length(d), ncol = 2)
for(i in 1:length(d)){
  test <- relic.pool(n.gamma = 1000, sd.gamma = 2, n = 20000, 
                     r = 0.1, m = 0.1, im = 2000, im.2 = 0,
                     d = d[i], t = 10^3, uniform.decay = TRUE,
                     a.decay = 2, b.decay = 2)
  uniform.test[i, 1] <- test$Prop
  uniform.test[i, 2] <- test$Rich_Ratio_S
}

imm.test <- matrix(NA, nrow = length(d), ncol = 2)
for(i in 1:length(d)){
  test <- relic.pool(n.gamma = 1000, sd.gamma = 2, n = 20000, 
                     r = 0.1, m = 0.1, im = 2000, im.2 = 200,
                     d = d[i], t = 10^3, uniform.decay = TRUE,
                     a.decay = 2, b.decay = 2)
  imm.test[i, 1] <- test$Prop
  imm.test[i, 2] <- test$Rich_Ratio_S
}


beta1 <- matrix(NA, nrow = length(d), ncol = 2)
for(i in 1:length(d)){
  test <- relic.pool(n.gamma = 1000, sd.gamma = 2, n = 20000, 
                     r = 0.1, m = 0.1, im = 2000, im.2 = 0,
                     d = d[i], t = 10^3, uniform.decay = FALSE,
                     a.decay = 1, b.decay = 1)
  beta1[i, 1] <- test$Prop
  beta1[i, 2] <- test$Rich_Ratio_S
}


beta2 <- matrix(NA, nrow = length(d), ncol = 2)
for(i in 1:length(d)){
  test <- relic.pool(n.gamma = 1000, sd.gamma = 2, n = 20000, 
                     r = 0.1, m = 0.1, im = 2000, im.2 = 0,
                     d = d[i], t = 10^3, uniform.decay = FALSE,
                     a.decay = 0.7, b.decay = 0.7)
  beta2[i, 1] <- test$Prop
  beta2[i, 2] <- test$Rich_Ratio_S
}


plot(uniform.test[, 1], uniform.test[, 2], pch = 23, bg = "white",
     ylim = c(0.5, 1.5), xlim = c(0.15, 0.95), las = 1,
     xlab = "Proportion Relic", ylab = "Richness Ratio")
abline(h = 1)
points(imm.test[, 1], imm.test[, 2], pch = 21, bg = "gray")
points(beta1[, 1], beta1[, 2], pch = 24, bg = "white")
points(beta2[, 1], beta2[, 2], pch = 22, bg = "gray20")
legend("bottomleft", c("Neutral", "Protection", "Relic Imm.", "Protection 2"), 
       pch = c(23, 22, 21, 24),
       pt.bg = c("white", "gray20", "gray", "white"), bty = "n")






