rm(list=ls())
setwd("~/GitHub/relicDNA/code")
library("vegan")
library("png")
library("grid")
set.seed(156)

# Evenness Functions
# Calculates Smith and Wilson's evenness index - E var
#   Smith & Wilson (1996) A consumer's guide to evenness indices. Oikos
e_var <- function(SAD = " "){
  SAD <- subset(SAD, SAD > 0)
  P <- log(SAD)
  S <- length(SAD)
  X <- rep(NA, S)
  for (i in 1:S){
    X[i] <- ((P[i] - mean(P))^2)/S
  }
  evar <- 1 - (2/(pi * atan(sum(X))))
  return(evar)
}

# Calculates Simpsons Evenness
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


gamma1 <- rlnorm(n=100000, meanlog = 1, sdlog = 0.98)
gamma1 <- gamma1[rev(order(gamma1))]

# Simple Plot
layout(1)
par(mar = c(5, 5, 3, 2.5))
plot(log10(gamma1), col = "forestgreen",
     main = "Regional Species Abundance Distribution",
     xlab = "Species Index", ylab = "Log Sampling Probability",
     las = 1)

# Define OTUs
otus <- paste("OTU", sprintf("%05d", seq(1:length(gamma1))), sep = "")

# Examples
N <- 100000
i <- 1000
r <- 0.001
m <- 0.001
d <- 0.0001
t <- 10^4


relic.pool <- function(names = otus, N = n, 
                       birth = r, mortality = m, immigration = i, 
                       decay = d, time = t){
  # Define Variables
  otus <- names; N = N; i <- immigration; r <- birth; m <- mortality 
  d <- decay; t <- time
  
  # Initiate Community
  site1 <- sample(otus, size = N, replace = T, prob = gamma1)

  # Run Simulation
  for (i in 1:t){
    i.comm <- sample(otus, size = i, replace = T, prob = gamma1)
    comm <- c(site1, i.comm)
    temp <- sample(c(1:length(comm)), size = length(comm) * r, replace = F)
    comm <- c(comm, comm[temp])
    if (!exists("relic")){
      relic <- vector(mode = "character", length = 0)
    }
    temp <- sample(c(1:length(comm)), size = length(comm) * m, replace = F)
    dead <- comm[temp]
    comm <- comm[-temp]
    relic <- c(relic, dead)
    temp <- sample(c(1:length(relic)), size = length(relic) * d, replace = F)
    if (length(temp) > 0){
      relic <- relic[-temp]
    }
    if (i %in% seq(0, t, 100)){
      print(paste("N_active = ", length(comm), sep = ""), quote = F)
      print(paste("S_active = ", length(unique(comm)), sep = ""), quote = F)
      print(paste("N_relic = ", length(relic), sep = ""), quote = F)
      print(paste("S_relic = ", length(unique(relic)), sep = ""), quote = F)
    }
    out <- list(comm, relic, 
                N_active = length(comm), S_active = length(unique(comm)),
                N_relic = length(relic), S_relic = length(unique(relic)))
  }
  return(out)
}

# Same Rates
# Examples
N <- 100000
i <- 1000
r <- 0.001
m <- 0.001
d <- 0.0001
t <- 10^4


same.rates <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                       birth = 0.01, mortality = 0.01, decay = 0.01, time = 10^4)

active1 <- table(same.rates[[1]])
relic1 <- table(same.rates[[2]])

active1 <- active1[rev(order(active1))]
relic1 <- relic1[rev(order(relic1))]

# Simple Plot
layout(1)
par(mar = c(3, 5, 3, 2.5))
plot(log10(active1), col = "forestgreen",
     main = "Local Species Abundance Distributions",
     xlab = "", ylab = "", 
     xlim = c(0, max(length(active1), length(relic1))),
     ylim = c(0, max(max(log10(active1)), max(log10(relic1)))),
     las = 1, type = "p", xaxt = "n", yaxt = "n", pch = 16)

points(log10(relic1), col = "darkblue", type = "p", pch = 16)
legend("topright", legend = c("Active", "Relic"), pch = 16, 
       bty = "n", col = c("forestgreen", "darkblue"))
axis(side = 1, labels = F, lwd = 1.5, tck = 0)
axis(side = 2, labels = T, lwd = 1.5, las = 1)
mtext("Log Abundance", side = 2, cex = 1.5, line = 2.5)
mtext("Species Index", side = 1, cex = 1.5, line = 0.75)
box(lwd = 1.5)


# Low Decay
low.decay <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                         birth = 0.01, mortality = 0.01, decay = 0.001, time = 10^4)

active1 <- table(low.decay[[1]])
relic1 <- table(low.decay[[2]])

active1 <- active1[rev(order(active1))]
relic1 <- relic1[rev(order(relic1))]

# Simple Plot
layout(1)
par(mar = c(3, 5, 3, 2.5))
plot(log10(active1), col = "forestgreen",
     main = "Local Species Abundance Distributions",
     xlab = "", ylab = "", 
     xlim = c(0, max(length(active1), length(relic1))),
     ylim = c(0, max(max(log10(active1)), max(log10(relic1)))),
     las = 1, type = "p", xaxt = "n", yaxt = "n", pch = 16)

points(log10(relic1), col = "darkblue", type = "p", pch = 16)
legend("topright", legend = c("Active", "Relic"), pch = 16, 
       bty = "n", col = c("forestgreen", "darkblue"))
axis(side = 1, labels = F, lwd = 1.5, tck = 0)
axis(side = 2, labels = T, lwd = 1.5, las = 1)
mtext("Log Abundance", side = 2, cex = 1.5, line = 2.5)
mtext("Species Index", side = 1, cex = 1.5, line = 0.75)
box(lwd = 1.5)

# High Decay
high.decay <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                         birth = 0.01, mortality = 0.01, decay = 0.1, time = 10^4)

active1 <- table(high.decay[[1]])
relic1 <- table(high.decay[[2]])

active1 <- active1[rev(order(active1))]
relic1 <- relic1[rev(order(relic1))]

# Simple Plot
layout(1)
par(mar = c(3, 5, 3, 2.5))
plot(log10(active1), col = "forestgreen",
     main = "Local Species Abundance Distributions",
     xlab = "", ylab = "", 
     xlim = c(0, max(length(active1), length(relic1))),
     ylim = c(0, max(max(log10(active1)), max(log10(relic1)))),
     las = 1, type = "p", xaxt = "n", yaxt = "n", pch = 16)

points(log10(relic1), col = "darkblue", type = "p", pch = 16)
legend("topright", legend = c("Active", "Relic"), pch = 16, 
       bty = "n", col = c("forestgreen", "darkblue"))
axis(side = 1, labels = F, lwd = 1.5, tck = 0)
axis(side = 2, labels = T, lwd = 1.5, las = 1)
mtext("Log Abundance", side = 2, cex = 1.5, line = 2.5)
mtext("Species Index", side = 1, cex = 1.5, line = 0.75)
box(lwd = 1.5)

# Low Decay
high.birth <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                        birth = 0.05, mortality = 0.01, decay = 0.01, time = 10^4)

active1 <- table(high.birth[[1]])
relic1 <- table(high.birth[[2]])

active1 <- active1[rev(order(active1))]
relic1 <- relic1[rev(order(relic1))]

# Simple Plot
layout(1)
par(mar = c(5, 5, 3, 2.5))
plot(log10(active1), col = "forestgreen",
     main = "Local Species Abundance Distributions",
     xlab = "", ylab = "", 
     xlim = c(0, max(length(active1), length(relic1))),
     ylim = c(0, max(max(log10(active1)), max(log10(relic1)))),
     las = 1, type = "p", xaxt = "n", yaxt = "n", pch = 16)

points(log10(relic1), col = "darkblue", type = "p", pch = 16)
legend("topright", legend = c("Active", "Relic"), pch = 16, 
       bty = "n", col = c("forestgreen", "darkblue"))
axis(side = 1, labels = F, lwd = 1.5, tck = 0)
axis(side = 2, labels = T, lwd = 1.5, las = 1)
mtext("Log Abundance", side = 2, cex = 1.5, line = 2.5)
mtext("Species Index", side = 1, cex = 1.5, line = 0.5)
box(lwd = 1.5)

# Save Outputs
save(same.rates, low.decay, high.decay, high.birth, file = "../data/simSAD.R")
