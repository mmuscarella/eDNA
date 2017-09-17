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
    out <- list(Comm = comm, Relic = relic, 
                N_active = length(comm), S_active = length(unique(comm)),
                Evar_active = e_var(table(comm)), 
                Esim_active = simp_even(table(comm)),
                N_relic = length(relic), S_relic = length(unique(relic)),
                Evar_relic = e_var(table(comm)), 
                Esim_active = simp_even(table(comm)))
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

# Goal: Keep birth and death the same and just alter decay
# Decay Range: 0.001 -> 0.01

D <- c(seq(0.001, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.01))
turnover <- vector("list", length(D))

for(i in 1:length(D)){
  print(paste("Simulation", i, " of ", length(D), sep = ""))
  turnover[[i]] <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                     birth = 0.01, mortality = 0.01, 
                     decay = D[i], time = 10^4)
}

save(turnover, file = "../data/Turnover.RData")
load(file = "../data/Turnover.RData")

str(turnover)
length(turnover)

for (i in 1:length(turnover)){
  turnover[[i]]$Esim_active <- simp_even(table(turnover[[i]][1]))
  turnover[[i]]$Esim_relic <- simp_even(table(turnover[[i]][2]))
}


prop.relic <- rep(NA, length(D))
for(i in 1:length(D)){
  prop.relic[i] <- turnover[[i]]$N_relic / turnover[[i]]$N_active
}


plot(log10(D), log10(prop.relic),
     xlab = "", ylab = "", las = 1)
mtext(expression(paste("Turnover Rate")),
      side = 1, line = 2.75, cex = 1.25)
mtext(expression(paste("(log" [10], ")")),
      side = 1, line = 3.75)
mtext(expression(paste("Proportion Relic DNA")),
      side = 2, line = 3.5, cex = 1.25)
mtext(expression(paste("(log" [10], ")")),
      side = 2, line = 2.5)


rich.ratio <- rep(NA, length(D))
rich.ratio2 <- rep(NA, length(D))
for(i in 1:length(D)){
  active <- turnover[[i]][1]
  relic <- turnover[[i]][2]
  total <- c(unlist(turnover[[i]][1]), unlist(turnover[[i]][2]))
  active2 <- sample(unlist(active), 5000)
  relic2 <- sample(unlist(relic), 5000)
  total2 <- sample(unlist(total), 5000)
  t.total <- table(total)
  t.active <- table(active)
  t.relic <- table(relic)
  t.total2 <- table(total2)
  t.active2 <- table(active2)
  t.relic2 <- table(relic2)
  rich.ratio[i] <- length(t.total) / length(t.active)
  rich.ratio2[i] <- length(t.total2) / length(t.active2)
}

plot(log10(prop.relic), rich.ratio, pch = 21, bg = "gray60", 
     xlab = "", ylab = "", las = 1, ylim = c(0.8, 1.4))
points(log10(prop.relic), rich.ratio2, pch = 22, bg = "gray40")
legend("topleft", c("No Sampling", "With Sampling"), 
       pch = c(21, 22), pt.bg = c("gray60", "gray40"), bty = "n")
mtext(expression(paste("Turnover Rate")),
      side = 1, line = 2.75, cex = 1.25)
mtext(expression(paste("(log" [10], ")")),
      side = 1, line = 3.75)
mtext(expression(paste("Richness Ratio")),
      side = 2, line = 3.5, cex = 1.25)



even.ratio <- rep(NA, length(D))
even.ratio2 <- rep(NA, length(D))
for(i in 1:length(D)){
  active <- turnover[[i]][1]
  relic <- turnover[[i]][2]
  total <- c(unlist(turnover[[i]][1]), unlist(turnover[[i]][2]))
  active2 <- sample(unlist(active), 5000)
  relic2 <- sample(unlist(relic), 5000)
  total2 <- sample(unlist(total), 5000)
  t.total <- table(total)
  t.active <- table(active)
  t.relic <- table(relic)
  t.total2 <- table(total2)
  t.active2 <- table(active2)
  t.relic2 <- table(relic2)
  even.ratio[i] <- simp_even(t.total) / simp_even(t.active)
  even.ratio2[i] <- simp_even(t.total2) / simp_even(t.active2)
}

plot(log10(prop.relic), even.ratio, pch = 21, bg = "gray60", 
     xlab = "", ylab = "", las = 1, ylim = c(0.5, 1.2))
points(log10(prop.relic), even.ratio2, pch = 22, bg = "gray40")
legend("topright", c("No Sampling", "With Sampling"), 
       pch = c(21, 22), pt.bg = c("gray60", "gray40"), bty = "n")
mtext(expression(paste("Turnover Rate")),
      side = 1, line = 2.75, cex = 1.25)
mtext(expression(paste("(log" [10], ")")),
      side = 1, line = 3.75)
mtext(expression(paste("Evenness Ratio")),
      side = 2, line = 3.5, cex = 1.25)


