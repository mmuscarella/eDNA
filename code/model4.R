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

dlnorm(20, meanlog = -1, sdlog = 2.5)
gamma1 <- rlnorm(n=1000000, meanlog = -1, sdlog = 2.5) # Changed sdlog param
gamma1 <- gamma1[rev(order(gamma1))]

# Simple Plot
layout(1)
par(mar = c(5, 5, 3, 2.5))
plot(log10(rel), col = "forestgreen",
     main = "Regional Species Abundance Distribution",
     xlab = "Species Index", ylab = "Log Sampling Probability",
     las = 1)

# Define OTUs and Rel Abundance in Regional Pool
otus <- paste("OTU", sprintf("%05d", seq(1:10000)), sep = "")
rel <- rep(NA, length(otus))
for(i in 1:length(otus)){
  rel[i] <- dlnorm(i, meanlog = -1, sdlog = 2)
}

# Plot Gamma 
layout(1)
par(mar = c(5, 5, 3, 2.5))
plot(log10(rel), col = "forestgreen",
     main = "Regional Species Abundance Distribution",
     xlab = "Species Index", ylab = "Log Sampling Probability",
     las = 1)

# Examples
N <- 100000
i <- 1000
r <- 0.001
m <- 0.001
d <- 0.0001
t <- 10^3


relic.pool <- function(names = otus, N = n, 
                       birth = r, mortality = m, immigration = i, 
                       decay = d, time = t, random.decay = T){
  # Define Variables
  otus <- names; N = N; i <- immigration; r <- birth; m <- mortality 
  d <- decay; t <- time
  
  # Define Birth Rate
  rel <- rep(NA, length(otus))
  for(i in 1:length(otus)){
    rel[i] <- dlnorm(i, meanlog = -1, sdlog = 2)
  }
  
  # Initiate Community
  #site1 <- sample(otus, size = N, replace = T, prob = rel)
  site1 <- table(otus)
  
  # Run Simulation
  for (j in 1:t){
    # Immigration
    i.comm <- sample(otus, size = i, replace = T, prob = rel)
    comm <- sample(c(site1, i.comm))
    
    # Birth
    temp <- sample(seq_along(otus), size = sum(site1) * r, replace = T, prob = rel)
    comm <- sample(c(comm, comm[temp]))
  }
    
    # Decay
    if (exists("relic")){
      if (random.decay == T){
        temp <- sample(c(1:length(relic)), size = length(relic) * d, replace = F)
      } else {
        temp <- c(1:(length(relic) * d))
      }
      if (length(temp) > 0){
        relic <- relic[-temp]
      }
    }
    
    # Death
    if (!exists("relic")){
      relic <- vector(mode = "character", length = 0)
    }
    temp <- sample(c(1:length(comm)), size = length(comm) * m, replace = F)
    dead <- comm[temp]
    comm <- sample(comm[-temp])
    relic <- c(relic, dead)
    
    # Calculations
    N_active  <- length(comm)
    S_active  <- length(unique(comm))
    N_relic   <- length(relic)
    S_relic   <- length(unique(relic))
    Prop_Relic <- round(length(relic) / (length(comm) + length(relic)), 2)
    Rich_Ratio <- round(length(unique(c(comm, relic))) / length(unique(comm)), 2)
    Sample_comm <- sample(comm, length(comm) / 5)
    Sample_total <- sample(c(relic, comm), (length(relic) + length(comm)) / 5)
    Rich_Ratio_S <- round(length(unique(Sample_total)) / 
                            length(unique(Sample_comm)), 2)
    #Sample_comm2 <- sample(comm, 1000, replace = T)
    #Sample_total2 <- sample(c(relic, comm), 1000, replace = T)
    #Rich_Ratio_S2 <- round(length(unique(Sample_total2)) / 
    #                        length(unique(Sample_comm2)), 2)
    
    # Print Statements
    if (j %in% seq(0, t, 1000)){
      print(paste("N_active = ", N_active, sep = ""), quote = F)
      print(paste("S_active = ", S_active, sep = ""), quote = F)
      print(paste("N_relic = ", N_relic, sep = ""), quote = F)
      print(paste("S_relic = ", S_relic, sep = ""), quote = F)
      print(paste("Prop = ", Prop_Relic, sep = ""), quote = F)
      print(paste("Rich_Ratio = ", Rich_Ratio, sep = ""), quote = F)
      print(paste("Rich_Ratio w/ Sampling = ", Rich_Ratio_S, sep = ""), quote = F)
      #print(paste("Rich_Ratio w/ Sampling #2 = ", Rich_Ratio_S2, sep = ""), quote = F)
    }
    
    # Output
    out <- list(Comm = comm, Relic = relic, 
                N_active = N_active, S_active = S_active,
                N_relic = N_relic, S_relic = S_relic,
                Prop = Prop_Relic, Rich_Ratio = Rich_Ratio, 
                Rich_Ratio_S = Rich_Ratio_S) #, Rich_Ratio_S2 = Rich_Ratio_S2)
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

# Goal: Explore Ways to Change SAD

test <- relic.pool(names = otus, N = 100000, immigration = 1000, 
           birth = 0.01, mortality = 0.01, 
           decay = 0.01, time = 10^4)


plot(log10(gamma1), col = "forestgreen",
     main = "Regional Species Abundance Distribution",
     xlab = "Species Index", ylab = "Log Sampling Probability",
     las = 1)
points



D <- c(seq(0.001, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.01))
turnover <- vector("list", length(D))

for(i in 1:length(D)){
  print(paste("Simulation", i, " of ", length(D), sep = ""))
  turnover[[i]] <- relic.pool(names = otus, N = 100000, immigration = 1000, 
                              birth = 0.01, mortality = 0.01, 
                              decay = D[i], time = 10^3)
}