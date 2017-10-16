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
  
  # Initiate Community
  site1 <- sample(otus, size = N, replace = T, prob = rel)
  
  # Run Simulation
  for (j in 1:t){
    # Immigration
    i.comm <- sample(otus, size = i, replace = T, prob = rel)
    comm <- sample(c(site1, i.comm))
    
    # Birth
    temp <- sample(c(1:length(comm)), size = length(comm) * r, replace = F)
    comm <- sample(c(comm, comm[temp]))
    
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

test <- relic.pool(names = otus, N = 1000, immigration = 100, 
                   birth = 0.01, mortality = 0.05, 
                   decay = 0.001, time = 10^3, random.decay = T)

test <- relic.pool(names = otus, N = 1000, immigration = 100, 
                   birth = 0.01, mortality = 0.02, 
                   decay = 0.005, time = 2 * 10^3, random.decay = F)

# Goal: Keep birth and the same and alter mortality and decay
# Mortality Range 0.001 -> 0.05
# Decay Range: 0.0001 -> 0.05
#D <- c(seq(0.0001, 0.001, by = 0.00025), seq(0.001, 0.01, by = 0.0025), 
#       seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.25))
#M <- c(seq(0.0001, 0.001, by = 0.00025), seq(0.001, 0.01, by = 0.0025), 
#       seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.25))


D <- c(seq(0.001, 0.01, by = 0.0025), 
       seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.25))
M <- c(seq(0.001, 0.01, by = 0.0025), 
       seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.25))

combo <- expand.grid(M, D) 
  
turnover <- vector("list", dim(combo)[1])
sim.summary <- as.data.frame(matrix(NA, nrow = dim(combo)[1], ncol = 5))
colnames(sim.summary) <- c("Mortality", "Decay", "Prop", "Rich_Ratio", "Rich_Ratio_S")

for(i in 1:dim(combo)[1]){
  print(paste("Simulation ", i, " of ", dim(combo)[1], sep = ""))
  turnover[[i]] <- relic.pool(names = otus, N = 1000, immigration = 100, 
                     birth = 0.01, mortality = combo[i, 1], 
                     decay = combo[i, 2], time = 2 * 10^3)
  # Print Statements
  print(paste("Mortality = ", combo[i, 1], sep = ""), quote = F)
  print(paste("Decay = ", combo[i, 2], sep = ""), quote = F)
  print(paste("Prop = ", turnover[[i]]$Prop, sep = ""), quote = F)
  print(paste("Rich_Ratio = ", turnover[[i]]$Rich_Ratio, sep = ""), quote = F)
  print(paste("Rich_Ratio w/ Sampling = ", turnover[[i]]$Rich_Ratio_S, sep = ""), quote = F)
  
  sim.summary[i, ] <- c(combo[i, 1], combo[i, 2], turnover[[i]]$Prop,
                        turnover[[i]]$Rich_Ratio, turnover[[i]]$Rich_Ratio_S)
}

turnover.2 <- vector("list", dim(combo)[1])
sim.summary.2 <- as.data.frame(matrix(NA, nrow = dim(combo)[1], ncol = 5))
colnames(sim.summary.2) <- c("Mortality", "Decay", "Prop", "Rich_Ratio", "Rich_Ratio_S")

for(i in 1:dim(combo)[1]){
  print(paste("Simulation ", i, " of ", dim(combo)[1], sep = ""))
  turnover.2[[i]] <- relic.pool(names = otus, N = 1000, immigration = 100, 
                              birth = 0.01, mortality = combo[i, 1], 
                              decay = combo[i, 2], time = 2 * 10^3, random.decay = F)
  # Print Statements
  print(paste("Mortality = ", combo[i, 1], sep = ""), quote = F)
  print(paste("Decay = ", combo[i, 2], sep = ""), quote = F)
  print(paste("Prop = ", turnover[[i]]$Prop, sep = ""), quote = F)
  print(paste("Rich_Ratio = ", turnover[[i]]$Rich_Ratio, sep = ""), quote = F)
  print(paste("Rich_Ratio w/ Sampling = ", turnover[[i]]$Rich_Ratio_S, sep = ""), quote = F)
  
  sim.summary.2[i, ] <- c(combo[i, 1], combo[i, 2], turnover[[i]]$Prop,
                        turnover[[i]]$Rich_Ratio, turnover[[i]]$Rich_Ratio_S)
}

sim.summary$MD <- sim.summary$Mortality / sim.summary$Decay
sim.summary[which(sim.summary$Rich_Ratio_S < 1.1), ]

plot(sim.summary$Rich_Ratio_S ~ sim.summary$Mortality)
plot(sim.summary$Rich_Ratio_S ~ sim.summary$Decay)
plot(sim.summary$Prop ~ sim.summary$Mortality)

mod <- lm(sim.summary$Prop ~ sim.summary$Decay * sim.summary$Mortality)
summary(mod)

mod.2 <- lm(sim.summary.2$Prop ~ sim.summary.2$Decay * sim.summary.2$Mortality)
summary(mod.2)

col.m <- rep(1, length(sim.summary$Mortality))
for(i in 1:length(col.m)){
  if(sim.summary$Mortality[i] < 0.1){col.m[i] = 2}
  if(sim.summary$Mortality[i] < 0.01){col.m[i] = 3}
  if(sim.summary$Mortality[i] < 0.001){col.m[i] = 4}
}
col.d <- rep(1, length(sim.summary$Decay))
for(i in 1:length(col.d)){
  if(sim.summary$Decay[i] < 0.1){col.d[i] = 2}
  if(sim.summary$Decay[i] < 0.01){col.d[i] = 3}
  if(sim.summary$Decay[i] < 0.001){col.d[i] = 4}
}

plot(sim.summary$Rich_Ratio ~ jitter(sim.summary$Prop), col = col.m, pch = 16, main = "color by mortality")

plot(sim.summary$Rich_Ratio_S ~ sim.summary$Prop, col = col.d, pch = 16, main = "color by decay")
abline(h = 1)

plot(sim.summary.2$Rich_Ratio_S ~ jitter(sim.summary.2$Prop), col = col.m, pch = 16, main = "color by mortality")
abline(h = 1)

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


