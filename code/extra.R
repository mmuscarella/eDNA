### OLD STUFF THAT NEEDS TO BE REMOVED ###


plot(rev(sort(test$Intact)))
plot(rev(sort(test$Relic)))
simp_even(test$Intact)
simp_even(test$Relic)

length(test$Comm)
length(test$Relic)



size.comm <- 1000
props <- seq(size.comm/10, size.comm, 50)


for(i in 1:length(props)){
  n <- props[i]
  # Initialize
  intact <- c(sample(test$Comm, size = size.comm, replace = F))
  total <- c(sample(test$Comm, size = size.comm, replace = F),
             sample(test$Relic, size = n, replace = F))
  ratio <- length(unique(total)) / length(unique(intact))
  print(paste("bias with proportion = ", props[i]/size.comm, " is ", ratio, sep = ""))
}






even.test.R <- matrix(NA, 100, 2)
for(i in 1:dim(even.test.R)[1]){
  test <- relic.pool(names = otus, N = 1000, immigration = 100, immigration2 = 0,
                     birth = 0.1, mortality = 0.1, 
                     decay = 0.1, time = 10^3, random.decay = T)
  even.test.R[i, 1] <- simp_even(table(test$Comm))
  even.test.R[i, 2] <- simp_even(table(test$Relic))
}

even.test.T <- matrix(NA, 100, 2)
colnames()
for(i in 1:dim(even.test.T)[1]){
  test <- relic.pool(names = otus, N = 1000, immigration = 100, immigration2 = 0,
                     birth = 0.1, mortality = 0.1, 
                     decay = 0.1, time = 10^3, random.decay = F)
  even.test.T[i, 1] <- simp_even(table(test$Comm))
  even.test.T[i, 2] <- simp_even(table(test$Relic))
}



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

D <- c(seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.2))
M <- c(seq(0.01, 0.1, by = 0.025), seq(0.1, 0.5, by = 0.2))
i <- c(10)
i2 <- c(0, 2)

combo <- expand.grid(M, D, i, i2) 

turnover <- vector("list", dim(combo)[1])
sim.summary <- as.data.frame(matrix(NA, nrow = dim(combo)[1], ncol = 7))
colnames(sim.summary) <- c("Mortality", "Decay", "Immigration", "Relic Immigration",
                           "Prop", "Rich_Ratio", "Rich_Ratio_S")

for(i in 1:dim(combo)[1]){
  print(paste("Simulation ", i, " of ", dim(combo)[1], sep = ""))
  turnover[[i]] <- relic.pool(names = otus, N = 100, 
                              immigration = combo[i, 3], 
                              immigration2 = combo[i, 4],
                              birth = 0.01, mortality = combo[i, 1], 
                              decay = combo[i, 2], time = 2 * 10^3)
  # Print Statements
  print(paste("Mortality = ", combo[i, 1], sep = ""), quote = F)
  print(paste("Decay = ", combo[i, 2], sep = ""), quote = F)
  print(paste("Prop = ", turnover[[i]]$Prop, sep = ""), quote = F)
  print(paste("Rich_Ratio = ", turnover[[i]]$Rich_Ratio, sep = ""), quote = F)
  print(paste("Rich_Ratio w/ Sampling = ", turnover[[i]]$Rich_Ratio_S, sep = ""), quote = F)
  
  sim.summary[i, ] <- c(combo[i, 1], combo[i, 2], combo[i, 3], combo[i, 4],
                        turnover[[i]]$Prop, turnover[[i]]$Rich_Ratio, 
                        turnover[[i]]$Rich_Ratio_S)
}



sim.summary$MD <- sim.summary$Mortality / sim.summary$Decay
sim.summary[which(sim.summary$Rich_Ratio_S < 1.1), ]

plot(sim.summary$Rich_Ratio_S ~ sim.summary$Mortality)
plot(sim.summary$Rich_Ratio_S ~ sim.summary$Decay)
plot(sim.summary$Prop ~ sim.summary$Mortality)

mod <- lm(sim.summary$Prop ~ sim.summary$Decay * sim.summary$Mortality)
summary(mod)

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



relic.pool <- function(names = otus, N = n, 
                       birth = r, mortality = m, immigration = i, immigration2 = i2, 
                       decay = d, time = t, random.decay = T){
  
  # Define Variables
  otus <- names; N = N; i <- immigration; r <- birth; m <- mortality 
  d <- decay; t <- time; i2 <- immigration2
  
  # Initiate Community
  comm <- sample(otus, size = N, replace = T, prob = rel)
  
  # Run Simulation
  for (j in 1:t){
    
    size.comm <- length(comm)
    
    # Immigration
    i.comm <- sample(otus, size = i, replace = T, prob = rel)
    comm2 <- c(comm, i.comm)
    
    # Birth
    temp.b <- sample(c(1:size.comm), size = (size.comm * r) - length(i.comm), replace = F)
    b.comm <- comm[temp.b]
    comm2 <- c(comm2, b.comm)
    
    # Death
    if (!exists("relic")){
      relic <- vector(mode = "character", length = 0)
      age.relic <- vector(mode = "numeric", length = 0)
    }
    temp.m <- sample(c(1:length(comm2)), size = size.comm * m, replace = F)
    r.comm <- comm2[temp.m]
    relic <- c(relic, r.comm)
    age.relic <- c(age.relic, rep(as.numeric(j), length(r.comm)))
    comm2 <- comm2[-temp.m]
    
    # Immigration to Relic
    if (i2 > 0){
      i.relic <- sample(otus, size = i2, replace = T, prob = rel)
      relic <- c(relic, i.relic)
      age.relic <- c(age.relic, rep(as.numeric(j), length(i.comm)))
    } 
    
    # Decay
    if (length(relic) > 0){
      if (random.decay == T){
        temp.d <- sample(c(1:length(relic)), 
                         size = ((length(relic) ) * d) , replace = F)
      } else {
        age <- j - age.relic
        #age <- seq(1:1000)
        c1 <- log(1/d - 1)
        c2 <- 1/750
        B <- 1
        log.decay <- 1 / (1 + exp(-c1 * (c2 * age - 1)))
        #c1 <- 150; c2 <- 500; B <- 
        log.decay2 <- exp(c1 * (c2  - age))/ (1 + exp(B * c1 * (c2 - age)))
        #plot(log.decay ~ age)
        #plot(log.decay2 ~ age)
        sp.decay.prop <- rep(NA, length(relic))
        for(k in 1:length(sp.decay.prop)){
          sp.decay.prop[k] <- sp.decay[which(names(sp.decay) == relic[k])]
        }
        temp.d <- sample(c(1:length(relic)), 
                         size = ((length(relic) ) * d) , replace = F, prob = sp.decay.prop)
      }
      if (length(temp.d) > 0){
        relic <- relic[-temp.d]
        age.relic <- age.relic[-temp.d]
      }
    }
    
    # Randomize Community
    comm <- sample(comm2)
    
    # Calculations
    N_active  <- length(comm)
    S_active  <- length(unique(comm))
    N_relic   <- length(relic)
    S_relic   <- length(unique(relic))
    Prop_Relic <- round(length(relic) / (length(comm) + length(relic)), 2)
    Rich_Ratio <- round(length(unique(c(comm, relic))) / length(unique(comm)), 2)
    Sample_comm <- sample(comm, length(comm) / 2)
    Sample_total <- sample(c(comm, relic), (length(comm)) / 2)
    Rich_Ratio_S <- round(length(unique(Sample_total)) / 
                            length(unique(Sample_comm)), 2)
    
    # Print Statements
    if (j %in% seq(0, t, 1000)){
      print(paste("N_active = ", N_active, sep = ""), quote = F)
      print(paste("S_active = ", S_active, sep = ""), quote = F)
      print(paste("N_relic = ", N_relic, sep = ""), quote = F)
      print(paste("S_relic = ", S_relic, sep = ""), quote = F)
      print(paste("Prop = ", Prop_Relic, sep = ""), quote = F)
      print(paste("Rich_Ratio = ", Rich_Ratio, sep = ""), quote = F)
      print(paste("Rich_Ratio w/ Sampling = ", Rich_Ratio_S, sep = ""), quote = F)
    }
    
    # Output
    out <- list(Comm = comm, Relic = relic, 
                N_active = N_active, S_active = S_active,
                N_relic = N_relic, S_relic = S_relic,
                Prop = Prop_Relic, Rich_Ratio = Rich_Ratio, 
                Rich_Ratio_S = Rich_Ratio_S)
  }
  return(out)
}

