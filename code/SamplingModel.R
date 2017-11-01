################################################################################
#
# Sampling Model To Test Bias From Relic DNA
#
#   The effects of mixing species distributions with different S/E
#
################################################################################
#
# Written by: Mario Muscarella
#
# Last update: 2017/10/31
#
################################################################################

# Initial Setup
rm(list=ls())
setwd("~/GitHub/relicDNA/code")
library("vegan")
library("png")
library("grid")
set.seed(22)

################################################################################
# Different Evenness
################################################################################

gamma1 <- rlnorm(n=10000, meanlog = 1, sdlog = 0.98)
gamma1 <- gamma1[rev(order(gamma1))]
gamma2 <- rlnorm(n=10000, meanlog = 1, sdlog = 1.8)
gamma2 <- gamma2[rev(order(gamma2))]
gamma3 <- rlnorm(n=10000, meanlog = 1, sdlog = 0.25)
gamma3 <- gamma3[rev(order(gamma3))]

# Define OTUs
otus <- paste("OTU", sprintf("%05d", seq(1:length(gamma1))), sep = "")

# Initiate Communities
site1 <- sample(otus, size = 100000, replace = T, prob = gamma1)
site2 <- sample(otus, size = 100000, replace = T, prob = gamma2)
site3 <- sample(otus, size = 100000, replace = T, prob = gamma3)

site1.t <- table(site1); site2.t <- table(site2); site3.t <- table(site3)

# Plot Species Abundance Distrubtions
layout(matrix(1:3, ncol = 3))
par(mar = c(4, 2, 4, 1), oma = c(0, 3, 1, 1))
plot(log10(site1.t[rev(order(site1.t))]),
     ylim = c(0, 4), xlim = c(0, 10000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 2, "Log Abundance", line = 3)
mtext(side = 3, "Intact\nCommunity", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")
plot(log10(site2.t[rev(order(site2.t))]), 
     ylim = c(0, 4), xlim = c(0, 10000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 1, "Species Index (OTU)", line = 1.5)
mtext(side = 3, "Species Abundance Distributions", line = 1.5)
mtext(side = 3, "Less Even\nRelic", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")
plot(log10(site3.t[rev(order(site3.t))]), 
     ylim = c(0, 4), xlim = c(0, 10000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 3, "More Even\nRelic", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")

# Basic Statistics
specnumber(site1.t)
specnumber(site2.t)
specnumber(site3.t)

# Sampling Simulation
size.comm <- 1000000
props <- seq(size.comm/100, size.comm, 50000)

ratio.tab <- matrix(NA, nrow = 5, ncol = length(props))
rownames(ratio.tab) <- c("Null", "Uneven", "Even", "Same", "Proportion")

bc.tab <- matrix(NA, nrow = 4, ncol = length(props))
rownames(bc.tab) <- c("1-2", "1-3", "1-4", "Proportion")

for(i in 1:length(props)){
  n <- props[i]
  # Initialize
  temp.null <- c(sample(otus, size = size.comm, replace = T, prob = gamma1))
  temp.uneven <- c(sample(otus, size = size.comm, replace = T, prob = gamma1),
                   sample(otus, size = n, replace = T, prob = gamma2))
  temp.even <- c(sample(otus, size = size.comm, replace = T, prob = gamma1),
                 sample(otus, size = n, replace = T, prob = gamma3))
  temp.same <- c(sample(otus, size = size.comm, replace = T, prob = gamma1),
                 sample(otus, size = n, replace = T, prob = gamma1))
  tab1 <- table(temp.null)
  tab2 <- table(temp.uneven)
  tab3 <- table(temp.even)
  tab4 <- table(temp.same)
  
  sbys <- matrix(NA, ncol = length(otus), nrow = 4)
  colnames(sbys) <- otus
  rownames(sbys) <- c("site1", "site2", "site3",
                      "site4")
  
  for(j in 1:dim(sbys)[2]){
    if (otus[j] %in% names(tab1)){
      sbys[1, j] <- tab1[which(names(tab1) == otus[j])]
    } else {sbys[1, j] <- 0}
    
    if (otus[j] %in% names(tab2)){
      sbys[2, j] <- tab2[which(names(tab2) == otus[j])]
    } else {sbys[2, j] <- 0}
    
    if (otus[j] %in% names(tab3)){
      sbys[3, j] <- tab3[which(names(tab3) == otus[j])]
    } else {sbys[3, j] <- 0}
    
    if (otus[j] %in% names(tab4)){
      sbys[4, j] <- tab4[which(names(tab4) == otus[j])]
    } else {sbys[4, j] <- 0}
    
  }
  
  sbys2 <- rrarefy(sbys, 10000)
  
  rich <- specnumber(sbys2)
  
  ratio.tab[1, i] <- round(rich[1]/rich[1], 2)
  ratio.tab[2, i] <- round(rich[2]/rich[1], 2)
  ratio.tab[3, i] <- round(rich[3]/rich[1], 2)
  ratio.tab[4, i] <- round(rich[4]/rich[1], 2)
  ratio.tab[5, i] <- props[i] / size.comm
  
  dis <- as.matrix(vegdist(sbys2, method = "bray"))
  bc.tab[1, i] <- dis[2, 1]
  bc.tab[2, i] <- dis[3, 1]
  bc.tab[3, i] <- dis[4, 1]
  bc.tab[4, i] <- props[i] / size.comm
  
  print(paste(i, " of ", dim(ratio.tab)[2], sep = ""))
}

# Rename Output
ratio.tab.even <- ratio.tab
bc.tab.even <- bc.tab


################################################################################
## Different Richness (S)
################################################################################

gamma1 <- rlnorm(n=10000, meanlog = 1, sdlog = 0.98)
gamma1 <- gamma1[rev(order(gamma1))]
gamma2 <- rlnorm(n=40000, meanlog = 1, sdlog = 0.98)
gamma2 <- gamma2[rev(order(gamma2))]
gamma3 <- rlnorm(n=7500, meanlog = 1, sdlog = 0.98)
gamma3 <- gamma3[rev(order(gamma3))]

# Define OTUs
otus1 <- paste("OTU", sprintf("%05d", seq(1:length(gamma1))), sep = "")
otus2 <- paste("OTU", sprintf("%05d", seq(1:length(gamma2))), sep = "")
otus3 <- paste("OTU", sprintf("%05d", seq(1:length(gamma3))), sep = "")

# Initiate Communities
site1 <- sample(otus1, size = 100000, replace = T, prob = gamma1)
site2 <- sample(otus2, size = 100000, replace = T, prob = gamma2)
site3 <- sample(otus3, size = 100000, replace = T, prob = gamma3)

site1.t <- table(site1); site2.t <- table(site2); site3.t <- table(site3)

# Plot Species Abundance Distrubtions
layout(matrix(1:3, ncol = 3))
par(mar = c(4, 2, 4, 1), oma = c(0, 3, 1, 1))
plot(log10(site1.t[rev(order(site1.t))]),
     ylim = c(0, 4), xlim = c(0, 20000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 2, "Log Abundance", line = 3)
mtext(side = 3, "Intact\nCommunity", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")
plot(log10(site2.t[rev(order(site2.t))]), 
     ylim = c(0, 4), xlim = c(0, 20000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 1, "Species Index (OTU)", line = 1.5)
mtext(side = 3, "Species Abundance Distributions", line = 1.5)
mtext(side = 3, "Larger S\nRelic", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")
plot(log10(site3.t[rev(order(site3.t))]), 
     ylim = c(0, 4), xlim = c(0, 20000),
     main = "", axes = F,
     xlab = "", ylab = "",
     las = 1, type = "p", 
     bg=rgb(0, 0, 1, 0.5), col = "gray40", pch=21, lwd = 0.05, cex = 1)
axis(side = 1, labels = F, tick = F)
axis(side = 2, labels = T, las = 1)
mtext(side = 3, "Smaller S\nRelic", line = -2, cex = 0.8)
box(lwd = 1.5, bty = "l")

# Basic Statistics
specnumber(site1.t)
specnumber(site2.t)
specnumber(site3.t)

# Sampling Simulation
size.comm <- 1000000
props <- seq(size.comm/100, size.comm, 50000)

ratio.tab <- matrix(NA, nrow = 5, ncol = length(props))
rownames(ratio.tab) <- c("Null", "Larger", "Smaller", "Same", "Proportion")

bc.tab <- matrix(NA, nrow = 4, ncol = length(props))
rownames(bc.tab) <- c("1-2", "1-3", "1-4", "Proportion")

for(i in 1:length(props)){
  n <- props[i]
  # Initialize
  temp.null <- c(sample(otus1, size = size.comm, replace = T, prob = gamma1))
  temp.uneven <- c(sample(otus1, size = size.comm - n, replace = T, prob = gamma1),
                   sample(otus2, size = n, replace = T, prob = gamma2))
  temp.even <- c(sample(otus1, size = size.comm - n, replace = T, prob = gamma1),
                 sample(otus3, size = n, replace = T, prob = gamma3))
  temp.same <- c(sample(otus1, size = size.comm - n, replace = T, prob = gamma1),
                 sample(otus1, size = n, replace = T, prob = gamma1))
  tab1 <- table(temp.null)
  tab2 <- table(temp.uneven)
  tab3 <- table(temp.even)
  tab4 <-table(temp.same)
  
  sbys <- matrix(NA, ncol = length(otus), nrow = 4)
  colnames(sbys) <- otus
  rownames(sbys) <- c("site1", "site2", "site3",
                      "site4")
  
  for(j in 1:dim(sbys)[2]){
    if (otus[j] %in% names(tab1)){
      sbys[1, j] <- tab1[which(names(tab1) == otus[j])]
    } else {sbys[1, j] <- 0}
    
    if (otus[j] %in% names(tab2)){
      sbys[2, j] <- tab2[which(names(tab2) == otus[j])]
    } else {sbys[2, j] <- 0}
    
    if (otus[j] %in% names(tab3)){
      sbys[3, j] <- tab3[which(names(tab3) == otus[j])]
    } else {sbys[3, j] <- 0}
    
    if (otus[j] %in% names(tab4)){
      sbys[4, j] <- tab4[which(names(tab4) == otus[j])]
    } else {sbys[4, j] <- 0}
    
  }
  
  sbys2 <- rrarefy(sbys, 10000)
  
  rich <- specnumber(sbys2)
  
  ratio.tab[1, i] <- round(rich[1]/rich[1], 2)
  ratio.tab[2, i] <- round(rich[2]/rich[1], 2)
  ratio.tab[3, i] <- round(rich[3]/rich[1], 2)
  ratio.tab[4, i] <- round(rich[4]/rich[1], 2)
  ratio.tab[5, i] <- props[i] / size.comm
  
  dis <- as.matrix(vegdist(sbys2, method = "bray"))
  bc.tab[1, i] <- dis[2, 1]
  bc.tab[2, i] <- dis[3, 1]
  bc.tab[3, i] <- dis[4, 1]
  bc.tab[4, i] <- props[i] / size.comm
  
  print(paste(i, " of ", dim(ratio.tab)[2], sep = ""))
}

# Rename Output
ratio.tab.s <- ratio.tab
bc.tab.s <- bc.tab


################################################################################
# Test Plots
################################################################################

# Layout
layout(1)
par(mar = c(4, 4, 1, 5) + 0.5)

# Evenness Plot
ratio.tab <- ratio.tab.even
plot(ratio.tab[5, ], ratio.tab[1, ], type = 'n',
     pch = 23, col = "black", bg = "white", axes = F,
     xlab = "",
     ylab = "",
     xlim = c(0, 1), ylim = c(0.78,1.22), las = 1)

points(ratio.tab[5, ], ratio.tab[2, ], pch = 22,bg = "gray40")
points(ratio.tab[5, ], ratio.tab[3, ], pch = 21, bg = "gray60")
points(ratio.tab[5, ], ratio.tab[4, ], pch = 23,bg = "white")

axis(side = 1, labels = T, tck = -0.05, lwd = 1.5, las = 1)
axis(side = 2, labels = T, tck = -0.05, lwd = 1.5, las = 1,
     at = c(seq(0.8, 1.4, 0.2)))
axis(side = 1, labels = F, tck = 0.02, lwd = 1.5)
axis(side = 2, labels = F, tck = 0.02, lwd = 1.5,
     at = c(seq(0.8, 1.4, 0.2)))
box(lwd = 1.5)

mtext(side = 1, "Proportion of Relic DNA", line = 2.5, cex = 1.25)
mtext(side = 2, "Richness Ratio", line = 2.5, cex = 1.25)

text(x = 1.02, y = 1.15, expression(paste(italic("E")["relic"],
     " > ", italic("E")["intact"])), pos = 4, cex = 1, xpd = T)
text(x = 1.02, y = 1, expression(paste(italic("E")["relic"],
     " = ", italic("E")["intact"])), pos = 4, cex = 1, xpd = T)
text(x = 1.02, y = 0.85, expression(paste(italic("E")["relic"],
    " < ", italic("E")["intact"])), pos = 4, cex = 1, xpd = T)

# Richness Plot
ratio.tab <- ratio.tab.s
plot(ratio.tab[5, ], ratio.tab[1, ], type = 'n',
     pch = 23, col = "black", bg = "white", axes = F,
     xlab = "",
     ylab = "",
     xlim = c(0, 1), ylim = c(0.78,1.22), las = 1)

points(ratio.tab[5, ], ratio.tab[2, ], pch = 22,bg = "gray40")
points(ratio.tab[5, ], ratio.tab[3, ], pch = 21, bg = "gray60")
points(ratio.tab[5, ], ratio.tab[4, ], pch = 23,bg = "white")

axis(side = 1, labels = T, tck = -0.05, lwd = 1.5, las = 1)
axis(side = 2, labels = T, tck = -0.05, lwd = 1.5, las = 1,
     at = c(seq(0.8, 1.4, 0.2)))
axis(side = 1, labels = F, tck = 0.02, lwd = 1.5)
axis(side = 2, labels = F, tck = 0.02, lwd = 1.5,
     at = c(seq(0.8, 1.4, 0.2)))
box(lwd = 1.5)

mtext(side = 1, "Proportion of Relic DNA", line = 2.5, cex = 1.25)
mtext(side = 2, "Richness Ratio", line = 2.5, cex = 1.25)

text(x = 1.02, y = 1.17, expression(paste(italic("S")["relic"],
     " > ", italic("S")["intact"])), pos = 4, cex = 1, xpd = T)
text(x = 1.02, y = 1, expression(paste(italic("S")["relic"],
     " = ", italic("S")["intact"])), pos = 4, cex = 1, xpd = T)
text(x = 1.02, y = 0.87, expression(paste(italic("S")["relic"],
     " < ", italic("S")["intact"])), pos = 4, cex = 1, xpd = T)

## Save Outputs
save(ratio.tab.even, bc.tab.even, ratio.tab.s, bc.tab.s, 
     file = "../data/SampleModel.RData")
#load("../data/SampleModel.RData")
