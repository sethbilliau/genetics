### See how long it takes a neutral mutation to go to fixity

# Step 0
library(shinythemes)
library(shiny)
library(ggplot2)
library(plyr)
library(shinyLP)
library(plotly)
source("NEUTRAL.R")


# Vectorize fixity simulator -> uncurry the function
fixitysimulatorvec <- function(x) {
  if ((fAA + fAa > 1) || fAA > 1 || fAA < 0 || fAa >1 || fAa < 0) {
    stop('Invalid Frequencies for genotypes.')
  }
  # Recursive helper to check if gex$Neration reaches fixity
  gex$Nerations <- function(num, gen0) {
    gen1 <- recombix$Nepop2(gen0, x$Ne)
    if (checkfixity(gen1) == 2) {
      return(num)
    } else if (num == x$gens || checkfixity(gen1) == 0) {
      return("Failure")
    } else {
      gex$Nerations(num+1, gen1)
    }
  }
  
  # Gex$Nerate a dataframe
  df <- data.frame()
  
  # Initialize a population
  gen0 <- allelefpop(x$Ne,x$fAA,x$fAa)
  
  if (checkfixity(gen0) == 2) {
    return(data.frame(genid = seq(from=1, to=x$NTrials, 1),
                      ReachFix=rep(1,x$NTrials), 
                      GensToFixity=rep(0,x$NTrials)))
  } else {
    for (i in 1:x$NTrials) {
      out <- gex$Nerations(1, gen0)
      if (out == "Failure") {
        df <- rbind(df,c(i,F,NA))
      } else {
        df <- rbind(df,c(i,T,out))
      }
    }
    names(df) <- c("genid","ReachFix","GensToFixity")
    return(df)
  }
}


# Sim Conditions
NTrials = 1000
Ne = seq(from=10, to=50, by =10)
Ngens = 10^6
fAA = seq(from=.05, to=.85, by=0.1)
fAa = rep(.1,8)

# simulator 
simulator <- function(fAA,Ne) {
  for (j in fAA) {
    p = allelef(j,.1)
    for (i in Ne) {
      start.time <- Sys.time()
      x <- fixitysimulator(1000,i,10^6, j, .1)
      end.time <- Sys.time()
      time.taken <- difftime(end.time,start.time, units="secs")
      
      result = data.frame(NTrials=NTrials, 
                          Ne=i, 
                          Ngens = Ngens, 
                          p=p, 
                          propfix = mean(x$ReachFix), 
                          gens = mean(x$GensToFixity, na.rm = T),
                          time = time.taken)
      final = rbind(final,result);final
      write.csv(final, file = "~/kimuradat.csv")
    }
  }
}


analyticalgens <- function(Ne, p) {
  z <- data.frame()
  for (j in 1:length(Ne)) {
    for (i in 1:length(p)) {
    result<-((-1 / p[i]) * (
             4*Ne[j] * (1-p[i]) * log(1-p[i], base = exp(1))
            ))
    row <- cbind(Ne=Ne[j], p=p[i], Aresult = result)
    z <- rbind(z,row)
    }
  }
  return(z)
}

p = seq(.1,.9, by=.1)
z<-analyticalgens(Ne,p); z

dat <- read.csv("kimuradat.csv")
dattrim <- dat[-c(seq(6,10,1)),]
datsmallp <- dat[-c(seq(11,55,1)),]

lin.model <- lm(gens ~ Ne, datsmallp)
lin.model$coefficients[1]
# Gens to Fixity vs Ne
allele10pop <- ggplot(data=datsmallp, aes(x=Ne, y=gens))+  
  geom_abline(intercept = 0, slope = 4, color="red") +
  geom_point() + 
  stat_smooth(method = lm) +
  scale_x_continuous(name="Effective Population", limits=c(0, 100), breaks = seq(0,100,10)) +
  scale_y_continuous(name="Generations to Fixity", limits=c(0, 400)) +
  labs(title="Generations to Fixity vs. Ne",
       subtitle="Number of Trials: 1000")

# Time vs Ne
allele10time <- ggplot(data=datsmallp, aes(x=Ne, y=time)) +
  geom_point() +
  scale_x_continuous(name="Effective Population", limits=c(0, 100), breaks = seq(0,100,10)) +
  stat_smooth(fullrange=TRUE, col = "red") +
  scale_y_continuous(name="Simulation Time") +
  labs(title="Time vs. Ne",
       subtitle=paste("Number of Trials: 1000; Total time: 5.08 hours"))

# Allelefreqfns for plotting
allelefreq10 <- function(p) {
  (-1 / p) * (40 * (1-p) * log(1-p, base = exp(1)))
}
allelefreq20 <- function(p) {
  (-1 / p) * (80 * (1-p) * log(1-p, base = exp(1)))
}
allelefreq30 <- function(p) {
  (-1 / p) * (120 * (1-p) * log(1-p, base = exp(1)))
}
allelefreq40 <- function(p) {
  (-1 / p) * (160 * (1-p) * log(1-p, base = exp(1)))
}
allelefreq50 <- function(p) {
  (-1 / p) * (200 * (1-p) * log(1-p, base = exp(1)))
}


aggplot <- ggplot(data=dattrim, aes(x=p, y=gens)) +
  geom_point(aes(colour = factor(Ne))) + 
  scale_color_manual("Ne", values=c( "#C44D58","#556270","#4ECDC4", "#C7F464", "#FF6B6B" )) + 
  scale_x_continuous(name="Initial Allele Frequency",
                     limits=c(0, 1),
                     breaks = seq(0,1,.1)) +
  scale_y_continuous(name="Generations to Fixity", limits=c(0, 200)) +
  stat_function(fun = allelefreq10, colour = "#C44D58") + 
  stat_function(fun = allelefreq20, colour = "#556270") + 
  stat_function(fun = allelefreq30, colour = "#4ECDC4") +
  stat_function(fun = allelefreq40, colour = "#C7F464") + 
  stat_function(fun = allelefreq50, colour = "#FF6B6B") +
  labs(title="Generations to Fixity vs. Inital Allele Frequency",
       subtitle="Number of Trials per dot: 1000") +
  theme(legend.position="bottom")
  