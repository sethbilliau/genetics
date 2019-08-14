# Assume that there are 4 genes with a dominant and a recessive. 
library(shinythemes)
library(shiny)
library(ggplot2)
library(plyr)
library(shinyLP)

# Create Random people
randpeople <- function(N) {
  df <- data.frame()
  for (i in 1:N) {
    person <- c(id=i,A=sample(c(2,1,0),1),
      B=sample(c(2,1,0),1),
      C=sample(c(2,1,0),1),
      D=sample(c(2,1,0),1))
    df <- rbind(df, cbind(t(person),height = sum(person[2:5])))
  }
  return(df)
}

# Testing constants
one<-randpeople(1)
thousand<-randpeople(1000)

# take in checked dataframe and return plot
plotheight <- function(pop) {
  ggplot(pop, aes(height)) +
    geom_histogram(color="darkblue", fill="white", binwidth=1) +
    labs(title="Distribution of Height",
         subtitle=paste("N = ", length(pop$id)),
         x="Height", y= "Frequency") +
    geom_vline(data=pop, aes(xintercept=mean(pop$height)),
               linetype="dashed") + coord_cartesian(xlim = c(0,8))
}

# Take in a checked population and a threshhold and eliminate losers
naturalselection <- function(pop,threshhold,strength) {
  fit <- subset(pop, height > threshhold)
  unfit <- subset(pop, height <= threshhold)
  unfitsurvive <- unfit[sample(nrow(unfit), (1-strength) * length(unfit$id), replace = F),]
  newpop <- rbind(fit, unfitsurvive)
  # Rescramble
  newpop <- newpop[sample(nrow(newpop)),]
  # Re-id
  newpop[,1] <- seq(from = 1, to = length(newpop$id), by = 1 )
  return(newpop)
}

# Helper to recombine an individual pair
recombinepair <- function(m,f) {
  child <- data.frame(id=0,A=0,B=0,C=0,D=0, height=0)
  gamete <- function(a){
    if (a == 2){
      1
    } else if (a == 0) {
      0
    } else {
      sample(c(1,0),1)
    }
  }
  child$A <-sum(gamete(m$A),gamete(f$A))
  child$B <-sum(gamete(m$B),gamete(f$B))
  child$C <-sum(gamete(m$C),gamete(f$C))
  child$D <-sum(gamete(m$D),gamete(f$D))
  child$height <- sum(child$A,child$B,child$C,child$D)
  child
}

recombinepop <- function(pop, N) {
  newpop <-data.frame(id=rep(NA,N),A=rep(NA,N),
                      B=rep(NA,N),C=rep(NA,N),D=rep(NA,N),height=rep(NA,N))
  males <- pop[sample(nrow(pop), N, replace = T),]
  females <- pop[sample(nrow(pop), N, replace = T),]
  for (i in 1:N) {
    newpop[i,]<-recombinepair(males[i,],females[i,])
    newpop[i,1] <- i
  }
  return(newpop)
}

countfreq <- function(pop) {
  As <- count(pop$A)
  Bs <- count(pop$B)
  Cs <- count(pop$C)
  Ds <- count(pop$D)
  tbl <- data.frame(Freq=c(0,1,2),A=As[,2],B=Bs[,2],C=Cs[,2],D=Ds[,2])
  tbl
}

cat("GENETICS.R Loaded...\n")
