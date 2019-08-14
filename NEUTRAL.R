### See how long it takes a neutral mutation to go to fixity

# Step 0
library(shinythemes)
library(shiny)
library(ggplot2)
library(plyr)
library(shinyLP)

# Assume no selection  

# Fun of effective population size, f of AA and F of heterozygotes
allelefpop <- function(N, fAA, fAa) {
  # Initialize df and genotypes
  pop <- data.frame()
  aa <- data.frame(id=0, A=0); Aa <- data.frame(id=0, A=1); AA <- data.frame(id=0, A=2)
  
  # Add appropriate number of genotypes
  if (round(1-fAA-fAa) > 0) {
    pop <- rbind(pop,aa[rep(seq_len(nrow(aa)), each=floor(N*(1-fAA-fAa))),])
  }
  if (fAA > 0) {
    pop <- rbind(pop,AA[rep(seq_len(nrow(AA)), each=floor(N*fAA)),])
  }
  if (fAa > 0) {
    pop <- rbind(pop,Aa[rep(seq_len(nrow(Aa)), each=floor(N*fAa)),])
  }
  
  # Randomly add a genotype or two if there's a funky remainder
  while (length(pop$id) < N) {
    pop <- rbind(pop,data.frame(id=0,A=sample(c(0,1,2),1)))
  }
  
  # Scramble, Correct id #'s and rownames and return population
  pop$A<-sample(pop$A); pop$id <- seq(1,N,by=1); rownames(pop) <- seq(1,N,by=1)
  return(pop)
}

# Recombine pair of individuals
recombinepair2 <- function(m,f) {
  child <- data.frame(id=0,A=0)
  # Check gamete
  gamete <- function(a){
    if (a == 2){
      1
    } else if (a == 0) {
      0
    } else {
      sample(c(1,0),1)
    }
  }
  # Create child from parent gametes
  child$A <-sum(gamete(m$A),gamete(f$A))
  return(child)
}

# Recombine populations calling the helper recombinepair
recombinepop2 <- function(pop, N) {
  # Init new pop
  newpop <- data.frame(id=rep(NA,N),A=rep(NA,N))
  
  # Randomly assign males and females for mating. Sample with replacement. 
  # NOTE: This makes it possible, though unlikely, for a species to self
  males <- pop[sample(nrow(pop), N, replace = T),]
  females <- pop[sample(nrow(pop), N, replace = T),]
  
  # create N people from random permutations of parents. 
  for (i in 1:N) {
    newpop[i,]<-recombinepair2(males[i,],females[i,])
    newpop$id[i] <- i
  }
  return(newpop)
}

# Return True if fixed, false else
checkfixity <- function(pop) {
  if (length(pop$A) == length(which(pop$A == 2))) {
    2
  } else if (length(pop$A) == length(which(pop$A == 0))) {
    0
  } else {1}
}

# Trials, Ne, Ngens, fAA, fAa
fixitysimulator <- function(NTrials, Ne, Ngens, fAA, fAa) {
  if ((fAA + fAa > 1) || fAA > 1 || fAA < 0 || fAa >1 || fAa < 0) {
    stop('Invalid Frequencies for genotypes.')
  }
  # Recursive helper to check if generation reaches fixity
  generations <- function(num, gen0) {
    gen1 <- recombinepop2(gen0, Ne)
    if (checkfixity(gen1) == 2) {
      return(num)
    } else if (num == Ngens || checkfixity(gen1) == 0) {
      return("Failure")
    } else {
      generations(num+1, gen1)
    }
  }
  
  # Generate a dataframe
  df <- data.frame()
  
  # Initialize a population
  gen0 <- allelefpop(Ne,fAA,fAa)
  
  if (checkfixity(gen0) == 2) {
    return(data.frame(genid = seq(from=1, to=NTrials, 1),
                      ReachFix=rep(1,NTrials), 
                      GensToFixity=rep(0,NTrials)))
  } else {
    for (i in 1:NTrials) {
      out <- generations(1, gen0)
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
cat("NEURTRAL.R Loaded...\n")
# Usage: fixitysimulator(Trials, Ne, Ngens, fAA, fAa)