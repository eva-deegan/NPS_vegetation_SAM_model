model{
  
  for(i in 1:n.data){
    for(s in 1:n.species){
      
      # Likelihood --------------------------------------------------------------
      
      #observed count data are assumed to follow a binomial distribution
      #with presence probability, p, and number of trials, N:
      y[i,s] ~ dbin(p[i,s], N[i])
      # Replicated data:
      y.rep[i,s] ~ dbin(p[i,s], N[i])
      # Compute squared difference between observed and replicate data:
      sqdiff[i,s] <- pow(y.rep[i,s]-y[i,s], 2)
      #logistic regression for the presence probability:
      logit(p[i,s]) <- b0.species[s] + #species-specific intercept
        b0.plot[Plot.ID[i],s] + #plot random effects
        b[1,s]*AntWintP[i] + #antecedent winter precip
        b[2,s]*AntSummP[i] +#antecedent summer precip
        b[3,s]*AntWintT[i] +#antecedent winter temp
        b[4,s]*AntSummT[i] +#antecedent summer temp
        b[5,s]*AntWintP[i]*AntSummP[i] +
        b[6,s]*AntWintP[i]*AntWintT[i] +
        b[7,s]*AntWintP[i]*AntSummT[i] +
        b[8,s]*AntSummP[i]*AntWintT[i] +
        b[9,s]*AntSummP[i]*AntSummT[i] +
        b[10,s]*AntWintT[i]*AntSummT[i] 
      
      
    }
    # SAM weighting -----------------------------------------------------------
    
    #antecedent effects are the sum of the weighted effects for all lags
    ## computing the non-standardized antecedent climate variables:
    AntWintTx[i] <- sum(WintTemp[i,])
    AntSummTx[i] <- sum(SummTemp[i,])
    AntWintPx[i] <- sum(WintPrecip[i,])
    AntSummPx[i] <- sum(SummPrecip[i,])

    ## computing the standardized climate variables to use in the above logistic regression:
    AntWintT[i] <- (AntWintTx[i]-mean.AntWintT)/sd.AntWintT
    AntSummT[i] <- (AntSummTx[i]-mean.AntSummT)/sd.AntSummT
    AntWintP[i] <- (AntWintPx[i]-mean.AntWintP)/sd.AntWintP
    AntSummP[i] <- (AntSummPx[i]-mean.AntSummP)/sd.AntSummP
    
    #Multiply each lag climate variable by the re-normalized
    #yearly weights for each of those variables
    for(l in 1:n.lag){
      #multiply climate variables by re-normalized weights
      #these are based on sampling our wW.bar and wS.bar
      #matrices based on the survey month for each observation
      #which is a number 1:12 (oct = 1, nov = 2, ...., sept = 12)
      WintTemp[i,l] <- WintT[i,l]*wWT.bar[survey.month[i],l]
      SummTemp[i,l] <- SummT[i,l]*wST.bar[survey.month[i],l]
      WintPrecip[i,l] <- WintP[i,l]*wWP.bar[survey.month[i],l]
      SummPrecip[i,l] <- SummP[i,l]*wSP.bar[survey.month[i],l]
      
    }
    
  } #likelihood n.plots

  for(s in 1:n.species){
    # sum of squared differences for each species:
    Dsum1[s] <- sum(sqdiff[,s])
    }

  #Total sum of squared differences
  Dsum <- sum(Dsum1[])

  ## computing mean and st dev of non-standardized SAM-climate variables as 
  ## we will likely want to monitor these and compute their posterior stats:
  mean.AntWintT <- mean(AntWintTx[])
  sd.AntWintT <- sd(AntWintTx[])
  mean.AntSummT <- mean(AntSummTx[])
  sd.AntSummT <- sd(AntSummTx[])
  mean.AntWintP <- mean(AntWintPx[])
  sd.AntWintP <- sd(AntWintPx[])
  mean.AntSummP <- mean(AntSummPx[])
  sd.AntSummP <- sd(AntSummPx[])
  
  # Renormalize SAM weights -------------------------------------------------
  
  #due to variable sampling months, need to correct OG SAM 
  #weights for the first month in any season to 
  #account for how many months we actually sampled
  #we're penalizing that weight for any locations where
  #we don't have a full survey of the months in that
  #season and then renormalizing the SAM weights to equal 1
  
  # Summer renormalized lag matrix ------------------------------------------
  
  #first, set the first six months for lag 1 (current year) equal to zero, 
  #since in these cases, no first summer weight should be incorporated
  #based on sampling dates that occur in the winter (r=1,..,6):
  for(r in 1:6){
    wST.bar[r,1] <- 0
    sumS2T[r] = sum(wST[2:n.lag]) + wS1T.hat[r]
    wSP.bar[r,1] <- 0
    sumS2P[r] = sum(wSP[2:n.lag]) + wS1P.hat[r]
  }
  #for the first lags for each month following the first 6, 
  #it will be the reweighted lag weight wS1.hat that we 
  #are then normalizing again
  for(r in 7:12){ #the remainder of the months
    # concurrent summer:
    wST.bar[r,1] <- wS1T.hat[r-6]/(sumS2T[r-6])
    wSP.bar[r,1] <- wS1P.hat[r-6]/(sumS2P[r-6])
  }
  
  #for the remainder of the first half of the summer lags,
  #the renormalized weight is just the weight divided by the 
  #sum of the three non-zero weights
  for(l in 2:n.lag){
    for(r in 1:6){ # sampled in winter, so only prior summer periods:
      wST.bar[r, l] <- wST[l]/sum(wST[2:n.lag])
      wSP.bar[r, l] <- wSP[l]/sum(wSP[2:n.lag])
    }
    #for the rest of the summer lags, the renormalized weights 
    #incorporate the wS1.hat values, which are un-normalized
    # sampled in the summer
    for(r in 7:12){ #the remainder of the months (sampled in summer)
      ## divided by sumS2
      wST.bar[r,l] <- wST[l]/sumS2T[r-6]
      wSP.bar[r,l] <- wSP[l]/sumS2P[r-6]
    }
  }
  
  # Winter renormalized lag matrix ------------------------------------------
  
  #renormalized weights
  ## need separate "special" calculations for lag = 1 (current year)
  for(r in 1:6){ #6 months for this weight
    sumW2T[r] = sum(wWT[2:n.lag]) + wW1T.hat[r]
    sumW2P[r] = sum(wWP[2:n.lag]) + wW1P.hat[r]
    # First (concurrent year) weight:
    wWT.bar[r,1] <- wW1T.hat[r]/(sumW2T[r])
    wWP.bar[r,1] <- wW1P.hat[r]/(sumW2P[r])
    # Below (r+6) corresponds to case when sampling occured in the summer
    # (so, full set of winter weights).
    wWT.bar[r+6,1] <- wWT[1]
    wWP.bar[r+6,1] <- wWP[1]
    for(l in 2:n.lag){
      #winter re-normalized weights
      wWT.bar[r,l] <- wWT[l]/(sumW2T[r])
      wWP.bar[r,l] <- wWP[l]/(sumW2P[r])
      #the rest of the winter weights are just their
      #original (when sampling occured in summer)- 
      # but kept this matrix of 12 x 4
      #dimensions so we wouldn't have to read in two
      #indexes for summer and winter seasons
      wWT.bar[r+6,l] <- wWT[l]
      wWP.bar[r+6,l] <- wWP[l]
    }
  }
  
  # First lag re-weighted lag weights ---------------------------------------
  
  #these are re-weighted based on the number of months
  #that were sampled out of the season
  #(gives less weight to lags when fewer months were 
  #sampled)
  for(r in 1:n.months){ #this will be a number 1-6 for the months in that season
    #r is just a vector 1:6
    #read in as "data"
    wW1T.hat[r] <- wWT[1]*(r/6)
    wS1T.hat[r] <- wST[1]*(r/6)
    wW1P.hat[r] <- wWP[1]*(r/6)
    wS1P.hat[r] <- wSP[1]*(r/6)
  }
  
  # Original SAM weighting --------------------------------------------------
  wWT[1:4]  ~ ddirch(alpha[1:4])
  wWP[1:4]~ ddirch(alpha[1:4])
  wST[1:4]  ~ ddirch(alpha[1:4])
  wSP[1:4] ~ ddirch(alpha[1:4])

  
  # Beta priors -------------------------------------------------------------
  
  
  for(c in 1:n.covs){
    # priors for the community-level effect of covariate c
    mu.b[c] ~ dnorm(0, 0.01)  # Mean effect across all species
    tau.b[c] <- 1/pow(sig.b[c], 2)
    ## Since there are likely relatively few species (e.g., only a few C3 or C4),
    ## using folded-t priors for the standard deviations.
    ## need to provide inits for sig.t.b (can use former sig.b inits initially).
    sig.t.b[c] ~ dt(0,0.1,2)   # Jags parameterization: dt(mu,tau,k), k = degrees of freedom.
    sig.b[c] <- abs(sig.t.b[c])
    
    # Priors for species-level effects
    for(s in 1:n.species){
      b[c,s] ~ dnorm(mu.b[c], tau.b[c]) # Each species' effect is drawn from the common distribution
    }
  }

  ## Prior for species-level intercepts, similar to the prior for b[c,s]; assuming you have at least 3
  ## species. If fewer than 3 species, then a hierarchical prior doesn't make sense and would likely cause
  ## convergence/mixing issues:
  for(s in 1:n.species){
      b0.species[s] ~ dnorm(mu.species, tau.species)
    ## back-transform species-level intercept to get the species-level baseline presence probability at
    ## average climate:
    p0.species[s] <- ilogit(b0.species[s])
    }

  ## Priors for the overall intercept and standard deviation describing variation among species-
  ## level intercepts:
  mu.species ~ dnorm(0,0.01)
  ## back-transform overall intercept to get the overall baseline presence probability
  p.overall <- ilogit(mu.species)

  tau.species <- 1/pow(sig.species,2)
  # Folded t prior for standard deviation (need inits for sig.t.species):
  sig.t.species ~ dt(0,0.1,2)
  sig.species <- abs(sig.t.species)  
  
  # Hierarchical prior for plot-level random effects.---------------------------------------------
  
  # Sum-to-zero constraint for plot random effects since there are many (>>10)
  # plots, and they are not nested in a higher level (e.g., site).

  for(s in 1:n.species){
    for(p in 1:(n.plots-1)){
      b0.plot[p, s] ~ dnorm(0, tau.plot)
    } 
    b0.plot[n.plots, s] <- -sum(b0.plot[1:(n.plots-1), s])
  }
  
  # Priors for precisions and standard deviations for random effects.

  tau.plot <- 1/pow(sig.plot,2)
  ## To be consistent with the folded-t prior for sig.b, and to account for some 
  ## ecosites having relatively few plots, also use a folded-t prior for sig.plot:
  ## Will need to provide inits for sig.t.plot (can use former sig.plot inits initially).
  sig.t.plot ~ dt(0,0.01,2)   # Jags parameterization: dt(mu,tau,k), k = degrees of freedom.
  sig.plot <- abs(sig.t.plot) 
}
