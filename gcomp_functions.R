##########################################################################
## G-Computation Functions
## Creators: Daisy Shepherd, Margarita Moreno-Betancur, Anthony Charsley
## Affiliations: The Murdoch Children's Research Institute & The University of Melbourne

## Code originally developed for ViCBiostat summer school (Feb 2020)
## Last updated: June 2024
##########################################################################

## Function: gcompstat(data,ind,out,exp,formula)
## Arguments:
##      data = datset
##      ind = indices to use of dataset
##      out = name of outcome variable
##      exp = name of exposure variable (factor, takes values "Yes" and "No")
##      formula = to pass to fitted model
##
## Notes: Function works for continuous or binary outcome
##        Continuous outcome - ACE is the difference in expected outcomes
##        Binary outcome - ACE is the causal risk ratio

gcompstat <- function(data,ind,out,exp,formula){
  dat <- data[ind,] 
  
  # Original dataset
  dat$index <- -1
  
  # Dataset untreated: Treatment set to "No"
  dat0 <- dat
  dat0$index <- 0
  dat0[[exp]] <- "No"
  dat0[[out]] <- NA
  
  # Dataset treated: Treatment set to "Yes"
  dat1 <- dat
  dat1$index <- 1
  dat1[[exp]] <- "Yes"
  dat1[[out]] <- NA
  
  # Combining into one sample
  onesample <- rbind(dat,dat0,dat1)
  
  ## Continuous Outcome
  if(is.factor(dat[[out]])==FALSE){
    if(is.numeric(dat[[out]])==TRUE){
      
      # Fitting the linear model
      fit <- lm(as.formula(formula), data=onesample)
      
      ## Predicted values under fitted model
      onesample$predicted_meanY <- predict(fit, onesample)
      
      ## Calculating average causal effect (ACE)
      EYA0 <- mean(onesample[which(onesample$index==0),]$predicted_meanY) 
      EYA1 <- mean(onesample[which(onesample$index==1),]$predicted_meanY) 
      ACE_diff <- EYA1 - EYA0
      return(ACE_diff)
    }
    else{
      print("Error: Outcome must be continuous or binary")
    }  
  }
  
  ## Binary outcome
  if(is.factor(dat[[out]])==TRUE){
    
    # Fitting the glm
    fit <- glm(as.formula(formula), data=onesample, family=binomial(link=logit))
    
    ## Predicted values under fitted model
    onesample$predicted_meanY <- predict(fit, onesample, type="response")
    
    ## Calculating average causal effect (ACE)
    EYA0 <- mean(onesample[which(onesample$index==0),]$predicted_meanY) 
    EYA1 <- mean(onesample[which(onesample$index==1),]$predicted_meanY) 
    ACE_risk_ratio <- EYA1/EYA0
    return(ACE_risk_ratio)
  }
}



##########################################################################

## Function: gcompboot(runboot_dat,R,stat,out,exp,formula,seed)
## Bootstrap function to generate CIs, SEs, P-vals
## Arguments:
##      runboot_dat = datset to perform bootstrapping on
##      R = number of bootstrap replicates
##      stat = function implementing estimator
##      out = name of outcome variable
##      exp = name of exposure variable
##      formula = to pass to fitted model
##      seed = numeric value to set seed for reproducibility

gcompboot <- function(runboot_dat,R=1000,stat,out,exp,formula,seed){
  set.seed(seed)
  
  ## Performing bootstrap
  bstrap <- boot(data=runboot_dat, statistic=stat, stype="i", R=R,
                 out=out, exp=exp, formula=formula)
  
  ## Obtaining summary statistics
  bt <- boot.ci(bstrap, index=1, type=c("perc"))  
  est <- bstrap$t0
  se <- apply(bstrap$t, 2, sd, na.rm=T)
  ci_low <- bt$percent[4]
  ci_upp <- bt$percent[5]
  pval <- 2*(1-pnorm(q=abs(est/se)))
  return(c(est=est,se=se,ci.low=ci_low,ci.upp=ci_upp,pval=pval))
}


##############################################################################



