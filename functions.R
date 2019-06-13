library(broom.mixed)
library(lme4)


# Pseudo R square LME4
r.lme4 <- function(Model){
  
  Info <-   data.frame( broom.mixed::augment(Model) )
  response      <- Info[,2]
  mean.response <- mean(response,na.rm = T)
  fitted        <- Info[,".fitted"]
  SS.fitted     <- sum( (response-fitted)^2 , na.rm = T)
  SS.response   <- sum( (response-mean.response)^2 , na.rm = T )
  
  R <- 1- SS.fitted/SS.response
  
  names(R) <- "r.square"
  return(round(R,3))
}

# Number of estimates
npars.lme4 <-  function(Model){
  p <- length(fixef(Model))
  q <- nrow(data.frame(VarCorr(Model)))
  return(p+q)
}
