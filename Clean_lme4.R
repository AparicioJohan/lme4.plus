


# Function for clean_datasets using lme4

library(lme4)
library(tidyverse)
library(data.table)
source("https://raw.githubusercontent.com/AparicioJohan/lme4.plus/master/functions.R")

Clean_lme4 <- function(Response, Geno , Num_desv=3, Show_results=TRUE, data=NULL , name='Exp', rep="rep", block="block"){
  
  Datos <- droplevels(data)
  
  w <- 1                                         # Starter 
  k <- Num_desv                                  # Number of standard deviations to consider extreme outliers 

  remFinal <- data.frame(Response=as.numeric(),  Genotype = as.character(), rep=as.numeric(),block=as.numeric() )

  m_alpha <- paste(Response, "~" , "1 + (1|" , Genotype, ") + " , rep , "+ (1|", rep,":", block, ")")
  m_rcbd <- paste(Response, "~" , "1 + (1|" , Genotype, ") + " , rep )
  Model <-  c(m_alpha, m_rcbd)
  cond=  length(unique(Datos$block))
  mvOn = F
  ix = c()
  
  reps <- sum(is.na(Datos[,"Response"] ))
  
  if(cond == nrow(Datos)) return()
  if(reps >= 0.6*nrow(Datos)) return()
  if(length(unique(Datos[,rep]))==1) return()
  if( sum(is.na(Datos[,Genotype]))>0.5*length(Datos[,Genotype]) )  return()
  if( sum(is.na(Datos[,Response]))>0.5*length(Datos[,Response])) return()
  
  if (cond>1) {
    form=as.formula(Model[1])
  } else {form=as.formula(Model[2])}
  
  cat(name,"\n")
  message("\nRemoving outliers from ", name, "\n")
  
  while (mvOn == F) {
    
    Mo <- lmer(form, data = Datos, 
               control = lmerControl(optimizer ="Nelder_Mead"),
               na.action = na.omit, REML = T) 
    Mo %>% 
      residuals(., scaled=TRUE) -> res
    
    if ( sum( abs(res) > 3 ) ){
      
      ix = c(ix, names( which( abs(res)>3 ) ))  # which.max( abs(res) )
      message("There are ", sum(abs(res) > 3), " NA values remaining")
      
      remTMP <- data.frame(Response=Datos[ix,Response],  Genotype = Datos[ix,Geno], col=Datos[ix,rep],row=Datos[ix,block] )
      remFinal <- rbind(remTMP,remFinal)
      
      Datos[ix,"YDHA"] = NA
      
    } else {
      
      mvOn = T
    }
  }
  
  d <- broom.mixed::augment(ranef(Mo))
  d <- d[d$grp==Genotype,c("level","estimate","std.error")]
  d <- data.frame(level=d[,1],estimate=round(d[,2],2),std.error=round(d[,3],2))
  BLUPs<-d
  
  # names(BLUPs)[ncol(BLUPs)] <- "Trial"
  BLUPs$file <- paste0(name)
  
  
  Modelo.vars   <- as.data.table(VarCorr(Mo))
  VarG <- as.numeric(Modelo.vars[grp==Genotype,'vcov'])
  VarE <- as.numeric(Modelo.vars[grp=='Residual','vcov'])
  replicate <- length(unique(Datos[,rep]))
  r2 <- as.numeric(r.lme4(Mo))
  
  Sum <- data.frame(file= paste0(name) , VarE=VarE, VarG=VarG, rep=replicate , r2=r2)
  Sum <- Sum %>% mutate(H=round(VarG/(VarG+VarE/replicate),2))
  
  
  k=list(BLUPs=BLUPs, data_clean=Datos, Model=Mo, Remove=remFinal, Summ=Sum)
  k

  
}  
  
  
#################################
##          Example 
#################################


# temp <- Clean_lme4(Response = "Response",name = alp,
#            Geno = "Genotype",
#            Num_desv = 3 , Show_results = T,
#            data = tmpdata,
#            rep = "rep",
#            block = "block"  )
