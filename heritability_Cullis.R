# Function Heritability cullis 

Heri.cullis <- function(Model,Nom.gen){
  require("data.table")
  
  Gen_levels=levels(Model@frame[,Nom.gen])
  
  # to obtain var-cov-matrix for BLUPs of gen effect.
  vc   <- as.data.table(VarCorr(Model)) # extract estimated variance components (vc)
  # Number of random effects
  n.ran <- nrow(vc)
  
  # R = varcov-matrix for error term
  n    <- length(summary(Model)$residuals) # numer of observations
  vc.e <- vc[grp=="Residual", vcov]        # error vc
  R    <- diag(n)*vc.e                     # R matrix = I * vc.e
  
  # names of random effects
  nomb <- names(summary(Model)$ngrps)
  
  # Genotype
  n.g <- summary(Model)$ngrps[which(nomb==Nom.gen)]
  vc.g <- vc[grp==Nom.gen, vcov]   
  
  # G matrix of random effects
  G.tmp <- list()
  if (n.ran>=2) {
    for (i in 1:(n.ran-1)) {   # remove the residual variance
      n.tmp <- summary(Model)$ngrps[which(nomb==nomb[i])]
      vc.tmp <- vc[grp==nomb[i], vcov]   
      G.tmp[[i]] <- diag(n.tmp)*vc.tmp   
    }
  }
  
  G <- bdiag(G.tmp) # G is blockdiagonal with G.g and G.b
  
  # Design Matrices
  X <- as.matrix(getME(Model, "X")) # Design matrix fixed effects
  Z <- as.matrix(getME(Model, "Z")) # Design matrix random effects
  
  # Mixed Model Equation 
  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 
  
  C <- as.matrix(rbind(cbind(C11, C12),  # Combine components into one matrix C
                       cbind(C21, C22)))
  
  # Mixed Model Equation Solutions 
  C.inv <- solve(C)                                # Inverse of C
  C22.g <- C.inv[Gen_levels, Gen_levels] # subset of C.inv that refers to genotypic BLUPs
  
  
  # Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
  vdBLUP.sum<- n.g*sum(diag(C22.g))-sum(C22.g)
  vdBLUP.avg <- vdBLUP.sum * (2/(n.g*(n.g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs
  
  #############
  # H2 Cullis #
  #############
  H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
  H2Cullis
}
