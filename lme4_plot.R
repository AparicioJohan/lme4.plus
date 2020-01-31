##########################
#  Spatial plot for lme4 
##########################

library(data.table)
library(dplyr)

# first fill the col/row in your dataset
"coords" <- function(data,col,row) {
  require(data.table)
  deparse(substitute(data))
  x.coord <- data[,col]
  y.coord <- data[,row]
  
  columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
  rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
  
  xy.coord <- data.table(expand.grid(col = columns, row = rows))
  h <- merge(xy.coord, data, by.x = c( "row","col"), by.y = c( "row","col"), 
             all.x = TRUE, all.y = TRUE)
  data <- arrange(h,col,row)
  
  return(data)
}

"dup" <- function(data,col="col",row="row"){
  dup <- data[,c(col,row)][duplicated(data[,c(col,row)]), ] 
  data  <-data[ which(!duplicated(data[,c(col,row)])) , ]
  return(data)
}


# Spatial plot 
lme4.plot <-
  function(x,Datos,col="col",row="row",gen="line",all.in.one = TRUE, main = NULL, annotated = FALSE, depict.missing = FALSE, ...) {
    xlab <- col
    ylab <- row
    x.coord <- Datos[,col]
    y.coord <- Datos[,row]
    response <- Datos[, names(x@frame)[1]]  # !is.na(Datos[,names(x@frame)[1]]) 
    Datos$rowname <- rownames(Datos)
    
    residuals <- data.frame(rowname=names(residuals(x)), residuals= as.numeric(residuals(x)))
    residuals <- merge(residuals,Datos,by="rowname", all = T)
    residuals <- arrange(residuals,col,row)$residuals
    
    fitted <- data.frame(rowname=names(fitted.values(x)), fitted=as.numeric(fitted.values(x)))
    fitted <- merge(fitted,Datos,by="rowname", all = T)
    fitted <- arrange(fitted,col,row)$fitted
    
    geno.pred <- data.frame(genotype=rownames(ranef(x)[[gen]])   , predicted.value=ranef(x)[[gen]][,1])
    names(geno.pred)[1] <- gen
    columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
    rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
    xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
    v <- merge(geno.pred,Datos,by = gen,all.y = TRUE)
    xx <- arrange(v,col,row)
    geno.pred <- as.vector(xx$predicted.value)
    
    environment <- fitted-geno.pred-fixef(x)["(Intercept)"]
    
    setNumericRounding(2)
    
    if(is.null(main)) main = paste("Trait: ",names(x@frame)[1] , sep = "")
    
    setkeyv(xy.coord, c("rows", "columns"))
    ONE <- rep(1, length(x.coord))    
    df <- data.table(columns = x.coord, rows = y.coord, 
                     response = response, fitted = fitted,environment,
                     residuals = residuals, geno.pred = geno.pred, ONE = ONE)
    setkeyv(df, c("rows", "columns"))
    df <- df[xy.coord]
    df <- df[order(df$columns, df$rows),]
    
    colors = topo.colors(100)
    
    main.legends <- c('Raw data', 'Fitted data', 'Residuals',"Effect Design"  ,"Genotypic BLUPs", 'Histogram')
    if(all.in.one) {
      op <- par(mfrow = c(2,3), oma = c(ifelse(annotated, 12, 2), 1, 3, 2), mar = c(2.7, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))                
    } else {
      if(!is.null(main))
        main.legends <- rep(main, length(main.legends))
    }
    
    range <- range(c(response, fitted), na.rm = TRUE)
    fields::image.plot(columns, rows, t(matrix(df$response, ncol = length(columns), nrow = length(rows))), main = main.legends[1], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
    if(!all.in.one)
      readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$fitted, ncol = length(columns), nrow = length(rows))), main = main.legends[2], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
    if(!all.in.one)
      readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$residuals, ncol = length(columns), nrow = length(rows))), main = main.legends[3], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
      readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$environment, ncol = length(columns), nrow = length(rows))), main = main.legends[4], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
      readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$geno.pred, ncol = length(columns), nrow = length(rows))), main = main.legends[5], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
      readline("Press return for next page....")
    # 
    suppressWarnings(hist(unique(geno.pred), main = main.legends[6], xlab = main.legends[6], ...))        
    title("")
    mtext(main, cex = 1.5, outer = TRUE, side = 3)
    invisible(df)
  }


##########################
#         Example
##########################


# library(SpATS)
# library(lme4)
# data("wheatdata")
# summary(wheatdata)
# 
# wheatdata <- coords(col = "col", row="row", data = wheatdata)
# 
# Mo <- lmer(formula = yield ~ 1+(1|geno) + rep, data = wheatdata,REML = T )
# 
# # Spatial plot for lme4 
# lme4.plot(Mo, wheatdata, col="col", row="row", gen="geno")

