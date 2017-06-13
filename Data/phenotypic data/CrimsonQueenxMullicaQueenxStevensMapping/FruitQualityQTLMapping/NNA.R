NNA <- function(datar, plantid, repetition, covar, response){
  
  v1 <- which(names(datar) %in% repetition)
  v2 <- which(names(datar) %in% covar)
  v3 <- which(names(datar) %in% plantid)
  v4 <- which(names(datar) %in% response)
  
  just <- c("Row", "Col")
  v5 <- which(just %in% names(datar))
  
  if(length(v1) == 0 | length(v2) == 0 | length(v3) == 0 |length(v4) == 0){
    warning("Please check the names of your dataframe match the names you used as arguments we did not find at least one match")
    stop
  }
  if(length(v5) != 2){
    warning("'Row' and 'Col' variables indicating row and column position in the field cannot be excempted, please add tem to your dataframe with uppercase in the first letter.")
    stop
  }

  
  repi <- unique(datar[,v1]) 
  covariate <- unique(datar[,v2]) 
  
  names(datar)[v3] <- "Genotype"
  names(datar)[v2] <- paste("Covariate",covar,sep="_")
  
  library(lme4)
  library(plyr)
  mo1 <- lmer(datar[,v4] ~(1|Genotype), data=datar)
 
  pro <- as.numeric(names(resid(mo1))) 
  datar$res <- NA 
  datar[pro,"res"] <- resid(mo1) 
  
  ###########################################
  # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
  ###########################################
  last.list <- list(NA)
  for(l in 1:length(covariate)){ 
    
    data <- datar[which(datar[,v2] == paste(covariate[l])),]
    rownames(data) <- NULL
    
    data$NS <- NA
    data$EW <- NA
    r <- unique(data$Row)
    c <- unique(data$Col)
    a <- repi 
    resp <- response
    
    seas.list <- list(NA)
    for(k in 1:length(a)){
      season <- a[k]
      data2 <- data[which(data[,v1] == paste(season)),]
      for(i in r){ 
        for(j in c){ 
         
          n <- data2[intersect(which(data2$Row == i) ,  which(data2$Col == j-1)),resp] #resp
          s <- data2[intersect(which(data2$Row == i) ,  which(data2$Col == j+1)),resp]
          e <- data2[intersect(which(data2$Row == i+1) ,  which(data2$Col == j)),resp]
          w <- data2[intersect(which(data2$Row == i-1) ,  which(data2$Col == j)),resp]
          ## getting north south estimate
          if(length(n) > 0 & length(s) > 0){
            ns <- mean(c(n,s),na.rm=T)
          }
          if(length(n) == 0 & length(s) > 0){
            ns <- s
          }
          if(length(n) > 0 & length(s) == 0){
            ns <- n
          }
          ## getting east west estimate
          if(length(e) > 0 & length(w) > 0){
            ew <- mean(c(e,w),na.rm=T)
          }
          if(length(e) == 0 & length(w) > 0){
            ew <- w
          }
          if(length(e) > 0 & length(w) == 0){
            ew <- e
          }
          data2[intersect(which(data2$Row == i) ,  which(data2$Col == j)),"NS"] <- ns
          data2[intersect(which(data2$Row == i) ,  which(data2$Col == j)),"EW"] <- ew
        }
      }
     
      seas.list[[k]] <- data2
    }
    newm <- seas.list[[1]]
    
   
    if(length(seas.list) > 1){
      for(t in 2:length(seas.list)){
        newm <- rbind(newm,seas.list[[t]])
      }
    }
   
    last.list[[l]] <- newm
    
  }
  
  
  all <- last.list[[1]]
  if(length(last.list) > 1){
    for(s in 2:length(last.list)){
      all <- rbind(all,last.list[[s]])
    }
  }
  return(all)
}