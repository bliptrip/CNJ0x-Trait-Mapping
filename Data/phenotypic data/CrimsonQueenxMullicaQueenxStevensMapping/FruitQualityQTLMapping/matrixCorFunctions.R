CorMatrixMonthYear<-function(x,y,z,l,n){
      
      #did the user didn't supply an additional factor vector
      if(l=="nope"){
            print(l)
            l<-1
            #create the matrix to be filled
            mat<-matrix(1,length(which(x[,z]==unique(x[,z])[2])),length(y)*length(unique(x[,z])))
            #number of factors in column
            a<-length(unique(x[,z]))
            b<-1
            
            counter<-1
            
            for(i in 1:length(y)){
                  for(j in 1:a){
                        mat[,counter]<-x[(x[,z]==unique(x[,z])[j]),y[i]]   
                        counter<-counter+1
                        
                  }
                  
            }
            #Vector of Z names
            colzvec<-c(rep(as.vector(unique(x[,z])),length(y)*b))
            #Get Vector of Trait Column Names
            colyvec<-c(rep(NA,length(y)*a*b))  
            count<-1
            for(i in 1:length(y)){
                  for(j in 1:(a*b)){
                        colyvec[count]<-colnames(x)[y[i]]
                        count<-count+1
                  }
            }
            #row names matrix
            rownames(mat)<-x[(x[,z]==unique(x[,z])[1]),n]
            #column names matrix
            colnames(mat)<-paste(colyvec,colzvec,sep=".")
            
            #the user supplied an additional factor vector      
      } else {
            print("l was not nope")
            #create data matrix
            mat<-matrix(1,length(x[(x[,z]==unique(x[,z])[2])&(x[,l]==unique(x[,l])[2]),1]),length(y)*length(unique(x[,z]))*length(unique(x[,l]))) 
            #number of factors in column z
            a<-length(unique(x[,z]))
            #number of factors in column l
            b<-length(unique(x[,l]))
            
            counter<-1
            # fill the matrix
            for(i in 1:length(y)){
                  for(j in 1:a){
                        for(k in 1:b){
                              mat[,counter]<-x[(x[,z]==unique(x[,z])[j])&(x[,l]==unique(x[,l])[k]),y[i]]   
                              counter<-counter+1
                              
                        }
                  }
                  
            } 
            #vector of names from factors in column z and l
            counted<-1
            colzvec<-c(rep(as.vector(unique(x[,z])),length(y)*b))
            for(i in 1:length(y)){
                  for(j in 1:a){
                        for(k in 1:b){
                              colzvec[counted]<-as.character(unique(x[,z])[j])
                              counted<-counted+1
                        }
                  }
                        
            }
                       
            collvec<-c(rep(as.vector(unique(x[,l])),length(y)*a))
            #Get Vector of Trait Column Names
            colyvec<-c(rep(NA,length(y)*a*b))  
            count<-1
            for(i in 1:length(y)){
                  for(j in 1:(a*b)){
                        colyvec[count]<-colnames(x)[y[i]]
                        count<-count+1
                  }
            }
            rownames(mat)<-x[(x[,z]==unique(x[,z])[1])&(x[,l]==unique(x[,l])[1]),n]
            colnames(mat)<-paste(colyvec,colzvec,collvec,sep=".")
      }
      
      spcor<-round(cor(mat,method="spearman",use="pairwise.complete.obs"),3)
      
      theList<-list(theMatrix=mat,Spearman=spcor)
      return(theList)                  
}  

#x is the dataframe
Heatplots1<-function(x){
      z<-colnames(x)[6]
      # heatplots to asses spatial variation
      all2011<-levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 All residuals distribution"), data=x[which(x[,4] == "y2011"),])
      all2012<-levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 All residuals distribution"), data=x[which(x[,4] == "y2012"),])
      all2013<-levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2013 All residuals distribution"), data=x[which(x[,4] == "y2013"),])
      #September
      sept2011<-levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 Sept. residuals distribution"), data=x[(x[,4] == "y2011")&(x[,5]=="September"),])
      sept2012=levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Sept. residuals distribution"), data=x[(x[,4] == "y2012")&(x[,5]=="September"),])
      sept2013=levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2013 Sept. residuals distribution"), data=x[(x[,4] == "y2013")&(x[,5]=="September"),])
      #October
      oct2011=levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 Oct. residuals distribution"), data=x[(x[,4] == "y2011")&(x[,5]=="October"),])
      oct2012=levelplot(res ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Oct. residuals distribution"), data=x[(x[,4] == "y2012")&(x[,5]=="October"),])
      oct2013=levelplot(res~ Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Oct. residuals distribution"), data=x[(x[,4] == "y2013")&(x[,5]=="October"),])
      return(list(all2011=all2011,all2012=all2012,all2013=all2013,
                  sept2011=sept2011,sept2012=sept2012,sept2013=sept2013,
                  oct2011=oct2011,oct2012=oct2012,oct2013=oct2013))
}

#x is the dataframe
Heatplots2<-function(x){
      z<-colnames(x)[6]
      # heatplots to asses spatial variation
      all2011<-levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 All residuals distribution"), data=x[which(x[,4] == "y2011"),])
      all2012<-levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 All residuals distribution"), data=x[which(x[,4] == "y2012"),])
      all2013<-levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2013 All residuals distribution"), data=x[which(x[,4] == "y2013"),])
      #September
      sept2011<-levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 Sept. residuals distribution"), data=x[(x[,4] == "y2011")&(x[,5]=="September"),])
      sept2012=levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Sept. residuals distribution"), data=x[(x[,4] == "y2012")&(x[,5]=="September"),])
      sept2013=levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2013 Sept. residuals distribution"), data=x[(x[,4] == "y2013")&(x[,5]=="September"),])
      #October
      oct2011=levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2011 Oct. residuals distribution"), data=x[(x[,4] == "y2011")&(x[,5]=="October"),])
      oct2012=levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Oct. residuals distribution"), data=x[(x[,4] == "y2012")&(x[,5]=="October"),])
      oct2013=levelplot(res2 ~Col*Row,col.regions = heat.colors(10, alpha = 1), main=c(z," 2012 Oct. residuals distribution"), data=x[(x[,4] == "y2013")&(x[,5]=="October"),])
      return(list(all2011=all2011,all2012=all2012,all2013=all2013,
                  sept2011=sept2011,sept2012=sept2012,sept2013=sept2013,
                  oct2011=oct2011,oct2012=oct2012,oct2013=oct2013))
}

#x is a summary object and y is a rand object 
byYear<-function(x,y){
      
      a<-c(round(x$varcor$Genotype[1],digits=3),
           round(x$varcor$Covariate_Year[1],digits=3),
           round(x$varcor[[1]][[1]],digits=3),
           round(x$varcor$NS[1],digits=3),
           round(x$varcor$EW[1],digits=3),
           round(x$sigma[[1]]^2, digits=3))
      
      b<-c(y[[1]][[3]])
      d<-c(rep("*",5))
      for(i in 1:5){
            if(b[i]<=0.001){d[i]<-"***"
            } else if(b[i]<=0.01){d[i]<-"**"
            } else if(b[i]<=0.05){d[i]<-"*"
            } else if(b[i]<=0.1){d[i]<-"."
            } else {d[i]<-""
            }
            
      }
      d<-c(d,"")
      e<-c(paste(a,d,sep=" "))
      return(e)
}

#x is summary object and y is rand object
byBothMonths<-function(x,y){
      a<-c(round(x$varcor$Genotype[1],digits=3),
           round(x$varcor$Month[1],digits=3),
           round(x$varcor[[1]][[1]],digits=3),
           round(x$varcor$NS[1],digits=3),
           round(x$varcor$EW[1],digits=3),
           round(x$sigma[[1]]^2, digits=3))
      
      b<-c(y[[1]][[3]])
      d<-c(rep("*",5))
      for(i in 1:5){
            if(b[i]<=0.001){d[i]<-"***"
            } else if(b[i]<=0.01){d[i]<-"**"
            } else if(b[i]<=0.05){d[i]<-"*"
            } else if(b[i]<=0.1){d[i]<-"."
            } else {d[i]<-""
            }
            
      }
      d<-c(d,"")
      e<-c(paste(a,d,sep=" "))
      return(e) 
}

#x is summary object and y is rand object
oneMonth<-function(x,y){
      a<-c(round(x$varcor$Genotype[1],digits=3),
           round(x$varcor$NS[1],digits=3),
           round(x$varcor$EW[1],digits=3),
           round(x$sigma[[1]]^2, digits=3))
      
      b<-c(y[[1]][[3]])
      d<-c(rep("*",3))
      for(i in 1:3){
            if(b[i]<=0.001){d[i]<-"***"
            } else if(b[i]<=0.01){d[i]<-"**"
            } else if(b[i]<=0.05){d[i]<-"*"
            } else if(b[i]<=0.1){d[i]<-"."
            } else {d[i]<-""
            }
            
      }
      d<-c(d,"")
      e<-c(paste(a,d,sep=" "))
      return(e) 
}