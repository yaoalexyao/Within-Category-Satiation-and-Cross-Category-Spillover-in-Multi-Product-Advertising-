########################################################################################
# This file does estimate the multi-choice model



###function to calculate probability
multichoice = function(xyvector) { 

  y = xyvector[(length(xyvector)/3+1):(length(xyvector)/3*2)]
  x = xyvector[1:(length(xyvector)/3)]
  avail = xyvector[(length(xyvector)/3*2+1):length(xyvector)]
  
  y=y[avail==1]
  x=x[avail==1]
  
  u0 = sum(x[y==0])
  u1 = x[y==1]
  for (i in 1:length(u1) ) {
    if (i == 1) {
     
      if (length(u1)==1) {
        prob0 = u1/c(u1+u0)
      } else {  
        matrixu0pluscombu1 = rbind( u0, combn(u1,i) ) 
        prob0 = 1 + (-1)^i*( sum(matrixu0pluscombu1[1,]/colSums(matrixu0pluscombu1) )  )      
      }
    } else {
      matrixu0pluscombu1 = rbind( u0, combn(u1,i) )
      prob0 = prob0 + (-1)^i*( sum(matrixu0pluscombu1[1,]/colSums(matrixu0pluscombu1) )  )      
    }  
  }
  return(prob0)
}



database = read.csv("cleaned data.csv",header=TRUE)

rprice = database$rprice
click = database$click

rpricema = cbind( t( matrix(rprice, nrow = 8) ), 0) 


###click choice data
clickma = t( matrix(click, nrow = 8) )
ymatrix = cbind( clickma, rowSums(clickma) == 0)


category1029 = t( matrix( database$subcat_id3, nrow = 8) )
clickcategory1029 = category1029*clickma

NROW = nrow(category1029)

clickcategory = category1029*clickma
avail = matrix(0, nrow = NROW, ncol = 8)


for (i in 1:NROW) {
  temp = clickcategory[i, ]
  if (sum(temp) != 0) {
    tempno0 = temp[temp!=0]
    tempunique = unique(tempno0)
    for (j in 1:length(tempunique)) {
      avail[i, category1029[i,]==tempunique[j] ] = 1
    }
  }
}

avail = cbind(avail, 1)

rpricemashort = rpricema[rowSums(avail[,1:8])>0, 1:8]
ymatrixshort = ymatrix[rowSums(avail[,1:8])>0, 1:8]
availshort = avail[rowSums(avail[,1:8])>0, 1:8]


####llh to maximize
multichoicellh <- function(par){
  ###xmatrix and ymatrix should be pre-determined, including the outside good. xmatrix should product attriabutes. ymatrix should be the choice
  xmatrix = exp( rpricemashort * par[1] )
  
  xymatrix = cbind(xmatrix,ymatrixshort,availshort)
  
  llh = -sum( log(apply(xymatrix, 1, multichoice ) ) )
  llh
  
}



set.seed(1011)
init=runif(1)

multichoice1=optim(init,multichoicellh,method="BFGS",hessian=T,control=list(maxit=1000 )) 
multichoice1
stderr=sqrt(abs(diag(solve(multichoice1$hessian))))
stderr

n = dim(ymatrix)[1]
BICmultichoice1=2*multichoice1$value+1*log( sum(rowSums(availshort)>0) ) 
BICmultichoice1

signif=abs(multichoice1$par)>=1.96*stderr
signif2=abs(multichoice1$par)>=2.58*stderr

writeout=cbind(multichoice1$par,stderr,signif2,signif,multichoice1$value,BICmultichoice1,multichoice1$convergence)
writeout


