## Creating variance- covariance matrix for compound symmetry and Auto-regressive of order 1 AR(1)

VCOV_mat = function(VAR=c(1,1,1), rho=0.5, type='CS'){
  VARCOV<-matrix(nrow=length(VAR),ncol=length(VAR),byrow=TRUE)
  for (i in 1:length(VAR)) {
    for (j in 1:length(VAR)){ 
      VARCOV[i,j] <- if (i == j) {VAR[i]} 
      else {if (type == 'CS') {rho*sqrt(VAR[i]*VAR[j])} else {(rho**(i+j-2))*sqrt(VAR[i]*VAR[j])}}
      
    }
  }
  return(VARCOV)
}

# sample call for 3x3(default) variance covariance matrix of type CS with variance=1 and correlation=0.5
# VCOV_mat()
# sample call for 5x5 variance covariance matrix of type AR1
# VCOV_mat(VAR=c(25,25,25,25,25), rho=0.5, type='AR1')


# test <- runif(1000)
# test1 <- data.frame(test)
# summary(test1$test <= 0.15)
# test1
# VCOV_mat(VAR=c(12.95,12.95, 12.95, 12.95), rho=0.8, type='AR1')

