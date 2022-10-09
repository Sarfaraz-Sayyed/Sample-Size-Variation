## Creating variance-covariance matrix for
## 1) compound symmetry and
## 2) Auto-regressive of order 1 AR(1)
## sample call for 3x3(default) variance covariance matrix
## for type CS with variance=1 and correlation=0.5 - vcov_mat()
## sample call for 5x5 variance covariance matrix
## for type AR1 - vcov_mat(VAR=c(25,25,25,25,25), rho=0.5, type="AR1")

vcov_mat <- function(varvec = c(1, 1, 1), rho = 0.5, type = "CS") {
  var_cov <- matrix(nrow = length(varvec), ncol = length(varvec), byrow = TRUE)
  for (i in 1:length(varvec)) {
    for (j in 1:length(varvec)) {
      var_cov[i, j] <- if (i == j) {
        varvec[i]
      } else {
        if (type == "CS") {
          rho * sqrt(varvec[i] * varvec[j])
        } else {
          (rho**(i + j - 2)) * sqrt(varvec[i] * varvec[j])
        }
      }
    }
  }
  return(var_cov)
}
