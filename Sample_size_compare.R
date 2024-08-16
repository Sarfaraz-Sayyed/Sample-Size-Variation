## Program for generating sample size for single time-point post-baseline vs
## Multiple time post-baseline

library(ggplot2)

## Creating variance-covariance matrix for
## 1) compound symmetry and
## 2) Auto-regressive of order 1 AR(1)
## sample call for 3x3(default) variance-covariance matrix
## for type CS with variance=1 and correlation=0.5 - vcov_mat()
## sample call for 5x5 variance-covariance matrix
## for type AR1 - vcov_mat(varvec=c(25,25,25,25,25), rho=0.5, type="AR1")

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


## Single post-baseline time-point case
power <- 0.85
alpha <- 0.05
## effect size = (mu1 - mu2)/sigma
effect_size <- (0.9) / sqrt(12.95)
npergrp <- (2 * ((qnorm(1 - (alpha / 2)) + qnorm(power))**2)) / (effect_size)**2
npergrp

## Multiple post-baseline time-point case
rm_sample <- function(alpha = 0.025,
                      power = 0.85,
                      meandiff = c(0, 0.393, 0.785, 1.2, 1.57),
                      contrast = c(-1, 0, 0, 0, 1),
                      var = vcov_mat(
                        varvec = c(26.52, 26.52, 26.52, 26.52, 26.52),
                        rho = 0.5,
                        type = "AR1"
                      ))
{
  {
    contdiff <- t(contrast) %*% meandiff
    contvar <- t(contrast) %*% var %*% contrast
    npergrp_rm <- (2 * ((qnorm(1 - alpha) +
      qnorm(power))**2) * contvar) / (contdiff**2)
  }
  return(round(npergrp_rm))
}

sample_size <- function(rho = seq(0, 1, by = 0.05),
                        alpha = 0.025, # 2 tailed  5% type I error
                        power = 0.85, # 85% power
                        meandiff = c(0, 0.5, 0.9, 0.9),
                        contrast = c(-1, 0, 0, 1),
                        var_mean = c(1, 9.43, 7.08, 12.95),
                        type = "AR1") {
  rm_dat <- matrix(nrow = 1, ncol = 2, byrow = TRUE)

  for (i in rho) {
    rm <- rm_sample(
      alpha = alpha,
      power = power,
      meandiff = meandiff,
      contrast = contrast,
      var = vcov_mat(varvec = var_mean, rho = i, type = type)
    )
    rm_dat <- rbind(rm_dat, c(i, rm))
  }
  return(rm_dat)
}

meandiff3 <- c(0, 0.9, 0.9)
meandiff4 <- c(0, 0.85, 0.95, 0.9)
meandiff5 <- c(0, 0.85, 0.98, 0.87, 0.9)
meandiff6 <- c(0, 0.83, 0.97, 0.88, 0.92, 0.9)
meandiff8 <- c(0, 0.82, 0.95, 0.87, 0.92, 0.89, 0.94, 0.9)
meandiff10 <- c(0, 0.82, 0.93, 0.92, 0.87, 0.92, 0.89, 0.93, 0.91, 0.9)
var_3 <- c(replicate(3, 12.95))
var_4 <- c(replicate(4, 12.95))
var_5 <- c(replicate(5, 12.95))
var_6 <- c(replicate(6, 12.95))
var_8 <- c(replicate(8, 12.95))
var_10 <- c(replicate(10, 12.95))

run4cas <- function(ntime = 3,
                    meanvector = meandiff3,
                    varcovar = VAR3) {
  sample1 <- data.frame(sample_size(
    rho = seq(0.1, 0.95, by = 0.05),
    meandiff = meanvector,
    contrast = c(-1, replicate(ntime - 2, 0), 1),
    var_mean = varcovar,
    type = "AR1"
  ))

  sample2 <- data.frame(sample_size(
    rho = seq(0.1, 0.95, by = 0.05),
    meandiff = meanvector,
    contrast = c(-1, replicate(ntime - 1, (1 / (ntime - 1)))),
    var_mean = varcovar,
    type = "AR1"
  ))

  sample3 <- data.frame(sample_size(
    rho = seq(0.1, 0.95, by = 0.05),
    meandiff = meanvector,
    contrast = c(-1, replicate(ntime - 2, 0), 1),
    var_mean = varcovar,
    type = "CS"
  ))

  sample4 <- data.frame(sample_size(
    rho = seq(0.1, 0.95, by = 0.05),
    meandiff = meanvector,
    contrast = c(-1, replicate(ntime - 1, (1 / (ntime - 1)))),
    var_mean = varcovar,
    type = "CS"
  ))
  return(list(sample1, sample2, sample3, sample4))
}

# 3 time points
dataout3 <- run4cas(
  ntime = 3,
  meanvector = meandiff3,
  varcovar = var_3
)
# 4 time points
dataout4 <- run4cas(
  ntime = 4,
  meanvector = meandiff4,
  varcovar = var_4
)
# 5 time points
dataout5 <- run4cas(
  ntime = 5,
  meanvector = meandiff5,
  varcovar = var_5
)
# 6 time points
dataout6 <- run4cas(
  ntime = 6,
  meanvector = meandiff6,
  varcovar = var_6
)
# 8 time points
dataout8 <- run4cas(
  ntime = 8,
  meanvector = meandiff8,
  varcovar = var_8
)
# 10 time points
dataout10 <- run4cas(
  ntime = 10,
  meanvector = meandiff10,
  varcovar = var_10
)

## combining all the AR(1) outputs into a single data frame
combined <- rbind(
  cbind(data.frame(dataout3[1]), type = "AR1_diff( 3)"),
  cbind(data.frame(dataout4[1]), type = "AR1_diff( 4)"),
  cbind(data.frame(dataout5[1]), type = "AR1_diff( 5)"),
  cbind(data.frame(dataout6[1]), type = "AR1_diff( 6)"),
  cbind(data.frame(dataout8[1]), type = "AR1_diff( 8)"),
  cbind(data.frame(dataout10[1]), type = "AR1_diff(10)"),
  cbind(data.frame(dataout3[2]), type = "AR1_mean( 3)"),
  cbind(data.frame(dataout4[2]), type = "AR1_mean( 4)"),
  cbind(data.frame(dataout5[2]), type = "AR1_mean( 5)"),
  cbind(data.frame(dataout6[2]), type = "AR1_mean( 6)"),
  cbind(data.frame(dataout8[2]), type = "AR1_mean( 8)"),
  cbind(data.frame(dataout10[2]), type = "AR1_mean(10)")
)

## combining all the CS outputs into a single data frame
combined1 <- rbind(
  cbind(data.frame(dataout3[3]), type = "CS_diff( 3)"),
  cbind(data.frame(dataout4[3]), type = "CS_diff( 4)"),
  cbind(data.frame(dataout5[3]), type = "CS_diff( 5)"),
  cbind(data.frame(dataout6[3]), type = "CS_diff( 6)"),
  cbind(data.frame(dataout8[3]), type = "CS_diff( 8)"),
  cbind(data.frame(dataout10[3]), type = "CS_diff(10)"),
  cbind(data.frame(dataout3[4]), type = "CS_mean( 3)"),
  cbind(data.frame(dataout4[4]), type = "CS_mean( 4)"),
  cbind(data.frame(dataout5[4]), type = "CS_mean( 5)"),
  cbind(data.frame(dataout6[4]), type = "CS_mean( 6)"),
  cbind(data.frame(dataout8[4]), type = "CS_mean( 8)"),
  cbind(data.frame(dataout10[4]), type = "CS_mean(10)")
)

# plotting the data for CS related outputs
ggplot(data = na.omit(combined1[1:3]), aes(x = X1, y = X2, group = type)) +
  scale_shape_manual(values = 1:12) +
  geom_point(aes(color = factor(type), shape = factor(type))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 200)) +
  labs(
    x = "Correlation(rho)",
    y = "Sample Size",
    title = "Relation between different contrast, correlation vs Sample size",
    subtitle = "Constant correlation over time(Compound Symmetry)",
    caption = c(
      "Contrast: _diff -> last visit - baseline",
      "          _mean -> average effect of all post baseline visits",
      "red line represent the sample with univariate case"
    )
  ) +
  geom_hline(
    yintercept = npergrp,
    linetype = "solid",
    color = "red",
    size = 0.7
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    strip.text.x = element_blank()
  )

# plotting the data for AR(1) related outputs
ggplot(data = na.omit(combined[1:3]), aes(x = X1, y = X2, group = type)) +
  scale_shape_manual(values = 1:12) +
  geom_point(aes(color = factor(type), shape = factor(type))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1200, by = 200)) +
  labs(
    x = "Correlation(rho)",
    y = "Sample Size",
    title = "Relation between different contrast, correlation vs Sample size",
    subtitle = "Decreasing correlation over time(Auto regressive of order 1)",
    caption = c(
      "Contrast: _diff -> last visit - baseline",
      "          _mean -> average effect of all post baseline visits",
      "red line represent the sample with univariate case"
    )
  ) +
  geom_hline(
    yintercept = npergrp,
    linetype = "solid",
    color = "red",
    size = 0.7
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    strip.text.x = element_blank()
  )
