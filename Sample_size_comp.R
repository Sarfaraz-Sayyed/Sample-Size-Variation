## Program for generating sample size for single time-point post-baseline vs multiple time post baseline

# library(tidyverse)
library(ggplot2)
source("/cloud/project/varcov.R")

#Year 9 relative to Year 6 between Z9 and Z6P3
#Approximately 195 patients per treatment Arm
#85% power to detect a difference between treatment groups 
#5% two-sided significance level 
#assuming difference between treatment groups with percentage change in BMD of the total hip at will be 1.57%. 
#with 20% dropout rate we would need 244 in each arm


# single post baseline time-point case
power <- 0.85 
alpha=0.05
effect_size <-  (0.9)/sqrt(12.95)  # (mu1 - mu2)/sigma

Npergrp <- (2*((qnorm(1 - (alpha/2)) + qnorm(power))**2))/(effect_size)**2
Npergrp 

## multiple post baseline time-point case
RMsample = function(alpha=0.025, 
                    power=0.85, 
                    meandiff = c(0, 0.393, 0.785, 1.2, 1.57), 
                    contrast = c(-1, 0, 0, 0, 1), 
                    var=VCOV_mat(VAR=c(26.52,26.52,26.52,26.52,26.52), rho=0.5, type='AR1')
) 
{
  {
    contdiff = t(contrast)%*%meandiff
    contvar = t(contrast)%*%var%*%contrast
    Npergrp_RM = (2*((qnorm(1 - alpha) + qnorm(power))**2) * contvar)/(contdiff**2)
  }
  return(round(Npergrp_RM))
}

sample.size = function(rho=seq(0, 1, by=0.05),
                       alpha=0.025,           #2 tailed  5% type I error
                       power=0.85,             #85% power
                       meandiff = c(0, 0.5, 0.9, 0.9),
                       contrast = c(-1, 0, 0, 1),
                       VAR=c(1,  9.43, 7.08, 12.95),
                       type='AR1'
)
{
  
  RM_dat<-matrix(nrow=1,ncol=2,byrow=TRUE)
  
  for(i in rho)
  {
    RM <- RMsample(alpha=alpha,              
                   power=power,              
                   meandiff = meandiff, 
                   contrast = contrast,
                   var=VCOV_mat(VAR=VAR, rho=i, type=type)
    ) 
    RM_dat <- rbind(RM_dat, c(i, RM))
  }
  return(RM_dat)
}

meandiff3 <- c(0, 0.9, 0.9)
meandiff4 <- c(0, 0.85, 0.95, 0.9)
meandiff5 <- c(0, 0.85, 0.98, 0.87, 0.9)
meandiff6 <- c(0, 0.83, 0.97, 0.88, 0.92, 0.9)
meandiff8 <- c(0, 0.82, 0.95, 0.87, 0.92, 0.89, 0.94, 0.9)
meandiff10 <- c(0, 0.82, 0.93, 0.92, 0.87, 0.92, 0.89, 0.93, 0.91, 0.9)
VAR3=c(replicate(3,12.95))
VAR4=c(replicate(4,12.95))
VAR5=c(replicate(5,12.95))
VAR6=c(replicate(6,12.95))
VAR8=c(replicate(8,12.95))
VAR10=c(replicate(10,12.95))

run4cas = function(ntime = 3,
                   meanvector=meandiff3,
                   RMCOVTYP = 'CS',
                   varcovar = VAR3)
{
  sample1 <- data.frame(sample.size(rho=seq(0, 1, by=0.05),
                                    #meandiff = c(0, .5, 1),
                                    meandiff = meanvector,
                                    contrast = c(-1, replicate(ntime-2,0), 1),
                                    VAR=varcovar,
                                    type='AR1'
  )) 
  
  sample2 <- data.frame(sample.size(rho=seq(0, 1, by=0.05),
                                    meandiff = meanvector,
                                    contrast = c(-1, replicate(ntime-1,(1/(ntime-1)))),
                                    VAR=varcovar,
                                    type='AR1'
  ))
  
  sample3<- data.frame(sample.size(rho=seq(0, 1, by=0.05),
                                   meandiff = meanvector,
                                   contrast = c(-1, replicate(ntime-2,0), 1),
                                   VAR=varcovar,
                                   type='CS'
  ))
  
  sample4 <- data.frame(sample.size(rho=seq(0, 1, by=0.05),
                                    meandiff = meanvector,
                                    contrast = c(-1, replicate(ntime-1,(1/(ntime-1)))),
                                    VAR=varcovar,
                                    type='CS'
  ))
  return(list(sample1, sample2, sample3, sample4))
}

# 3 time points
dataout3 <- run4cas(ntime = 3,
                    meanvector=meandiff3,
                    RMCOVTYP = 'CS',
                    varcovar = VAR3)
# 4 time points
dataout4 <- run4cas(ntime = 4,
                    meanvector=meandiff4,
                    RMCOVTYP = 'CS',
                    varcovar = VAR4)
# 5 time points
dataout5 <- run4cas(ntime = 5,
                    meanvector=meandiff5,
                    RMCOVTYP = 'CS',
                    varcovar = VAR5)
# 6 time points
dataout6 <- run4cas(ntime = 6,
                    meanvector=meandiff6,
                    RMCOVTYP = 'CS',
                    varcovar = VAR6)
# 8 time points
dataout8 <- run4cas(ntime = 8,
                    meanvector=meandiff8,
                    RMCOVTYP = 'CS',
                    varcovar = VAR8)
# 10 time points
dataout10 <- run4cas(ntime = 10,
                     meanvector=meandiff10,
                     RMCOVTYP = 'CS',
                     varcovar = VAR10)


combined <- rbind(cbind(data.frame(dataout3[1]),type='AR1_diff( 3)'),
                  cbind(data.frame(dataout4[1]),type='AR1_diff( 4)'),
                  cbind(data.frame(dataout5[1]),type='AR1_diff( 5)'),
                  cbind(data.frame(dataout6[1]),type='AR1_diff( 6)'),
                  cbind(data.frame(dataout8[1]),type='AR1_diff( 8)'),
                  cbind(data.frame(dataout10[1]),type='AR1_diff(10)'),
                  cbind(data.frame(dataout3[2]),type='AR1_mean( 3)'),
                  cbind(data.frame(dataout4[2]),type='AR1_mean( 4)'),
                  cbind(data.frame(dataout5[2]),type='AR1_mean( 5)'),
                  cbind(data.frame(dataout6[2]),type='AR1_mean( 6)'),
                  cbind(data.frame(dataout8[2]),type='AR1_mean( 8)'),
                  cbind(data.frame(dataout10[2]),type='AR1_mean(10)')
)

combined1 <- rbind(cbind(data.frame(dataout3[3]),type='CS_diff( 3)'),
                   cbind(data.frame(dataout4[3]),type='CS_diff( 4)'),
                   cbind(data.frame(dataout5[3]),type='CS_diff( 5)'),
                   cbind(data.frame(dataout6[3]),type='CS_diff( 6)'),
                   cbind(data.frame(dataout8[3]),type='CS_diff( 8)'),
                   cbind(data.frame(dataout10[3]),type='CS_diff(10)'),
                   cbind(data.frame(dataout3[4]),type='CS_mean( 3)'),
                   cbind(data.frame(dataout4[4]),type='CS_mean( 4)'),
                   cbind(data.frame(dataout5[4]),type='CS_mean( 5)'),
                   cbind(data.frame(dataout6[4]),type='CS_mean( 6)'),
                   cbind(data.frame(dataout8[4]),type='CS_mean( 8)'),
                   cbind(data.frame(dataout10[4]),type='CS_mean(10)')
)

#plotting the data
tiff("/cloud/project/Figure_1_1.tiff", units="in", width=6.1, height=4, res=310, compression = 'lzw')
ggplot(data = na.omit(combined[1:3]), aes(x = X1, y = X2 , group=type
)) + 
  scale_shape_manual(values = 1:12)+
  geom_point(aes(color = factor(type), shape=factor(type))) + 
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1200, by = 200)) +
  labs(
    x = "Correlation(rho)",
    y = "Sample Size",
    title = "Relation between different contrast, correlation vs Sample size",
    subtitle = "Decreasing correlation over time(Auto regressive of order 1)",
    caption = c("Contrast: _diff -> last visit - baseline; _mean -> average effect of all post baseline visits",
                "red line represent the sample with univariate case")
  ) +
  geom_hline (yintercept=287, linetype='solid', color='red', size=0.7) +
  theme(
    legend.title = element_blank(), 
    legend.text=element_text(size=7),
    strip.text.x = element_blank())
dev.off()
#________________________________________________________________________
tiff("/cloud/project/Figure_1_2.tiff", units="in", width=6.1, height=4, res=310, compression = 'lzw')
ggplot(data = na.omit(combined1[1:3]), aes(x = X1, y = X2 , group=type
)) + 
  scale_shape_manual(values = 1:12)+
  geom_point(aes(color = factor(type), shape=factor(type))) + 
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1600, by = 200)) +
  labs(
    x = "Correlation(rho)",
    y = "Sample Size",
    title = "Relation between different contrast, correlation vs Sample size",
    subtitle = "Constant correlation over time(Compound Symmetry)",
    caption = c("Contrast: _diff -> last visit - baseline; _mean -> average effect of all post baseline visits",
                "red line represent the sample with univariate case")
  ) +
  geom_hline (yintercept=287, linetype='solid', color='red', size=0.7)+
  theme(legend.title = element_blank(), 
        legend.text=element_text(size=7),
        strip.text.x = element_blank())
#___________________________________________________________________________________
dev.off()