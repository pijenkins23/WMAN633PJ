#1. Describe a sampling procedure that may have given rise to this dataset.

Randomly selecting 100 sites and doing 3 replicate surveys in which you collect a species is 
detected or not detected (1 or 0) at each site for each replicate.
Collecting site abundance covariates x1 and x2 once per site.
Record detection covariatates detcov1 and detcove 2 for each replicate survey.  



#2. Import data and fit an occupancy model that assumes detection probability is an additive function of
#obscovs1 and obscovs2; and that occupancy probability is an additive function of x1 and x2.

y <- read.csv('detect.csv')
head(y)

library(unmarked)
# converting count data into a matrix
# matrix format required by unmarkedFrameoccu
sosp_mat <- as.matrix(y)
nmix_data <- unmarkedFrameOccu(y = sosp_mat)


# Covariate s
##################

obs_covs1 <- read.csv('obscovs1.csv')
obs_covs2 <- read.csv('obscovs2.csv')
head(obs_covs1)
head(obs_covs2)

obs_covs <- list(
  obscovs1 = data.frame(obs_covs1),
  obscovs2 = data.frame(obs_covs2)
  
)


site_covs <- read.csv('sitecovs.csv')
head(site_covs)




nmix_data <- unmarkedFrameOccu(y = sosp_mat, # detection / non-detection
                                 siteCovs = site_covs, # site-level covs
                                 obsCovs = obs_covs) 


fit <- occu(formula = ~ obscovs1 + obscovs2 ~ x1 + x2,
              data = nmix_data)
fit

#3. Use contrasts to determine if occupancy probability different when x1 = 2 vs. when x1 = -2?
  

b <- coef(fit)


cm <- matrix(c(0, 2, 0, 0, -2, 0), nrow = 2, byrow = T ) #unsure? 
cnt <- linearComb(obj = fit, coefficients = cm, type = 'state')

pnorm(-1 * abs(coef(cnt) / SE(cnt))) * 2

new_psi <- data.frame(x1 = rep(mean(site_covs$x1), length.out = 1),
                      x2 = rep(c(2,-2), 1))
(predict <- predict(object = fit, newdata = new_psi, type = 'state'))


  
  4. Use model selection to compare the following 4 models. Which model is the "top" model? How do you
know?
  (a) ∼ obscovs1 + obscovs2 ∼ x1 + x2
(b) ∼ obscovs1 + obscovs2 ∼ x1
(c) ∼ obscovs1 + obscovs2 ∼ x2
(d) ∼ obscovs1 + obscovs2 ∼ 1


fit1 <- occu(formula = ~ obscovs1 + obscovs2 ~ x1 + x2,
               data = nmix_data)
fit2 <- occu(formula = ~ obscovs1 + obscovs2 ~ x1 ,
               data = nmix_data)
fit3 <- occu(formula = ~  obscovs1 + obscovs2 ~  x2,
               data = nmix_data)
fit4 <- occu(formula = ~ obscovs1 + obscovs2 ~ 1,
               data = nmix_data)

cand.set <- list(P1 = fit1, P2 = fit2, P3= fit3, P4 = fit4)

library(AICcmodavg)
mods <- aictab(cand.set = cand.set, second.ord = F)
head(mods)



#5. Obtain model-averaged estimates of x1. What conclusions do you draw regarding this variable?
  
avg_obscov <- modavgShrink(cand.set = cand.set,
                             parm = 'x1',
                             second.ord = F,
                             parm.type = 'psi')

avg_obscov
  
The model averaged confidence intervals overlap 0 so I conclude this coeffiecient
has zero influence on detection probability.
  
#  6. Plot model-averaged predictions of how detection probability changes across the observed range of
#obscovs2.


nd <- data.frame(
  obscovs1 = rep(0,100),
  obscovs2 = seq(min(obs_covs$obscovs2), max(obs_covs$obscovs2), length.out = 100))


prd <- modavgPred(cand.set = cand.set, newdata= nd, second.ord = F, parm.type = "detect")
plot(x = nd$obscovs2, y = prd$mod.avg.pred,  type = 'l', 
     ylim = c(min(prd$lower.CL), max(prd$upper.CL)))
lines(x = nd$obscovs2, y = prd$upper, lty = 2)
lines(x = nd$obscovs2, y = prd$lower, lty = 2)




#7. Evaluate the fit of the top model using the sum of squared Pearson’s residuals 
#as a test statistic. A function for evaluating this test statistic is 
#provided at the bottom of the exam.
mod <- fit
chisq <- function(mod){ # mod is fitted model
  obs <- getY(mod@data) # observed
  ex <- fitted(mod) # expected
  ts <- (ex - obs) ^ 2 / # chi-square statistic
    (ex * (1 - ex))
  return(sum(ts))
}
chisq(fit3)
chisq(fit4)
chisq(fit1)





