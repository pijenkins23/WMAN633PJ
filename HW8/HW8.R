
y <- read.csv('sosp_nmix (2).csv')
head(y)

library(unmarked)
# converting count data into a matrix
# matrix format required by unmarkedFramePCount
sosp_mat <- as.matrix(y)
nmix_data <- unmarkedFramePCount(y = sosp_mat)


# Covariate s
##################

p_covs <- read.csv('p_covs_nmix (2).csv')
head(det_covs)

det_covs <- list(
  time = data.frame(p_covs[, c('time.1', 'time.2')]),
  sky = data.frame(sky.1 = factor(p_covs$sky.1),
                   sky.2 = factor(p_covs$sky.2))
  
)


site_covs <- read.csv('n_covs_nmix (2).csv')
head(site_covs)

site_covs$type <- factor(site_covs$type)


nmix_data <- unmarkedFramePCount(y = sosp_mat, # detection / non-detection
                                 siteCovs = site_covs, # site-level covs
                                 obsCovs = det_covs) 










#Fit an N-mixture model that assumes abundance is a function of wetland size and type, and detection
#probability is a function of sky and time (5 points).


fit <- pcount(formula = ~ time + sky ~ size + type,
              data = nmix_data, K = 100)
fit







#2. Write a function that calculates the sum of squared Pearson residuals from a fitted model. This test
#statistic can be calculated as:
mod <- fit  

pear <- function(mod){ # mod is fitted model
  obs <- getY(mod@data) # observed count
   ex <- fitted(mod) #obs prob * expected count
   #exabun <- predict(mod, type = "state")
   det <- predict(mod, type = "det")
   detm <- matrix(det$Predicted, byrow = TRUE, ncol = 2) # 
  ts <- (obs - (ex)) ^ 2 / # chi-square statistic
    ((ex) * (1 - detm))
  return(sum(ts))
}
pear(fit)  # 1025.856 expected value

 exP <- exabun$Predicted*det$Predicted


#3. Use the parboot() function in R to simulate the distribution of this test statistic under the assumption
#that your fitted model is the data-generating model. Simulate 1000 values of the test statistic. Note
#that this may take several minutes (5 points).

sims <- parboot(fit, pear, nsim = 1000)
sims


#4. Plot the distribution of the simulated test statistic. Include in this plot the value of your test statistic
#calculated from your fitted model. 
#What is the null hypothesis you are testing when conducting model checking?
#Do you reject or fail to reject this null hypothesis? 
#What are the implications for how well you model fits the data (5 points)?
  
#plot(sims)

hist(sims@t.star[, 1], xlab = 'Pearson Residuals',
       main = 'distribution of test statistic',
       cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, xlim= c(0,1000))
lines(x = rep(pear(fit), 2),
      y = c(0, 1000),
      col = 'red', lwd = 3)

sum(sims@t.star[, 1] > pear(fit)) / 1000

# the null hypothesis is that the test statistic could have reasonably arrisen
# from the theoretical distribution associated with the null hypothesis. 

# we reject the null hypothesis given the extreme test statistic and the pvalue of 0
# Our test statistic does not fall within the distribuition of the generated data, 
meaning that our model does not fit well with our data. 
