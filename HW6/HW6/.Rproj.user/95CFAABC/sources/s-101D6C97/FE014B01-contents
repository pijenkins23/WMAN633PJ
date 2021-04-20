#1. Load data and place into an unmarkedFramePCount object

y <- read.csv('count.csv')
head(y)

library(unmarked)
# converting count data into a matrix
# matrix format required by unmarkedFramePCount
sosp_mat <- as.matrix(y)
nmix_data <- unmarkedFramePCount(y = sosp_mat)





obs_covs <- read.csv("obs_covs.csv")
head(obs_covs)

det_covs <- list(
  rep = obs_covs)


site_covs <- read.csv('site_covs.csv')
head(site_covs)
site_covs$x2 <- as.factor(site_covs$x2)

nmix_data <- unmarkedFramePCount(y = sosp_mat, # counts
                                 siteCovs = site_covs,# site cov
                                 obsCovs = det_covs  # detection covariates
                                 )



#2. Fit an N-mixture model that assumes conditional detection probability is a function of the detection
#covariate provided, and expected abundance is a additive function of variables x1 and x2.

fit <- pcount(formula = ~ rep ~ x1 + x2 ,
              data = nmix_data, K = 100)
fit



#3. Interpret the effect of x1 on the expected count at each site. Verity your interpretation in R.
summary(fit)
betas <- coef(fit)

(x1.0a <- exp(betas[1] + betas[2] * 0 ))
(x1.1a <- exp(betas[1] + betas[2] * 1 ))

newdata <- data.frame(rep = rep(0.005, 100),
                     x1 = seq(min(site_covs$x1), max(site_covs$x1), length.out = 100) ,
                     x2 = factor(c('a'),
                                 levels = c('a', 'b','c','d'))
)

prda <- predict(object = fit, newdata = newdata,
               
               type = 'state')


par(mfrow=c(2,2))

rep = seq(min(-2.4), max(3.38), length.out = 100)
plot(rep,prda$Predicted, xlab = "x1", ylab = "Predicted Count", main = "a")

p <- ggplot(prda, aes(rep, Predicted)) +
  geom_point() +
  xlab('x1') +# for the x axis label
  ylab('predicted count')

# 3. Add prediction intervals
p + geom_line(aes(y = lower), color = "red", linetype = "dashed")+
  geom_line(aes(y = upper), color = "red", linetype = "dashed")
###############b 
newdata <- data.frame(rep = rep(0.005, 100),
                      x1 = seq(min(site_covs$x1), max(site_covs$x1), length.out = 100) ,
                      x2 = factor(c('b'),
                                  levels = c('a', 'b','c','d'))
)

prdb <- predict(object = fit, newdata = newdata,
                
                type = 'state')



rep = seq(min(-2.4), max(3.38), length.out = 100)
plot(rep,prdb$Predicted, xlab = "x1", ylab = "Predicted Count", main = "b")
######################## c
newdata <- data.frame(rep = rep(0.005, 100),
                      x1 = seq(min(site_covs$x1), max(site_covs$x1), length.out = 100) ,
                      x2 = factor(c('c'),
                                  levels = c('a', 'b','c','d'))
)

prdc <- predict(object = fit, newdata = newdata,
                
                type = 'state')



rep = seq(min(-2.4), max(3.38), length.out = 100)
plot(rep,prdc$Predicted, xlab = "x1", ylab = "Predicted Count", main = "c")

######################## d
newdata <- data.frame(rep = rep(0.005, 100),
                      x1 = seq(min(site_covs$x1), max(site_covs$x1), length.out = 100) ,
                      x2 = factor(c('d'),
                                  levels = c('a', 'b','c','d'))
)

prdd <- predict(object = fit, newdata = newdata,
                
                type = 'state')


rep = seq(min(-2.4), max(3.38), length.out = 100)
plot(rep,prdd$Predicted, xlab = "x1", ylab = "Predicted Count", main = "d")




#The expected count at each site increases by 1.1 for each 1 unit increase in x1 for site a. 



#4. Predict and plot the effect of the supplied detection covariate. Do this over the range of this covariate.
apply(X = det_covs$rep, MARGIN = 2, FUN = max, na.rm = TRUE) # determing max number of people observed 
apply(X = det_covs$rep, MARGIN = 2, FUN = min, na.rm = TRUE) # determing max number of people observed 

newdat <- data.frame(rep = seq(min(-2.4), max(3.38), length.out = 100),
                     x1 = mean(site_covs$x1),
                     x2 = factor(c('a', 'b','c','d'),
                                 levels = c('a', 'b','c','d'))
)

prd <- predict(object = fit, newdata = newdat,
               
               type = 'state')

prd

library(ggplot2)

rep = seq(min(-2.4), max(3.38), length.out = 100)
par(mfrow = (c(1,1)))
plot(x = c(0, 0), y = prd[1, c('lower', 'upper')],
     ylim = c(0, prd[4, 'upper']), 
      type = 'l', xlim = c(-0.5, 3.5), xaxt = "n",
     xlab = '', ylab = 'Expected count',
     cex.axis = 1.5, cex.lab = 1.5, lwd = 3)
lines(x = c(1, 1), y = prd[2, c('lower', 'upper')], lwd = 3)
lines(x = c(2, 2), y = prd[3, c('lower', 'upper')], lwd = 3)
lines(x = c(3, 3), y = prd[4, c('lower', 'upper')], lwd = 3)
points(x = c(0, 1,2,3), y = prd[, 'Predicted'], pch = 16,
       cex = 3)
axis(side = 1, at = c(0, 1,2,3),
     labels = c('a', 'b','c','d'),
     cex.axis = 1.5)



#5. Use contrasts to compare expected abundance between all pairwise 
#levels of variable x2. Obtain p-values
#associated with each contrast and tell me whether you reject or fail to reject each null hypothesis tested.


x <- matrix(
  c(0, 0, 1, -1, 0,
    0, 0, 1, 0, -1,
    0, 0, 0, 1, -1),
  nrow = 3, byrow = T
)
x

lin_com <- linearComb(obj = fit, coefficients = x,
                      
                      type = 'state')

lin_com

# wald test
w <- coef(lin_com) / SE(lin_com)
w

#Calculating p-values
2 * pnorm(-1 * abs(w))


We reject each of the null hypothesis that there is no difference between
each of the groups, since the pvalues are all below the alpha of 0.05. 