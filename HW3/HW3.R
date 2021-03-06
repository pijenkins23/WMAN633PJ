# HW 3 Peter JEnkins 2.26.2021

dat <- read.csv('Homework 3 Data.csv')
head(dat)


#1. Fit a logistic regression model that assumes the probability of success is an additive function of variables
#x1 and x2.

fit <- glm(y ~ x1 + x2, family= binomial, data= dat)
summary(fit)

#2. Interpret the effect of variable x1 on the log odds of success. Verify your interpretation in R.

#Let’s make these quantities concrete in R. Let’s evaluate the log odds ratio
#associated with a 1-unit change in latitude. First, calculate the probability
#of presence at 1 and 2 x1 , both fixed at b:

betas <- coef(fit)
betas
summary(dat$x1)

# probability of occurrence at 1
p_1 <- plogis(betas[1] + betas[2] * 1 + betas[3] *1)
p_1

# probability of occurrence at 2
p_2 <- plogis(betas[1] + betas[2] * 2 + betas[3] * 1)
p_2

#Let’s re-create the slope coefficient from the equation in slide 44:
# log odds ratio
log((p_2 / (1 - p_2)) / (p_1 / (1 - p_1)))

betas[2]
# slope = the change in y associated with a 1-unit change in the predictor variable.


# odds 
exp(betas[2])
(p_2 / (1 - p_2)) / (p_1 / (1 - p_1)) # odds ratio by hand
 # this represents the log odds ratio of a point being in b category relative to a point being in the a category



#3. Interpret the effect of variable x2 on the log odds of success. Verify your interpretation in R.

betas <- coef(fit)
p_b <- plogis(betas[1] + betas[3]) # set x1 = 0
p_a <- plogis(betas[1])
log((p_b / (1 - p_b)) / (p_a / (1 - p_a)))
# Slope of category b 
betas[3]



#4. Duplicate the Wald Test and p-values for variables x1 and x2 performed by the glm() function. Do
#you reject or fail to reject your null hypothesis?

#########
# wald test 
betas <- coef(fit)
summary(fit)
#Let’s conduct a Wald Test for the effect of x1:
ts <-  betas[2] / summary(fit)[['coefficients']]['x1', 'Std. Error']

summary(fit)[['coefficients']]['x1', 'z value']

# p value 
2 * pnorm(-1 * abs(ts), mean= 0, sd= 1)
summary(fit)[['coefficients']]["x1",'Pr(>|z|)']

##### 
##wald test x2
#######
betas <- coef(fit)
summary(fit)
#Let’s conduct a Wald Test for the effect of x2:
ts <-  betas[3] / summary(fit)[['coefficients']]['x2b', 'Std. Error']

summary(fit)[['coefficients']]['x2b', 'z value']

# p value 
2 * pnorm(-1 * abs(ts), mean= 0, sd= 1)
summary(fit)[['coefficients']]["x2b",'Pr(>|z|)']




#5. Predict and plot the mean probability of success over the range of values of x1.


# log odds of use
y <- betas[1] + betas[2] * mean(dat$x1) + betas[3] 
# probability of use
exp(y) / (1 + exp(y))

plogis(y)

# predicting across x1 with category b 

x1 <- seq(from = min(dat$x1), to = max(dat$x1),
           length.out = 100)
y <- betas[1] + betas[2] * x1 + betas[3] 
plot(x = x1, y = plogis(y), ylab = 'Probability of occurance',
     xlab = 'x1', cex.axis = 1.5, cex.lab = 1.5, type = 'l')


