# trying to use kernel desnsity estimate of home range musky 
# using https://jamesepaterson.github.io/jamespatersonblog/04_trackingworkshop_kernels#:~:text=Kernel%20density%20estimators%2C%20which%20map,each%20pixel%20within%20a%20grid.&text=The%20choice%20of%20kernel%20is,typically%20return%20very%20similar%20results.
# peter jenkins 2.24.2021
#install.packages('sp')
library('sp')
library('sf')
#install.packages('adehabitatHR')
library(adehabitatHR)
#install.packages("rgdal")
library('rgdal')
library('maps')
library("tidyverse")
library(dplyr)
# loading data
dat <- read.csv("CompletedatawWQ.csv")
head(dat)

hr <- dat[c("Radio.Tag",'lat', 'lon')]
head(hr)

hr <- hr[!is.na(hr$lat) & !is.na(hr$lon),]

coordinates(hr) <- c('lat','lon')
str(hr)
head(hr)


# The sample data are lat long from WGS 84
proj4string(hr) <- CRS( "EPSG:4326" )
head(hr)


# converting to utm from lat lon
res <- spTransform(hr, CRS("+proj=utm +zone=17N ellps=WGS84"))
head(res)


#kernel.ref <- kernelUD(res, h = "href" , grid=1000, extent= 5)  # href = the reference bandwidth
#image(kernel.ref$`13`) # plot

#kernel.ref[[1]]@h # The smoothing factor is stored for each animal in the "h" slot

#kernel.lscv <- kernelUD(hr, h = "LSCV") # Least square cross validation

#image(kernel.lscv) # plot
#plotLSCV(kernel.lscv) # Look for a dip
# hectares by default.
#turtle.kernel.poly.95 <- getverticeshr(kernel.ref, percent = 95) 
#print(turtle.kernel.poly)  # returns the area of each polygon

#turtle.kernel.poly.50 <- getverticeshr(kernel.ref, percent = 50) 
#print(turtle.kernel.poly)  # returns the area of each polygon

#plot(turtle.kernel.poly.95 ) 
#par(mfrow=c(1,1))

musky <- data.frame(
  "Radio.Tag" = unique(dat$Radio.Tag))

head(musky)


# using MCP instead of KDE 
# MCP 
musky.mcp <- adehabitatHR::mcp(res, percent = 95)

#Calculate the MCP by including 50 to 100 percent of points
var.per <- seq(50, 100, by = 5)
hrs <- mcp.area(res, percent = var.per)
str(hrs)


musky$MCP95 <- musky.mcp$area
musky$MCP95ID <- musky.mcp$id


head(musky)
write.csv(musky, "MuskyHomeRange.csv")


tag <- read.csv('MuskyTaggingData.csv')
head(tag)

dat <- merge(tag,musky, by = "Radio.Tag")
datc <- dat[-c(3,36),]
datc

boxplot(MCP95~ Rel..Site, data=datc,
        par(mar = c(15, 5, 4, 2)+ 0.1), 
        las = 2, xlab="", ylab="MCP Home Range (ha)",
        cex.axis = 1.5,
        cex.lab = 1.5)

boxplot(MCP95~ Recaptured..summer.angling., data=datc, main="MCP",
        xlab="Capture Site", ylab="Home Range MCP ha")


# model to test if home range changes with sex, length or weight 
fit <- glm(MCP95 ~ Sex + TL..mm. + Wt..kg., data = datc , family= Gamma)
summary(fit)

# model to check if home range size varies based on tagging location 
fit <- glm(MCP95 ~ Cap..Site , data = datc , family= Gamma)
summary(fit)



# model selection AIC for Gamma

fit1 <- glm(MCP95 ~ TL..mm., data = datc , family= Gamma)
fit2 <- glm(MCP95 ~ TL..mm. + Sex, data = datc , family= Gamma)
fit3 <- glm(MCP95 ~ Rel..Site, data = datc , family= Gamma)
fit4 <- glm(MCP95 ~ TL..mm.+ Rel..Site, data = datc , family= Gamma)
fit5 <- glm(MCP95 ~ TL..mm.+ Rel..Site + Sex, data = datc , family= Gamma)


cand.set <- list(
  P1 = fit1, 
  P2 = fit2, 
  P3 = fit3,
  P4 = fit4,
  P5 = fit5
)

#install.packages("AICcmodavg")
library(AICcmodavg)
mods <- aictab(cand.set = cand.set, second.ord = F)
mods[,1] <-  c('TL', 'TL +Sex' , 'Release Site' , "TL + Release Site", "TL + Release Site + Sex" )

head(mods)

write.csv(mods,"G:\\My Drive\\QuantEcology\\HomeRangeProject\\AICtab.csv", row.names = FALSE)

# model selection AIC for lm

fit1 <- lm(MCP95 ~ TL..mm., data = datc )
fit2 <- lm(MCP95 ~ Sex, data = datc )
fit3 <- lm(MCP95 ~ Cap..Site, data = datc )
fit4 <- lm(MCP95 ~ TL..mm.+ Cap..Site, data = datc )
fit5 <- lm(MCP95 ~ TL..mm.+ Cap..Site+ Sex, data = datc)


cand.set <- list(
  P1 = fit1, 
  P2 = fit2, 
  P3 = fit3,
  P4 = fit4,
  P5 = fit5
  
)

mods <- aictab(cand.set = cand.set, second.ord = F)
head(mods)


betas <- coef(fit3)
betas[1] + betas[2] #Jacksonville Home Range in ha
betas[1] + betas[3] # LSC Home Range in ha
betas[1] + betas[4] #  SandFork Home Range in ha
betas[1] + betas[5] #  State Park Home Range in ha
betas[1] + betas[6] # Vandalia Home Range in ha

