## R-Script for the workshop 'Introduction to Structural Equation Modelling' ##
## Author: Moritz Link
###############################################################################

## load packages
library(piecewiseSEM)
library(lavaan)

## ----------------------------------------------------------- ##
## 1. Lets generate some random data and have a look at it 
## ----------------------------------------------------------- ##

## Generate random data
set.seed(1) # Generates identical values to the Powerpoint

data <- data.frame(x1 = rnorm(100))
data$x2 <- data$x1 + runif(100, 0, 3)
data$y1 <- data$x2 + runif(100, 0, 6)
data$y2 <- data$x2 + runif(100, 0, 9)

## Know your data, always have a look at it
data
par(mfrow=c(2,2))
apply(data, 2, hist)
par(mfrow=c(1,1))

apply(data, 2, mean)
apply(data, 2, sd)

## When we are at it, lets also check if data is normal distributed 
## ggplot function
qqnorm_plot <- function(x) {qqnorm(x);
  qqline(x)}

## get panels
par(mfrow=c(2,2))
apply(data, 2, qqnorm_plot)
par(mfrow = c(1,1))
## checking visually using qq plots
## plots the theoretical quantiles against the sample quantiles
## normal distribution gives a 1to1 slope

## if you are not sure if your plot indicates normality
library(DAAG)
qreference(data$x1, nrep = 8)
qreference(data$y2, nrep = 8)
## could you pick your data, without the color? No? GOOD!

library(MVN)
mt<- mvn(data)



# ----------------------------------------------------------- ##
## 2. Lets have a look at the variances and covariances
## ----------------------------------------------------------- ##

## get variances (diagonal entries) and covariances
cov(data)
## here you can't compare the covariances between
## x1y1 and x1y2! the values depend on the measuring units!

## get correlations
cor(data)
## here you can compare the data

## we can also scale data to 0 mean and 1 std
data.scaled <- as.data.frame(apply(data, 2, scale))
par(mfrow=c(2,2))
apply(data.scaled, 2, hist)
par(mfrow= c(1,1))
cov(data.scaled)
## now results from cov and cor are the same!
cor(data.scaled)


## build a simple linear model
## lets look at the model coefficients
# Standardized coefficients: Bxy * sd(x) / sd(y)
mod <- lm(y2 ~ y1, data)

# `coefs` returns the coefficient table (both standardized and unstandardized)
coefs(mod)
cov(data)

# `unstdCoefs` and `stdCoefs` return unrounded coefficients
stdCoefs(mod)

BetaStd <- stdCoefs(mod)$Estimate * sd(data$y1) / sd(data$y2)

# Compare manually standardized to automatically standardized output
BetaStd; stdCoefs(mod)$Std.Estimate
cor(data)

# The same as scaling the data beforehand and retrieving the raw coefficients
mod2 <- lm(y2 ~ y1, data.scaled)

stdCoefs(mod2) # Estimte and Std.Estimate are the same

# ----------------------------------------------------------- ##
## 3. Examples for 8 rules of path coefficients
## ----------------------------------------------------------- ##

### RULE 1: PATH COEFFICIENTS FOR UNANALYZED RELATIONSHIPS (BIDIRECTIONAL ARROWS) AMONG
### EXOGENOUS VARIABLES ARE THEIR BIVARIATE CORRELATION

cor(data[, -4]) ## check slide 22

### RULE 2: WHEN VARIABLES ARE CONNECTED BY A SINGLE PATH, THE (STANDARDIZED) PATH COEFFICIENT IS
### SIMPLY THE CORRELATION COEFFICIENT

cor(data[, -2])

# Path 1
mody1.x1 <- lm(y1 ~ x1, data)

stdCoefs(mody1.x1)$Std.Estimate; cor(data[, c("y1", "x1")])[2, 1]

# Path 2
mody2.y1 <- lm(y2 ~ y1, data)

stdCoefs(mody2.y1)$Std.Estimate; cor(data[, c("y2", "y1")])[2, 1]
## Slide 22

### RULE 3: STRENGTH OF A COMPOOUND PATH IS THE PRODUCT OF THE (STANDARDIZED) COEFFICIENTS ALONG
### THE PATH

stdCoefs(mody1.x1)$Std.Estimate * stdCoefs(mody2.y1)$Std.Estimate; cor(data[, c("y2", "x1")])[2, 1]
## Slide 23

# Wait a minute...!

### RULE 4: WHEN VARIABLES ARE CONNECTED BY MORE THAN ONE CAUSAL PATHWAY, THE PATH COEFFICIENTS
### ARE "PARTIAL" REGRESSION COEFFICIENTS

# Path 1
mody2.x1 <- lm(y2 ~ y1 + x1, data)

## calculation looks like this: 
Gamma.y2.x1 <- (cor(data$y2, data$x1) - (cor(data$y2, data$y1) * cor(data$y1, data$x1))) /
  (1 - cor(data$y1, data$x1) ^ 2)

stdCoefs(mody2.x1)[2, 8]; Gamma.y2.x1
## slide 24

# Path 2
Gamma.y2.y1 <- (cor(data$y2, data$y1) - (cor(data$y2, data$x1) * cor(data$y1, data$x1))) /
  (1 - cor(data$y1, data$x1) ^ 2)

stdCoefs(mody2.x1)[1, 8]; Gamma.y2.y1
## Slide 24


### RULE 5: PATHS FROM ERROR VARIABLES REPRESENT PREDICTION ERROR

# Path 1
(Zeta.y1 <- 1 - summary(mody1.x1)$r.squared)

sqrt(Zeta.y1)
## slide 25

# Path 2
(Zeta.y2 <- 1 - summary(mody2.x1)$r.squared)

sqrt(Zeta.y2)
## slide 25


### RULE 6: UNANALYZED RESIDUAL CORRELATIONS BETWEEN ENDOGENOUS VARIABLES ARE PARTIAL CORRELATIONS
mody2.x1 <- lm(y2 ~ x1, data)
mody1.x1 <- lm(y1 ~ x1, data)

# Get residuals from models and compute correlation
(pcor <- cor(
  # effect of x1 on y2
  resid(mody2.x1),
  # effect of x1 on y1
  resid(mody1.x1)
))

# Can also use function from piecewiseSEM
partialCorr(y1 %~~% y2, list(mody2.x1, mody1.x1))

# Get residual correlation

(Zeta.y2 = sqrt(1 - summary(mody2.x1)$r.squared))


(Zeta.y1 = sqrt(1 - summary(mody1.x1)$r.squared))

# Total correlation between y1 and y2 equals sum of direct and indirect paths (Also Rule 8)
stdCoefs(mody2.x1)$Std.Estimate * stdCoefs(mody1.x1)$Std.Estimate +
  Zeta.y2 *
  pcor *
  Zeta.y1; cor(data$y1, data$y2)
## Slide 26


### RULE 7: TOTAL EFFECT ONE VARIABLE HAS ON ANOTHER IS THE SUM OF ITS DIRECT AND INDIRECT EFFECTS

mody2 <- lm(y2 ~ x1 + y1, data)

mody1 <- lm(y1 ~ x1, data)

Gamma.y2.x1 <- stdCoefs(mody2)$Std.Estimate[1] + (stdCoefs(mody1)$Std.Estimate[1] * stdCoefs(mody2)$Std.Estimate[2])

Gamma.y2.x1
## Slide 27

### RULE 8: SUM OF ALL PATHWAYS BETWEEN TWO VARIABLES (DIRECT AND INDIRECT) EQUALS THE CORRELATION

Gamma.y2.x1; cor(data$y2, data$x1)

# Suppression effect: when the estimated coefficients deviates strongly from its correlation
mody1.1 <- lm(y1 ~ x1 + x2, data)

stdCoefs(mody1.1)$Std.Estimate[2]; cor(data$y1, data$x2)
## slide 27



# ----------------------------------------------------------- ##
## 4. Into into lavaan
## ----------------------------------------------------------- ##

#Load your Data File
keeley<-read.csv("/home/moritz/Nextcloud/PHD/Projects/Romania/SEM_Course/Dropbox/data/Keeley_rawdata_select4.csv")
head(keeley)

#regression syntax versus lavaan
aLM<-lm(cover ~ age, data=keeley)

aSEM<-sem('cover ~ age', data=keeley)

#output from both
summary(aSEM)## yeahh it converged!, very important
## no negative degrees of freedom!, model is identified!? Look this up
## gives no intercepts, thats normal because it is based on variance covariance matrices

summary(aLM)

#we want out intercept!
aMeanSEM<-sem('cover ~ age', data=keeley, meanstructure=T)
## some functions require the intercepts!
summary(aMeanSEM)


#standardized coefficients
standardizedSolution(aSEM)

summary(aSEM, standardized=T, rsq=T)

####################################
# SEM with multiple paths ####
####################################

partialMedModel<-' firesev ~ age
cover ~ firesev + age'
## feed it a string
## one relationship per line!
## create the string object for the model

## now put that string in the function
partialMedSEM<-sem(partialMedModel, 
                   data=keeley)

#look at coefficients
summary(partialMedSEM, rsquare=T, standardized=T)
## in a paper also present somewhere the whole table of results
## dropping paths from the model depends on your ultimate goal!
## it does influence other paths down the line, it is better to 
## design different models ahead and test/compare them
#################################
#Direct and Indirect effects ####
#################################

partialMedModelInd <-'

  #model
  firesev ~ af*age
  cover ~ fc*firesev + ac*age

  #Derived Calcuations
  indirect := af*fc 
  total := ac + (af*fc)
'

## := calculates the indirect and the total effect of age on cover
partialMedSEMInd<-sem(partialMedModelInd, 
                      data=keeley)

summary(partialMedSEMInd)

standardizedSolution(partialMedSEMInd)

# ----------------------------------------------------------- ##
## 5. Exercise
## ----------------------------------------------------------- ##

##The Richness Partial Mediation Model
distModel <- 'rich ~ distance + abiotic + hetero
hetero ~ distance
abiotic ~ distance'

distFit <- sem(distModel, data=keeley)

summary(distFit, std=T, rsquare=T)
standardizedSolution(distFit)

distModelEff <- '
rich ~ dr*distance + ar*abiotic + hr*hetero
hetero ~ dh*distance
abiotic ~ da*distance

#The effects
direct := dr
indirect := dh*hr + da*ar
total := direct + indirect
'

distFitEff <- sem(distModelEff, data=keeley)

standardizedSolution(distFitEff)
