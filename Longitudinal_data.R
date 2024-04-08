###################################################################################################################
#          Illustration of Multiple imputation methods for 
#          paper titled : "Multiple imputation for longitudinal data: A tutorial" 
#          Rushani Wijesuriya, Margarita Moreno-Betancur, John B Carlin, Ian R White,
#          Matteo Quartagno and Katherine J Lee (Senior Author)
#          Last updated : 03 April 2024
#          Author responsible for the code:  Rushani Wijesuriya    
###################################################################################################################

rm(list = ls())

library(jomo)
library(mitml)
library(lme4)
library(mice)
library(mitools)
library(micemd)

# set working directory
setwd("~/Multilevel MI- Tutorial paper/Longitudinal_Tutorial")

# load the data
CATS_long <- read.csv("CATS_dataL.csv", header = T)
summary(CATS_long)

# check missing data proportions
sum(complete.cases(CATS_long))

# fit analysis models

# 2-level
lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
    factor(ses) + (1 | id), data = CATS_long)

# 3-level
lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
    factor(ses) + (1 | id), data = CATS_long)

###------------------------------------Illustration of the reshape() function--------------------------------------------

# Reshaping CATS_long (in long format) to wide format
CATS_wide <- reshape(CATS_long, v.names = c("numeracy_score", "prev_dep",
    "prev_sdq"), timevar = "time", idvar = "id", direction = "wide")

# Reshaping CATS_wide (in wide format) to long format
CATS_long <- reshape(CATS_wide, varying = list(c("prev_dep.3", "prev_dep.5",
    "prev_dep.7"), c("numeracy_score.3", "numeracy_score.5", "numeracy_score.7"),
    c("prev_sdq.3", "prev_sdq.5", "prev_sdq.7")), idvar = "id", v.names = c("prev_dep",
    "numeracy_score", "prev_sdq"), times = c(3, 5, 7), direction = "long")

###-------------------------------------------Imputing in wide format-----------------------------------------------------

#-------------------------#
# Approach 1: JM-1L-wide  #
#-------------------------#

# To ensure correct function is used,set binary and categorical
# variables to impute as factors
cat.vars <- c("ses", "prev_dep.3", "prev_dep.5", "prev_dep.7")
CATS_wide[cat.vars] <- lapply(CATS_wide[cat.vars], factor)

######### Create imputations using jomo->jomo1->jomo1mix

# Step 1: Create a data frame with variables to be imputed
dataw_inc <- CATS_wide[, substr(names(CATS_wide), 1, 6) %in% c("prev_d",
    "numera", "ses")]

# Step 2: Create a data frame with complete variables to be used as
# predictors of the imputation model Note: A column of 1's must be
# included for the intercept Fully observed binary covariates(such as
# sex in this example)can be included as type numeric.  To include
# fully observed categorical variables with k (>2) categories, (k-1)
# dummy variables need to be created.
dataw_comp <- cbind(Intercept = rep(1, nrow(CATS_wide)), CATS_wide[, substr(names(CATS_wide),
    1, 6) %in% c("age", "sex", "prev_s")])

# Step 3: Perform imputations (set number of imputations using nimp,
# burn in iterations using nburn and between imputations using
# between options)
set.seed(2946)
imp1 <- jomo(Y = dataw_inc, X = dataw_comp, nimp = 66, nburn = 1000, nbetween = 1000)

# Step 4 : Check convergence of imputation procedure
set.seed(2946)
impCheck <- jomo.MCMCchain(Y = dataw_inc, X = dataw_comp, nburn = 1000)

# outputs of jomo.MCMCchain() are : (1)finimp:the final state of the
# data set, which would be the first imputation if we ran the jomo
# function with nburn burn-in iterations; (2)collectbeta: fixed
# effect parameter draws at each of the nburn iterations;
# (3)collectomega: level-1 covariance matrix draws at each of the
# nburn iterations; Convergence of the sampler can be assessed by
# looking at the trace plot for each parameter.  Below we plot for
# beta 0,1
plot(c(1:1000), impCheck$collectbeta[1, 1, 1:1000], type = "l", ylab = expression(beta["0,1"]),
    xlab = "Iteration number")

######### Analyze imputed data by fitting substantive model to each
######### imputed data set and combine using Rubins rules

# Step 1: Extract imputed datasets and store in a list
imp.list <- imputationList(split(imp1, imp1$Imputation)[-1])

# Step 2: Reshape all data back to long format for analysis
imp.long <- lapply(imp.list$imputations, function(d) {
    reshape(d, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"),
        c("numeracy_score.3", "numeracy_score.5", "numeracy_score.7"),
        c("prev_sdq.3", "prev_sdq.5", "prev_sdq.7")), idvar = "id", v.names = c("prev_dep",
        "numeracy_score", "prev_sdq"), times = c(3, 5, 7), direction = "long")
})

# Step 3: Fit the analysis of interest on the imputed data sets
fit.imp1 <- lapply(imp.long, function(d) {
    lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
        factor(ses) + (1 | id), data = d)
})

# Step 4: Pool the results
testEstimates(fit.imp1, extra.pars = TRUE)

#----------------------------#
# Approach 2: FCS-1L-wide    #
#----------------------------#

# Set categorical and binary variables to impute as factors
cat.vars <- c("ses", "prev_dep.3", "prev_dep.5", "prev_dep.7")
CATS_wide[cat.vars] <- lapply(CATS_wide[cat.vars], factor)

######### Create imputations

# Step 1: Set predictor matrix
# In the predictor matrix, a value of 0: Indicates that the column
# variable is not used as predictor for the row variable 1: Indicates
# that the column variable is used as a predictor with a fixed effect
# for the row variable 2: Indicates that the column variable is used
# as a predictor with a fixed and a random effect for the row
# variable
#-2: Indicates that column variable is the cluster/group variable (Only one variable is allowed)
pred1 <- make.predictorMatrix(CATS_wide)
pred1[, c("id", "school")] <- 0

# Step 2: Set imputation methods (logreg,norm,and polr for specifying
# a logistic regression model,linear regression model and a
# proportional odds model respectively)
meth1 <- make.method(CATS_wide)
meth1[substr(names(CATS_wide), 1, 6) %in% c("prev_d")] <- "logreg"
meth1[substr(names(CATS_wide), 1, 6) %in% c("numera")] <- "norm"
meth1[substr(names(CATS_wide), 1, 5) %in% c("ses")] <- "polr"

# Step 3: Perform imputations(set number of imputations using
# m,predictors using predictorMatrix imputation method using method
# and burn in iterations using maxit options)
set.seed(3726)
imp2 <- mice(data = CATS_wide, m = 66, predictorMatrix = pred1, method = meth1,
    maxit = 10)

# Step 4: Check convergence
plot(imp2, c("prev_dep.3", "prev_dep.5", "prev_dep.7"))

######### Analyze imputed data by fitting substantive model to each
######### imputed data set and combine using Rubins rules

# Step 1: Extract imputed datasets and store in a list
imp.list <- complete(imp2, "all")

# Step 2: Reshape all data back to long format for analysis
imp.long <- lapply(imp.list, function(d) {
    reshape(d, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"),
        c("numeracy_score.3", "numeracy_score.5", "numeracy_score.7"),
        c("prev_sdq.3", "prev_sdq.5", "prev_sdq.7")), idvar = "id", v.names = c("prev_dep",
        "numeracy_score", "prev_sdq"), times = c(3, 5, 7), direction = "long")
})

# Step 3: Fit the analysis of interest on the imputed data sets
fit.imp2 <- lapply(imp.long, function(d) {
    lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
        factor(ses) + (1 | id), data = d)
})

# Step 4: Pool the results
testEstimates(fit.imp2, extra.pars = TRUE)

#---------------------------------#
# Approach 3: FCS-1L-wide-MTW     #
#---------------------------------#

# Set categorical and binary variables to impute as factors
cat.vars <- c("ses", "prev_dep.3", "prev_dep.5", "prev_dep.7")
CATS_wide[cat.vars] <- lapply(CATS_wide[cat.vars], factor)

######### Create imputations

# Step 1: Set predictor matrix
# In the predictor matrix, a value of 0: Indicates that the column
# variable is not used as predictor for the row variable 1: Indicates
# that the column variable is used as a predictor with a fixed effect
# for the row variable 2: Indicates that the column variable is used
# as a predictor with a fixed and a random effect for the row
# variable
#-2: Indicates that column variable is the cluster/group variable (Only one variable is allowed)
# Here for a repeated measure measured at wave/time point t, only
# those measured at (t-1) and (t+1) are used as predictors, but in
# other practical applications using a wider time window might be
# more appropriate
pred2 <- make.predictorMatrix(CATS_wide)
pred2[, c("id", "school")] <- 0
pred2["numeracy_scoreW1", c(grep("5", names(CATS_wide)), grep("7", names(CATS_wide)))] <- 0
pred2[c("prev_dep.3", "prev_sdq.3", "numeracy_score.3"), grep("7", names(CATS_wide))] <- 0
pred2[c("prev_dep.5", "prev_sdq.5", "numeracy_score.5"), "numeracy_scoreW1"] <- 0
pred2[c("prev_dep.7", "prev_sdq.7", "numeracy_score.7"), grep("3", names(CATS_wide))] <- 0
pred2[c("prev_dep.7", "prev_sdq.7", "numeracy_score.7"), "numeracy_scoreW1"] <- 0

# Step 2: Set imputation methods
meth2 <- make.method(CATS_wide)
meth2[substr(names(CATS_wide), 1, 6) %in% c("prev_d")] <- "logreg"
meth2[substr(names(CATS_wide), 1, 6) %in% c("numera")] <- "norm"
meth2[substr(names(CATS_wide), 1, 5) %in% c("c_ses")] <- "polr"

# Step 3: Perform imputations(set number of imputations using
# m,predictors using predictorMatrix imputation method using method
# and burn in iterations using maxit options)
set.seed(8920)
imp3 <- mice(data = CATS_wide, m = 66, predictorMatrix = pred2, method = meth2,
    maxit = 10)

######### Analyze imputed data by fitting substantive model to each
######### imputed data set and combine using Rubins rules

# Step 1: Extract imputed datasets and store in a list
imp.list <- complete(imp3, "all")

# Step 2: Reshape all data back to long format for analysis
imp.long <- lapply(imp.list, function(d) {
    reshape(d, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"),
        c("numeracy_score.3", "numeracy_score.5", "numeracy_score.7"),
        c("prev_sdq.3", "prev_sdq.5", "prev_sdq.7")), idvar = "id", v.names = c("prev_dep",
        "numeracy_score", "prev_sdq"), times = c(3, 5, 7), direction = "long")
})

# Step 3: Fit the analysis of interest on the imputed data sets
fit.imp3 <- lapply(imp.long, function(d) {
    lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
        factor(ses) + (1 | id), data = d)
})

# Step 4: Pool the results
testEstimates(fit.imp3, extra.pars = TRUE)

###------------------------------------------------------Imputing in long format--------------------------------------------

#-------------------------#
# Approach 4: JM-2L       #
#-------------------------#

# To ensure correct function is used,set binary and categorical
# variables to impute as factors
CATS_long[, "ses"] <- as.factor(CATS_long[, "ses"])
CATS_long[, "prev_dep"] <- as.factor(CATS_long[, "prev_dep"])

######### Create imputations using jomo->jomo2->jomo2hr

# Step 1: Create a data frame with variables to be imputed

# Level 1 variables
dataL_inc1 <- data.frame(prev_dep = CATS_long[, c("prev_dep")], numeracy_score = CATS_long[,
    c("numeracy_score")])

# Level 2 variables
dataL_inc2 <- data.frame(ses = CATS_long[, c("ses")], numeracy_scoreW1 = CATS_long[,
    c("numeracy_scoreW1")])


# Step 2: Create a data frame with complete variables to be used as
# predictors of the imputation model, with only FIXED EFFECTS Note: A
# column of 1's must be included for the intercept Fully observed
# binary covariates(such as sex in this example)can be included as
# type numeric.  To include fully observed categorical variables with
# k (>2) categories, (k-1) dummy variables need to be created.

# Level 1 variables
dataL_compFE1 <- cbind(Intercept = rep(1, nrow(CATS_long)), CATS_long[,
    substr(names(CATS_long), 1, 6) %in% c("age", "sex", "prev_s", "time")])

# Level 2 variables
dataL_compFE2 <- cbind(Intercept = rep(1, nrow(CATS_long)), CATS_long[,
    substr(names(CATS_long), 1, 6) %in% c("age", "sex")])

# Step 3: Create a data frame with complete variables to be used as
# predictors of the imputation model, with RANDOM EFFECTS Note: A
# column of 1's must be included for the random intercept Fully
# observed binary covariates can be included as type numeric.  To
# include fully observed categorical variables with k (>2)
# categories, (k-1) dummy variables need to be created.
dataL_compRE <- cbind(Intercept = rep(1, nrow(CATS_long)), time = CATS_long[,
    "time"])

## Step 4: Perform imputations (set number of imputations using nimp,
## burn in iterations using nburn and between imputations using
## between options, meth='random' invokes heterogeneous covariance
## matrices)
set.seed(7251)
imp4 <- jomo(Y = dataL_inc1, Y2 = dataL_inc2, X = dataL_compFE1, X2 = dataL_compFE2,
    Z = dataL_compRE, clus = CATS_long$id, nimp = 66, nburn = 1000, nbetween = 1000)

# #check convergence set.seed(7251)
# impCheck<-jomo.MCMCchain(Y=dataL_miss1, Y2=dataL_miss2,
# X=dataL_compFE1, X2=dataL_compFE2,
# meth='random',Z=dataL_compRE,clus=CATS_dataL$c_id,nburn=1000)
# #outputs of jomo.MCMCchain() are : #(1)finimp:the final state of
# the data set, which would be the first imputation if we ran the
# jomo #function with nburn burn-in iterations; #(2)collectbeta:
# fixed effect parameter draws at each of the nburn iterations;
# #(3)collectomega: level-1 covariance matrix draws at each of the
# nburn iterations; #Convergence of the sampler can be assessed by
# looking at the trace plot for each parameter.  #Below we plot for
# beta 0,1 plot(c(1:1000),impCheck$collectbeta[1,1,1:1000],type='l',
# ylab=expression(beta['0,1']), xlab='Iteration number')


######### Analyze imputed data by fitting substantive model to each
######### imputed data set and combine using Rubins rules

# Step 1: Extract imputed datasets and store in a list
imp.list <- imputationList(split(imp4, imp4$Imputation)[-1])
imp.list <- imp.list$imputations

# Step 2: Fit the analysis of interest on the imputed data sets
fit.imp4 <- lapply(imp.list, function(d) {
    lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
        factor(ses) + (1 | clus), data = d)
})

# Step 3: Pool the results
testEstimates(fit.imp4, extra.pars = TRUE)


#-------------------------#
# Approach 5: FCS-2L      #
#-------------------------#

# Set categorical and binary variables to impute as factors
cat.vars <- c("ses", "prev_dep")
CATS_long[cat.vars] <- lapply(CATS_long[cat.vars], factor)

######### Create imputations

# Step 1: Set predictor matrix
# In the predictor matrix, a value of 0: Indicates that the column
# variable is not used as predictor for the row variable 1: Indicates
# that the column variable is used as a predictor with a fixed effect
# for the row variable 2: Indicates that the column variable is used
# as a predictor with a fixed and a random effect for the row
# variable
#-2: Indicates that column variable is the cluster/group variable (Only one variable is allowed)
# 3: Indicates that the cluster means of the covariates are added as
# a predictor to the imputation model.
pred3 <- make.predictorMatrix(CATS_long)
pred3[, c("school", "id")] <- 0
pred3[, c("id")] <- -2
pred3[c("ses", "numeracy_scoreW1"), c("time")] <- 0
pred3[c("prev_dep", "numeracy_score"), c("time")] <- 2
pred3[c("prev_dep"), c("numeracy_score","prev_sdq")] <- 3
pred3[c("numeracy_score"), c("prev_dep","prev_sdq")] <- 3

# Step 2: Set imputation methods
meth3 <- make.method(CATS_long)
meth3["prev_dep"] <- "2l.jomo"
meth3["numeracy_score"] <- "2l.pan"
meth3["ses"] <- "2lonly.pmm"
meth3["numeracy_scoreW1"] <- "2lonly.norm"

# Step 3: Perform imputations(set number of imputations using
# m,predictors using predictorMatrix imputation method using method
# and burn in iterations using maxit options)
set.seed(3528)
imp5 <- mice(data = CATS_long, m = 66, predictorMatrix = pred3, method = meth3,
    maxit = 10)

######### Analyze imputed data by fitting substantive model to each
######### imputed data set and combine using Rubins rules

# Step 1: Extract imputed datasets and store in a list
imp.list <- complete(imp5, "all")

# Step 2: Fit the analysis of interest on the imputed data sets
fit.imp5 <- lapply(imp.list, function(d) {
    lmer(numeracy_score ~ prev_dep + time + age + numeracy_scoreW1 + sex +
        factor(ses) + (1 | id), data = d)
})

# Step 3: Pool the results
testEstimates(fit.imp5, extra.pars = TRUE)






