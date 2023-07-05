###################################################################################################################
#          Illustration of Multiple imputation methods for longitudinal data with higher-level clustering for the
#          paper titled : "Multiple imputation for longitudinal data: A tutorial" 
#          08th of September 2022
#          Author:  Rushani Wijesuriya    
###################################################################################################################

rm(list=ls())

library(jomo)
library(mitml)
library(lme4)
library(mice)
library(mitools)
library(miceadds)
library('fastDummies')
library(micemd)

#set working directory
setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/Multilevel MI- Tutorial paper/Longitudinal_Tutorial")

#load the data 
CATS_dataL <- read.csv("CATS_dataL.csv", header=T)

#Reshape CATS_dataL (in long format) to wide format
CATS_dataw <-reshape(CATS_dataL,v.names=c("numeracy_score","prev_dep","prev_sdq"),timevar = "time",idvar="c_id",
                     direction= "wide") 

###-------------------------------------------Imputing in wide format-----------------------------------------------------


ST <- Sys.time()

#1. JM-1L-DI-wide

#To ensure correct function is used,set binary and categorical variables to impute as factors
cat.vars <- c("c_ses","prev_dep.3","prev_dep.5","prev_dep.7")
CATS_dataw[cat.vars]=lapply(CATS_dataw[cat.vars],factor)

######### Create imputations using jomo->jomo1->jomo1mix

##Step 1: Create a data frame with variables to be imputed
dataw_miss<-CATS_dataw[,substr(names(CATS_dataw),1,6)%in%c("prev_d","numera","c_ses")]

## Step 2: Create a data frame with complete variables to be used as predictors of the imputation model
# Note: A column of 1's must be included for the intercept
# Fully observed binary covariates(such as sex in this example)can be included as type numeric. 
# To include fully observed  categorical variables(such as the school cluster variable) with k (>2) categories, 
#(k-1) dummy variables need to be created.
dataw_comp <-  cbind(Intercept=rep(1,nrow(CATS_dataw)),
                     CATS_dataw[,substr(names(CATS_dataw),1,6)%in%c( "c_age","c_sex",
                                                                     "prev_s","school")])

#Create dummy indicators for school cluster variable
dataw_comp <- dummy_cols(dataw_comp, select_columns = 'school',remove_first_dummy = TRUE,
                         remove_selected_columns=TRUE)

## Step 3: Perform imputations (set number of imputations using nimp, burn in iterations using nburn and 
#between imputations using between options)
set.seed(37291)
imp1<-jomo(Y=dataw_miss, X=dataw_comp,nimp=66, nburn=1000, nbetween=100)
# 
# ## Step 4 : Check convergence of imputation procedure
# set.seed(37291)
# impCheck<-jomo.MCMCchain(Y=dataw_miss, X=dataw_comp,nburn=1000)
# #outputs of jomo.MCMCchain() are : 
# #(1)finimp:the final state of the data set, which would be the first imputation if we ran the jomo
# #function with nburn burn-in iterations;
# #(2)collectbeta: fixed effect parameter draws at each of the nburn iterations;
# #(3)collectomega: level-1 covariance matrix draws at each of the nburn iterations;
# #Convergence of the sampler can be assessed by looking at the trace plot for each parameter. 
# #Below we plot for beta 0,1 
# plot(c(1:1000),impCheck$collectbeta[1,1,1:1000],type="l", ylab=expression(beta["0,1"]),
#      xlab="Iteration number")
# 

######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed datasets and store in a list
imp.list <- imputationList(split(imp1,imp1$Imputation)[-1])

#Step 2: Attach the original school cluster variable and reshape each data set to long format
imp.long <- lapply(imp.list$imputations,
                   function(d) {
                     d$school <- CATS_dataw$school
                     reshape(d,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                        c("numeracy_score.3","numeracy_score.5","numeracy_score.7"),
                                         c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),
                                        idvar="id",v.names=c("prev_dep","numeracy_score","prev_sdq"), 
                                        times=c(3,5,7), direction= "long") } )

#Step 3: Fit the analysis of interest on the imputed data sets 
fit.imp1 <- lapply(imp.long,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+
                                                factor(c_ses)+(1|school/id), data=d) } )
#Step 4: Pool the results
testEstimates(fit.imp1, extra.pars =TRUE)

ET <- Sys.time()

ET-ST
#2. FCS-1L-DI-wide

#Step 1: Set categorical and binary variables to impute as factors
cat.vars <- c("c_ses","prev_dep.3","prev_dep.5","prev_dep.7","school")
CATS_dataw[cat.vars]=lapply(CATS_dataw[cat.vars],factor)

#Step 2: Set predictor matrix
pred1 <- make.predictorMatrix(CATS_dataw)
pred1[,c("c_id")]=0

#Step 3: Set imputation methods 
meth1 <- make.method(CATS_dataw)
meth1[substr(names(CATS_dataw),1,6)%in%c("prev_d")]="logreg"
meth1[substr(names(CATS_dataw),1,6)%in%c("numera")]="norm"
meth1[substr(names(CATS_dataw),1,5)%in%c("c_ses")]="polr"

#Step 4: Perform imputations
set.seed(8392)
imp2=mice(data=CATS_dataw,m=66,predictorMatrix = pred1,method=meth1,maxit=10)

#Step 5: Extract imputed data sets and store in a list
imp.list <- complete(imp2, "all")

#Step 6: Reshape all data back to long format for analysis 
imp.long <- lapply(imp.list,
                   function(d) {reshape(d,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                                        c("numeracy_score.3","numeracy_score.5",
                                                          "numeracy_score.7"),
                                                        c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),
                                        idvar="id",v.names=c("prev_dep","numeracy_score","prev_sdq"), 
                                        times=c(3,5,7), direction= "long") } )

#Step 7: Fit the analysis of interest on the imputed data sets 
fit.imp2 <- lapply(imp.long,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+ 
                                                factor(c_ses)+(1|school/c_id), data=d) } )
#Step 8: Pool the results
testEstimates(fit.imp2, extra.pars =TRUE)

#Step 4: Check convergence
plot(imp2,c("prev_dep.3","prev_dep.5","prev_dep.7"))


#3. JM-2L-wide

#To ensure correct function is used,set binary and categorical variables to impute as factors
cat.vars <- c("c_ses","prev_dep.3","prev_dep.5","prev_dep.7")
CATS_dataw[cat.vars]=lapply(CATS_dataw[cat.vars],factor)

######### Create imputations using jomo->jomo1-> jomo1ranmixhr

##Step 1: Create a data frame with variables to be imputed

#Level 1 variables (individual-specific variables)
dataw_miss<-CATS_dataw[,substr(names(CATS_dataw),1,6)%in%c("prev_d","numera","c_ses")]

#Need to create a separate level 2 variable data frame if there are incomplete cluster(school) specific
#variables

## Step 2: Create a data frame with complete variables to be used as predictors of the imputation model
# Note: A column of 1's must be included for the intercept
# Fully observed binary covariates(such as sex in this example)can be included as type numeric. 
# To include fully observed  categorical variables with k (>2) categories,(k-1) dummy variables need to be created.

# Level 1 variables
dataw_compFE <-  cbind(Intercept=rep(1,nrow(CATS_dataw)),
                     CATS_dataw[,substr(names(CATS_dataw),1,6)%in%c( "c_age","c_sex","prev_s")])


#Need to create a separate level 2 variable data frame if there were incomplete cluster(school) specific
#variables in step 1


## Step 3: Create a data frame with complete variables to be used as predictors of the imputation model, 
# with  RANDOM EFFECTS 
# Note: A column of 1's must be included for the random intercept 
# Fully observed binary covariates can be included as type numeric. 
# To include fully observed  categorical variables with k (>2) categories, (k-1) dummy variables 
#need to be created.
dataw_compRE <- cbind(Intercept=rep(1,nrow(CATS_dataw)))

## Step 4: Perform imputations (set number of imputations using nimp, burn in iterations using nburn and 
#between imputations using between options, meth="random" invokes heterogeneous covariance matrices)
set.seed(29102)
imp3 <- jomo(Y=dataw_miss, X=dataw_compFE, Z=dataw_compRE, clus=CATS_dataw$school,
             meth="random", nimp=66, nburn=1000, nbetween=100)

######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed datasets and store in a list
imp.list <- imputationList(split(imp3,imp3$Imputation)[-1])

#Step 2: Reshape all data back to long format for analysis 
imp.long <- lapply(imp.list$imputations,
                   function(d) {
                     reshape(d,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                             c("numeracy_score.3","numeracy_score.5","numeracy_score.7"),
                                             c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),
                             idvar="id",v.names=c("prev_dep","numeracy_score","prev_sdq"), 
                             times=c(3,5,7), direction= "long") } )

#Step 3: Attach the original school cluster variable and fit the analysis of interest on the imputed data sets 
fit.imp <- lapply(imp.long,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+
                                                factor(c_ses)+(1|clus/id), data=d) } )
#Step 4: Pool the results
testEstimates(fit.imp, extra.pars =TRUE)

ST=Sys.time()
#.4 FCS-2L-wide

#Set categorical and binary variables to impute as factors
cat.vars <- c("c_ses","prev_dep.3","prev_dep.5","prev_dep.7")
CATS_dataw[cat.vars]=lapply(CATS_dataw[cat.vars],factor)

#Set cluster variable to integer
CATS_dataw$school <- as.integer(CATS_dataw$school)

######### Create imputations 

#Step 1: Set predictor matrix

# In the predictor matrix, a value of
#0: Indicates that the column variable is not used as predictor for the row variable
#1: Indicates that the column variable is used as a predictor with a fixed effect for the row variable
#2: Indicates that the column variable is used as a predictor with a fixed and a random effect for 
#   the row variable
#-2: Indicates that column variable is the cluster/group variable (Only one variable is allowed)
pred2 <- make.predictorMatrix(CATS_dataw)
pred2[,c("c_id")]=0
pred2[,c("school")]=-2

#Step 2: Set imputation methods 
meth2 <- make.method(CATS_dataw)
meth2[substr(names(CATS_dataw),1,6)%in%c("prev_d")]="2l.jomo"
meth2[substr(names(CATS_dataw),1,6)%in%c("numera")]="2l.pan"
meth2[substr(names(CATS_dataw),1,5)%in%c("c_ses")]="2l.pmm"

#Step 3: Perform imputations(set number of imputations using m,predictors using predictorMatrix
#imputation method using method and burn in iterations using maxit options)
set.seed(937202)
imp4=mice(data=CATS_dataw,m=66,predictorMatrix = pred2,method=meth2,maxit=10)

######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed datasets and store in a list
imp.list <- complete(imp4, "all")

#Step 2: Reshape all data back to long format for analysis 
imp.long <- lapply(imp.list,
                   function(d) {reshape(d,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                                        c("numeracy_score.3","numeracy_score.5",
                                                          "numeracy_score.7"),
                                                        c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),
                                        idvar="id",v.names=c("prev_dep","numeracy_score","prev_sdq"), 
                                        times=c(3,5,7), direction= "long") } )

#Step 3: Fit the analysis of interest on the imputed data sets 
fit.imp4 <- lapply(imp.long,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+ 
                                                factor(c_ses)+(1|school/c_id), data=d) } )
#Step 4: Pool the results
testEstimates(fit.imp4, extra.pars =TRUE)
ET=Sys.time()
###------------------------------------------------------Imputing in long format--------------------------------------------

#5. JM-2L-DI

#To ensure correct function is used,set binary and categorical variables to impute as factors
CATS_dataL[,"c_ses"]<-as.factor(CATS_dataL[,"c_ses"])
CATS_dataL[,"prev_dep"]<-as.factor(CATS_dataL[,"prev_dep"])

######### Create imputations using jomo->jomo2->jomo2hr

##Step 1: Create a data frame with variables to be imputed

# Level 1 variables
dataL_miss1<-data.frame(prev_dep=CATS_dataL[,c("prev_dep")],numeracy_score=CATS_dataL[,c("numeracy_score")])

# Level 2 variables
dataL_miss2<-data.frame(c_ses=CATS_dataL[,c("c_ses")],numeracy_scoreW1=CATS_dataL[,c("numeracy_scoreW1")])


## Step 2: Create a data frame with complete variables to be used as predictors of the imputation model, 
# with only FIXED EFFECTS
# Note: A column of 1's must be included for the intercept
# Fully observed binary covariates(such as sex in this example)can be included as type numeric. 
# To include fully observed  categorical variables with k (>2) categories, (k-1) dummy variables 
#need to be created.

# Level 1 variables
dataL_compFE1 <-  cbind(Intercept=rep(1,nrow(CATS_dataL)),
                        CATS_dataL[,substr(names(CATS_dataL),1,6)%in%c( "c_age","c_sex","prev_s","time","school")])

#Create dummy indicators for school cluster variable
dataL_compFE1 <- dummy_cols(dataL_compFE1, select_columns = 'school',remove_first_dummy = TRUE,
                         remove_selected_columns=TRUE)
# Level 2 variables
dataL_compFE2 <-  cbind(Intercept=rep(1,nrow(CATS_dataL)),
                        CATS_dataL[,substr(names(CATS_dataL),1,6)%in%c( "c_age","c_sex","school")])

#Create dummy indicators for school cluster variable
dataL_compFE2 <- dummy_cols(dataL_compFE2, select_columns = 'school',remove_first_dummy = TRUE,
                         remove_selected_columns=TRUE)

## Step 3: Create a data frame with complete variables to be used as predictors of the imputation model, 
# with  RANDOM EFFECTS 
# Note: A column of 1's must be included for the random intercept 
# Fully observed binary covariates can be included as type numeric. 
# To include fully observed  categorical variables with k (>2) categories, (k-1) dummy variables 
#need to be created.
dataL_compRE <- cbind(Intercept=rep(1,nrow(CATS_dataL)),time=CATS_dataL[,"time"])

## Step 4: Perform imputations (set number of imputations using nimp, burn in iterations using nburn and 
#between imputations using between options, meth="random" invokes heterogeneous covariance matrices)
set.seed(398326)
imp5 <- jomo(Y=dataL_miss1, Y2=dataL_miss2, X=dataL_compFE1, X2=dataL_compFE2, Z=dataL_compRE, clus=CATS_dataL$c_id,
          nimp=66, nburn=1000, nbetween=100)



set.seed(398326)
impCheck<-jomo.MCMCchain(Y=dataL_miss1, Y2=dataL_miss2, X=dataL_compFE1, X2=dataL_compFE2, Z=dataL_compRE, clus=CATS_dataL$c_id,
                         meth="random", nburn=1000)
#outputs of jomo.MCMCchain() are : 
#(1)finimp:the final state of the data set, which would be the first imputation if we ran the jomo
#function with nburn burn-in iterations;
#(2)collectbeta: fixed effect parameter draws at each of the nburn iterations;
#(3)collectomega: level-1 covariance matrix draws at each of the nburn iterations;
#Convergence of the sampler can be assessed by looking at the trace plot for each parameter. 
#Below we plot for beta 0,1 
plot(c(1:1000),impCheck$collectbeta[1,1,1:1000],type="l", ylab=expression(beta["0,1"]),
     xlab="Iteration number")

######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed data sets and store in a list
imp.list <- imputationList(split(imp5,imp5$Imputation)[-1])
imp.list <- imp.list$imputations

#Step 2: Fit the analysis of interest on the imputed data sets 
fit.imp5 <- lapply(imp.list,function(d) {d$school <- CATS_dataL$school
                          lmer(numeracy_score ~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+
                                                factor(c_ses)+ (1|school/clus), data=d) } )
#Step 3: Pool the results
testEstimates(fit.imp5, extra.pars =TRUE)


#6. FCS-2L-DI

#Set categorical and binary variables to impute as factors
cat.vars <- c("c_ses","prev_dep","school")
CATS_dataL[cat.vars]=lapply(CATS_dataL[cat.vars],factor)

CATS_dataL$c_ses=NULL
CATS_dataL$numeracy_scoreW1=NULL

######### Create imputations 

#Step 1: Set predictor matrix

# In the predictor matrix, a value of
#0: Indicates that the column variable is not used as predictor for the row variable
#1: Indicates that the column variable is used as a predictor with a fixed effect for the row variable
#2: Indicates that the column variable is used as a predictor with a fixed and a random effect for 
#   the row variable
#-2: Indicates that column variable is the cluster/group variable (Only one variable is allowed)
#3: Indicates that the cluster means of the covariates are added as a predictor to the imputation model.

pred3 <- make.predictorMatrix(CATS_dataL)
pred3[,c("c_id")]=(-2)
pred3[c("c_ses","numeracy_scoreW1"),c("time")]=0
pred3[c("prev_dep","numeracy_score"),c("time")]=1


#Step 2: Set imputation methods
meth3 <- make.method(CATS_dataL)
meth3["prev_dep"]="2l.bin"
meth3["numeracy_score"]="2l.pan"
meth3["c_ses"]="2lonly.pmm"
meth3["numeracy_scoreW1"]="2lonly.norm"

#Step 3: Perform imputations(set number of imputations using m,predictors using predictorMatrix
#imputation method using method and burn in iterations using maxit options)
set.seed(3528)
imp5=mice(data=CATS_dataL,m=66,predictorMatrix = pred3,method=meth3,maxit=10)

imp5=mice(data=CATS_dataL,m=6,predictorMatrix = pred3,method=meth3,maxit=1)


######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed datasets and store in a list
imp.list <- complete(imp5, "all")

#Step 2: Fit the analysis of interest on the imputed data sets 
fit.imp5 <- lapply(imp.list,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+ factor(c_ses)+
                                                (1|c_id), data=d) } )
#Step 3: Pool the results
testEstimates(fit.imp5, extra.pars =TRUE)


#7. FCS-3L

#Set categorical and binary variables as numeric (0,1,2,..) 
cat.vars <- c("c_ses","prev_dep")
CATS_dataL[cat.vars]=lapply(CATS_dataL[cat.vars],as.numeric)
CATS_dataL$c_ses <- CATS_dataL$c_ses-1
CATS_dataL$prev_dep <- CATS_dataL$prev_dep-1

######### Create imputations 

# For variables imputed with ml.lmer, the setup  requires a few extra arguments.
# Specifically, we need to specify 
#(a) the level at which the higher-level variables (level 3 and 2) are measured and
#(b) the cluster variables that define the clustered structure in the imputation model of y.

#Step1 :Specify levels of higher-level variables (only relevant for variables with missing data) 
#Note: leave variables at lowest level blank (i.e., "")
variables_levels <- miceadds:::mice_imputation_create_type_vector( colnames(CATS_dataL), value="")
variables_levels["c_ses"] <- "c_id"    
variables_levels["numeracy_scoreW1"] <- "c_id" 
#if school-level variables with missing values were present:
#variables_levels["variable_name"] <- "school"    

#Step 2: Specify predictor matrix
pred4 <- mice::make.predictorMatrix(data=CATS_dataL)
pred4[ ,c("c_id", "school") ] <- 0
pred4[("numeracy_scoreW1") ,c("time") ] <- 0
pred4[("c_ses"),c("time") ] <- 0

# set -2 for cluster identifier for school level incomplete variables (if any and 2lonly function
#is used in step 3)
#predmat["variable_name" , "school" ] <- -2

#Step3: Specify imputation methods AND model 
# method for lower-level variables (time-varying and time-fixed)
meth4 <- mice::make.method(data=CATS_dataL)
meth4[c("prev_dep","numeracy_score","c_ses","numeracy_scoreW1")] <- "ml.lmer"
#If incomplete variables at school level were there, can use 2lonly. functions 

#specify models (impute c_ses and prev_dep with pmm and others with normally distributed residuals)
model <- list(prev_dep="pmm",c_ses="pmm", numeracy_score="continuous",numeracy_scoreW1="continuous")

#specify random slopes
random_slopes <- list()
random_slopes[["prev_dep"]] <- list( "c_id"=c("time"))
random_slopes[["numeracy_score"]] <- list( "c_id"=c("time"))

# Specify cluster indicators (as list)
cluster <- list()
cluster[["prev_dep"]] <-  c("c_id","school")  
cluster[["numeracy_score"]] <-  c("c_id","school")  
cluster[["c_ses"]] <-  c("school")  
cluster[["numeracy_scoreW1"]] <-  c("school")  
set.seed(73980)
imp7 <- mice(CATS_dataL, method = meth4, predictorMatrix = pred4, maxit = 10, random_slopes=random_slopes,
            m = 66, levels_id = cluster, variables_levels = variables_levels,model=model,aggregate_automatically=FALSE)

######### Analyze imputed data by fitting substantive model to each imputed data set and combine using 
#          Rubins rules

#Step 1: Extract imputed data sets and store in a list
imp.list <- complete(imp7, "all")

#Step 2: Fit the analysis of interest on the imputed data sets 
fit.imp7 <- lapply(imp.list,function(d) {lmer(numeracy_score~ prev_dep+time+c_age+numeracy_scoreW1+c_sex+ factor(c_ses)+
                                                (1|c_id)+(1|school), data=d) } )
#Step 3: Pool the results
testEstimates(fit.imp7, extra.pars =TRUE)



