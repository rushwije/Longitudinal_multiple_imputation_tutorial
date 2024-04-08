###################################################################################################################
#          Simulation of example data mimicking Childhood to Adolescence Transition Study (CATS) for 
#          paper titled : "Multiple imputation for longitudinal data: A tutorial" 
#          Rushani Wijesuriya, Margarita Moreno-Betancur, John B Carlin, Ian R White,
#          Matteo Quartagno and Katherine J Lee (Senior Author) 
#          Last updated: 08th of April 2024
#          Author responsible for the code:  Rushani Wijesuriya    
###################################################################################################################

rm(list=ls())

library(ReIns) #for generating random numbers from log-normal distribution 
library(splitstackshape) #to exapand the rows
library(boot)#for inverting the logit function 
library(dplyr)
library(DataCombine) #for the slide function 
library(fastDummies) #for creating dummy indicators for the school clusters
library(mice)
library(rstudioapi)

#set working directory to working folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#set seed
set.seed(937254)


#the number of school clusters at wave 1
I<-40

#the number of students in each school at wave 1(do we fix this or vary)
J<-30

#parameters for log school frequencies 
log.mu <- 3.27
log.sigma <- 0.57
log.min <- 2.10
log.max <- 4.20


##parameters for the random effects

#exposure-depressive symptoms
sd_depL3 <- 0.25  #this was very small in the data (sd=0.0002)
sd_depL2 <- 1.50

#outcome-NAPLAN 
sd_NAPLANL3 <- sqrt(0.05)
sd_NAPLANL2 <- sqrt(0.25)
sd_NAPLANL1 <-sqrt(0.25)

#Auxiliary-SDQ
sd_sdqL3 <- 0.8
sd_sdqL2 <- 4.0
sd_sdqL1 <- 3.0


  ##generate the school clusters at wave 1
  D <- data.frame(school=c(1:I))
  
  #generate school level (Level 3) RE 
  D$L3_RE_dep <- rnorm(40,0,sd_depL3)
  D$L3_RE_NAPLAN <- rnorm(40,0,sd_NAPLANL3)
  D$L3_RE_sdq <- rnorm(40,0,sd_sdqL3)
  
  #Generate the class sizes
  D$freq <-round(exp(rtlnorm(40,log.mu,log.sigma,endpoint = log.max)))
  
  #obtain the total number of students
  total <- sum(D$freq)
  
  #recale to 1200
  N <- I*J
  D$class_size <- round(D$freq*(N/total))
  
  #add/remove students from the last cluster to total the no of students to 1200
  if ((sum(D$class_size)-N)<0){D$class_size[I] <- D$class_size[I]+(N-sum(D$class_size))} 
  if((sum(D$class_size)-N)>0){D$class_size[I] <- D$class_size[I]-(sum(D$class_size)-N)} 
  
  #replicate the rows for adding students
  D <- expandRows(D, "class_size")
  
  D$class_size <- NULL
  
  #generate the individuals
  D$c_id <- c(1:N)
  
  D <- D[order(D$school),]
  
  D$freq=NULL
  
  D <- D[order(D$c_id),]
  
  #child's age at wave 1
  D$age <- runif((I*J),7,10)
  
  #Simulate sex (M,F) groups (males=1, females=0)
  D$uran=runif((I*J),0,1)
  D$sex=ifelse(D$uran<=0.5,1,0)
  D$uran=NULL
  
  ##simulate SES values
  D$ses <- sample( 1:5,1200,replace=TRUE, prob=c(0.1, 0.1, 0.2, 0.3 ,0.3) )
  D$ses <- as.factor(D$ses)
  
  #Simulate NAPLAN scores at wave 1
  eta0 <- -1.20
  eta1 <- 0.22
  eta2 <- 0.08
  eta31 <- 0.01
  eta32 <-  0.37
  eta33 <- 0.33
  eta34 <- 0.65
 
  
  e_teacherw1 <- rnorm(N,0,1)
  
  #generate indictor variables for SES and school sector variables
  D=fastDummies::dummy_cols(D, select_columns =c("ses"))
  
  
  #generate NAPLAN scores at wave 1
  D$NAPLAN_w1=eta0+eta1*D$sex+eta2*D$age+eta31*D$ses_2 + eta32*D$ses_3 + eta33*D$ses_4 +
    eta34*D$ses_5 +e_teacherw1
  
  #generate individual level (Level 2) REs
  D$L2_RE_dep <- rnorm(N,0,sd_depL2)
  D$L2_RE_NAPLAN <- rnorm(N,0,sd_NAPLANL2)
  D$L2_RE_sdq <- rnorm(N,0,sd_sdqL2)
  
  #expand each row for repeated measures
  D <- D[rep(seq_len(nrow(D)), each = 3), ]
  D$wave <- rep(c(2,4,6))
  
  #generate the exposure (depressive symptoms at waves 2,4 and 6)
  
  #freqs- No:2652,Yes:532,No:533
  delta0 <- -4.0 #increase this to increase the prob of a negative resp
  delta1 <- 0.31
  delta2 <- -0.52
  delta3 <- -0.05
  delta4 <- 0.08
  delta51 <- -0.30
  delta52 <- -0.40
  delta53 <- -0.57
  delta54 <- -0.86
  
  uran <- runif(3600,0,1)
  D$dep <- ifelse(uran<inv.logit(delta0+delta1*D$age+delta2*D$sex+delta3*D$NAPLAN_w1+
                         delta4*D$wave+delta51*D$ses_2+delta52*D$ses_3+delta53*D$ses_4+
                           delta54*D$ses_5+D$L3_RE_dep+D$L2_RE_dep),1,0)
  
  #check
  D$dep <- as.factor(D$dep)
  summary(D$dep)

  colnames(D)[colnames(D)=="dep"] <- "prev_dep"
  
  #generate the auxiliary variable (SDQ at waves 2,4 and 6)
  beta0 <- 16.0
  beta1 <-1.60
  beta2 <- -0.1
  
  e_sdq <- rnorm(3600,0,sd_sdqL1)
  D$sdq <- beta0+beta1*(as.numeric(D$prev_dep)-1)+beta2*D$wave+D$L3_RE_sdq+ D$L2_RE_sdq+e_sdq
  
  colnames(D)[colnames(D)=="sdq"] <- "prev_sdq"
  
  D$wave <- D$wave+1
  
 
  #generate the outcome (NAPLAN scores at waves 3,5 and 7)
  gamma0 <- 2.0
  gamma1 <- -0.02
  gamma2 <- -0.20
  gamma3 <- 0.15
  gamma4 <- 0.70
  gamma51 <- -0.02
  gamma52 <- -0.10
  gamma53 <- 0.02
  gamma54 <- -0.02
  gamma6 <- -0.01
  gamma7 <- -0.01
  e_NAPLAN <- rnorm(3600,0,sd_NAPLANL1)
  D$NAPLAN <- gamma0+gamma1*(as.numeric(D$prev_dep)-1)+gamma2*D$age+gamma3*D$sex+gamma4*D$NAPLAN_w1+gamma51*D$ses_2+
    gamma52*D$ses_3+gamma53*D$ses_4+gamma54*D$ses_5+
    gamma6*D$wave+gamma7*D$prev_sdq+D$L3_RE_NAPLAN+ D$L2_RE_NAPLAN+e_NAPLAN
  
  #drop variables
  D <- D %>% dplyr::select(!contains("L3")& !contains("L2"))

  ##########MIssing data generation######################
  

  D <- D %>% group_by(school) %>% 
    mutate(L3_RE_r_dep = rnorm(1,0,0.01),L3_RE_r_NAPLAN = rnorm(1,0,0.4),L3_RE_r_sdq=rnorm(1,0,0.8) )
  
  D <- D %>% group_by(c_id) %>% 
    mutate(L2_RE_r_dep = rnorm(1,0,0.05),L2_RE_r_NAPLAN = rnorm(1,0,2.0),L2_RE_r_sdq = rnorm(1,0,1.2) )
  
  #missing data generation in prev_dep
  D$r_prev_dep <- as.numeric(runif(3600,0,1)<inv.logit(-9.8+0.72*D$age+0.16*D$sex+(-0.17)*D$NAPLAN_w1+
                                   (0.20)*(D$wave-1)+(-0.39)*D$ses_2+(0.27)*D$ses_3+(0.19)*D$ses_4+
                                   (-0.03)*D$ses_5+(-0.13)*D$NAPLAN+(0.04)*D$prev_sdq+D$L3_RE_r_dep+D$L2_RE_r_dep))
  
  D$r_prev_dep <- as.factor(D$r_prev_dep)

  #missing data generation in NAPLAN
  D$r_NAPLAN <- as.numeric(runif(3600,0,1)<inv.logit(-22.8+1.77*D$age+0.01*D$sex+(-0.70)*D$NAPLAN_w1+
                                                         (0.7)*D$wave+(-4.9)*D$ses_2+(-1.9)*D$ses_3+(2.19)*D$ses_4+
                                                         (-2.35)*D$ses_5+(-0.25)*(as.numeric(D$prev_dep)-1)+(0.11)*D$prev_sdq+
                                                       D$L3_RE_r_NAPLAN+D$L2_RE_r_NAPLAN))

  D$r_NAPLAN <- as.factor(D$r_NAPLAN)
  
  
  summary(D)
  
  #assign NA values 
  D$prev_dep <- ifelse(D$r_prev_dep==1,NA,D$prev_dep)
  #D$prev_sdq <- ifelse(D$r_prev_sdq==1,NA,D$prev_sdq)
  D$NAPLAN <- ifelse(D$r_NAPLAN==1,NA,D$NAPLAN)
  
  
  #check missingness (intermittent,drop out)
  
  D <- D %>% dplyr::select(!contains("L3")& !contains("L2") & !contains("c_ses_"))
  
  D$r_prev_dep=NULL
  D$r_NAPLAN=NULL
 # D$r_prev_sdq=NULL
  
  
  md.pattern(D)
  
  #check the drop out and intermittent missingness
  #NAPLAN
  drop_out <- aggregate(NAPLAN ~ c_id, data=D, function(x) {sum(is.na(x))}, na.action = NULL)
  
  table(drop_out$NAPLAN)
  
  #prev_dep
  drop_out2 <- aggregate(prev_dep ~ c_id, data=D, function(x) {sum(is.na(x))}, na.action = NULL)
  table(drop_out2$prev_dep)
  
  
  # #prev_sdq
  # drop_out3 <- aggregate(prev_sdq ~ c_id, data=D, function(x) {sum(is.na(x))}, na.action = NULL)
  # table(drop_out3$prev_sdq)
  # 
  
  
  D <- as.data.frame(D)
  
  D_wide <- reshape(D,v.names=c("NAPLAN","prev_sdq","prev_dep"),timevar = "wave",idvar="c_id",direction= "wide")
  
  #Missing values in NAPLAN- 15% of NAPLAN at wave 1 data missing 
  D_wide$r_NAPLANW1 <- as.numeric(runif(1200,0,1)<inv.logit(-2.1+0.05*D_wide$age+0.02*D_wide$sex))
  
  
  #Missing values in SES- 10% of SES missing 
  D_wide$r_SES <- as.numeric(runif(1200,0,1)<inv.logit(-2.5+0.03*D_wide$age+0.01*D_wide$sex))
  
  
  #assign NA values 
  D_wide$ses <- ifelse(D_wide$r_SES==1,NA,D_wide$ses)
  D_wide$NAPLAN_w1 <- ifelse(D_wide$r_NAPLANW1==1,NA,D_wide$NAPLAN_w1)
  
  
  D_wide$r_NAPLANW1=NULL
  D_wide$r_SES=NULL
  
  summary(D_wide)
  
 #reshape to long
  D <- reshape(D_wide,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                                c("NAPLAN.3","NAPLAN.5","NAPLAN.7"),
                                    c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),idvar="c_id", 
                        v.names=c("prev_dep","NAPLAN","prev_sdq"), times=c(3,5,7),direction= "long")
  
  
  D <- D[order(D$c_id),]
  
  
  #rename NAPLAN score to a generic name
  
  D <- dplyr::rename(D, numeracy_score="NAPLAN",numeracy_scoreW1="NAPLAN_w1",id="c_id")
  D <- D %>% select(!contains("es_"))
#Save the data set
  
  write.csv(D, "CATS_dataL.csv", row.names = F)
  
  
  
  

  
  
  
  