********************************************************************************
* Illustration of Multiple imputation methods for paper titled : 
*  "Multiple imputation for longitudinal and clustered data: A tutorial"       *
* 30th of June 2023                                                                *
* Author: Rushani Wijesuriya                                                   *
********************************************************************************

version 18.0

clear all
set more off

*cd "your_working_directory_path"

cd "C:\Users\rushani.wijesuriya\OneDrive - Murdoch Children's Research Institute\Multilevel MI- Tutorial paper\Longitudinal_Tutorial"


import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

/*------------------------------------------------------------------------------
                  Illustration of the reshape function
-------------------------------------------------------------------------------*/

*reshape data from wide to long format 
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)


*reshape data from long to wide format 
reshape long prev_dep numeracy_score prev_sdq, i(id) j(time)


/*------------------------------------------------------------------------------
 Illustration of MI approaches for longitudinal data 
-------------------------------------------------------------------------------*/

timer on 1

*****************1. JM-1L-wide***********

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: Change data to wide format
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* ses 
mi register regular age sex prev_sdq*

*optional step if using projected distance based rounding approach for categirical variables-
*create inidictor variables-

tabulate ses, generate(ind)


*Step 3: Imputation
set seed 2937405
mi impute mvn prev_dep3 prev_dep5 prev_dep7 numeracy_score3 numeracy_score5 numeracy_score7 numeracy_scorew1 ses = ///
	age sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) burnin(1000) saveptrace(trace, replace)
	
	
*Step 4: Asess convergence (Trace plot)
preserve
mi ptrace use trace, clear 
tsset iter
tsline b_y1x1
restore 

*Step 5: Analysis 

* (i) Round binary varaiables using adaptive rounding 
foreach var of varlist prev_dep* {
	sum `var' if _mi_m==0
	local omegabar  = `r(mean)'
	local threshold = `omegabar' - invnorm(`omegabar') * sqrt(`omegabar'*(1 - `omegabar'))
	replace `var' = 0 if `var' < `threshold' & _mi_m!=0
	replace `var' = 1 if `var' > `threshold' & _mi_m!=0
	}

* (ii) Reshape to long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(time)

* (iii) Fit  analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 age i.sex ses || id:

timer off 1


******************1. JM-1L-wide using continous calibartion for categorical variables **********


clear all


import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: Change data to wide format
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)



*Missing indicator
 gen missing=1 if ses==.
 
 tabulate ses, generate(ses)
 drop ses1
 
 * Proportion in each category
	forvalues i=2(1)5 {
		sum ses`i'
		local prop`i'=r(mean)
		}
	local prop1=1-`prop2'-`prop3'-`prop4'-`prop5'
	noi di "prop1=`prop1', prop2=`prop2', prop3=`prop3', prop4=`prop4', prop5=`prop5'"
	
	
	* Number of categories
	noi sum ses
	local ncat=r(max)
	local ncat1=`ncat'-1
	di "ncat=`ncat', ncat1=`ncat1'"
	
	
	* Duplicate data
	gen expand=2
	expand expand
	bysort id: gen dup=1 if _n==_N
	foreach var of varlist ses ses2 ses3 ses4 ses5 {
		replace `var'=. if dup==1
		}
		
		* Impute in combined dataset
	mi set flong
	mi register imputed prev_dep* numeracy_score* ses

*Step 4: Imputation
mi impute mvn prev_dep3 prev_dep5 prev_dep7 numeracy_score3 numeracy_score5 numeracy_score7 numeracy_scorew1 ses= ///
	age sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) rseed(2937405) burnin(1000)
		
*Step 5: Analysis 


*(i) continous calibration for ordinal categorical variables 
	*** Find cut-points - proportions in the observed data applied to the imputed observed data
	preserve
	keep if missing~=1 & dup==1 & _mi_m~=0
	local nobs=_N
	sort ses
	gen ses_cat=1	if _n<=`nobs'*`prop1'
	replace ses_cat=2	if _n<=`nobs'*(`prop1'+`prop2') & ses_cat==.
	replace ses_cat=3	if _n<=`nobs'*(`prop1'+`prop2'+`prop3') & ses_cat==.
	replace ses_cat=4	if _n<=`nobs'*(`prop1'+`prop2'+`prop3'+`prop4') & ses_cat==.
	replace ses_cat=5	if _n> `nobs'*(`prop1'+`prop2'+`prop3'+`prop4') & ses_cat==.
	
	
	* Take cut-point from a uniform discribution betwen values
	gen unif=uniform()
	forvalues i=1(1)`ncat1' {
		gen cut`i'_ind=1 if ses_cat==`i' & ses_cat[_n+1]==`i'+1
		gen cut`i'_val=ses+(unif*(ses-ses[_n+1])) if cut`i'_ind==1
		sum cut`i'_val
		local cut`i'=r(mean)
		di "cut`i'=`cut`i''"
	}
	if "`ncat1'"=="3" {
		local cut4=1000000
		}
	restore
	
	
	
	*** Apply to imputed data
	drop if dup==1
	gen     ses_cal=ses if _mi_m==0
	replace ses_cal=ses if _mi_m~=0 & missing~=1
	replace ses_cal=1 if ses<=`cut1' & missing==1 & _mi_m~=0
	replace ses_cal=2 if ses>`cut1' & ses<=`cut2' & missing==1 & _mi_m~=0
	replace ses_cal=3 if ses>`cut2' & ses<=`cut3' & missing==1 & _mi_m~=0
	replace ses_cal=4 if ses>`cut3' & ses<=`cut4' & missing==1 & _mi_m~=0
	replace ses_cal=5 if ses>`cut4' & ses~=. & missing==1 & _mi_m~=0
	assert ses_cal~=. if _mi_m~=0
	
	
	
	* Generate indicators
	forvalues i=2(1)5 {
		gen ses_cal`i'=1 if ses_cal==`i'
		replace ses_cal`i'=0 if ses_cal~=`i' & ses_cal~=.
		}
		
	*** Proportion in each category
	count if _mi_m~=0
	local n_tot=r(N)	
	forvalues i=1(1)5 {
		count if ses_cal==`i' & _mi_m~=0
		local n`i'=r(N)
		local p`i'_cal= `n`i''/`n_tot'
	}
	matrix p_cal = J(1,5,.)
	matrix p_cal[1,1] = `p1_cal'
	matrix p_cal[1,2] = `p2_cal'
	matrix p_cal[1,3] = `p3_cal'
	matrix p_cal[1,4] = `p4_cal'
	matrix p_cal[1,5] = `p5_cal'


* (i) Round binary varaiables using adaptive rounding 
foreach var of varlist prev_dep* {
	sum `var' if _mi_m==0
	local omegabar  = `r(mean)'
	local threshold = `omegabar' - invnorm(`omegabar') * sqrt(`omegabar'*(1 - `omegabar'))
	replace `var' = 0 if `var' < `threshold' & _mi_m!=0
	replace `var' = 1 if `var' > `threshold' & _mi_m!=0
	}

* (ii) Reshape to long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(time)

* (iii) Fit  analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 age i.sex ses_cal || id:

mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 age i.sex ses_cal2 ses_cal3 ses_cal4 ses_cal5  || id:

******************2. FCS-1L-wide***********

timer on 2

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* ses 
mi register regular age sex prev_sdq*

*Step 3: Imputation
set seed 238130
mi impute chained (regress) numeracy_score3 (regress) numeracy_score5 (regress) numeracy_score7 (regress) numeracy_scorew1 ///
	(logit) prev_dep3 (logit) prev_dep5 (logit) prev_dep7 (ologit) ses= ///
	age sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) savetrace(trace, replace)
	

* Step 4: Asess convergence (Trace plot)
preserve
use trace, clear 
keep if m==1
tsset iter
tsline numeracy_score3_mean
restore

*Step 5: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 age i.sex i.ses || id:

timer off 2

******************3. FCS-1L-MTW*************


timer on 3

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* ses 
mi register regular age sex prev_sdq*

*Step 3: Imputation
set seed 83928354

mi impute chained ///
	(regress, omit( numeracy_score7 i.prev_dep7 prev_sdq7 )) numeracy_score3 /// 
	(regress, omit(numeracy_scorew1)) numeracy_score5 /// 
	(regress, omit( numeracy_score3  i.prev_dep3 prev_sdq3 )) numeracy_score7 /// 
	(regress, omit(numeracy_score5 i.prev_dep5 prev_sdq5 )) numeracy_scorew1 /// 
	(logit, omit(numeracy_score7 i.prev_dep7 prev_sdq7 )) prev_dep3 /// 
	(logit, omit( numeracy_scorew1)) prev_dep5 /// 
	(logit, omit(numeracy_score3 i.prev_dep3 prev_sdq3 )) prev_dep7 /// 
	(ologit) ses ///
	= age sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) 

*Step 4: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 age i.sex i.ses || id:

timer off 3

/*------------------------------------------------------------------------------
 Illustration of MI approaches for longitudinal data with additional higher 
 level clustering
-------------------------------------------------------------------------------*/

timer on 4

******************1.JM-1L-DI-wide***************

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: Change data to wide format
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)

*Step 2: Create dummy indicators for schools 

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* ses 
mi register regular age sex prev_sdq* school /*add school cluster variable*/


*Step 3: Imputation
*note binary and categorical variables are imputed as continous 
set seed 848462
mi impute mvn prev_dep3 prev_dep5 prev_dep7 numeracy_score3 numeracy_score5 numeracy_score7 numeracy_scorew1 ses = ///
	age sex prev_sdq3 prev_sdq5 prev_sdq7 i.school, add(66) burnin(1000) saveptrace(trace, replace)
	
	

*Step 4: Analysis 

* (i) Round binary varaiables using adaptive rounding
foreach var of varlist prev_dep* {
	sum `var' if _mi_m==0
	local omegabar  = `r(mean)'
	local threshold = `omegabar' - invnorm(`omegabar') * sqrt(`omegabar'*(1 - `omegabar'))
	replace `var' = 0 if `var' < `threshold' & _mi_m!=0
	replace `var' = 1 if `var' > `threshold' & _mi_m!=0
	}


* (ii) Reshape to long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(time)

* (iii) Fit  analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 age i.sex ses || school: || id: 
	


timer off 4

******************2.FCS-1L-DI-wide***************

timer on 5

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* ses 
mi register regular age sex prev_sdq* school   /*add school cluster variable*/

*Step 3: Imputation
set seed 9274560
mi impute chained (regress) numeracy_score3 (regress) numeracy_score5 (regress) numeracy_score7 (regress) numeracy_scorew1 ///
	(logit) prev_dep3 (logit) prev_dep5 (logit) prev_dep7 (ologit) ses= ///
	age sex prev_sdq3 prev_sdq5 prev_sdq7 i.school, add(66) augment
	

*Step 4: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 age i.sex i.ses || school: || id: 


timer off 5
