********************************************************************************
* Illustration of Multiple imputation methods for paper titled :                *
* "Multiple imputation for longitudinal data: A tutorial"                       *
* Rushani Wijesuriya, Margarita Moreno-Betancur, John B Carlin, Ian R White,    *
* Matteo Quartagno and Katherine J Lee (Senior Author)                          *
*  Last updated: 30th of June 2023                                              *
*  Author responsible for the code: Rushani Wijesuriya                          *
********************************************************************************

version 17.0

clear all
set more off

*cd "your_working_directory_path"

cd "C:\Users\rushani.wijesuriya\OneDrive - Murdoch Children's Research Institute\Multilevel MI- Tutorial paper\Longitudinal_Tutorial"


import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

/*------------------------------------------------------------------------------
                  Illustration of the reshape function
-------------------------------------------------------------------------------*/

*reshape data from wide to long format 
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)


*reshape data from long to wide format 
reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(time)


/*------------------------------------------------------------------------------
 Illustration of MI approaches for longitudinal data 
-------------------------------------------------------------------------------*/

timer on 1

*****************1. JM-1L-wide***********

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: Change data to wide format
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* c_ses 
mi register regular c_age c_sex prev_sdq*

*Step 3: Imputation
set seed 2937405
mi impute mvn prev_dep3 prev_dep5 prev_dep7 numeracy_score3 numeracy_score5 numeracy_score7 numeracy_scorew1 c_ses = ///
	c_age c_sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) saveptrace(trace, replace)
	
	
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
mi reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(time)

* (iii) Fit  analysis model to each of the mi datasets
mi xtset c_id
mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 c_age i.c_sex c_ses || c_id:

timer off 1



******************2. FCS-1L-wide***********

timer on 2

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* c_ses 
mi register regular c_age c_sex prev_sdq*

*Step 3: Imputation
set seed 238130
mi impute chained (regress) numeracy_score3 (regress) numeracy_score5 (regress) numeracy_score7 (regress) numeracy_scorew1 ///
	(logit) prev_dep3 (logit) prev_dep5 (logit) prev_dep7 (ologit) c_ses= ///
	c_age c_sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) savetrace(trace, replace)
	

* Step 4: Asess convergence (Trace plot)
preserve
use trace, clear 
keep if m==1
tsset iter
tsline numeracy_score3_mean
restore

*Step 5: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset c_id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 c_age i.c_sex i.c_ses || c_id:

timer off 2
******************3. FCS-1L-MTW*************


timer on 3

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* c_ses 
mi register regular c_age c_sex prev_sdq*

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
	(ologit) c_ses ///
	= c_age c_sex prev_sdq3 prev_sdq5 prev_sdq7, add(66) 

*Step 4: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset c_id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 c_age i.c_sex i.c_ses || c_id:

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
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)

*Step 2: Create dummy indicators for schools 

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* c_ses 
mi register regular c_age c_sex prev_sdq* school /*add school cluster variable*/


*Step 3: Imputation
*note binary and categorical variables are imputed as continous 
set seed 848462
mi impute mvn prev_dep3 prev_dep5 prev_dep7 numeracy_score3 numeracy_score5 numeracy_score7 numeracy_scorew1 c_ses = ///
	c_age c_sex prev_sdq3 prev_sdq5 prev_sdq7 i.school, add(66) saveptrace(trace, replace)
	
	

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
mi reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(time)

* (iii) Fit  analysis model to each of the mi datasets
mi xtset c_id
mi estimate: xtmixed numeracy_score i.prev_dep time numeracy_scorew1 c_age i.c_sex c_ses || school: || c_id: 
	


timer off 4

******************2.FCS-1L-DI-wide***************

timer on 5

import delimited using CATS_dataL.csv,numericcols(5,6,8,9) clear

*recode binary to 0,1 
recode prev_dep (2=1) (1=0) 

*Step 1: reshape data from long to wide format 
reshape wide prev_dep numeracy_score prev_sdq, i(c_id) j(time)

*Step 2: Set up data for imputation
mi set flong
mi register imputed prev_dep* numeracy_score* c_ses 
mi register regular c_age c_sex prev_sdq* school   /*add school cluster variable*/

*Step 3: Imputation
set seed 9274560
mi impute chained (regress) numeracy_score3 (regress) numeracy_score5 (regress) numeracy_score7 (regress) numeracy_scorew1 ///
	(logit) prev_dep3 (logit) prev_dep5 (logit) prev_dep7 (ologit) c_ses= ///
	c_age c_sex prev_sdq3 prev_sdq5 prev_sdq7 i.school, add(66) augment
	

*Step 4: Analysis 

* (i) First need to reshape long
mi reshape long prev_dep numeracy_score prev_sdq, i(c_id) j(wave)


* (ii) Fit standard analysis model to each of the mi datasets
mi xtset c_id
mi estimate: xtmixed numeracy_score i.prev_dep wave numeracy_scorew1 c_age i.c_sex i.c_ses || school: || c_id: 


timer off 5
