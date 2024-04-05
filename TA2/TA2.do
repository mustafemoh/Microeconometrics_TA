** =============================================================================
** TA Session 2 - CHAPTER 2: PANEL DATA
** Microeconometrics with Hanna Wang, IDEA, FALL 2022
** TA: Conghan Zheng

cls, clear all
cap set more off
	
cd "specify path here"	

use "TA2.dta", clear


/*  CONTENTS:

	I) Manipulating Panel Data
	
	II) Static Models
		
	III) Dynamic Models
*/
	
	
** PART I: MANIPULATING PANEL DATA ---------------------------------------------

** Reshape: wide -> long
reshape long n w k y, i(firm) j(year)
	
** Sort
sort firm year
	
** Overview
summarize
	
** Drop missing (in all variables)
drop if n==. | w==. | k==. | y==.
	
** An alternative way
* drop if missing(n, w, k, y)
	
** One more way
* foreach var of varlist n w k y { drop if missing(`var') } 

** Exclude single observations (firms that only appear once in the panel)
** Generate a variable that counts the number of times each firm appears
by firm: gen count = _N
** or
* bysort firm: gen count = _N

** Check
tab count

** Drop single obserations (in this example there are no single observations) 
drop if count == 1
drop count
	
** Exclude duplicates
** For duplicates, keep first appearances
bysort firm year: keep if _n==1
	
bysort firm year: gen order = _n


** PART II: STATIC MODELS ------------------------------------------------------


** II.0 Panel Setting

** Set the panel data structure (which allows us to use panel data operators)
xtset firm year
	
** We can further specify whether the data is yearly, quarterly, monthly...
* xtset firm year, yearly
	
* Take a look at the structure of the panel 
xtdescribe 
	
** Another way of studying patterns in panel data
* ssc install xtpatternvar
xtpatternvar, gen(pattern)
tab pattern
	
** Panel summary stats: within and between variation
xtsum n k w y
	 

** II.1 FE
	
xtreg n k w y, fe vce(cluster firm)
eststo FE

xtreg n k w y i.year, fe vce(robust)
eststo FE_twoway

** II.2 Least Squares Dummy Variables (LDSV) Estimator
	
** Generate firm dummies
xi, noomit prefix("di") gen i.firm
	
** Estimate the model (leave firm 1 as base category)
reg n k w y difirm_2-difirm_140
eststo OLS
	
** Another way of including dummies
* reg n k w y i.firm

** Comparison
esttab OLS FE, se ar2 nogaps compress mtitle keep(k w y)


** II.3 Large number of fixed effects

/* When you have a large number of individuals in your sample, and you want to 
   add individual fixed effects, the 'xtreg' command might not able to do it. 
   
   - 'areg': when fixed effects are only associated with one variable. 
   - 'reghdfe': two or more fixed effects
*/
	   
** One fixed effects: firm
areg n k w y, absorb(firm)
	
** Extract the estimated fixed effects
predict fe_firm, d
	
** Two fixed effects: firm and year
* cap ssc install reghdfe
* ssc install ftools
reghdfe n k w y, absorb(firm year)
** ... is equivalent to command
* xtreg n k w y i.year, fe
	
** Save fixed effects
reghdfe n k w y, absorb(firm year, savefe)

	
** II.4 First-Differenced Least Squares (FDLS)

/* Stata does not have any command to automatically implement FDLS. We first 
   need to generate the first differences. */
	
sort firm year
	
reg D.(n k w y), vce(cluster firm)
eststo FD
	
reg D.(n k w y), vce(cluster firm) noconstant
eststo FD_nocons
	
** Compare with previous approaches
esttab FE FD FD_nocons, se ar2 nogaps compress mtitle
	
	
** II.5 Random Effects Model (FGLS estimator)
	
xtreg n k w y, re vce(cluster firm)
eststo RE

** Compare with WG estimator
esttab FE RE, se ar2 nogaps compress mtitle


** II.6 FE or RE

** Comparison
esttab FE FE_twoway RE, keep(k w y) se mtitle compress

** Specification test: Hausman (1978) Test

xtreg n k w y, fe
eststo fe

xtreg n k w y, re
eststo re

** Notice that the order of estimations in the command matters
hausman fe re, constant sigmamore
hausman re fe, constant sigmamore

** II.7 Panel IV

** FE
xtivreg y k (n = w), fe
eststo IV_FE

** FD
xtivreg y k (n = w), fd
eststo IV_FD
	
** RE
xtivreg y k (n = w), re
eststo IV_RE
	
esttab IV_FE IV_FD IV_RE, se nogaps compress mtitle


** PART III: DYNAMIC MODELS ----------------------------------------------------

** III.1 Anderson and Hsiao (1981, 1982)
	
** Generate year dummies 
xi, noomit prefix("di") gen i.year
	
** Generate lags
sort firm year
foreach var of varlist n k w y {
	gen `var'L1 = L1.`var'
	gen `var'L2 = L2.`var'
}
	
/* An alternative way of generating lags:
   
   by firm: gen nL1 = n[_n-1]		//Lag 1
   by firm: gen nL2 = n[_n-2]		//Lag 2

   foreach var of varlist n k w y {
	  by firm: gen `var'L1 = `var'[_n-1]
	  by firm: gen `var'L2 = `var'[_n-2]
	}

*/
	
/* Instrument the lagged differenced variable (D.nL1) with the twice-lagged 
   level (nL2).
*/
ivregress 2sls D.n (D.nL1 = nL2) D.(nL2 w wL1 k kL1 kL2 y yL1 yL2 ///
          diyear_1979 diyear_1980 diyear_1981 diyear_1982 diyear_1983)
eststo AH81
				   
	
** III.2 Arellano and Bond (1991): Difference GMM

** Install command 'xtabond2'
* ssc install xtabond2 
	
/* Useful Link: https://journals.sagepub.com/doi/pdf/10.1177/1536867X0900900106

   - Option 'gmmstyle' specifies the endogeneous variables (L.n in this case)

   - Option 'ivstyle' specifies variables to serve as instruments (no need to
	 put L2.n in this option since the command already put it there)

   - Option 'noleveleq' is specified to exclude level equation, yielding
	 Difference GMM instead of System GMM.
	   
   - Option 'h' controls the form of the covariance matrix of the idiosyncratic 
     errors.
	 
   - Option 'small' requests t-statistics instead of z-statistics and an F-test 
     instead of a Wald test for overall model fit. 
*/

/* Notation L(2/.).n indicates that GMM-type instruments were created using 
   lag-2 of n from back. 
   
   For example, L(2/4).n would indicate that GMM-type instruments were created 
   using lags 2, 3, and 4 of n. 
*/
xtabond2 n L.n L2.n w L.w L(0/2).(k y) diyear_*, ///
		 gmmstyle(L.n) ///
		 ivstyle(w L.w L(0/2).(k y) diyear_*) ///
		 noleveleq robust
eststo D_GMM

** Two-step estimation
xtabond2 n L.n L2.n w L.w L(0/2).(k y) diyear_*, ///
		 gmmstyle(L.n) ///
		 ivstyle(w L.w L(0/2).(k y) diyear_*) ///
		 noleveleq twostep
eststo D_GMM_2S

** Include lags of more variables: is w and k are not strictly exogenous
xtabond2 n L.n L2.n w L.w L(0/2).(k y) diyear_*, ///
		 gmmstyle(L.(n w k)) ///
		 ivstyle(L(0/2).y diyear_*) ///
		 noleveleq robust small
	
** With all lags available (no need for your PS, and may take a while to run)
/* 
xtabond2 n L(1/7).n L(0/7).(w k y) diyear_*, ///
		 gmmstyle(L.n) ///
		 ivstyle(L(0/7).(w k y) diyear_*) ///
		 nolevel robust small 
*/

	
** III.3 Blundell and Bond (1998): System GMM
	
/* We use 'xtabond2', although this estimator can also be implemented by 
   'xtdpdsys'
*/
	
xtabond2 n L.n L2.n w L.w L(0/2).(k y) diyear_*, ///
		 gmmstyle(L.n) ///
		 ivstyle(w L.w L(0/2).(k y) diyear_*) ///
		 robust
eststo S_GMM
	
** Two-step estimator
xtabond2 n L.n L2.n w L.w L(0/2).(k y) diyear_*, ///
		 gmmstyle(L.n) ///
		 ivstyle(w L.w L(0/2).(k y) diyear_*) ///
		 robust twostep
eststo S_GMM_2S
	
** Compare
esttab D_GMM D_GMM_2S S_GMM S_GMM_2S, se ar2 nogaps compress mtitle
