** =============================================================================
** TA Session 3 - CHAPTER 3: DISCRETE CHOICE
** Microeconometrics with Hanna Wang, IDEA, FALL 2022
** TA: Conghan Zheng

cls, clear all
set seed 13
cap set more off

cd "..."

** Packages to be installed ----
** - estout: Tools for making regression tables
* cap ssc install estout
** - margeff: Partial effects after estimation
* cap ssc install margeff
** - fitstat: Compute fit statistics for single equation regression models
* cap ssc install fitstat 


/*  Contents
	
I)   Binary Outcome Models

II)  Multinomial Models

III) Remaining topics:
	 - Endogenous Variables
	 - Binary Models in Panel Data	
		
Inputs:
	- TA3_1.dta: I.1-I.6, II.1, II.5-II.6, III.1
	- TA3_2.dta: II.2-II.4
	- TA3_3.dta: III.2
	
Outputs:
    - Graphs
*/
	
	
** PART I: BINARY OUTCOME MODELS -----------------------------------------------
	
use "TA3_1.dta", clear
	
describe 
summarize


** I.0 LPM

reg lfp age age2 married educ black nchild citiz

		
** I.1 Logit

** partical effects
logit lfp age age2 married educ black nchild citiz

** odds ratio
logit lfp age age2 married educ black nchild citiz, or


** I.2 Marginal effects

** I.2.1 At-means marginal effects (MEM):

quietly logit lfp age age2 married educ black nchild citiz
margins, dydx(*) atmeans

** Stata version <= 11
quietly logit lfp age age2 married educ black nchild citiz
mfx

	
** I.2.2 Marginal effects at a representative value (MER):

quietly logit lfp age age2 married educ black nchild citiz
margins, dydx(*) at(age=20 age2=400 married=1 educ=73 black=1 nchild=2 citiz=1)
	
** Stata version <= 11
quietly logit lfp age age2 married educ black nchild citiz
mfx, at(20 400 1 73 1 2 1)
	
	
** I.2.3 Average marginal effect (AME):

quietly logit lfp age age2 married educ black nchild citiz
margins, dydx(*)

** Stata version <= 11
quietly logit lfp age age2 married educ black nchild citiz
margeff

** I.2.4 Example: margins by education 

** At-means marginal effects, factorize educ
logit lfp age age2 i.married i.educ black nchild citiz
margins educ, atmeans
	
** Generate another variable to simplify the education measure
gen educ_1 = 0
replace educ_1 = 1 if educ == 73 
replace educ_1 = 2 if educ == 81
replace educ_1 = 3 if educ>=91 & educ <= 111
replace educ_1 = 4 if educ > 111

label define labeleduc 0 "Less than high-school" ///
	1 "High-school diploma or equivalent" ///
	2 "Some college but not degree" ///
	3 "College degree" ///
	4 "Master, professional school, doctorate degree", replace
label values educ_1 labeleduc

** At-means marginal effects, new education measure: educ_1
logit lfp age age2 i.married i.educ_1 black nchild citiz
margins educ_1, atmeans
		
/* Margins by education level:

   - Option 'vsquish' specifies that blank spaces must be suppressed. 
	
   - Option 'post' causes 'margins' to save the vector of estimated margins 
     along with the  estimated variance-covariance matrix to 'e()', so that you 
	 can use them later. 
	   
   Check the stata manual of 'margins' for more info. 
*/
margins, at(educ_1=(2(1)4)) atmeans vsquish post
		
** Plot the margins
quietly {
	logit lfp age age2 i.married i.educ_1 black nchild citiz
	margins educ_1, atmeans
	}
marginsplot
marginsplot, recast(line) recastci(rarea) nolabels // Plot confidence intervals
graph export "graph_edu_margins", as(png) replace 


** I.3 Goodness of fit

/* The classification is 'correct' if the predicted y_hat is above the threshold 
   (0.5 by default) and the actual y = 1 or if y_hat is below the threshold and 
   y = 0. 
	   
   'Sensitivity' is the fraction of y = 1 observations that are correctly 
   specified. 
   
   'Specificity' is the percentage of y= 0 observations that are correctly 
   classified. 
*/
	   
quietly logit lfp $x_list
estat classification
	
quietly probit lfp $x_list
estat classification
	
** Comparing fitted probabilities
/* Option 'pr' predicts the probability of a positive outcome (y = 1) using the 
   estimates. 
*/
quietly logit lfp $x_list
predict plogit, pr

quietly probit lfp $x_list
predict pprobit, pr

quietly reg lfp $x_list
predict pols, xb

summarize lfp plogit pprobit pols

** Plotting predicted probabilities against observed choices by age
	
sort age
	
egen m_lfp = mean(lfp), by(age)
egen m_logit = mean(plogit), by(age)
egen m_probit = mean(pprobit), by(age)
egen m_ols = mean(pols), by(age)

graph twoway (line m_lfp age, lcolor(blue)) ///
	  (line m_logit age, lcolor(red) lpattern (shortdash)) ///
	  (line m_probit age, lcolor(green) lpattern (dash_dot)) ///
      (line m_ols age, lcolor(magenta) lpattern(longdash)), ///
	  scale(1) plotregion(style(none)) ///
	  title("Actual vs. Predicted Choice Probabilities across Models", color(black)) ///
	  xtitle("Age") ///
	  ylabel(,angle(horizontal)) ///
	  subtitle("Choice probability, p(lfp = 1)", position(11) justification(left)) ///
	  ytitle("") ///
	  legend(pos(5) ring(0) col(1)) legend(size(small)) ///
	  legend(label(1 "Actual Data") label(2 "Logit")) ///
	  legend(label(3 "Probit") label(4 "OLS")) 
graph export "graph_effects_by_age", as(png) replace 


** I.4 Model specification tests
	
** I.4.1 Wald test:
	
/* Let's test for the presence interaction between marriage and age.
   
   Test whether the coefficients of the interaction terms. Only the estimation 
   of the restricted model (with interaction terms) is required.
   
   H0: coeffients are jointly 0.
*/
	   
gen age_married = age*married
gen age2_married = age2*married
	
global age_eff age_married age2_married
	
quietly logit lfp age age2 married i.educ_1 black nchild citiz $age_eff
test $age_eff
   
   
** I.4.2 Likelihood-ratio test:
	
/* Equivalent to the Wald test if the model is correctly specified. 
   
   Implementation requires to estimate both the general and the restricted 
   models (with and without interactions).
   
   H0: L = L(+)
*/
	   
global x_list age age2 married i.educ_1 black nchild citiz
	
quietly logit lfp $x_list $age_eff
estimates store full
	
quietly logit lfp $x_list
lrtest full

	
** I.4.3 Lagrange multiplier test:

/* Let's test whether adding the generated regressor (x'betahat)^2 improves the 
   fitness of the model.
   
   H0: coefficient of the new regressor is 0
*/
	   
quietly logit lfp $x_list
predict xbhat, xb
	
gen xbhat2 = xbhat^2
	
quietly logit lfp $x_list xbhat2
test xbhat2


** I.5 Maximum likelihood estimation programs
	
capture program drop lf_logit
	program lf_logit
	  version 15
	  args lnf xb
	  local y "$ML_y1"
	  quietly replace `lnf' = ln(  invlogit(`xb')) if `y'==1
	  quietly replace `lnf' = ln(1-invlogit(`xb')) if `y'==0
end	

ml model lf lf_logit (lfp = age age2 married educ black nchild citiz)
ml maximize

** Compare manually MLE with the logit command (supposed to deliver the same results)

quietly {
	logit lfp age age2 married educ black nchild citiz
	estimates store logit_res

	ml model lf lf_logit (lfp = age age2 married educ black nchild citiz)
	estimates store logit_ml
	}
estimates table logit_res logit_ml, star stats(N r2 r2_a)
	

** I.6 Probit: 

probit lfp age age2 married educ black nchild citiz

** Marginal effects
quietly probit lfp age age2 i.married educ black nchild citiz
margins married, atmeans
	
/* Measures of fit (also works with other discrete choice models)
   
   Package 'fitstat' is required. 
*/ 
fitstat

/* Heteroskedastic probit
	
   Option 'het' tells Stata which regressor to use for modelling the variance.
*/
hetprob lfp $x_list, het(black) nolog
	
	
** I.7 Compare logit, probit and OLS models
	
** A comparison of estimates
quietly eststo logit : logit lfp $x_list
quietly eststo logit_r: logit lfp $x_list, vce(robust)
quietly eststo probit : probit lfp $x_list
quietly eststo probit_r: probit lfp $x_list, vce(robust)
quietly eststo ols : reg lfp $x_list
quietly eststo ols_r: reg lfp $x_list, vce(robust)

esttab logit logit_r probit probit_r ols ols_r, ///
       b se mtitle("Logit" "Logit r" "Probit" "Probit r" "OLS" "OLS r")
		

** PART II: MULTINOMIAL MODELS -------------------------------------------------
	
** Sector variables

** 0 = not participating 
gen sector = 0 if lfp==0
** 1 = self-employed
replace sector = 1 if inlist(classwkr,13,14)
** 2 = private sector employee
replace sector = 2 if inlist(classwkr,22,23)
	
label define labelsec 0 "Not participating" 1 "Self-employed" 2 "Private sector employee"
label values sector labelsec
	
** keep if 'sector' is not missing
keep if sector!=.
	

** II.1 Multinomial logit

** All regressors are alternative-invariant
mlogit sector age age2 married i.educ_1 black nchild citiz, nolog
	
** Specify the base category ('not participating')
mlogit sector age age2 married i.educ_1 black nchild citiz, baseoutcome(0) nolog
	
** Report estimates in terms of relative-risk ratios (RRR)
mlogit sector age age2 married i.educ_1 black nchild citiz, rrr baseoutcome(0) nolog
		
** Predicted probabilities 	
xi, noomit prefix("i_") gen i.sector 
	
quietly mlogit sector age age2 married i.educ_1 black nchild citiz, baseoutcome(0) nolog
predict pmlogit1 pmlogit2 pmlogit3, pr
	
summarize pmlogit* i_sector_0-i_sector_2, separator(3)
	
** Marginal effects
quietly mlogit sector age age2 married i.educ_1 black nchild citiz, baseoutcome(0) nolog
	
margins, dydx(*) predict(outcome(0))
margins, dydx(*) predict(outcome(1))
margins, dydx(*) predict(outcome(2))

	
** II.2 Conditional logit

use "TA3_2.dta", clear

/* Now we use data on individual choice of whether to fish using one of four 
   possible alternatives: 
   - from the beach, 
   - the pier, 
   - a private boat, or 
   - a charter boat. 
   
   One explanatory variable is case specific (income) and all the others are 
   alternative specific.
*/
	   
/* Reshape tha data: wide -> long. 
   
   The parameters of conditional logit models are estimated with commands that 
   require the data to be in long form, with one observation for one 
   alternative*individual (including alternatives not chose by the individual). 
*/

** Obtain four observations for each individual 
reshape long d p c, i(id) j(fishmode beach pier boat charter) string
	
label var d "outcome"
label var p "price"
label var c "catch rate"
	
list in 1/8, table sepby(id) 
	
drop mode // [mode] is misleading after reshaping and we are not going to use it

/* Conditional logit with alternative-varying regressors price 'p' and catch 
   rate 'c'. 
	   
   'case' specifies the numeric variable that identifies each case (each 
   individual).
	   
   'alternatives' specifies the variable that identifies alternatives for each 
   case/individual. 
*/
asclogit d p c, case(id) alternatives(fishmode) noconstant

** Marginal effects
estat mfx, varlist(p)
estat mfx, varlist(c)

** A more genral condtional logit with alternative-invariant regressor income. 
asclogit d p c, case(id) alternatives(fishmode) casevars(income)
	 
** Marginal effects
estat mfx, varlist(p)
estat mfx, varlist(c)
estat mfx, varlist(income)

	 
/* The conditional logit model can also be fitted by the 'clogit' command, which
   does not have an option for case-specific variables (yes, designed for pure 
   conditional logit, not for mixed). 
   
   Instead, a case-specific variable can be added to the regression using 
   dummies for m-1 alternatives, and these m-1 variables enter as regressors. 
*/

	   
** II.3 Nested logit: 

/* Define the nested structure: 'coast' and 'water' are the two limbs of the 
   tree, and both have two branches (alternatives). 
*/
nlogitgen type = fishmode(coast: beach|pier, water: charter|boat)

/* Display the tree structure in the results window. 'choice()' lists sample 
   frequencies for each alternative */
nlogittree fishmode type, choice(d)
	 
/* Estimate the model. 'p' is estimated by limbs (coast/water) and 'c' by 
   branches (fishing modes).
*/
nlogit d || type: p || fishmode: c, case(id) noconst nolog notree
	
/* or try: depvar, d; 
           alternative-specifc regressor, p; 
		   level 1, no regressor;
           level 2, income and an intercept;
		   notree: suppresses the tree.
*/
* nlogit d p || type: , base(coast) || fishmode: income, case(id) nolog notree

	
/* Predicted probabilities: 
	
   The 'predict, pr' command provides predicted probabilities for level 1, 
   level 2, and so on. 
   
   In this example there are two levels. The first-level probabilities are for 
   'coast' or 'water'. The second-level probabilities are for each of the four 
   fishing modes. 
*/
predict plevel1 plevel2, pr
tab fishmode, summarize(plevel2)

/* Marginal effects: 
   
   Neither the 'mfx' command nor the 'margeff' command are available after 
   'nlogit'. We need to compute AMEs manually.
*/
	
preserve 
	quietly summarize p
	gen delta=r(sd)/1000
	quietly replace p = p + delta if fishmode == "beach"
		
	* Predict at the new values
	predict pnew1 pnew2, pr
		
	* Subtract the two predictions and divide by the amount of the change
	gen dpdbeach = (pnew2-plevel2)/delta
		
	* The AME is the average of the previous quantity
	tab fishmode, summarize(dpdbeach)
restore 
	
** Compare
quietly eststo CL : asclogit d p c, case(id) alternatives(fishmode) noconstant
quietly eststo NL : nlogit d || type: p || fishmode: c, noconst case(id) nolog
	
esttab CL NL, keep(p c) stats(N aic bic) eq(1) b(%7.3f) mtitle se
* esttab CL NL, keep(#1:p) stats(N aic bic) eq(1) b(%7.3f) mtitle se


** II.4 Random Parameters Logit (Mixed Logit, optional)
/* Random Parameters Logit allows the coefficients (eta) on the 
   alternative-varying regressors to be random across individuals.
   
   For example, 
   
   U = X*beta_j + Z_j*eta + epsilon
   
   eta ~ N(gamma, Sigma)
   
   This is called mixed logit by some researchers. 
*/

/* The 'mixlogit' command does not give any option to declare case-specific 
   regressors. We need to manually create regressors for the intercepts and
   income. 
   
   For case-specific regressors, a normalization is needed. We set one 
   fishing mode as base category and construct three intercepts and interactions
   with income. 
*/

gen dbeach = (fishmode=="pier")
gen dboat = (fishmode=="boat")
gen dcharter = (fishmode=="charter")

gen ybeach = dbeach*income
gen yboat = dboat*income
gen ycharter = dcharter*income

drop if fishmode=="charter" | mode==4

** You may need to install package 'mixlogit'
* ssc install mixlogit, replace
	
/* Estimate the model. 'group()' is used to identify each case or individual. 
   Regressors with random coefficients are listed in 'rand()'. The rest are 
   regresors with nonrandom coefficients. */
mixlogit d c dbeach dboat ybeach yboat, group(id) rand(p) 


** II.5 Multinomial Probit
	
use "TA3_1.dta", clear

gen sector = 0 if lfp==0
replace sector = 1 if inlist(classwkr,13,14)
replace sector = 2 if inlist(classwkr,22,23)

** drop missings
keep if sector!=.

/* Multinomial probit: 

   The output is qualitatively similar to that from mlogit, though parameters 
   estimates are scaled differently, as in the binary case. The fitted log is 
   very close to that for the multinomial logit. 
*/
mprobit sector age age2 married educ black nchild citiz, baseoutcome(0)


** II.6 Ordered outcome models

/* In some cases, categorical data are naturally ordered. For example, education
   level, consider using it as our dependent variable. 
*/

gen educ_1 = 0
replace educ_1 = 1 if educ == 73
replace educ_1 = 2 if educ == 81
replace educ_1 = 3 if educ>=91 & educ <= 111
replace educ_1 = 4 if educ > 111
	
label define labeleduc 0 "Less than high-school" ///
	1 "High-school diploma or equivalent" ///
	2 "Some college but not degree" ///
	3 "College degree" ///
	4 "Master, professional school, doctorate degree", replace
label values educ_1 labeleduc
	
/* Ordered Logit model:

   The threshold parameters (four \texttt{/cut}s) appear to be statistically 
   significantly different from each other. So the four categories should not be
   collapsed into three or even less.
*/
ologit educ_1 age age2 educ_sp black citiz married if educ_sp!=., nolog

** Predicted probabilities
predict p1ologit p2ologit p3ologit p4ologit p5ologit, pr
summarize i.educ_1 p2ologit-p5ologit, separator(5)


** PART III: REMAINING TOPICS --------------------------------------------------

** III.1 Endogenous variables

** Probit with endogenous regressor (FIML estimates) 
ivprobit lfp age age2 (educ = educ_sp) married black nchild citiz if educ_sp!=., vce(robust)

** Two-step sequential estimates (LIML) 
ivprobit lfp age age2 (educ = educ_sp) married black nchild citiz if educ_sp!=., twostep first


** III.2 Binary models for panel data

use "TA3_3.dta", clear

** Random-effects probit model
xtset cpsidp year
xtprobit lfp age age2 educ married black nchild citiz

** compare probit and xtprobit
quietly eststo r_prob : probit lfp age age2 educ married black nchild citiz
quietly eststo r_xtprob : xtprobit lfp age age2 educ married black nchild citiz
esttab r_prob r_xtprob, star stats(N)

** Random-effects logit model
xtlogit lfp age age2 educ married black nchild citiz, i(cpsidp)

** compare logit and xtlogit 
quietly eststo r_logit : logit lfp age age2 educ married black nchild citiz
quietly eststo r_xtlog : xtlogit lfp age age2 educ married black nchild citiz, i(cpsidp)
esttab r_logit r_xtlog, star stats(N)
