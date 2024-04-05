** =============================================================================
** TA Session 4 - CHAPTER 4: CENSORING, TRUNCATION, AND SELECTION
** Microeconometrics with Hanna Wang, IDEA, FALL 2022
** TA: Conghan Zheng

cls, clear all
set seed 13
cap set more off

cd "/Users/zheng/Documents/conghan/UAB/IDEA_PhD/Teaching/TA_Microeconometrics_Fall_IDEA/03_TA_sessions/TA4/TA4_Stata"

/* Contents
   
   I)   Censoring 
   II)  Truncation
   III) Selection
	
*/

/* Dataset: 

   Cross-section from 2001 Medical Expenditure Panel Survey (MEPS). We are 
   interested in studying ambulatory expenditure ('ambexp') 
*/

use "TA4.dta", clear


** PART I: Censoring -----------------------------------------------------------

summarize ambexp age female educ blhisp totchr ins

global xlist age female educ blhisp totchr ins

** I.1 Tobit Regression and post-estimation
	
/* Tobit model: 'll()' indicates the left-censoring limit. The coefficients are 
   interpreted as the partial derivative of the latent variable y* w.r.t. x. 
   
   MLE for tobit regression in Stata has the following basic syntax: 
   
   . tobit depvar [indepvars] ..., ll[(#)] ul[(#)] [options]
   
   The specification ll[(#)] and ul[(#)] refer to the lower limit 
   (left-censoring point) and the upper limit (right censoring point), 
   respectively.
*/
tobit ambexp $xlist , ll(0) vce(robust)
	
** Tobit prediction: in-sample fitted values of the latent variable y*
predict yhatlin
summarize yhatlin, detail
summarize ambexp, detail

** Plot: data vs. fitted values
twoway kdensity ambexp || kdensity yhatlin, ///
	title("Ambulatory Expenditure: Data vs. Tobit Prediction") ///
	ytitle(" ") ///
	xlabel(,format(%9.0gc)) ///
	xtitle("Ambulatory expenditure") ///
	legend(label(1 "Data") label(2 "Prediction")) ///
	subtitle("Frequency density", position(11) justification(left)) ///	
	lwidth(medthick)
graph export "graph_tobit_fit", as(png) replace

** Marginal effects
	
** At-means ME for the left-truncated mean: E(y|x, y > 0)
** e(a,b) is for the interval a<y<b
mfx compute, predict(e(0, .))

** At-means ME for right-truncated (at the median) mean: E(y|x, 0<y<535)
mfx compute, predict(e(0, 535))
	
** At-means ME for the censored mean E(y|x)
mfx compute, predict(ystar(0, .))
* estout ME, cells("b p") margin style(fixed)
* esttab ME, cells("b") margin

	
** Marginal impact on probabilities
** (275, 1913) = (the 25th percentile, the 75th percentile)
quietly tobit ambexp $xlist , ll(0) vce(robust)
mfx compute, predict(pr(113, 1618)) 


** I.2 Tobit with log-normal data
	
** Summarize amexp with the detail optio
summarize ambexp, detail

** Summarize ignoring zero values
summarize ambexp if ambexp > 0, detail
	
** Let's see what happens if we condition on regressors
reg ambexp $xlist , vce(robust)
predict olsres, resid
	
** Summarize OLS residuals
summarize olsres, detail
	
** Let's try OLS on observations with positive values of ambexp 
reg ambexp age $xlist if ambexp > 0, vce(robust)
predict olsres2, resid
	
** Same conclusion from the residuals
summarize olsres2, detail
drop olsres olsres2

** Summarize the dependent variable in logs
summarize lambexp, detail
		
/* Notice that the smallest positive value of 'ambexp' is 1, in which case 
   'lambexp' is 0. Stata's ll(0) option mistakenly treats these observations as 
   censored rather than as zero, leading to shrinkage in the sample size for 
   noncensored observations. 
   
   To avoid this loss, we set all censored observations of lambexp to an amount 
   slightly smaller than the minimum noncensored value of lambexp. 
*/

generate y = ambexp
generate dy = (ambexp > 0)
gen lny = ln(y)
quietly summarize lny 
scalar gamma = r(min)
display "gamma = " gamma
scalar gamma01 = gamma - 0.0000001 
replace lny = gamma01 if lny ==.
	
** 526 observations have expenditure equal to zero. In logs, these zeros become -1.00e-07
tabulate y if y < 0.02
tabulate lny if lny < gamma + 0.02
	
** Estimate the tobit model for the new lognormal variable
tobit lny $xlist , ll(gamma01) vce(robust)
eststo tobit_log

** Compare Tobit with OLS
reg lny $xlist , vce(robust)
eststo ols_log

esttab tobit_log ols_log, se compress mtitle wide
	
/* Two-limit Tobit
   
   In less than 1.5% of the sample, ambexp exceeds 10,000 USD. Suppose that we 
   want to exclude these high values that contribute to the nonnormal kurtosis. 
   
   Choosing 10,000 as the upper censoring point, let's estimate a two-limit 
   Tobit.
*/
scalar upper = log(10000)
tobit lny $xlist, ll(gamma01) ul(upper) vce(robust)

/* Prediction: 
   We have been dealing with a log-transformed dependent variable. If we want 
   predictions in levels, we need to do some transformations.
*/

qui tobit ambexp age female educ blhisp totchr ins, ll(0)
** Censoring point
scalar gamma = 0
** Estimates
matrix btobit = e(b)
/* Estimated standard error of the regression, e(df_m): model degrees of freedom,
   btobit[1,e(df_m)+2]: var(e.ambexp)*/
scalar sigma = sqrt(btobit[1,e(df_m)+2])
** Linear prediction
predict xb
** Standardized censoring point
gen threshold = (gamma-xb)/sigma
	
** Predict y_hat (notice that our dependent variable is not in logs)
gen yhat = exp(xb + 0.5*sigma^2)*(1-normal((gamma-xb-sigma^2)/sigma))
gen ytrunchat = yhat/(1 - normal(threshold)) if dy == 1
summarize y yhat ytrunchat if dy == 1


** I.3 Two-part model

/* Part 1- A binary outcome equation that models Pr(ambexp > 0) using probit or 
	   logit.
   Part 2- Linear regression to model E(ln(ambexp)| ambexp > 0).
		
   In two-part models, unlike Tobit, neither homoskedasticity nor normality are 
   necessary for the estimator to be consistent. The key assumption is that
   E(ln(y)|d = 1, x) is linear in x.
*/
		
** Part 1
probit dy $xlist , nolog vce(robust)
eststo part1_probit
** Save the value of the log-likelihood
scalar llprobit = e(ll) 
		
** Part 2
reg lny $xlist if dy==1, vce(robust)
eststo part2_ols
scalar lllognormal = e(ll)
predict rlambexp, residuals

esttab part1_probit part2_ols, se compress mtitle wide
		
** Joint log-likelihood: sum of the two log-likelihoods
scalar lltwopart = llprobit + lllognormal
display "lltwopart = " lltwopart
		
/* If we compare with the tobit model used before, the two-part model fits the 
   data considerably better. */
quietly tobit $xlist, ll(0) vce(robust)
display e(ll)

** Prediction 
qui probit dy age female educ blhisp totchr ins
predict dyhat, pr
qui reg lny age female educ blhisp totchr ins if dy==1
predict xbpos, xb
	
** Fitted log values from the second part of the two-part model
gen yhatpos = exp(xbpos+0.5*e(rmse)^2)
	
/* Estimate the unconditional values by multiplying by the fitted probability of
   the positive expenditure from the probit regression.
*/
gen yhat2step = dyhat*yhatpos
summarize y yhat2step
summarize yhatpos if dy==1

** PART II: TRUNCATION ---------------------------------------------------------

** Truncation regression
truncreg lny $xlist , ll(gamma01) vce(robust)
eststo truncate_log

tobit lny $xlist , ll(gamma01) nolog vce(robust)
eststo censor_log

esttab truncate_log censor_log, se compress mtitle wide


** PART III: SELECTION ---------------------------------------------------------
	 
** 1. FIML (Stata's default) without exclusion restriction (x1 = x2)

heckman lny $xlist, ///
		select(dy = $xlist) nolog
eststo HKM_FIML	

** 2. LIML without exclusion restriction
heckman lny $xlist , ///
		select(dy = $xlist) ///
		nolog twostep
eststo HKM_LIML	

** 3. LIML with exlusion restriction (x1 != x2)
heckman lny $xlist, ///
		select(dy = $xlist income) ///
		nolog twostep
eststo HKM_LIML_ex
	
** Compare the three sets of estimates
esttab HKM_FIML HKM_LIML HKM_LIML_ex, se wide mtitle

** Prediction

qui heckman lny age female educ blhisp totchr ins, select(dy = age female educ blhisp totchr ins)
* Predicted probability of observing the outcome 
predict probpos, psel
* Linear prediction for selection equation 
predict x1b1, xbsel
* Linear prediction for main equation
predict x2b2, xb
* Estimated variance of the main equation
scalar sig2sq = e(sigma)^2
* Covariance of the errors 
scalar sig12sq = (e(rho)*e(sigma))^2
display "sigma1sq = 1" "sigma12sq = " sig12sq "sigma2sq = " sig2sq
	
** Predict y_hat
gen yhatheck = exp(x2b2 + 0.5*(sig2sq))*(1-normal(-x1b1-sig12sq))
gen yhatposheck = yhatheck/probpos
summarize yhatheck y probpos dy
summarize yhatposheck y probpos dy if dy == 1


** PART IV: IV TOBIT ----------------------------------------------------------
	
** Use the number of chronic conditions as instrument for having insurance
ivtobit ambexp age female educ blhisp (ins = totchr), ll(0)
ivtobit ambexp age female educ blhisp (ins = totchr), ll
	
** Two-step estimation 
ivtobit ambexp age female educ blhisp (ins = totchr), ll twostep
	
** At-means ME for the left-truncated mean E(y|x, y > 0)	
qui ivtobit ambexp age female educ blhisp (ins = totchr), ll
margins, dydx(*) predict(e(0, .))
