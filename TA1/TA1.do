** =============================================================================
** TA Session 1 - CHAPTER 1: GMM and MLE
** Microeconometrics with Hanna Wang, IDEA, FALL 2022
** TA: Conghan Zheng

cls, clear all
cap set more off

cd "..."

	
** PART I: GMM =================================================================

use "TA1.dta", clear

/*
- Commands: 'ivregress', 'gmm', or 'ivreg2' (or a moment evaluator program), 
please check their manuals: 

	- iveregress: single-equation IV
	manual: https://www.stata.com/manuals/rivregress.pdf
	
	- gmm: multiple equation GMM (3SLS)
	manual: https://www.stata.com/manuals13/rgmm.pdf
	
	- ivreg2: community-contributed, IV/2SLS/GMM/LIML...
	manual: http://www.repec.org/bocode/i/ivreg2.html
	
For the same model, these commands deliver you the same results, with a bit 
different syntax. It doesn't matter which command or even which language you 
use.
*/

* est sto clear

global xlist totchr age female blhisp linc

global olsmodel "ldrugexp hi_empunion $xlist"

global just_identified_model "ldrugexp (hi_empunion = ssiratio) $xlist"

global over_identified_model "ldrugexp (hi_empunion = ssiratio multlc) $xlist"


** I.1: JUST-IDENTIFIED MODEL ==================================================
/*
In just-identified case, the choice of weighting matrix has no impact on the 
solution to the minimization problem, IV and 2SLS are equivalent. 
*/

** just-identified model, OLS
qui ivregress 2sls $olsmodel
est sto ols1

** just-identified model, IV/2SLS
qui ivregress 2sls $just_identified_model
est sto iv1

** DWH test
/* we can't use classical Hausman test since it's valid only under homoskedasticity.
 if you cannot run the following commands, run "update all" in your console. */
* estat endogenous 
/* p < 0.05, hi_empunion is endogenous. */

** just-identified model, one-step GMM
/* option wmatrix(unadjusted): do not adjust for heteroskedasticity, GMM reduces
 to 2SLS*/
qui ivregress gmm $just_identified_model, wmatrix(unadjusted) 
est sto gmm_1s1

** just-identified model, two-step optimal GMM
/* option wmatrix(robust): the optimal weighting matrix, and is the default */
qui ivregress gmm $just_identified_model, wmatrix(robust)
est sto gmm_optimal1
	  
** just-identified model, two-step iterated optimal GMM
qui ivregress gmm $just_identified_model, wmatrix(robust) igmm
est sto gmm_iterated1

** export regression table
/* You can also write it to latex file with extra options, see the guide at
http://repec.org/bocode/e/estout/esttab.html#esttab005
*/
esttab ols1 iv1 gmm_1s1 gmm_optimal1 gmm_iterated1, se nostar compress mtitles("OLS" "IV" "GMM_1S" "GMM_Optimal" "GMM_iterated")


** I.2: OVER-IDENTIFIED MODEL ==================================================

** over-identified model, IV / 2SLS
qui ivregress 2sls $over_identified_model
est sto iv2

** over-identified model, one-step GMM
qui ivregress gmm $over_identified_model, wmatrix(unadjusted) 
est sto gmm_1s2

** over-identified model, two-step optimal GMM
qui ivregress gmm $over_identified_model, wmatrix(robust)
est sto gmm_optimal2
	  
** over-identified model, two-step iterated optimal GMM
qui ivregress gmm $over_identified_model, wmatrix(robust) igmm
est sto gmm_iterated2

** export regression table
esttab iv2 gmm_1s2 gmm_optimal2 gmm_iterated2, se nostar compress mtitles("IV" "GMM_1S" "GMM_Optimal" "GMM_iterated")


** I.3: multiple equation GMM ==================================================

/*
simultanous equations: 
(eq1)
(eq2)

Stata commands:
	- 2SLS: ivregress, ivreg2, gmm
	- 3SLS: gmm, reg3

*/

** use stata built-in data
use http://www.stata-press.com/data/r13/klein, clear

** y: consump
** endogenous x: wagepriv
** available instruments z: wagegovt, govt, capital1 
global eq1 "consump - {eq1: _cons wagepriv wagegovt}"
global eq2 "wagepriv - {eq2: _cons consump govt capital1}"
global eq2_yhat "wagepriv - {eq2: _cons consump_hat govt capital1}"

** Separate estimation ====
/* Two GMM are done separetely, the predicted value from the first estimation 
is pluged into the second estimation. */

** GMM estimation of eq1 
gmm ($eq1), instruments(wagegovt govt capital1) onestep
est sto gmm_eq1

** calculate fitted value from eq1
gen consump_hat =  _b[_cons] + _b[wagepriv]*wagepriv + _b[wagegovt]*wagegovt

** GMM estimation of eq2, with fitted value of consump_hat plugged in 
gmm ($eq2_yhat), instruments(wagegovt govt capital1) onestep
est sto gmm_eq2

** Joint estimation ====
/* Equation-by-equation GMM, the minimization of distance is done jointly. 
If independence of errors is assumed across equations, this becomes 3SLS. */

gmm (eq1: $eq1) (eq2: $eq2), ///
instruments(eq1: wagegovt govt capital1) ///
instruments(eq2: wagegovt govt capital1) ///
winitial(unadjusted, independent) wmatrix(robust) twostep
est sto gmm_joint

esttab gmm_eq1 gmm_eq2 gmm_joint, se nostar compress


** PART II: MLE ================================================================

use "TA1.dta", clear
	
/* Command: 'ml'
	- manual: https://www.stata.com/manuals/rml.pdf
*/

**  II.1: Example 1: Normal distribution =======================================

** drop previous program in environment called "lf_norm", if any
cap program drop lf_norm 

/* 
- program lf_norm ... end: define the evaluator program, calling it 'if_norm' 
- version: specify your version of Stata
- args: arguments
	- lnf: local macro, individual contribution of each observation to log-likelihood 
	- mu: conditional mean, in the current case is xb = b0 + b1*x1 + b2*x2 + ...
	- sigma: standard deviation (better to write sigma in log space)
- $ML_y1: the first and only dependent variable here
*/

program lf_norm
    version 17
	args lnf mu sigma 
	quietly replace `lnf' = ln(normalden($ML_y1, `mu', `sigma'))
end 

/* - ml model: define the context of the problem to solve
     - method 'lf': linear form
   
   - ml maximize: maximizes the likelihood function and reports results
*/

** We specify a homoskedastic model, sigma is not a function of other variables
ml model lf lf_norm (mu: ldrugexp = $xlist) (sigma: ) 
ml maximize

** We specify a heteroskedastic model, sigma is modelled as a linear function of age
ml model lf lf_norm (mu: ldrugexp = $xlist) (sigma: age) 
ml maximize

	
** II.2: Example 2: Logistic regression ========================================

cap program drop lf_logit

program lf_logit
	version 17
	args lnf xb
	quietly replace `lnf' = ln(invlogit(`xb')) if $ML_y1==1
	quietly replace `lnf' = ln(1-invlogit(`xb')) if $ML_y1==0
end

** only one parameter: xb, you just need to specify one equation
ml model lf lf_logit (beta: good = age age2 marry educyr msa)
ml maximize
