** =============================================================================
** TA Session 6 - CHAPTER 6: TREATMENT EFFECTS
** Microeconometrics with Hanna Wang, IDEA, FALL 2022
** TA: Conghan Zheng

/* Contents

I) Regression Adjustment
	
	I.1 Homogenous response to treatment
	I.2 Heterogenous response to treatment
		
II) Matching
	
	II.1 Covariates matching
	
	II.2 Propensity score matching
	
III) Instrumental Variables
	
IV) Regression Discontinuity

	- Sharp design 
	- Fuzzy design 
	
V) Difference in Differences 	
	
	V.1 Repeated cross-sections
	V.2 Panel data 
*/	

/* Datasets

TA6_1.dta: Ham and LaLonde (1996), link: https://www.jstor.org/stable/2171928.
		   Source: National Supported Work (NSW) labor market experiment.
		   Treatment: on-the-job training lasting between 9 months and one year (1976-1977).
		   
		   y = real earnings in 1978; d = on-the-job training.
		   
		   Used in Part I: RA and Part II: Matching.

TA6_2.dta: Cross-sectional data on 4,361 women of childbearing age in Botswana.

TA6_3.dta: An artificial cross section of 20 observations.

TA6_4.dta: An artificial panel of 12 individuals and two periods.
*/	

** Packages to be installed ----

** 1. estout - Making regression tables
* ssc install estout 
	
** 2. psmatch2 - Perform propensity score matching 
* ssc install psmatch2
	
** 3. pscore - Estimation of average treatment effects by Becker and Ichino
* help st0026_2
	
** 4. Module to perform coarsened exact matching 
* ssc install cem

** 5. (Binary) treatment models
* ssc install ivtreatreg
	
** 6. RD
* ssc install rd 
	
** 7. DID
* ssc install diff


cls, clear all
cap set more off

cd "/Users/zheng/Documents/conghan/UAB/IDEA_PhD/Teaching/TA_Microeconometrics_Fall_IDEA/03_TA_sessions/TA6/TA6_Stata"
	

** PART I: REGRESSION ADJUSTMENT -----------------------------------------------
** Stata command: teffects ra

** TA6_1.dta: y = real earnings in 1978; d = on-the-job training
use "TA6_1.dta", clear
	
describe 
summarize 
	
gen d = train
tab d

gen y = re78

** Compare y of treated and untreated individuals
tabulate d, summarize(y) means standard 
	
** Controls 
global x re74 re75 age agesq nodegree married black hisp
	
** Two groups are not very different in observed controls
bysort d: tabstat y $x, columns(statistics) s(mean sd)

** I.1 Homogenous response to treatment: difference in mean ----

reg y d
eststo ols_y_d

reg y d $x
eststo ols_y_d_x
	
** Compare the estimates 
esttab ols_y_d ols_y_d_x, se compress keep(d) stats(r2 N) mtitles("OLS w/o X" "OLS w X")
	
** I.2 Heterogeneous response to treatment ----

/* 'teffects ra' estimates treatment effects via regression adjustment (RA). 
   
   It contrasts averages of treatment-specific predicted outcomes. It offers the
   choice of several functional forms to deal with continuous, binary, count and
   nonnegative outcomes. 
   
   Here our outcomes y is linear.
*/
		
** ATE 	
teffects ra (y $x, linear) (d)
eststo ate
	
** ATE on the treated
teffects ra (y $x, linear) (d), atet
eststo atet

esttab ols_y_d ols_y_d_x ate atet, se compress keep(d r1vs0.d) stats(N) mtitles("Hom:w/o X" "Hom:w X" "Het:ATE" "Het:ATT")
	
predict ATE_x, te
sum ATE 
sum ATE_x if d==1 // ATE on the treated
sum ATE_x if d==0 // ATE on the non-treated

** Cofficients for controls of the two goups
teffects ra (y $x, linear) (d), aequations
	
/* Treatment effect as a percentage of the untreated potential outcome mean (the
   mean earnings that would occur if nobody was trained)
   
   Training program increases earnings by 33.8%, significantly.
*/
qui teffects ra (y $x, linear) (d), coeflegend
nlcom _b[ATE:r1vs0.d]/_b[POmean:0.d]					
	
** A nonlinear example: binary outcome (getting unemployed)
	
** ATE 
teffects ra (unem78 $x, probit) (d)

** Potential outcome means
teffects ra (unem78 $x, probit) (d), pomeans


** PART II: MATCHING -----------------------------------------------------------

gen re74sq = re74^2
gen u74hisp = unem74*hisp
gen educsq = educ^2

global x_ps age agesq educ marr nodegree black hisp re74 re75  

** II.1 Covariates matching ----
	
** Nearest-neighbor (NN) matching 
teffects nnmatch (y $x_ps) (d)
	
** Specify the minimum number of matches per observation (default is 1)
teffects nnmatch (y $x_ps) (d), nneighbor(3)

** NN matching with exact matching on some binary variables 
teffects nnmatch (y $x_ps) (d), nneighbor(3) ematch(hisp black)	

** Bias adjustment for large samples
teffects nnmatch (y $x_ps) (d), biasadj(re74 re75)	


** II.2 Propensity score matching ----

/* Stata commands: 
   
   - estimating propensity score: pscore
   
   - matching propensity score: psmatch2
   
   - estimating treatment effect: teffects psmatch; att*
*/
	
** NN matching, propensity scores are predicted from a logit of D on X
teffects psmatch (y) (d $x_ps)
	
** Assess the overalap graphically 
qui teffects psmatch (y) (d $x_ps), generate(near_obs)
** Nearest neighbors 
sum near_*
teffects overlap
graph export "ta6_psoverlap.png", as(png) replace


/* Command 'pscore': estimates the propensity score and tests the Balancing 
   Hypothesis.

   Link: https://www.sscnet.ucla.edu/soc/faculty/mason/readings/becker_ichino_pscore_sj_2002.pdf
   
   - Option 'pscore' is a compulsory option and specifies the variable name for 
   the estimated propensity score.
   
   - Option 'comsup' restricts the analysis of the balancing property to all 
   treated plus those controls in the region of common support. A dummy variable
   named comsup is added to the dataset to identify the observations in the 
   common support.
   
   - Option 'blockid' specifies the variable name for the block number of the 
   estimated propensity score.
   
   - What is a block? 
   Blocks defined over intervals of propensity score. Stratification or interval
   matching is based on the idea of dividing the range of variation of the 
   propensity score in intervals such that within each interval, the treated and
   control units have, on the average, the same propensity score.
   
   The overall treatment effect is a weighted average of block-specific 
   treatment effects, where the weight for each block is given by the 
   corresponding fraction of treated units.
   
   - Option 'level(#)' specifies the significance level of the tests of the 
   balancing property.
   
   - Option 'logit' specifies that a logit model to estimate the propensity 
   score be used instead of the default probit model
*/
pscore d $x_ps, pscore(myscore) comsup blockid(myblock) level(0.005) logit


/** 'psmatch2': matching the propensity score

   Link: http://repec.org/bocode/p/psmatch2.html
*/

** k-Nearest neighbors matching: neighbor(k)
psmatch2 d, out(y) pscore(myscore) neighbor(3) common 
	
** Test of the balancing property 
pstest myscore
pstest $x_ps
	
** pstest results for treated and untreated 
summarize myscore [aweight=_weight] if d==0
summarize myscore [aweight=_weight] if d==1 


** Graph: quality of our matching
label define tstatus 0 "Control sample" 1 "Treated sample"
label values d tstatus
label variable d "Treatment Status"
	
** Before matching 
qui graph twoway (kdensity myscore if d==1, msize(small)) ///
	(kdensity myscore if d==0, msize(small) ///
	lpattern(shortdash_dot)), ///
	subtitle(, bfcolor(none)) ///
	xtitle("propensity–score (Before)", size(medlarge)) xscale(r(0.1 0.7)) /// 
	ytitle("Density", size(medlarge)) yscale(r(0 9)) ///
	legend(pos(12) ring(0) col(1)) ///
	legend( label(1 "Treated") label(2 "Untreated")) saving(BEFORE, replace)

** After matching 
qui graph twoway (kdensity myscore [aweight=_weight] if d==1,  msize(small)) ///
	(kdensity myscore [aweight=_weight] if d==0, msize(small) ///
	lpattern(shortdash_dot)), ///
	subtitle(, bfcolor(none)) ///
	xtitle(" propensity–score (After)", size(medlarge)) xscale(r(0.1 0.7)) ///
	ytitle("Density", size(medlarge)) yscale(r(0 9)) ///
	legend(pos(12) ring(0) col(1)) ///
	legend( label(1 "Treated") label(2 "Untreated")) saving(AFTER, replace)
		
** Combine graphs
graph combine BEFORE.gph AFTER.gph
graph export "ta6_match.png", as(png) replace


/* Commands 'att*': estimate the ATT based on the propensity score

   - 'attnd' and 'attnw': NN matching, standard errors are obtained analytically
   or by bootstrapping
   
   Difference: 
   'attnw' gives equal weight to the groups of forward and backward 
   matches; 
   'attnd' randomly draws either the forward or backward matches.
   
   
   - 'attr': radius matching
   
   - 'attk': kernel matching, default: Gaussian kernel, standard errors are 
   obtained by bootstrapping.
   
   - 'atts': stratification method

*/

** ATT estimate based on NN matching
attnd y d $x_ps, comsup

** ATT estimate based on radius matching
attr y d $x_ps, comsup logit
	
attr y d $x_ps, comsup logit radius(0.001)

attr y d $x_ps, comsup logit radius(0.0001)

** ATT estimate based on statification matching 
atts y d, pscore(myscore) blockid(myblock) comsup

** ATT estimate based on Kernel matching, bootstrap repetitions: 50
attk y d $x_ps, comsup boot reps(50) dots logit

	
** III) INSTRUMENTAL VARIABLES ---------------------------------------------------
** Link: http://fmwww.bc.edu/repec/bocode/i/ivtreatreg_cerulli.pdf

** TA6_2.dta: y = number of children; d = having more than 7 years of education
use "TA6_2.dta", replace

** Difference in Means: y = alpha*D + epsilon
reg children educ7 
eststo DIM

** Control-Function OLS: Homogeneous treatment effects
ivtreatreg children educ7 age agesq evermarr urban electric tv, iv(frsthalf) model(cf-ols)
	
** Control-Function OLS: Heterogeneous treatment effects
ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(cf-ols) graphic
graph save ta6_cf_ols, replace
eststo CF_OLS
	
/* Direct 2SLS: Homogeneous treatment effects: 
	
   Instrument: proposed by Angrist and Krueger (1991), individuals born at the 
   beginning of the year (frsthalf) start school at a slightly older age.
*/
ivtreatreg children educ7 age agesq evermarr urban electric tv, iv(frsthalf) model(direct-2sls)

/* Direct 2SLS: Heterogeneous treatment effects: 

   We can identify LATE, the treatment effect on the group of compliers (those 
   who were born in the first half of the year and receive more than seven years
   of education, but would not otherwise receive that level of education if they
   were born in the second half of the year). We assume monotonicity (no 
   defiers).
   
   Option 'hetero' specifies the variables over which to calculate the 
   idyosincratic Average Treatment Effect ATE(x), ATT(x) and ATNT(x).
*/
ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(direct-2sls) graphic
graph save ta6_direct_2sls, replace
eststo DIRECT_2SLS
			
/* Probit-2SLS: Homogeneous treatment effects 

Probit: Estimate the probability of receiving more than seven years of education
using a probit model on the covariates and the instrument. 
	  
OLS of the treatment on the covariates and the predicted probability of 
receiving the treatment.
*/
ivtreatreg children educ7 age agesq evermarr urban electric tv, iv(frsthalf) model(probit-2sls)
	
** Probit-2SLS: Heterogeneous treatment effects 
ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(probit-2sls) graphic
graph save ta6_probit_2sls, replace
eststo PROBIT_2SLS

** Probit-OLS: Homogeneous treatment effects
ivtreatreg children educ7 age agesq evermarr urban electric tv, iv(frsthalf) model(probit-ols)
	
** Probit-OLS: Heterogeneous treatment effects 
ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(probit-ols) graphic
graph save ta6_probit_ols, replace
eststo PROBIT_OLS

** Selection on unobservables (Heckit)
ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(heckit) graphic
graph save ta6_heckit, replace
eststo HECKIT	

** Obtain standard errors for ATT and ATNT via bootstrapping
bootstrap atet=e(atet) atent=e(atent), rep(100): ivtreatreg children educ7 age agesq evermarr urban electric tv, hetero(age agesq evermarr urban) iv(frsthalf) model(heckit)
	
esttab DIM CF_OLS DIRECT_2SLS PROBIT_2SLS PROBIT_OLS HECKIT, mtitle se compress stats(N) keep(educ7 G_fv)

** The distribution of ATE(x), ATET(x) and ATENT(x), x: number of children
graph combine ta6_cf_ols.gph ta6_direct_2sls.gph ta6_probit_2sls.gph ta6_probit_ols.gph ta6_heckit.gph
graph export "ta6_IV.png", as(png) replace

	
** IV) REGRESSION DISCONTINUITY ------------------------------------------------

** IV.1 Sharp design ----

** Simulate some data (N = 1000)
clear all
set obs 1000
gen s = 10 + 5*invnormal(uniform())
	
** Generate X
global s_star = 10
gen x = s - $s_star
	
** Generate D (sharp design)
gen d = 1 if s > $s_star
replace d = 0 if s <= $s_star
	
** Define y. The true value of the treatment effect will be E[Y1] - E[Y0] = 400
gen y1 = 600 + 6.5*x - 2*x^2 + 0.001*x^3 + 300*invnorm(uniform()) 
gen y0 = 200 + 6.5*x - 0.20*x^2 + 0.01*x^3 + 300*invnorm(uniform()) 
	
** Generate the observable outcome 
gen y = y0 + d*(y1 - y0)
	
** Visualize outcome discontinuity 
twoway histogram y if s > $s_star, barw(60) color(green%50) || ///
	hist y if s < $s_star, barw(60) color(red%30) ///
	legend(order(1 "Right side" 2 "Left side") pos(11) col(1) ring(0)) ///
	xtitle() ytitle(Frequency) ylabel()
	
** Polynomials for the control function
gen dx = d*x
gen dx2 = d*(x^2)
gen dx3 = d*x^3
gen x2 = x^2
gen x3 = x^3
	
** Estimation without high order polynomials 
reg y d x dx 
cap drop y_hat_1
predict y_hat_1, xb
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(line y_hat_1 s if s>=$s_star , msymbol(o)) ///
	(line y_hat_1 s if s<=$s_star , msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Sharp–RDD – Parametric linear regression") ///
	legend( label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right Prediction") label(4 "Left Prediction"))
		
graph export "ta6_sharprd1.png", as(png) replace
	
** Improve the fit by including higher order polynomials 
reg y d x x2 x3 dx dx2 dx3
cap drop y_hat_2
predict y_hat_2, xb
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(scatter y_hat_2 s if s>=$s_star , msymbol(o)) ///
	(scatter y_hat_2 s if s<=$s_star , msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Sharp–RDD – Parametric Polynomial Regression") ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right Prediction") label(4 "Left Prediction"))

graph export "ta6_sharprd2.png", as(png) replace
	
/* Local polynomial functions f0 and f1: We can select different kernels, 
   bandwiths and degrees. 
*/
global bandw 5 
cap drop f0 f1
lpoly y s if s <= $s_star, gen(f0) at(s) kernel(triangle) bwidth($bandw) degree(3) nograph
lpoly y s if s > $s_star, gen(f1) at(s) kernel(triangle) bwidth($bandw) degree(3) nograph

graph twoway (scatter y s if s >= $s_star , clstyle(p1)) ///
	(scatter y s if s <= $s_star , clstyle(p1)) ///
	(scatter f0 s if s < $s_star, msize(medsmall) msymbol(o)) ///
	(scatter f1 s if s >= $s_star, msize(medsmall) msymbol(o)), ///
	xline($s_star, lpattern(dash)) ///
	title("Sharp RDD – Local polynomial regression (LPR)") ///
	legend(label(1 "Right actual data") label(2 "Left actual data") ///
	label(3 "Right LPR prediction") label(4 "Left LPR prediction")) ///
	note(Bandwidth = $bandw)

graph export "ta6_sharprd3.png", as(png) replace

** Calculate the treatment effect 
gen z=$s_star
cap drop f0 f1 
qui lpoly y s if s<=$s_star, gen(f0) at(z) k(tri) bw($bandw) deg(3) nogr
qui lpoly y s if s>$s_star, gen(f1) at(z) k(tri) bw($bandw) deg(3) nogr
scalar rdef=f1[1]-f0[1]
display rdef

** Obtain standard errors for the local polynomial estimates 
cap program drop rdd_s
program rdd_s, rclass
	version 16 
	args deg ker band cut
	cap drop f0 f1 z
	gen z = `cut'
	qui lpoly y s if s < `cut', gen(f0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly y s if s >= `cut', gen(f1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	return scalar rdef = f1[1]-f0[1]
end
	
/* 3-degree polynomial with triangular kernel, bandwith of 5, and discontinuity 
   cut equal to 10 
*/ 
rdd_s 3 tri 5 10
return list 
bootstrap r(rdef), reps(50): rdd_s 3 tri 5 10
	
/* Regression discontinuity estimates using the 'rd' command 

   'rd' estimates local linear or kernel regression models on both sides of the
   cutoff.  Estimates are sensitive to the choice of bandwidth. 
*/ 
rd y s, z0($s_star) bwidth($bandw)
	
** Various bandwiths
rd y s, z0($s_star) bwidth($bandw) mbw	


** b) Fuzzy design ----
	
** Simulate some data 
clear all
set obs 1000
	
gen s = -1 + 2*runiform()
gen Z = (s >= 0)
gen v = rnormal(0,1)
gen d = (-0.5 + Z + s + v >= 0)
gen y1 = 2 + s + s^2 + 3*s^3 + invnorm(uniform())
gen y0 = 1 + s + s^2 + 3*s^3 + invnorm(uniform())
gen y = y0 + d*(y1 - y0)

gen s2 = s^2
gen s3 = s^3

** Under the fuzzy design, the treatment is endogeneous
ivregress 2sls y s s2 s3 (d=Z), first 
	
** Nonparametric estimation using a 3rd-degree polynomial
global s_star = 0
global bandw = 5 
cap drop f0 f1 
lpoly y s if s<$s_star, gen(f0) at(s) k(tri) bw($bandw) deg(3) nogr
lpoly y s if s>=$s_star, gen(f1) at(s) k(tri) bw($bandw) deg(3) nogr

** discountinuity in y
graph twoway (scatter y s if s>=$s_star , clstyle(p1)) ///
	(scatter y s if s<=$s_star , clstyle(p1)) ///
	(scatter f0 s if s<$s_star, msize(medsmall) msymbol(o)) ///
	(scatter f1 s if s>=$s_star, msize(medsmall) msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Fuzzy–RDD – Outcome Non–parametric Local Linear Regression", size(medlarge)) ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right LLR Prediction") label(4 "Left LLR Prediction")) ///
	note(Bandwidth = $bandw)
graph export "ta6_fuzzyrd1.png", as(png) replace
	
capture drop g0 g1
lpoly d s if s<$s_star, gen(g0) at(s) k(tri) bw($bandw) deg(3) nogr
lpoly d s if s>=$s_star, gen(g1) at(s) k(tri) bw($bandw) deg(3) nogr

** discontinuity in P(D=1)
graph twoway ///
	(scatter d s if s>=$s_star & d==1, clstyle(p1)) ///
	(scatter d s if s<=$s_star & d==0, clstyle(p1)) ///
	(scatter g0 s if s<$s_star, msize(medsmall) msymbol(o)) ///
	(scatter g1 s if s>=$s_star, msize(medsmall) msymbol(o)), xline($s_star, lpattern(dash)) ///
	title("Fuzzy–RDD – Probability Non–parametric Local Linear Regression") ///
	legend(label(1 "Right Actual Data") label(2 "Left Actual Data") ///
	label(3 "Right LLR Prediction") label(4 "Left LLR Prediction")) ///
	note(Bandwidth = $bandw)
graph export "ta6_fuzzyrd2.png", as(png) replace

** Obtain standard errors for the local polynomial estimates by bootstrapping
cap program drop rdd_f
program rdd_f, rclass
	version 16
	args deg ker band cut
	* Outcome discontinuity
	cap drop z f0 f1 g0 g1
	gen z = `cut'
	qui lpoly y s if s<`cut', gen(f0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly y s if s>=`cut', gen(f1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	scalar disc_y = f1[1]-f0[1]
	* Probability discontinuity
	cap drop g0 g1
	qui lpoly d s if s<`cut', gen(g0) at(z) k(`ker') bw(`band') deg(`deg') nogr
	qui lpoly d s if s>=`cut', gen(g1) at(z) k(`ker') bw(`band') deg(`deg') nogr
	scalar disc_d = g1[1]-g0[1]
	return scalar rddef = disc_y/disc_d
end

cap drop z 
cap drop f0 f1 g0 g1
bootstrap r(rddef), reps(10): rdd_f 3 tri 5 0
	
** Or we can use the 'rd' command 
rd y d s, z0(0) bwidth(5)
	
	
** V) DIFFERENCE IN DIFFERENCES ----------------------------------------------
	
** V.1 Repeated cross-sections

** We assume that the sample composition does not vary over time
use "TA6_3.dta", clear
	
** Method 1: DID estimate for ATE is the coefficient estimate of DT
reg Y D T DT 
	
** Method 2: step by step: 

** Average over treated units in T = 0
qui summarize Y if D==1 & T==0
scalar treat0=r(mean)
	
** Average over treated units in T = 1
qui summarize Y if D==1 & T==1 
scalar treat1=r(mean)
	
** Average over untreated units in T = 0 
qui summarize Y if D==0 & T==0
scalar untreat0=r(mean)
	
** Average over untreated units in T = 1
qui summarize Y if D==0 & T==1
scalar untreat1=r(mean)	
	
** Diff-in-diff coefficient
scalar didbeta = treat1 - treat0 - (untreat1 - untreat0)
di "The DID coefficient is " didbeta
	
** Method 3: the 'diff' command
diff Y, treated(D) period(T)
	

** V.2 Panel data

** By construction, the sample composition is the same across time
use "TA6_4.dta", clear

xtset id year 
	
** Lags
gen y_1 = L.y 
gen d_1 = L.d 
gen x1_1 = L.x1 
gen x2_1 = L.x2 

** Differences
gen delta_y = D.y
gen delta_x1 = D.x1 
gen delta_x2 = D.x2 
gen delta_d = D.d 

** Method 1: estimate for ATE is the coefficient estimate for d
reg delta_y d if d_1==0

** Method 2: step by step
qui sum delta_y if d == 1 & d_1==0 
scalar mean_t = r(mean)
qui sum delta_y if d == 0 & d_1==0
scalar mean_c = r(mean)
scalar did = mean_t - mean_c 
di "The DID coefficient is " did 
	
** If fixed effects are considered, the first-differenced model can be used  
reg delta_y delta_d, noconst
	
** WG estimator: the coefficient estimate for d*dyear_2000 is the DID estimate
xi, noomit prefix("d") gen i.year
xtreg y d##dyear_2000, fe		

** DID with covariates 
reg delta_y d delta_x1 delta_x2
