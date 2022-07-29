/*********************************************************************************************************************************************************************************************************************************************************************
EXAMPLE 1: CARD DATASET
**************************************************************************************************************************************************************************************************************************************************/

***************Simple examples for ivqte***************

use http://fmwww.bc.edu/ec-p/data/wooldridge/CARD.dta

*CASE 1: CONDITIONAL EXOGENOUS QUANTILE TREATMENT EFFECTS (Koenker and Bassett, 1978)
*using the official command qreg 
qreg lwage college age black fatheduc motheduc reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669, q(0.1)
*the same point estimates but different standard errors (consistent in case of heteroscedasticity) are obtained with ivqte
ivqte lwage age black fatheduc motheduc reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669 (college), q(0.1) variance

*CASE 2: CONDITIONAL IV QUANTILE TREATMENT EFFECTS (Abadie, Angrist and Imbens, 2002)
*point estimates and variance at the first decile
ivqte lwage (college=nearc4), q(0.1) variance dummy(black) continuous(age fatheduc motheduc) unordered(region) aai

*CASE 3: UNCONDITIONAL EXOGENOUS QUANTILE TREATMENT EFFECTS (Firpo, 2007)
*point estimates and variance at the 9 deciles
ivqte lwage (college), variance dummy(black) continuous(age fatheduc motheduc) unordered(region)

*CASE 4: UNCONDITIONAL ENDOGENOUS QUANTILE TREATMENT EFFECTS (Froelich and Melly, 2007)
*point estimates and variance at the 9 deciles
ivqte lwage (college = nearc4), variance dummy(black) continuous(age fatheduc motheduc) unordered(region)
*using estimated positive weights
ivqte lwage (college = nearc4), variance dummy(black) continuous(age fatheduc motheduc) unordered(region) positive

***************Advanced examples for ivqte***************

set seed 123
sample 500, count

* UNCONDITIONAL EXOGENOUS QUANTILE TREATMENT EFFECTS (Firpo, 2007)
*the bandwidth is set to 2, the lambda to 0.8, the bandwidth used to estimate positive weights is set to 0.5, the lambda used to estimate the positive weights is set to 1
ivqte lwage (college), q(0.5) dummy(black) continuous(exper) unordered(region) b(2) l(0.8) variance vb(.) vl(1)

*estimation of the optimal smoothing parameters by cross-validation
locreg college, dummy(black) continuous(exper) unordered(region) b(1 2) l(0.8 1) logit
*the optimal smoothing parameters are 1 and 1
*estimation with these optimal parameters
ivqte lwage (college), q(0.5) dummy(black) continuous(exper) unordered(region) b(1) l(1)
*robustness check with undersmoothing
ivqte lwage (college), q(0.5) dummy(black) continuous(exper) unordered(region) b(0.5) l(0.5)

* CONDITIONAL IV QUANTILE TREATMENT EFFECTS (Abadie, Angrist and Imbens, 2002)
*estimation of the first set of optimal smoothing parameters by cross-validation
locreg nearc4, dummy(black) continuous(exper) unordered(region) b(0.5 0.8) l(0.8 1) generate(ps) logit
*generate the positive and negative weights
generate waai=1-college*(1-nearc4)/(1-ps)-(1-college)*nearc4/ps
*estimation of the second set of optimal smoothing parameters by cross-validation
locreg waai, dummy(black) continuous(exper lwage) unordered(region) b(0.5 0.8) l(0.8 1)
*estimation with these optimal parameters
ivqte lwage (college=nearc4), aai q(0.5) dummy(black) continuous(age fatheduc motheduc) unordered(region) b(0.8) l(0.8) v vb(0.8) vl(1) pb(0.8) pl(0.8)

/*********************************************************************************************************************************************************************************************************************************************************************
EXAMPLE 2: ABADIE, ANGRIST AND IMBENS (2002) DATASET
**************************************************************************************************************************************************************************************************************************************************/

* Data set of Abadie, Angrist and Imbens (2002). Can be found at http://econ-www.mit.edu/faculty/angrist/data/abangim02
use "C:/ado/personal/jtpa.dta", clear

global reg "highschool black hispanic married part_time classroom OJT_JSA age5 age4 age3 age2 age1 second_follow"

*CASE 1: CONDITIONAL QUANTILE TREATMENT EFFECTS UNDER EXOGENEITY

*replicate the point estimate of the median regression in table II of AAI using the official command qreg 
qreg earnings treatment $reg, q(0.5)

*replicate the same point estimate using ivqte
ivqte earnings $reg (treatment), q(0.5) variance
*Note that the variance estimates are still different from the AAI standard errors because of slightly different kernel and bandwidth.
*without the option aai ivqte use the kernel and bandwidth proposed by Koenker (2005)

*replicate the same point estimate and standard errors using ivqte and the same bandwidth choice as AAI
ivqte earnings (treatment=treatment), q(0.5) dummy($reg) variance aai

*CASE 2: CONDITIONAL QUANTILE TREATMENT EFFECTS UNDER ENDOGENEITY

*estimation using a linear estimator for the weights
ivqte earnings (treatment=assignment), q(0.5) dummy($reg) variance aai pb(2)

*In order to replicate the results of AAI table III, we estimate the propensity score and the weights by a series estimator
sum assignment
gen pi=r(mean)
generate kappa=1-treatment*(1-assignment)/(1-pi)-(1-treatment)*assignment/pi
forvalues i=1/5{
gen e`i'=earnings^`i'
gen de`i'=e`i'*treatment
}
reg kappa earnings treatment e2 e3 e4 e5 de1 de2 de3 de4 de5
predict kappa_pos
*finally, we give the estimated propensity scores and weights as input for ivqte
ivqte earnings (treatment=assignment) if kappa_pos>0, d($reg) q(0.5) variance aai what(kappa_pos) phat(pi)
