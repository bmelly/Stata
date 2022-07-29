***************Examples for locreg***************

use http://fmwww.bc.edu/ec-p/data/wooldridge/CARD.dta
set seed 123
sample 200, count
*Estimation of the probabilities by local linear regression
locreg nearc4, generate(fitted1) bandwidth(0.5) lambda(0.5) continuous(exper fatheduc) dummy(black) unordered(region)
*Some of the fitted probabilities are negative:
sum fitted1
*Estimation of the probabilities by local logistic regression
locreg nearc4, generate(fitted2) bandwidth(0.5) lambda(0.5) continuous(exper fatheduc) dummy(black) unordered(region) logit
*by definition, all fitted probabilities are between 0 and 1
sum fitted2
*Estimation of the optimal smoothing parameters by cross-validation
locreg nearc4, bandwidth(0.2 0.5) lambda(0.5 0.8) continuous(exper fatheduc) dummy(black) unordered(region)
*the estimated optimal parameters are h = 0.2 and lambda = 0.8:
return list
*Simulataneous cross-validation of the parameters and estimation of the probabilities with the optimal parameters
locreg nearc4, generate(fitted3) bandwidth(0.2 0.5) lambda(0.5 0.8) continuous(exper fatheduc) dummy(black) unordered(region)

