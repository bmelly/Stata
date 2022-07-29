*Examples: Engel curves, counterfactual option
*Setup
use "C:\ado\Personal\engel.dta", clear
summarize income
generate counter_income = r(mean) + 0.75*(income-r(mean))

*Quantile regression estimator, no inference
counterfactual foodexp income, counterfactual(counter_income) noboot

*Location shift estimator, no inference
counterfactual foodexp income, counterfactual(counter_income) method(loc) noboot

*Logit estimator, no inference
counterfactual foodexp income, counterfactual(counter_income) method(logit) noboot

*Location shift estimator, inference
counterfactual foodexp income, counterfactual(counter_income) method(loc)
counterfactual foodexp income, counterfactual(counter_income) method(loc) first(0.05) last(0.95) cons_test(-50 50)


*Examples: unions, group option

*Setup
webuse nlsw88, clear
gen lwage=log(wage)

*Quantile regression estimator, no inference
counterfactual lwage tenure ttl_exp grade, group(union) noboot

*Location estimator, no inference
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(loc)

*Location scale estimator, no inference
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(locsca) scale(tenure ttl_exp grade)

*Logit estimator, no inference
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(logit)

*Probit estimator, no inference
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(probit) 

*Linear probability model estimator, no inference 
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(lpm) 

*Cox estimator, no inference 
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(cox) 

*Location estimator, inference 
counterfactual lwage tenure ttl_exp grade, group(union) method(loc)  

*Censored quantile regression estimator, no inference  
replace lwage=1.2 if lwage<1.2  
generate censored=1.2  
counterfactual lwage tenure ttl_exp grade, group(union) noboot method(cqr) censoring(censored)  

