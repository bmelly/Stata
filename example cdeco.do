*Examples: unions
*Setup
webuse nlsw88, clear
gen lwage=log(wage)

*Point estimation, different conditional models
*Quantile regression estimator, no inference
cdeco lwage tenure ttl_exp grade, group(union) noboot

*Location estimator, no inference
cdeco lwage tenure ttl_exp grade, group(union) noboot method(loc)

*Location scale estimator, no inference
cdeco lwage tenure ttl_exp grade, group(union) noboot method(locsca) scale(tenure ttl_exp grade)

*Logit estimator, no inference
cdeco lwage tenure ttl_exp grade, group(union) noboot method(logit)

*Probit estimator, no inference
cdeco lwage tenure ttl_exp grade, group(union) noboot method(probit) 

*Linear probability model estimator, no inference 
cdeco lwage tenure ttl_exp grade, group(union) noboot method(lpm) 

*Cox estimator, no inference 
cdeco lwage tenure ttl_exp grade, group(union) noboot method(cox) 

*Censored quantile regression estimator, no inference  
preserve
replace lwage=1.2 if lwage<1.2  
generate censored=1.2  
cdeco lwage tenure ttl_exp grade, group(union) noboot method(cqr) censoring(censored)  
restore


*Inference
*Location estimator, inference 
cdeco lwage tenure ttl_exp grade, group(union) method(loc)  

*It's better to estimate the decomposition at more quantiles.
*To avoid overloading the screen, we do not print the results of the decomposition but only the tests
cdeco lwage tenure ttl_exp grade, group(union) method(logit) quantiles(0.01(0.01)0.99) noprintdeco

*The results of the decomposition can be found in the matrices returned by cdeco
ereturn list
*We plot the results
*First all elements of the decomposition
mat tot=e(total_difference)
mat char=e(characteristics)
mat coef=e(coefficients)
mat quant=e(quantiles)
svmat tot 
svmat char 
svmat coef
svmat quant
twoway (line tot1 quant1) (line char1 quant1) (line coef1 quant1), ytitle(Quantile Effect) xtitle(Quantile) legend(order(1 "Total difference" 2 "Effects of characteristics" 3 "Effects of coefficients"))
*We plot also the effects of coefficients ("discrimination") with a 95% pointwise confidence interval and a 95% functional confidence interval
gen coef_point_lb=coef1-1.96*coef2
gen coef_point_ub=coef1+1.96*coef2
twoway (rarea coef3 coef4 quant1, bcolor(gs5)) (rarea coef_point_lb coef_point_ub quant1, bcolor(gs10)) (line coef1 quant, lcolor(black)), xtitle("Quantile") ytitle("Quantile Effect") title(Effects of coefficients) legend(order(3 "Point estimates" 1 "Uniform 95% confidence bands" 2 "Pointwise 95% confidence intervals") rows(3))





