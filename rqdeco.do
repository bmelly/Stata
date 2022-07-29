*Set the number of observations and the seed
set obs 1000
set seed 1

*Generate female, experience and lwage
generate female=(uniform()<0.5)
generate experience=4*invchi2(5,uniform())*(1-0.2*female)
generate lwage=2+experience*0.03 ///
              -0.1*female+invnormal(uniform())*(0.6-0.2*female)

*Decomposition of the median difference in lwage between men and women.
*We use 100 quantile regression in the first step and we don't estimate
*the standard errors (both are default values):
rqdeco lwage experience, by(female) quantile(0.5)
*Interpretation: the observed median gender gap is 31%.
*About 17% is explained by gender differences in the distribution of 
*experience. About 14% is due to differing coefficients between men 
*and women and can be interpreted as discrimination.

*Decomposition of the 99 percentile differences in lwage between men 
*and women. We estimate 100 quantile regression in the first step and 
*estimate the standard errors by bootstraping the results 100 times. 
*We don't require a print of the results:
rqdeco lwage experience, by(female) qlow(0.01) qhigh(0.99) ///
	 qstep(0.01) vce(boot) reps(100) noprint

*We can find the point estimates in the matrix r(results):
matrix list r(results)
*We can find the standard errors in the matrix r(se):
matrix list r(se)

*We prepare the data to plot the results:
matrix results=r(results)
matrix se=r(se)
svmat results, names(col)
svmat se, names(col)

*We plot the decomposition as a function of the quantile:
twoway (line total_differential quantile)                             ///
       (line characteristics quantile) (line coefficients quantile),  ///        
       title(Decomposition of differences in distribution)      	    ///
       ytitle(Log wage effects) xtitle(Quantile)            	    ///
	 legend(order(1 "Total differential"                            ///
	 2 "Effects of characteristics" 3 "Effects of coefficients"))
*Interpretation: the observed gap is increasing (in absolute value) 
*when we move up on the wage distribution. Actually, women are     
*positively discriminated at the bottom of the distribution. Both  
*the experience distribution and the coefficients are responsible  
*for this fact. The experience distribution is less dispersed for  
*women than for men. The residuals also are less dispersed for     
*women than for men. Quantitatively, the second effect is more     
*important than the first one. Looking at these results, we can    
*write that there is a glass ceiling effect for women: the         
*discrimination increases as we move up on the wage distribution.

*We prepare the data to plot a 95% confidence interval for the 
*effects of coefficients (discrimination):
generate lo_coef=coefficients-1.96*se_coefficients
generate hi_coef=coefficients+1.96*se_coefficients

*We plot the effects of coefficients with a 95% confidence interval:
twoway (rarea hi_coef lo_coef quantile, bcolor(gs13) legend(off))   ///
       (line coefficients quantile,                                 ///
       title(Effects of coefficients (discrimination))		  ///
       ytitle(Log wage effects) xtitle(Quantile))             


/********************************************************************
How to calculate the interdecile range (difference between the 
the 9th decile and the 1st decile) and its standard error
*********************************************************************/
*Decomposition of the 99 percentile differences in lwage between men 
*and women. We estimate 100 quantile regression in the first step and 
*estimate the standard errors by bootstraping the results 100 times.
*We save the bootstrap results in the file "C:/ado/personal/rqdeco_boot".
*We don't require a print of the results:
rqdeco lwage experience, by(female) qlow(0.01) qhigh(0.99) ///
	 qstep(0.01) vce(boot) reps(100) noprint 		     ///
	 saving("C:/ado/personal/rqdeco_boot", replace)

*The results are saved in the matrix r(results)
mat results=r(results)

*The 90-10 ranges can be calculated using this matrix:
*Q90-q10 for the fitted_1 distribution
sca q90q10_fitted1=results[90,2]-results[10,2]
*Q90-q10 for the counterfactual distribution
sca q90q10_counter=results[90,3]-results[10,3]
*Q90-q10 for the fitted_0 distribution
sca q90q10_fitted0=results[90,4]-results[10,4]

*Decomposition of the q90-q10 range:
*total difference
sca q90q10_tot=q90q10_fitted1-q90q10_fitted0
*characteristics
sca q90q10_coef=q90q10_fitted1-q90q10_counter
*coefficients
sca q90q10_char=q90q10_counter-q90q10_fitted0
sca dir

*Standard errors
*open the bootstrap results
preserve
use "C:\ado\Personal\rqdeco_boot.dta", clear
*generate the q90q10 ranges for each bootstrap draw
generate q90_q10_fitted1= X1_C1_90- X1_C1_10
generate q90_q10_fitted0= X0_C0_90- X0_C0_10
generate q90_q10_counter= X1_C0_90- X1_C0_10
*generate the effects for each bootstrap draw
generate q90_q10_tot=q90_q10_fitted1-q90_q10_fitted0
generate q90_q10_coef=q90_q10_fitted1-q90_q10_counter
generate q90_q10_char=q90_q10_counter-q90_q10_fitted0
*The simplest way to obtain the standard errors consists in taking
*the standard errors over the bootstrap draws. (Alternatively, confidence
*intervals may be obtain by the percentiles of the effects over the 
*draws)
sum q90_q10_fitted1 q90_q10_fitted0 q90_q10_counter q90_q10_tot q90_q10_coef q90_q10_char

