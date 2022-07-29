*Set the number of observations and the seed
set obs 1000
set seed 1

*Generate female, experience and lwage
generate female=(uniform()<0.5)
generate experience=4*invchi2(5,uniform())*(1-0.4*female)
generate lwage=2+experience*0.03 ///
              -0.1*female+invnormal(uniform())*(0.6-0.2*female)

*Decomposition of the median difference in lwage between men and women.
*We use 100 quantile regression in the first step and we don't estimate
*the standard errors (both are default values):
rqdeco3 lwage experience, by(female) quantile(0.5)
*Interpretation: the observed median gender gap is about 41%.
*About 29% is explained by gender differences in the distribution of experience. 
*About 12% is due to differing median coefficients between men and women. 
*The part due to the residuals is negligible.

*Decomposition of the 99 percentile differences in lwage between men 
*and women. We estimate 100 quantile regression in the first step and 
*estimate the standard errors by bootstraping the results 100 times. 
*We don't require a print of the results:
rqdeco3 lwage experience, by(female) qlow(0.01) qhigh(0.99) 	///
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
twoway (line total_differential quantile)                            ///
	  (line residuals quantile) (line median quantile)             ///
        (line characteristics quantile), 					   ///
	  legend(order(1 "Total differential"                          ///
	  2 "Effects of residuals" 3 "Effects of median coefficients"  ///
	  4 "Effects of characteristics"))					   ///
        title(Decomposition of differences in distribution)          ///
        ytitle(Log wage effects) xtitle(Quantile)                  

*Interpretation: the observed gap is increasing (in absolute value) 
*when we move up on the wage distribution. Actually, women are      
*positively discriminated at the bottom of the distribution. Both   
*the experience distribution, the median coefficients and the       
*residuals are responsible for this result.   				  
*The experience distribution is less dispersed for   			  
*women than for men. The residuals also are less dispersed for     
*women than for men. Quantitatively, the second effect is more     
*important than the first one. The difference between the median   
*coefficients explain the different location but not the dispersion
*of the distributions. 

*We prepare the data to plot a 95% confidence interval for the 
*effects of median coefficients:
generate lo_residuals=residuals-1.96*se_residuals
generate hi_residuals=residuals+1.96*se_residuals

*We plot the effects of coefficients with a 95% confidence interval:
twoway (rarea hi_residuals lo_residuals quantile, bcolor(gs13) legend(off))   ///
       (line residuals quantile,                                              ///
             title(Effects of residuals)       					 	///
             ytitle(Log wage effects) xtitle(Quantile))             
