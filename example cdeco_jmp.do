*Set the number of observations and the seed
clear
set obs 1000
set seed 1

*Generate female, experience and lwage
generate female=(uniform()<0.5)
generate experience=4*invchi2(5,uniform())*(1-0.4*female)
generate lwage=2+experience*0.03 ///
              -0.1*female+invnormal(uniform())*0.5

*Decomposition of the median difference in lwage between men and women.
*First using quantile regression
cdeco_jmp lwage experience, group(female) quantile(0.5) noboot
*Second using the location model (OLS)
cdeco_jmp lwage experience, group(female) quantile(0.5) noboot method(loc)
*Interpretation: the results are similar using either method because 
*the residuals are independent from experience.
*The observed median gender gap is about 40%.
*About 28% is explained by gender differences in the distribution of experience. 
*About 11% is due to differing median coefficients between men and women. 
*The part due to the residuals is negligible.

*Plot of the decomposition results
cdeco_jmp lwage experience, group(female) quantile(0.01(0.01)0.99) noboot noprint
matrix total=e(total_difference)
matrix char=e(characteristics)
matrix coef=e(coefficients)
matrix resid=e(residuals)
matrix quant=e(quantiles)
mat total=total[1..99,1]
mat coef=coef[1..99,1]
mat char=char[1..99,1]
mat resid=resid[1..99,1]
svmat total
svmat coef
svmat resid
svmat char
svmat quant

*We plot the decomposition as a function of the quantile:
twoway (line total1 quant1) (line resid1 quant1) (line coef1 quant1) (line char1 quant1), legend(order(1 "Total differential" 2 "Effects of residuals" 3 "Effects of median coefficients" 4 "Effects of characteristics")) title(Decomposition of differences in distribution) ytitle(Log wage effects) xtitle(Quantile)                  
