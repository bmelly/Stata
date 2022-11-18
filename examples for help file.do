clear
set seed 1234
*100 groups
set obs 100
gen group = _n
*binary treatment variable
gen treatment = rbinomial(1,0.5)
*group-level variable
gen g_var = rnormal()
*30 individuals in each group
expand 30
*individual level variable
gen i_var = rnormal()
*outcome
gen y = 1 + treatment + g_var + i_var + rnormal()*(1 + 0.2 * treatment)
*Estimate treatment effect
mdqr y i_var g_var treatment, group(group)
*True treatment effect
mata: 1:+invnormal((0.1, 0.25, 0.5, 0.75, 0.9)):*0.2
*Now estimate the results at 19 quantiles and plot the treatment effect
quietly mdqr y i_var g_var treatment, group(group) quantile(0.05(0.05)0.95)
plotprocess treatment

set seed 2345
clear
*100 groups
set obs 100
gen group = _n
*treatment, instrumental variable, and (unobservable) group effects 
matrix C = (1, .5, 0\ .5, 1, 0.5 \ 0, 0.5, 1)
drawnorm g_effect treatment instrument, corr(C)
*group-level variable
gen g_var = rnormal()
*30 individuals in each group
expand 30
*individual level variable
gen i_var = rnormal()
*outcome
gen y = 1 + treatment + g_var + i_var + g_effect + rnormal()*(1 + 0.2 * treatment)
*Estimate treatment effect
*Inconsistent without instrument:
mdqr y i_var g_var treatment, group(group)
*Consistent with the instrument
mdqr y i_var g_var (treatment = instrument), group(group)
*True treatment effect
mata: 1:+invnormal((0.1, 0.25, 0.5, 0.75, 0.9)):*0.2

*Example with many fixed effects in the second stage
***generate data
clear
set seed 12345
***50 states
set obs 50
generate state = _n
generate state_fe = 0.2 * rnormal()
generate first_treated = ceil(runiform()*20) 
***20 years
expand 20
bysort state: generate year = _n
generate year_fe = year * 0.001
generate treated = year >= first_treated
generate state_char = rnormal()
***30 individuals per state
expand 30
generate ind_char = rnormal()
generate y = 1 + treated + 0.5 * state_char + 0.5 * ind_char + state_fe + year_fe + rnormal() * (1 + treated * 0.2)
***True quantile treatment effects at the 0.1, 0.25, 0.5, 0.75 and 0.9 quantiles:
mata: 1:+invnormal((0.1, 0.25, 0.5, 0.75, 0.9)):*0.2
***first possibility: include manually the fixed effects in the regression (many coefficients!)
mdqr y state_char ind_char i.year i.state treated, group(state year)
***second possibility: use areg to absorb one category of fixed effects
mdqr y state_char ind_char i.year treated, group(state year) est_command(areg) est_opts(absorb(state))
***third possibility: use reghdfe (or another user-written command) to absorb all fixed effects
***to install this command type "ssc install reghdfe"
mdqr y state_char ind_char treated, group(state year) est_command(reghdfe) est_opts(absorb(state year))
*** parallel processing
parallel initialize
mdqr y state_char ind_char treated, group(state year) parallel est_command(reghdfe) est_opts(absorb(state year))

*Panel data example
**Generate data
clear
set seed 12345
set obs 100
generate gid = _n
generate alpha = invnormal(uniform())
expand 20
generate x = invnormal(uniform())
generate y = x + alpha + invnormal(uniform()) * (1 + x * 0.2)
*** True values
mata: 1:+invnormal((0.1, 0.25, 0.5, 0.75, 0.9)):*0.2
*** Estimation
xtset gid
xtmdqr y x, fe
xtmdqr y x, be
xtmdqr y x
