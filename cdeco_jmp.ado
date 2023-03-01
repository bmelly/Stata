*Version 1.0.3 01mar2023

program cdeco_jmp, eclass
	version 9.2
	capt findfile lmoremata.mlib
	if _rc {
      	di as error "-moremata- is required; type {stata ssc install moremata} and restart Stata."
		error 499
	}
	syntax varlist [if] [in] [pweight/], Group(varname) [Method(string) Quantiles(numlist >0 <1 sort) NReg(real 100) Reps(integer 100) Level(cilevel) First(real 0.1) Last(real 0.9) noboot noprint noprintdeco noprinttest SCale(varlist) SAVing(string) CONS_test(string) beta(real 0.9995) small(real 0.00001) max_it(real 50) Censoring(varname) Firstc(real 0.1) Secondc(real 0.05) NSteps(integer 3) RIght est_opts(string)]
	local nreg=round(`nreg')
	if `nreg'<1{
		dis as error "The option nreg must be a strictly positive integer."
		exit
	}
	if "`method'"==""{
		local method "qr"
	}
	if "`method'"!="qr" & "`method'"!="loc" & "`method'"!="locsca" & "`method'"!="cqr"{
		dis as error "The selected method has not been implemented"
		exit
	}
	if "`method'"=="cqr"{
		if `nsteps'<3{
			dis in red "The options nsteps must be at least 3"
			exit
		}
		if "`censoring'"==""{
			dis as error "The option censoring must be provided to use the censored quantile regression estimator."
			exit
		}
	}
	if "`method'" == "qr"{
		*check that qrprocess is installed
		capt findfile qrprocess.ado
		if _rc {
			di as error `"-qrprocess- is required; type {stata "net install qrprocess, from(https://raw.githubusercontent.com/bmelly/Stata/main/)"}"'
			error 499
		}
	}
	if "`method'"=="logit" | "`method'"=="probit" | "`method'"=="lpm" | "`method'"=="cloglog"{
		*check that drprocess is installed
		capt findfile drprocess.ado
		if _rc {
			di as error `"-drprocess- is required; type {stata "net install drprocess, from(https://raw.githubusercontent.com/bmelly/Stata/main/)"}"'
			error 499
		}
	}
	marksample touse
	markout `touse' `group' `scale' `censoring' 
	gettoken dep varlist : varlist
	if "`exp'"==""{
		tempvar exp
		gen `exp'=1
	}
	quietly _rmcollright `varlist' [pw=`exp'] if `touse'
	local varlist `r(varlist)'
	if "`method'"=="cox"{
		quietly sum `dep' if `touse'
		if r(min)<0{
			dis in red "Only positive dependent variables allowed when the cox method is used."
			exit
		}
	}
	quietly tab `group' if `touse'
	if r(r)!=2{
		di in red "The group variable, `group', is not binary"
		exit
	}
	quietly sum `group'  if `touse'
	if r(min)!=0 | r(max)!=1{
		di in red "The group variable, `group', is not a 0/1 variable"
		exit
	}
	cdeco_int `dep' `varlist' `if' `in' [pweight=`exp'], group(`group') method(`method') quantiles(`quantiles') nreg(`nreg') scale(`scale') beta(`beta') small(`small') max_it(`max_it') censoring(`censoring') firstc(`firstc') secondc(`secondc') nsteps(`nsteps') `right' est_opts(`est_opts')
	tempname results coef0 coef1 covariance tot_dif char coef resid miss0 miss1 test_tot test_char test_coef test_resid test_miss0 test_miss1 disto0 disto1 dist0 dist1 distc1 distc2 obs0 obs1 quants constest ref nreg0 nreg1
	sca `nreg0'=r(nreg0)
	sca `nreg1'=r(nreg1)
	matrix `quants'=r(quants)
	local obs=r(obs)
	sca `obs0'=r(obs0)
	sca `obs1'=r(obs1)
	mata: st_matrix("`coef0'", qte_cov_coef0)
	mata: st_matrix("`coef1'", qte_cov_coef1)
	if "`boot'"==""{
		if "`cons_test'"=="" | "`constest'"=="0" {
			mat `constest'=0
		}
		else{
			tokenize `cons_test'
			local i=1
			while "``i''"!=""{
				if ``i''!=0{
					mat `constest'=(nullmat(`constest'))\(``i'')
				}
				local i=`i'+1
			}
		}
		if `"`saving'"'!="" {
			_prefix_saving `saving'
			local saving `"`s(filename)'"'
			if "`double'" == "" {
				local double `"`s(double)'"'
			}
			if "`double'" == "" {
				local double "double"
			}
			local every	`"`s(every)'"'
			if "`every'"==""{
				local every=1
			}
			local replace `"`s(replace)'"'
		}
		else local every=999999
		mata: qte_cov_defo0=qte_cov_obs0
		mata: qte_cov_defo1=qte_cov_obs1
		mata: qte_cov_def0=qte_cov_fitted0
		mata: qte_cov_def1=qte_cov_fitted1
		mata: qte_cov_defc1=qte_cov_counter1
		mata: qte_cov_defc2=qte_cov_counter2
		mata: qte_cov_boot=J(0,rows(qte_cov_quant)*6,.)
		preserve
		local actual_more=c(more)
		set more off
		di in gr "(bootstrapping " _c
		forvalues i=1/`reps'{
			bsample
			capture cdeco_int `dep' `varlist' `if' `in' [pweight=`exp'], group(`group') method(`method') quantiles(`quantiles') nreg(`nreg') scale(`scale') beta(`beta') small(`small') max_it(`max_it') censoring(`censoring') firstc(`firstc') secondc(`secondc') nsteps(`nsteps') `right'  est_opts(`est_opts')
			if _rc == 0 {
				mata: qte_cov_boot=qte_cov_boot\(qte_cov_obs0,qte_cov_obs1,qte_cov_fitted0,qte_cov_fitted1,qte_cov_counter1,qte_cov_counter2)
				di in gr "." _c
			}
			else{
				dis in red "x" _continue
			}
			if round(`i'/`every')==(`i'/`every'){
				drop _all
				mata: st_addobs(rows(qte_cov_boot))
				mata: idx = st_addvar(st_local("double"), st_tempname(cols(qte_cov_boot)))
				mata: st_store(.,idx,qte_cov_boot)
				if `i'==1{	
					quietly save `saving', `replace'
				}
				else{ 
					quietly save `saving', replace
				}
			}
			restore, preserve
		}
		set more `actual_more'
		di in gr ")"
		mata: ev_boot(qte_cov_boot, qte_cov_defo0, qte_cov_defo1, qte_cov_def0, qte_cov_def1, qte_cov_defc1, qte_cov_defc2, qte_cov_quant, "`constest'", `last', `first', `level', "`results'", "`covariance'", "`tot_dif'", "`char'", "`coef'", "`resid'", "`miss0'", "`miss1'", "`test_tot'", "`test_char'", "`test_coef'", "`test_resid'", "`test_miss0'", "`test_miss1'", "`disto0'", "`disto1'", "`dist0'", "`dist1'", "`distc1'", "`distc2'")
		mata: mata drop qte_cov_defo0 qte_cov_defo1 qte_cov_boot qte_cov_obs0 qte_cov_obs1 qte_cov_def0 qte_cov_def1 qte_cov_defc1 qte_cov_defc2 qte_cov_coef0 qte_cov_coef1 qte_cov_counter1 qte_cov_counter2 qte_cov_fitted0 qte_cov_fitted1 qte_cov_quant
		mat colnames `test_tot'=KS CMS
		mat colnames `test_char'=KS CMS
		mat colnames `test_coef'=KS CMS
		mat colnames `test_resid'=KS CMS
		mat colnames `test_miss0'=KS CMS
		mat colnames `test_miss1'=KS CMS
		local conam "no_effect"
		if `constest'[1,1]!=0{
			local nct=rowsof(`constest')
			forvalues i=1/`nct'{
				local conam "`conam' constant_`i'"
			}
		}
		local conam "`conam' constant_m stoch_dom_pos stoch_dom_neg"
		mat rownames `test_tot'=`conam'
		mat rownames `test_char'=`conam'
		mat rownames `test_coef'=`conam'
		mat rownames `test_resid'=`conam'
		mat `test_miss0'=`test_miss0'[1,1..2]
		mat `test_miss1'=`test_miss1'[1,1..2]
		local conam "no_difference"
		mat rownames `test_miss0'=`conam'
		mat rownames `test_miss1'=`conam'
	}
	else{
		mata: st_matrix("`results'",(qte_cov_fitted0:-qte_cov_fitted1,qte_cov_fitted0:-qte_cov_counter1,qte_cov_counter1:-qte_cov_counter2,qte_cov_counter2:-qte_cov_fitted1))
		mata: st_matrix("`covariance'",J(4*rows(qte_cov_quant),4*rows(qte_cov_quant),0))
		mata: st_matrix("`tot_dif'",(qte_cov_fitted0:-qte_cov_fitted1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`char'",(qte_cov_fitted0:-qte_cov_counter1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`coef'",(qte_cov_counter1:-qte_cov_counter2\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`resid'",(qte_cov_counter2:-qte_cov_fitted1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`miss0'",(qte_cov_obs0:-qte_cov_fitted0\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`miss1'",(qte_cov_obs1:-qte_cov_fitted1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`disto0'",(qte_cov_obs0\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`disto1'",(qte_cov_obs1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`dist0'",(qte_cov_fitted0\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`dist1'",(qte_cov_fitted1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`distc1'",(qte_cov_counter1\J(3,rows(qte_cov_quant),.))')
		mata: st_matrix("`distc2'",(qte_cov_counter2\J(3,rows(qte_cov_quant),.))')
		mata: mata drop qte_cov_coef0 qte_cov_coef1 qte_cov_counter1 qte_cov_counter2 qte_cov_fitted0 qte_cov_fitted1 qte_cov_obs0 qte_cov_obs1 qte_cov_quant
	}
	local conam ""
	local nq=rowsof(`quants')
	forvalues i=1/`nq'{
		local conam "`conam' t_q`i'"
	}
	forvalues i=1/`nq'{
		local conam "`conam' x_q`i'"
	}
	forvalues i=1/`nq'{
		local conam "`conam' b_q`i'"
	}
	forvalues i=1/`nq'{
		local conam "`conam' r_q`i'"
	}
	mat colnames `results'=`conam'
	mat rownames `results'=`dep'
	mat colnames `covariance'=`conam'
	mat rownames `covariance'=`conam'
	local conam ""
	forvalues i=1/`nq'{
		local conam "`conam' q`i'"
	}
	mat rownames `tot_dif'=`conam'
	mat colnames `tot_dif'=total_difference point_se uniform_lb uniform_ub
	mat rownames `char'=`conam'
	mat colnames `char'=characteristics point_se uniform_lb uniform_ub
	mat rownames `coef'=`conam'
	mat colnames `coef'=coefficients point_se uniform_lb uniform_ub
	mat rownames `resid'=`conam'
	mat colnames `resid'=residuals point_se uniform_lb uniform_ub
	mat rownames `miss0'=`conam'
	mat colnames `miss0'=misspecification_0 point_se uniform_lb uniform_ub
	mat rownames `miss1'=`conam'
	mat colnames `miss1'=misspecification_1 point_se uniform_lb uniform_ub
	mat rownames `disto0'=`conam'
	mat colnames `disto0'=observed_0 point_se uniform_lb uniform_ub
	mat rownames `disto1'=`conam'
	mat colnames `disto1'=observed_1 point_se uniform_lb uniform_ub
	mat rownames `dist0'=`conam'
	mat colnames `dist0'=fitted_0 point_se uniform_lb uniform_ub
	mat rownames `dist1'=`conam'
	mat colnames `dist1'=fitted_0 point_se uniform_lb uniform_ub
	mat rownames `distc1'=`conam'
	mat colnames `distc1'=counterfactual1 point_se uniform_lb uniform_ub
	mat rownames `distc2'=`conam'
	mat colnames `distc2'=counterfactual2 point_se uniform_lb uniform_ub
	mat rownames `quants'=`conam'
	mat colnames `quants'=quantile
*Display the results
*header
	if "`print'"=="" & "`printdeco'"==""{
		dis
		dis as text _column(0) "Conditional model" _c
		if "`method'"=="qr"{
			di as result _column(43) "linear quantile regression"
		}
		if "`method'"=="cqr"{
			di as result _column(43) %-8.0g "linear censored quantile regression"
		}
		if "`method'"=="logit"{
			di as result _column(43) "logit"
		}
		if "`method'"=="probit"{
			di as result _column(43) "probit"
		}
		if "`method'"=="lpm"{
			di as result _column(43) %-8.0g "linear probability model"
		}
		if "`method'"=="loc"{
			di as result _column(43) "location model"
		}
		if "`method'"=="locsca"{
			di as result _column(43) "location scale model"
		}
		if "`method'"=="cox"{
			di as result _column(43) "cox duration model"
		}
		if "`method'"!="cox"{
			dis as text _column(0) "Number of regressions estimated" _c
			di as result _column(43) %-8.0f `nreg0'
		}
		di
		if "`boot'"==""{
			di as text "The variance has been estimated by bootstraping the results " as result `reps'  as text " times."
		}
		else{
			di as text "The variance has not been computed." _new "Do not turn the option boot off if you want to compute it."
		}
		dis
		dis as text _column(0) "No. of obs. in the reference group" _c
		dis as result _column(43) %-8.0f `obs0'
		dis as text _column(0) "No. of obs. in the counterfactual group" _c
		dis as result _column(43) %-8.0f `obs1'
		dis
		display as text _n "Differences between the observable distributions (based on the conditional model)" 
		dis as text "{hline 12}" "{c TT}" "{hline 68}"
		dis as text _column(13) "{c |}"  _column(15) %~10s "Quantile" _column(28) %~10s "Pointwise" _column(46) %~8s "Pointwise" _column(62) %~21s "Functional"
		dis as text _column(0) "Quantile" _column(13) "{c |}" _column(15) %~10s "effect" _column(28) %~10s "Std. Err." _column(41) %~8s "[`level'% Conf. Interval]" _column(62) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 12}" "{c +}" "{hline 68}"
		local k=1
		local nq=rowsof(`quants')
*Results for each quantile 
		while `k'<=`nq'{
			dis as text _column(0) %17s `quants'[`k',1] as text _column(13) "{c |}" _column(16) as result %8.0g (`tot_dif'[`k',1]) _column (28) as result %8.0g (`tot_dif'[`k',2]) _column (40) as result %8.0g (`tot_dif'[`k',1]-invnormal(0.5*(1+`level'/100))*`tot_dif'[`k',2]) _column(51) as result %8.0g (`tot_dif'[`k',1]+invnormal(0.5*(1+`level'/100))*`tot_dif'[`k',2]) _column(63) as result %8.0g (`tot_dif'[`k',3]) _column(73) as result %8.0g (`tot_dif'[`k',4]) 
			if `k'==`nq' {
				dis as text "{hline 12}" "{c BT}" "{hline 68}" 
			}
			local k=`k'+1
		}
		dis
		display as text _n "Effects of characteristics" 
		dis as text "{hline 12}" "{c TT}" "{hline 68}"
		dis as text _column(13) "{c |}"  _column(15) %~10s "Quantile" _column(28) %~10s "Pointwise" _column(46) %~8s "Pointwise" _column(62) %~21s "Functional"
		dis as text _column(0) "Quantile" _column(13) "{c |}" _column(15) %~10s "effect" _column(28) %~10s "Std. Err." _column(41) %~8s "[`level'% Conf. Interval]" _column(62) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 12}" "{c +}" "{hline 68}"
		local k=1
		local nq=rowsof(`quants')
*Results for each quantile 
		while `k'<=`nq'{
			dis as text _column(0) %17s `quants'[`k',1] as text _column(13) "{c |}" _column(16) as result %8.0g (`char'[`k',1]) _column (28) as result %8.0g (`char'[`k',2]) _column (40) as result %8.0g (`char'[`k',1]-invnormal(0.5*(1+`level'/100))*`char'[`k',2]) _column(51) as result %8.0g (`char'[`k',1]+invnormal(0.5*(1+`level'/100))*`char'[`k',2]) _column(63) as result %8.0g (`char'[`k',3]) _column(73) as result %8.0g (`char'[`k',4]) 
			if `k'==`nq' {
				dis as text "{hline 12}" "{c BT}" "{hline 68}" 
			}
			local k=`k'+1
		}
		dis
		display as text _n "Effects of coefficients" 
		dis as text "{hline 12}" "{c TT}" "{hline 68}"
		dis as text _column(13) "{c |}"  _column(15) %~10s "Quantile" _column(28) %~10s "Pointwise" _column(46) %~8s "Pointwise" _column(62) %~21s "Functional"
		dis as text _column(0) "Quantile" _column(13) "{c |}" _column(15) %~10s "effect" _column(28) %~10s "Std. Err." _column(41) %~8s "[`level'% Conf. Interval]" _column(62) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 12}" "{c +}" "{hline 68}"
		local k=1
		local nq=rowsof(`quants')
*Results for each quantile 
		while `k'<=`nq'{
			dis as text _column(0) %17s `quants'[`k',1] as text _column(13) "{c |}" _column(16) as result %8.0g (`coef'[`k',1]) _column (28) as result %8.0g (`coef'[`k',2]) _column (40) as result %8.0g (`coef'[`k',1]-invnormal(0.5*(1+`level'/100))*`coef'[`k',2]) _column(51) as result %8.0g (`coef'[`k',1]+invnormal(0.5*(1+`level'/100))*`coef'[`k',2]) _column(63) as result %8.0g (`coef'[`k',3]) _column(73) as result %8.0g (`coef'[`k',4]) 
			if `k'==`nq' {
				dis as text "{hline 12}" "{c BT}" "{hline 68}" 
			}
			local k=`k'+1
		}
		dis
		display as text _n "Effects of residuals" 
		dis as text "{hline 12}" "{c TT}" "{hline 68}"
		dis as text _column(13) "{c |}"  _column(15) %~10s "Quantile" _column(28) %~10s "Pointwise" _column(46) %~8s "Pointwise" _column(62) %~21s "Functional"
		dis as text _column(0) "Quantile" _column(13) "{c |}" _column(15) %~10s "effect" _column(28) %~10s "Std. Err." _column(41) %~8s "[`level'% Conf. Interval]" _column(62) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 12}" "{c +}" "{hline 68}"
		local k=1
		local nq=rowsof(`quants')
*Results for each quantile 
		while `k'<=`nq'{
			dis as text _column(0) %17s `quants'[`k',1] as text _column(13) "{c |}" _column(16) as result %8.0g (`resid'[`k',1]) _column (28) as result %8.0g (`resid'[`k',2]) _column (40) as result %8.0g (`resid'[`k',1]-invnormal(0.5*(1+`level'/100))*`resid'[`k',2]) _column(51) as result %8.0g (`resid'[`k',1]+invnormal(0.5*(1+`level'/100))*`resid'[`k',2]) _column(63) as result %8.0g (`resid'[`k',3]) _column(73) as result %8.0g (`resid'[`k',4]) 
			if `k'==`nq' {
				dis as text "{hline 12}" "{c BT}" "{hline 68}" 
			}
			local k=`k'+1
		}
	}
*Tests
	if "`boot'"=="" & "`print'"=="" & "`printtest'"==""{
		dis
		dis as text _n "Bootstrap inference on the counterfactual quantile processes"
		dis as text "{hline 51}" "{c TT}" "{hline 29}"
		dis as text _column(52) "{c |}"  _column(62) %~10s "P-values"
		dis as text _column(0) "Null-hypothesis" _column(52) "{c |}" _column(54) %~10s "KS-statistic" _column(69) %~10s "CMS-statistic"
		dis as text "{hline 51}" "{c +}" "{hline 29}"
		dis as text _column(0) %17s "Correct specification of the parametric model 0" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_miss0'[1,1]) _column(69) as result %8.0g (`test_miss0'[1,2]) 
		dis as text _column(0) %17s "Correct specification of the parametric model 1" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_miss1'[1,1]) _column(69) as result %8.0g (`test_miss1'[1,2]) 
		dis as text _column(0) %17s "Differences between the observable distributions" _column(52) "{c |}" 
		dis as text _column(5) %17s "No effect: QE(tau)=0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_tot'[1,1]) _column(69) as result %8.0g (`test_tot'[1,2]) 
		if `constest'[1,1]!=0{
			local nct=rowsof(`constest')
			forvalues i=1/`nct'{
				local temp=`constest'[`i',1]
				dis as text _column(5) %17s "Constant effect: QE(tau)=`temp' for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_tot'[1+`i',1]) _column(69) as result %8.0g (`test_tot'[1+`i',2]) 
			}
		}
		dis as text _column(5) %17s "Constant effect: QE(tau)=QE(0.5) for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_tot'[rowsof(`test_tot')-2,1]) _column(69) as result %8.0g (`test_tot'[rowsof(`test_tot')-2,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)>0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_tot'[rowsof(`test_tot')-1,1]) _column(69) as result %8.0g (`test_tot'[rowsof(`test_tot')-1,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)<0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_tot'[rowsof(`test_tot'),1]) _column(69) as result %8.0g (`test_tot'[rowsof(`test_tot'),2]) 
		dis as text _column(0) %17s "Effects of characteristics" _column(52) "{c |}" 
		dis as text _column(5) %17s "No effect: QTE(tau)=0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_char'[1,1]) _column(69) as result %8.0g (`test_char'[1,2]) 
		if `constest'[1,1]!=0{
			local nct=rowsof(`constest')
			forvalues i=1/`nct'{
				local temp=`constest'[`i',1]
				dis as text _column(5) %17s "Constant effect: QE(tau)=`temp' for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_char'[1+`i',1]) _column(69) as result %8.0g (`test_char'[1+`i',2]) 
			}
		}
		dis as text _column(5) %17s "Constant effect: QE(tau)=QE(0.5) for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_char'[rowsof(`test_tot')-2,1]) _column(69) as result %8.0g (`test_char'[rowsof(`test_tot')-2,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)>0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_char'[rowsof(`test_tot')-1,1]) _column(69) as result %8.0g (`test_char'[rowsof(`test_tot')-1,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)<0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_char'[rowsof(`test_tot'),1]) _column(69) as result %8.0g (`test_char'[rowsof(`test_tot'),2]) 
		dis as text _column(0) %17s "Effects of coefficients" _column(52) "{c |}" 
		dis as text _column(5) %17s "No effect: QE(tau)=0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_coef'[1,1]) _column(69) as result %8.0g (`test_coef'[1,2]) 
		if `constest'[1,1]!=0{
			local nct=rowsof(`constest')
			forvalues i=1/`nct'{
				local temp=`constest'[`i',1]
				dis as text _column(5) %17s "Constant effect: QTE(tau)=`temp' for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_coef'[1+`i',1]) _column(69) as result %8.0g (`test_coef'[1+`i',2]) 
			}
		}
		dis as text _column(5) %17s "Constant effect: QE(tau)=QE(0.5) for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_coef'[rowsof(`test_tot')-2,1]) _column(69) as result %8.0g (`test_coef'[rowsof(`test_tot')-2,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)>0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_coef'[rowsof(`test_tot')-1,1]) _column(69) as result %8.0g (`test_coef'[rowsof(`test_tot')-1,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)<0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_coef'[rowsof(`test_tot'),1]) _column(69) as result %8.0g (`test_coef'[rowsof(`test_tot'),2]) 
		dis as text _column(0) %17s "Effects of residuals" _column(52) "{c |}" 
		dis as text _column(5) %17s "No effect: QE(tau)=0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_resid'[1,1]) _column(69) as result %8.0g (`test_resid'[1,2]) 
		if `constest'[1,1]!=0{
			local nct=rowsof(`constest')
			forvalues i=1/`nct'{
				local temp=`constest'[`i',1]
				dis as text _column(5) %17s "Constant effect: QTE(tau)=`temp' for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_resid'[1+`i',1]) _column(69) as result %8.0g (`test_resid'[1+`i',2]) 
			}
		}
		dis as text _column(5) %17s "Constant effect: QE(tau)=QE(0.5) for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_resid'[rowsof(`test_tot')-2,1]) _column(69) as result %8.0g (`test_resid'[rowsof(`test_tot')-2,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)>0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_resid'[rowsof(`test_tot')-1,1]) _column(69) as result %8.0g (`test_resid'[rowsof(`test_tot')-1,2]) 
		dis as text _column(5) %17s "Stochastic dominance: QE(tau)<0 for all taus" as text _column(52) "{c |}" _column(54) as result %8.0g (`test_resid'[rowsof(`test_tot'),1]) _column(69) as result %8.0g (`test_resid'[rowsof(`test_tot'),2]) 
		dis as text "{hline 51}" "{c BT}" "{hline 29}"
	}
	if "`boot'" == "" {
		ereturn post `results' `covariance', dep(`dep') obs(`obs') esample(`touse')
		ereturn local vce "bootstrap"
	}
	else{
		ereturn post `results', dep(`dep') obs(`obs') esample(`touse')
		ereturn local vce "novar"
	}
	ereturn matrix total_difference=`tot_dif'
	ereturn matrix characteristics=`char'
	ereturn matrix coefficients=`coef'
	ereturn matrix residuals=`resid'
	ereturn matrix misspecification_0 =`miss0'
	ereturn matrix misspecification_1 =`miss1'
	ereturn matrix observed_0=`disto0'
	ereturn matrix observed_1 =`disto1'
	ereturn matrix fitted_0 =`dist0'
	ereturn matrix fitted_1 =`dist1'
	ereturn matrix counterfactual1 =`distc1'
	ereturn matrix counterfactual2 =`distc2'
	ereturn local cmd "counterfactual"
	ereturn local plotdeco "total_difference characteristics coefficients residuals"
	ereturn scalar obs0=`obs0'
	ereturn scalar obs1=`obs1'
	ereturn scalar nreg0=`nreg0'
	ereturn scalar nreg1=`nreg1'
	ereturn matrix quantiles=`quants'
	ereturn matrix coef0=`coef0'
	ereturn matrix coef1=`coef1'
	if "`boot'"==""{
		ereturn matrix test_tot=`test_tot'
		ereturn matrix test_char=`test_char'
		ereturn matrix test_coef=`test_coef'
		ereturn matrix test_resid=`test_resid'
		ereturn matrix test_miss0=`test_miss0'
		ereturn matrix test_miss1=`test_miss1'
	}
end

*Mata function doing the evaluation of the bootstrap results
version 9.2
mata void ev_boot(real matrix qte_cov_booti, qte_cov_obs0, qte_cov_obs1, qte_cov_def0, real rowvector qte_cov_def1, real rowvector qte_cov_defc1, real rowvector qte_cov_defc2, real colvector qte_cov_quant, string scalar contest, real scalar max, real scalar min, real scalar level, string scalar results, string scalar covariance, string scalar tot_dif, string scalar char, string scalar coef, string scalar resid, string scalar miss0, string scalar miss1, string scalar test_tot, string scalar test_char, string scalar test_coef, string scalar test_resid, string scalar test_miss0, string scalar test_miss1, string scalar disto0, string scalar disto1, string scalar dist0, string scalar dist1, string scalar distc1, string scalar distc2)
{
	nq=rows(qte_cov_quant)
	qte_cov_boot_o0=qte_cov_booti[.,1..nq]
	qte_cov_boot_o1=qte_cov_booti[.,(1..nq):+nq]
	qte_cov_boot_0=qte_cov_booti[.,(1..nq):+2*nq]
	qte_cov_boot_1=qte_cov_booti[.,(1..nq):+3*nq]
	qte_cov_boot_c1=qte_cov_booti[.,(1..nq):+4*nq]
	qte_cov_boot_c2=qte_cov_booti[.,(1..nq):+5*nq]
	st_matrix(results,(qte_cov_def0:-qte_cov_def1,qte_cov_def0:-qte_cov_defc1,qte_cov_defc1:-qte_cov_defc2,qte_cov_defc2:-qte_cov_def1))
	st_matrix(covariance,variance((qte_cov_boot_0:-qte_cov_boot_1,qte_cov_boot_0:-qte_cov_boot_c1,qte_cov_boot_c1:-qte_cov_boot_c2,qte_cov_boot_c2:-qte_cov_boot_1)))
	constest=st_matrix(contest)
	test_boot((qte_cov_boot_0:-qte_cov_boot_1),(qte_cov_def0:-qte_cov_def1),qte_cov_quant,constest,max,min,level,tot_dif,test_tot)
	test_boot((qte_cov_boot_0:-qte_cov_boot_c1),(qte_cov_def0:-qte_cov_defc1),qte_cov_quant,constest,max,min,level,char,test_char)
	test_boot((qte_cov_boot_c1:-qte_cov_boot_c2),(qte_cov_defc1:-qte_cov_defc2),qte_cov_quant,constest,max,min,level,coef,test_coef)
	test_boot((qte_cov_boot_c2:-qte_cov_boot_1),(qte_cov_defc2:-qte_cov_def1),qte_cov_quant,constest,max,min,level,resid,test_resid)
	test_boot((qte_cov_boot_o0:-qte_cov_boot_0),(qte_cov_obs0:-qte_cov_def0),qte_cov_quant,a=0,max,min,level,miss0,test_miss0)
	test_boot((qte_cov_boot_o1:-qte_cov_boot_1),(qte_cov_obs1:-qte_cov_def1),qte_cov_quant,a=0,max,min,level,miss1,test_miss1)
	test_boot(qte_cov_boot_o0,qte_cov_obs0,qte_cov_quant,a=.,max,min,level,disto0,"")
	test_boot(qte_cov_boot_o1,qte_cov_obs1,qte_cov_quant,a=.,max,min,level,disto1,"")
	test_boot(qte_cov_boot_0,qte_cov_def0,qte_cov_quant,a=.,max,min,level,dist0,"")
	test_boot(qte_cov_boot_1,qte_cov_def1,qte_cov_quant,a=.,max,min,level,dist1,"")
	test_boot(qte_cov_boot_c1,qte_cov_defc1,qte_cov_quant,a=.,max,min,level,distc1,"")
	test_boot(qte_cov_boot_c2,qte_cov_defc2,qte_cov_quant,a=.,max,min,level,distc2,"")
}

mata void test_boot(real matrix qte_cov_boot, real rowvector qte_cov_def, real colvector qte_cov_quant, real matrix constest, real scalar max, real scalar min, real scalar level, string scalar results, string scalar tests)
{
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-qte_cov_def#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	sel=(sum(qte_cov_quant:<min)+1)..sum(qte_cov_quant:<=max)
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	lb=qte_cov_def-seuqf*Kalpha
	ub=qte_cov_def+seuqf*Kalpha
	st_matrix(results,(qte_cov_def\seuqf\lb\ub)')
	if(missing(constest)==0){
		//test of no effect, KS
		if(constest==J(1,1,0)) nc=0 
		else nc=rows(constest)
		KSstat=max((qte_cov_def[sel]:^2:/Vuqf[sel]):^0.5)
		test=mean(Kmaxuqf:>KSstat)
		if(constest[1]!=0){
			for(i=1;i<=nc;i++){
				KSstat=max(((qte_cov_def[sel]:-constest[i]):^2:/Vuqf[sel]):^0.5)
				test=test\mean(Kmaxuqf:>KSstat)
			}
		}
		//test of no effects, CMS
		test=test,test
		Kuqf=((qte_cov_boot-qte_cov_def#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1)))
		Kmeanuqf=mean(Kuqf[.,sel]')
		CMSstat=mean((qte_cov_def[sel]:^2:/Vuqf[sel])')
		test[1,2]=mean(Kmeanuqf':>CMSstat)
		if(constest[1]!=0){
			for(i=1;i<=nc;i++){
				CMSstat=mean(((qte_cov_def[sel]:-constest[i]):^2:/Vuqf[sel])')
				test[i+1,2]=mean(Kmeanuqf':>CMSstat)
			}
		}
		//test of constant effect, KS
		test=test\J(3,2,.)
		median=ceil(length(qte_cov_def)/2)
		sel=((sum(qte_cov_quant:<min)+1)..(median-1)),((median+1)..sum(qte_cov_quant:<=max))
		qte_cov_def1=qte_cov_def[sel]:-qte_cov_def[median]
		qte_cov_boot1=qte_cov_boot[.,sel]:-qte_cov_boot[.,median]
		Vuqf=diagonal(variance(qte_cov_boot1))'
		Kuqf=((qte_cov_boot1-qte_cov_def1#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
		Kmaxuqf=rowmax(Kuqf)
		KSstat=max((qte_cov_def1:^2:/Vuqf):^0.5)
		test[nc+2,1]=mean(Kmaxuqf:>KSstat)
		//test of no effects, CMS
		Kuqf=((qte_cov_boot1-qte_cov_def1#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1)))
		Kmeanuqf=mean(Kuqf')
		CMSstat=mean((qte_cov_def1:^2:/Vuqf)')
		test[nc+2,2]=mean(Kmeanuqf':>CMSstat)
		//test of stochastic dominance, KS
		sel=(sum(qte_cov_quant:<min)+1)..sum(qte_cov_quant:<=max)
		qte_cov_boot1=qte_cov_boot[.,sel]-qte_cov_def[sel]#J(rows(qte_cov_boot),1,1)
		qte_cov_boot1=qte_cov_boot1:*(qte_cov_boot1:<=0)
		qte_cov_def1=qte_cov_def[sel]:*(qte_cov_def[sel]:<=0)
		Vuqf=diagonal(variance(qte_cov_boot))'[sel]
		Kuqf=(qte_cov_boot1:^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
		Kmaxuqf=rowmax(Kuqf)
		KSstat=max((qte_cov_def1:^2:/Vuqf):^0.5)
		test[nc+3,1]=mean(Kmaxuqf:>KSstat)
		//test of stochastic dominance, CMS
		Kuqf=((qte_cov_boot1-qte_cov_def1#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1)))
		Kmeanuqf=mean(Kuqf')
		CMSstat=mean((qte_cov_def1:^2:/Vuqf)')
		test[nc+3,2]=mean(Kmeanuqf':>CMSstat)
		//test of being stochastically dominated, KS
		qte_cov_boot1=qte_cov_boot[.,sel]-qte_cov_def[sel]#J(rows(qte_cov_boot),1,1)
		qte_cov_boot1=qte_cov_boot1:*(qte_cov_boot1:>=0)
		qte_cov_def1=qte_cov_def[sel]:*(qte_cov_def[sel]:>=0)
		Kuqf=(qte_cov_boot1:^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
		Kmaxuqf=rowmax(Kuqf)
		KSstat=max((qte_cov_def1:^2:/Vuqf):^0.5)
		test[nc+4,1]=mean(Kmaxuqf:>KSstat)
		//test of being stochastically dominated, CMS
		Kuqf=qte_cov_boot1:^2:/(Vuqf#J(rows(qte_cov_boot),1,1))
		Kmeanuqf=mean(Kuqf')
		CMSstat=mean((qte_cov_def1:^2:/Vuqf)')
		test[nc+4,2]=mean(Kmeanuqf':>CMSstat)
		st_matrix(tests,test)
	}
}

*Generic function, counterfactual distribution using group==0 to estimate and group==1 to predict
program cdeco_int, rclas
	version 9.2
	syntax varlist [if] [in] [pweight/] , Group(varname) [ Method(string) Quantiles(numlist >0 <1 sort) NReg(real 100) scale(varlist) beta(real 0.9995) small(real 0.00001) max_it(real 50) Censoring(varname) Firstc(real 0.1) Secondc(real 0.05) NSteps(integer 3) right est_opts(string)]
		marksample touse
		markout `touse' `counterfactual' `group'
		tempname quants results obs0 obs1 obs nreg0 nreg1
		if "`quantiles'"==""{
			local quantiles "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
		}
		tokenize "`quantiles'", parse(" ")
		local i=1
		while "`1'" != "" {
			matrix `quants'=nullmat(`quants')\(`1')
			mac shift 
			local i=`i'+1
		}
		gettoken dep varlist : varlist
		tempvar ref counter 
		sca `nreg0'=`nreg'
		sca `nreg1'=`nreg'
		quietly gen `ref'=`touse'
		quietly replace `ref'=0 if `group'==1
		quietly gen `counter'=`touse'
		quietly replace `counter'=0 if `group'==0
		quietly sum `ref'
		sca `obs0'=r(sum)
		quietly sum `counter'
		sca `obs1'=r(sum)
		sca `obs'=`obs0'+`obs1'
		mata: qte_cov_quant=st_matrix("`quants'")
		mata: qte_cov_obs0=mm_quantile(st_data(.,"`dep'","`ref'"),st_data(.,"`exp'","`ref'"),qte_cov_quant)'
		mata: qte_cov_obs1=mm_quantile(st_data(.,"`dep'","`counter'"),st_data(.,"`exp'","`counter'"),qte_cov_quant)'
		if "`method'"=="qr"{
			local qlow = 0.5/`nreg0'
			local qhigh = 1-0.5/`nreg0'
			local qstep = 1/`nreg0'
			qrprocess `dep' `varlist' if `ref' == 1 [pw=`exp'], vce(novar) qlow(`qlow') qhigh(`qhigh') qstep(`qstep') noprint `est_opts'
			mata: qte_cov_coef0 = st_matrix("e(quantiles)")' \ st_matrix("e(coefmat)")
			qrprocess `dep' `varlist' if `ref' == 1 [pw=`exp'], vce(novar) quantile(0.5) noprint `est_opts'
			mata: qte_cov_coef0 = (0.5\st_matrix("e(coefmat)")), qte_cov_coef0
			local qlow = 0.5/`nreg0'
			local qhigh = 1-0.5/`nreg0'
			local qstep = 1/`nreg0'
			qrprocess `dep' `varlist' if `counter' == 1 [pw=`exp'], vce(novar) qlow(`qlow') qhigh(`qhigh') qstep(`qstep') noprint `est_opts'
			mata: qte_cov_coef1 = st_matrix("e(quantiles)")' \ st_matrix("e(coefmat)")
			qrprocess `dep' `varlist' if `counter' == 1 [pw=`exp'], vce(novar) quantile(0.5) noprint `est_opts'
			mata: qte_cov_coef1 = (0.5\st_matrix("e(coefmat)")), qte_cov_coef1
			mata: qte_cov_fitted0=rqpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
			mata: qte_cov_fitted1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef1[.,2..cols(qte_cov_coef1)])
			mata: qte_cov_counter2=rqpred("`quants'","`varlist'","`exp'","`counter'",(qte_cov_coef1[1,2..cols(qte_cov_coef1)]\(qte_cov_coef1[2..rows(qte_cov_coef1),1]:+qte_cov_coef0[2..rows(qte_cov_coef0),2..cols(qte_cov_coef0)]:-qte_cov_coef0[2..rows(qte_cov_coef0),1])))
			mata: qte_cov_counter1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
		}
		if "`method'"=="cqr"{
			if "`right'"==""{	
				local right=0
			}
			else{
				local right=1
			}
			mata: qte_cov_coef0=est_cqr("`dep'", "`censoring'", "`varlist'", "`exp'", "`ref'", "`nreg0'", `firstc', `secondc', `nsteps', `right', `beta',`small',`max_it')
			mata: qte_cov_coef1=est_cqr("`dep'", "`censoring'", "`varlist'", "`exp'", "`counter'", "`nreg1'", `firstc', `secondc', `nsteps', `right', `beta',`small',`max_it')
			mata: qte_cov_fitted0=rqpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
			mata: qte_cov_fitted1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef1[.,2..cols(qte_cov_coef1)])
			mata: qte_cov_counter2=rqpred("`quants'","`varlist'","`exp'","`counter'",(qte_cov_coef1[1,2..cols(qte_cov_coef1)]\(qte_cov_coef1[2..rows(qte_cov_coef1),1]:+qte_cov_coef0[2..rows(qte_cov_coef0),2..cols(qte_cov_coef0)]:-qte_cov_coef0[2..rows(qte_cov_coef0),1])))
			mata: qte_cov_counter1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
		}
		if "`method'"=="loc"{
			mata: qte_cov_coef0=est_loc("`dep'","`varlist'","`exp'","`ref'","`nreg0'")
			mata: qte_cov_coef1=est_loc("`dep'","`varlist'","`exp'","`counter'","`nreg1'")
			mata: qte_cov_fitted0=rqpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
			mata: qte_cov_fitted1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef1[.,2..cols(qte_cov_coef1)])
			mata: qte_cov_counter2=rqpred("`quants'","`varlist'","`exp'","`counter'",(qte_cov_coef1[1,2..cols(qte_cov_coef1)]\(qte_cov_coef1[2..rows(qte_cov_coef1),1]:+qte_cov_coef0[2..rows(qte_cov_coef0),2..cols(qte_cov_coef0)]:-qte_cov_coef0[2..rows(qte_cov_coef0),1])))
			mata: qte_cov_counter1=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef0[.,2..cols(qte_cov_coef0)])
		}
		if "`method'"=="locsca"{
			if "`scale'"==""{
				local scale "`varlist'"
			}
			mata: qte_cov_coef0=est_locsca("`dep'","`varlist'","`scale'","`exp'","`ref'","`nreg0'")
			mata: qte_cov_coef1=est_locsca("`dep'","`varlist'","`scale'","`exp'","`counter'","`nreg1'")
			mata: qte_cov_fitted0=lspred("`quants'","`varlist'","`scale'","`exp'","`ref'",qte_cov_coef0)
			mata: qte_cov_fitted1=lspred("`quants'","`varlist'","`scale'","`exp'","`counter'",qte_cov_coef1)
			mata: qte_cov_counter2=lspred("`quants'","`varlist'","`scale'","`exp'","`counter'",(qte_cov_coef1[1..temp_lreg]\qte_cov_coef0[(temp_lreg+1)..rows(qte_cov_coef0)]))
			mata: qte_cov_counter1=lspred("`quants'","`varlist'","`scale'","`exp'","`counter'",qte_cov_coef0)
			mata: mata drop temp_lreg
		}
		return scalar obs=`obs'
		return scalar obs0=`obs0'
		return scalar obs1=`obs1'
		return matrix quants `quants'
		return scalar nreg0=`nreg0'
		return scalar nreg1=`nreg1'
end

*Mata function doing the unconditional estimation using locsca
version 9.2
mata real matrix lspred(string scalar quants, string scalar varlist, string scalar scale1, string scalar weights, string scalar touse, real matrix Coef1)
{
	Quants=st_matrix(quants)
	Reg=st_data(.,tokens(varlist),touse)
	Reg=Reg,J(rows(Reg),1,1)
	scale=st_data(.,tokens(scale1),touse)
	scale=scale,J(rows(scale),1,1)
	Wei=st_data(.,weights,touse)
	beta=Coef1[1..cols(Reg)]
	betas=Coef1[(cols(Reg)+1)..(cols(Reg)+cols(scale))]
	resid=Coef1[(cols(Reg)+cols(scale)+1)..length(Coef1)]
	pred=Reg*beta
	predsca=exp(scale*betas):^0.5
	pred=pred#J(1,length(resid),1)+predsca#resid'
	RQ_deco_ReT=mm_quantile(vec(pred),vec(Wei#J(1,length(resid),1)),Quants)'
	return(RQ_deco_ReT)
}

*Mata function doing the unconditional estimation using qr
version 9.2
mata real matrix rqpred(string scalar quants, string scalar varlist, string scalar weights, string scalar touse, real matrix Coef1)
{
	Quants=st_matrix(quants)
	Reg=st_data(.,tokens(varlist),touse)
	Wei=st_data(.,weights,touse)
	wq=Coef1[1,.]
	wq=(0,(wq[1..(cols(wq)-1)]+wq[2..cols(wq)]):/2,1)
	wq=wq[2..cols(wq)]-wq[1..(cols(wq)-1)]
	Coef=Coef1[2..rows(Coef1),.]
	Pred=cross(Reg'\J(1,rows(Reg),1),Coef)
	RQ_deco_ReT=mm_quantile(vec(Pred),vec(Wei*wq),Quants)'
	return(RQ_deco_ReT)
}

*Mata function doing the conditional estimation using qr
mata real matrix est_qr(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar nquant1, real scalar beta, real scalar small, real scalar max_it)
{
	y=st_data(.,dep,touse)
	x=st_data(.,tokens(reg),touse)
	n=rows(x)
	x=x,J(n,1,1)
	k=cols(x)
	w=st_data(.,weight,touse)
	nquant=st_numscalar(nquant1)
	coef=J(cols(x)+1,nquant+1,.)
	coef[1,.]=0.5,(0.5/nquant:+(0..(nquant-1)):/nquant)
	for (i=1; i<=(nquant+1); i++) {
		coef[2..(k+1),i]=rq_fnm(x, y, w, coef[1,i], beta, small, max_it, dep, reg, weight, touse)
	}
	return(coef)
}

*interior QR
mata real vector rq_fnm(real matrix X, real colvector dep, real colvector weight, real scalar p, real scalar beta, real scalar small, real scalar max_it, string scalar depo, string scalar rego, string scalar weighto, string scalar touse)
{
	weight=weight:/mean(weight)
	n=rows(X)
	u=J(n,1,1)
	a=(1-p):*u
	A=X'
	c=-dep'
	b=X'a
	x=a
	s = u - x
	s=s:*weight
	x=x:*weight
	y = (invsym(cross(A', A'))*cross(A',c'))'
	r = c - y * A
	r = r :+ 0.001 * (r :== 0)
	z = r :* (r :> 0)
	w = z - r
	gap = c * x - y * b + w * u
	it = 0
	while( gap > small & it < max_it){
 	   	it = it + 1
   	 	q=1:/(z':/x+w':/s)
   	 	r = z - w
		AQ = (A':*sqrt(q))'
		rhs=r':*sqrt(q)
		dy = (invsym(cross(AQ', AQ'))*cross(AQ',rhs))'
		dx = q :* (dy * A - r)'    
		ds = -dx
		dz = -z :* (1 :+ dx :/ x)'
		dw = -w :* (1 :+ ds :/ s)'
		fx = bound(x', dx')'
		fs = bound(s', ds')'
		fw = bound(w, dw)
		fz = bound(z, dz)
		fp = rowmin((fx, fs))
		fd = colmin((fw, fz))
		fp = min((beta * fp\ 1))
		fd = min((beta * fd, 1))
		if(min((fp, fd)) < 1){
 			mu = z * x + w * s
			g = (z + fd * dz) * (x + fp * dx) + (w + fd * dw) * (s + fp * ds)
			mu = mu * (g / mu) ^3 / ( 2* n)
			dxdz = dx :* dz'
			dsdw = ds :* dw'
			xinv = 1 :/ x
			sinv = 1 :/ s
			xi = mu * (xinv - sinv)
			rhs = rhs + sqrt(q):*(dxdz - dsdw - xi)
			dy = (invsym(cross(AQ', AQ'))*cross(AQ',rhs))'
			dx = q :* (A' dy' + xi - r' -dxdz + dsdw)
			ds = -dx
			dz = mu * xinv' - z - xinv' :* z :* dx' - dxdz'
			dw = mu * sinv' - w - sinv' :* w :* ds' - dsdw'		
			fx = bound(x', dx')'
			fs = bound(s', ds')'
			fw = bound(w, dw)
			fz = bound(z, dz)
			fp = rowmin((fx, fs))
			fd = colmin((fw, fz))		
			fp = min((beta * fp\ 1))
			fd = min((beta * fd, 1))
		}    
		x = x + fp * dx
		s = s + fp * ds
		y = y + fd * dy
		w = w + fd * dw
		z = z + fd * dz
		gap = c * x - y * b + w * u
	}
	y=-y'
	if(missing(y)>0){
		stata("qreg "+depo+" "+rego+" [pweight="+weighto+"] if "+touse+" ,quantile("+strofreal(p)+")",1)
		y=st_matrix("e(b)")'
	}
	return(y)
}

*internal function for interior QR
mata real vector bound(real vector x, real vector dx)
{
	b = J(1,length(x),1e20)
	f=select(1..length(x),dx:<0)
	b[f] = -x[f] :/ dx[f]
	return(b)
}

*Mata function doing the conditional estimation using location model
version 9.2
mata real matrix est_loc(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	coef=J((cols(reg)+1),(nreg+1),.)
	coef[1,.]=0.5,(0.5/nreg:+(0..(nreg-1)):/nreg)
	beta=invsym(cross(reg,wei,reg))*cross(reg,wei,dep)
	resid=dep:-reg*beta
	coef[2..rows(coef),.]=beta#J(1,(nreg+1),1)
	coef[(cols(reg)+1),2..(nreg+1)]=coef[(cols(reg)+1),2..(nreg+1)]+mm_quantile(resid,wei,coef[1,2..(nreg+1)])
	return(coef)
}

*Mata function doing the conditional estimation using location model
version 9.2
mata real matrix est_locsca(string scalar dep1, string scalar reg1, string scalar scale1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	external real scalar temp_lreg
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	scale=st_data(.,tokens(scale1),touse)
	scale=scale,J(rows(scale),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	coef=J(cols(reg)+cols(scale),nreg,.)
	beta=invsym(cross(reg,wei,reg))*cross(reg,wei,dep)
	resid=dep:-reg*beta
	betas=invsym(cross(scale,wei,scale))*cross(scale,wei,log(resid:^2))
	predsca=exp(scale*betas)
	coef=beta\betas\mm_quantile(resid:/(predsca:^0.5),wei,(0.5/nreg:+(0..(nreg-1)):/nreg))'
	temp_lreg=rows(beta)
	return(coef)
}

*Mata function doing the conditional estimation using cqr
mata real matrix est_cqr(string scalar dep, string scalar censoring, string scalar reg, string scalar weight, string scalar touse, string scalar nquant1, real scalar firstc, real scalar secondc, real scalar nsteps, real scalar right, real scalar beta, real scalar small, real scalar max_it)
{
	nquant=st_numscalar(nquant1)
	c=st_data(.,censoring,touse)
	y=st_data(.,dep,touse)
	x=st_data(.,tokens(reg),touse)
	x=x,J(rows(x),1,1)
	w=st_data(.,weight,touse)
	quants=0.5\(0.5/nquant:+(0..(nquant-1)):/nquant)'
	if(right==0) coef=est_cqrl(y, c, x, w, quants, firstc, secondc, nsteps, beta, small, max_it, dep, reg, weight, touse)
	if(right==1){
		y=-y
		c=-c	
		coef=-est_cqrl(y, c, x, w, 1:-quants, firstc, secondc, nsteps, beta, small, max_it, dep, reg, weight, touse)
		coef[1,.]=quants
	}
	return(coef)
}

mata real matrix est_cqrl(real colvector y, real colvector c, real matrix x, real colvector w, real colvector quants, real scalar c1, real scalar c2, real scalar nsteps, real scalar beta, real scalar small, real scalar max_it, string scalar dep, string scalar reg, string scalar weight, string scalar touse)
{
	coef=J(cols(x)+1,rows(quants),.)
	coef[1,.]=quants'
	ncensored=(y:>c)
	idx = st_addvar("double", st_tempname())
	st_store(.,idx,touse,ncensored)
	stata("logit "+st_varname(idx)+" "+reg+" [iweight="+weight+"] if "+touse+"==1, asis",1)
	temp=st_matrix("e(b)")'
	pred=invlogit(x*temp)
	nq=rows(quants)
	for(i=1; i<=nq; i++) {
		delta1=mm_quantile(select(pred,pred:>(1-quants[i])),select(w,pred:>(1-quants[i])),c1)
		select1=(pred:>=delta1)
		temp=rq_fnm(select(x,select1),select(y,select1),select(w,select1),quants[i],beta,small,max_it, dep, reg, weight, touse)
		step=3
		while(step<=nsteps){
			pred1=x*temp
			delta2=mm_quantile(select(pred1:-c,pred1:>c),select(w,pred1:>c),c2)
			select2=(pred1:>=delta2)
			temp=rq_fnm(select(x,select2),select(y,select2),select(w,select2),quants[i],beta,small,max_it, dep, reg, weight, touse)
			step=step+1
			select1=select2
		}
		coef[2..rows(coef),i]=temp
	}
	return(coef)
}
