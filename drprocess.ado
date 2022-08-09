*drprocess: distribution regression process

program drprocess, eclass byable(recall) sortpreserve
	version 9.2
*check that moremata is installed
	capt findfile lmoremata.mlib
	if _rc {
      	di as error "-moremata- is required; type {stata ssc install moremata}"
		error 499
	}
	if replay() {
		if "`e(cmd)'"!="drprocess" { 
			error 301 
		} 
		if _by() {
			error 190 
		}
        syntax [, Level(cilevel)]
		ereturn display, level(`level')
 	}
	else {
		if _caller() >= 11{
			version 11.1: syntax varlist(numeric fv) [if] [in] [pweight fweight/] [, Method(string) NDreg(integer -1234) Thresholds(numlist sort) Functional Vce(string) Level(cilevel) noprint]
		}
		else{
			syntax varlist(numeric) [if] [in] [pweight fweight/] [, Method(string) NDreg(integer -1234) Thresholds(numlist sort) Functional Vce(string) Level(cilevel) noprint]
		}
		local fvops = "`s(fvops)'" == "true" & _caller() >= 11
		if `fvops' {
			local vv: di "version " string(max(11,_caller())) ", missing: " 
		}
*sample definition
		marksample touse
		markout `touse' `cluster' `strata'
*weights
		local isweight=("`exp'"!="")
		if `isweight'==0{
			tempvar exp
			quietly gen byte `exp'=1 if `touse'
			local weight "pweight"
		}
*separate dependent from regressors
		gettoken dep varlist : varlist
		if `fvops' {
			`vv' _fv_check_depvar `dep'
			`vv'  quiet _rmcoll `varlist' [pw=`exp'] if `touse'
			local varlist "`r(varlist)'"
			`vv' fvexpand `varlist' if `touse'
			local varlistt "`r(varlist)'"
			local varlist ""
			local i=1
			foreach vn of local varlistt{						
				`vv' _ms_parse_parts `vn'
				if ~`r(omit)'{
					if r(type)=="variable"{
						local varlist `varlist' `vn'
						local variable_names `variable_names' `vn'
					}
					else{
						tempvar factor`i'
						quiet gen `factor`i''=`vn' if `touse'
						local varlist `varlist' `factor`i''
						local variable_names `variable_names' `vn'
					}
					local i=`i'+1
				}
			}
		}
		else{
			*check multicolinearity
			quiet _rmcoll `varlist' [pw=`exp'] if `touse'
			local varlist `r(varlist)'
			local varlistt `r(varlist)'
			local variable_names `varlist'
		}

*number of regressors
		if "`varlist'"!=""{
			local k=wordcount("`varlist'")+1
		} 
		else{
			local k=1
		}
*number of observations
		quietly sum `dep' if `touse'
		local obs=r(N)
		local min=r(min)
		local max=r(max)
	*temporary names
		tempname quants covariance convergence coefficients rawdev mindev coefmat nc ns
			
*Estimation of the variance, cleaning
		if strpos("`vce'",",")==0 dr_VCEParse `vce', weight(`isweight') funcint(`functional')
		else dr_VCEParse `vce' weight(`isweight') funcint(`functional')
		local variance `r(var)'
		local reps `r(reps)'
		local boot_method `r(boot_method)'
		local cluster `r(cluster)'
		local strata `r(strata)'
		local sub_size `r(sub_size)'
		local noreplacement `r(noreplacement)'
		local HC1 `r(HC1)'
		local nodfadjust `r(nodfadjust)'
	*cleaning of the values at which the dr will be estimated
		if "`thresholds'"==""{
			if `ndreg'==-1234 {
				if "`functional'"=="" {
					local ndreg = 1
				}
				else{
					local ndreg = 100
				}
			}
			if `ndreg'<1{
				dis as error "The option ndreg must contain a strictly positive integer."
				error 400
			}
			else if `ndreg' == 1 {
				mata: st_matrix("`quants'", mm_median(st_data(.,"`dep'","`touse'")))
			}
			else if `ndreg'>=`obs'{
				mata: st_matrix("`quants'", uniqrows(st_data(.,"`dep'","`touse'")))
			}
			else{			
				local trim=min(15*`k'/`obs',0.3)
				mata: st_matrix("`quants'",uniqrows(mm_quantile(st_data(.,"`dep'","`touse'"),st_data(.,"`exp'","`touse'"),range(`trim',1-`trim',(1-2*`trim')/(`ndreg'-1)))))
			}
		}
		else{
			tokenize "`thresholds'", parse(" ")
			local i=1
			while "`1'" != "" {
				if `1'<`min' | `1'>`max'{
					dis as error "DR cannot be estimated for thresholds below the minimum or above the maximum values taken by the dependent variable."
					error 400
				}
				matrix `quants'=nullmat(`quants')\(`1')
				mac shift 
				local i=`i'+1
			}
		}
*Algorithm, cleaning		
		local nq=rowsof(`quants')
		dr_methodParse `method'
		local method `r(method)'
		local max_it `r(max_it)'
		local onestep `r(onestep)'
		if "`functional'"=="functional"{
			if `nq'<2{
				dis as error "It is not possible to perform functional inference with only 1 distribution regression."
				error 400
			}
			else if `nq'<5{
				dis as error "Warning: we recommend to use many distribution regressions to perform functional inference."
			}
		}
*Names of the columns of the matrices, for later. Taken from sqreg.
		capture{
			tokenize "`variable_names'"
			local nq=rowsof(`quants')
			forvalues j=1/`nq' {
				local i=1
				while "``i''" != "" {
					local conams "`conams' ``i''"
					local eqnams "`eqnams' d`j'"
					local i = `i' + 1
				}
				local conams "`conams' _cons"
				local eqnams "`eqnams' d`j'"
				local ueqnams "`ueqnams' d`j'"
 			}
		}
		if _rc{
			dis in red "Due to the high number of regressors and thresholds the rows and columns of the matrices cannot be named correctly."
			tokenize "`variable_names'"
			local nq=rowsof(`quants')
			local conams ""
			local eqnams ""
			forvalues j=1/`nq' {
				local i=1
				while "``i''" != "" {
					local conams "`conams' x`i'"
					local eqnams "`eqnams' d`j'"
					local i = `i' + 1
				}
				local conams "`conams' c"
				local eqnams "`eqnams' d`j'"
				local ueqnams "`ueqnams' d`j'"
 			}
		}
*Estimation
		preserve
		quiet keep if `touse'
		tempvar tempdep tempbweight
		quietly generate `tempdep'=.
		quietly generate `tempbweight'=.
		if "`method'"=="lpm"{
			mata: dr_lpm("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", `ndreg'>=`obs', "`nc'", "`ns'")
		}
		else if "`method'"=="logit" & "`onestep'"==""{
			mata: dr_logit("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		else if "`method'"=="logit" & "`onestep'"=="onestep"{
			mata: dr_logit1("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		else if "`method'"=="probit" & "`onestep'"==""{
			mata: dr_probit("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		else if "`method'"=="probit" & "`onestep'"=="onestep"{
			mata: dr_probit1("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		else if "`method'"=="cloglog" & "`onestep'"==""{
			mata: dr_cloglog("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		else if "`method'"=="cloglog" & "`onestep'"=="onestep"{
			mata: dr_cloglog1("`dep'", "`varlist'", "`weight'", "`exp'", "`touse'", "`quants'", "`variance'", "`cluster'", "`strata'", "`coefmat'", "`covariance'", "`nc'", "`ns'", "`tempdep'")
		}
		mata: st_matrix("`coefficients'",vec(st_matrix("`coefmat'"))')
*degrees of freedom
		if "`cluster'"!=""{
			cap confirm scalar `nc'
			if _rc{
				quiet tab `cluster' if `touse'
				sca `nc'=r(r)
			}
		}
		if "`strata'"!=""{
			cap confirm scalar `ns'
			if _rc{
				quiet tab `strata' if `touse'
				sca `ns'=r(r)
			}
		}
		if "`nodfadjust'"=="nodfadjust" | ("`cluster'"=="" & "`strata'"==""){
			local df=`obs'-colsof(`coefficients')/`nq'
		}
		else if "`nodfadjust'"=="" & "`cluster'"!="" & "`strata'"==""{
			local df=`nc'-1
		}
		else if "`nodfadjust'"=="" & "`cluster'"=="" & "`strata'"!=""{
			local df=`obs'-`ns'-colsof(`coefficients')/`nq'+1
		}
		else if "`nodfadjust'"=="" & "`cluster'"!="" & "`strata'"!=""{
			local df=`nc'-`ns'-colsof(`coefficients')/`nq'+1
		}

*Bootstrap
		if "`variance'"=="bootstrap" | "`variance'"=="multiplier"{
			if "`boot_method'"!="subsampling"{
				local sub_size `obs'
			}
			tempname pointwise uniform tests
			local actual_more=c(more)
			set more off
			if "`strata'"!="" & "`cluster'"!=""{
				sort `strata' `cluster', stable
				mata: st_view(strata=., ., "`strata'","`touse'")
				mata: st_view(cluster=., ., "`cluster'","`touse'")
				mata: mm_panels(strata, Sinfo=., cluster, Cinfo=.)
			}
			else if "`strata'"!="" & "`cluster'"==""{
				sort `strata', stable
				mata: st_view(strata=., ., "`strata'","`touse'")
				mata: mm_panels(strata, Sinfo=.)
				mata: Cinfo=.
			}
			else if "`strata'"=="" & "`cluster'"!=""{
				sort `cluster', stable
				mata: st_view(cluster=., ., "`cluster'","`touse'")
				mata: mm_panels(cluster, Cinfo=.)
				mata: Sinfo=.
			}
			else{
				mata: Sinfo=.
				mata: Cinfo=.
			}
			di in gr "(bootstrapping " _c
			mata: Qreg2_Results=st_matrix("`coefficients'")
			if "`variance'"=="bootstrap"{
				if "`onestep'"=="" | "`method'"=="lpm"{
					mata: Qreg2_Boot=dr_boot("`dep'", "`varlist'", "`exp'", "`touse'", "`quants'", st_matrix("`coefmat'"), `reps', "`boot_method'", "`cluster'", "`strata'", `noreplacement', `sub_size', Sinfo, Cinfo, "`method'", `ndreg'>=`obs', "`tempdep'", "`tempbweight'")
				}
				else{
					mata: Qreg2_Boot=dr_boot_1step("`dep'", "`varlist'","`exp'", "`touse'", "`quants'", `reps', "`boot_method'", "`cluster'", "`strata'", `noreplacement', `sub_size', Sinfo, Cinfo, "`method'", st_matrix("`coefmat'"), "`tempdep'", "`tempbweight'")
				}
			}
			else{
				mata: Qreg2_Boot=dr_multiplier("`dep'", "`varlist'","`exp'", "`touse'", "`quants'", st_matrix("`coefmat'"), `reps', "`boot_method'", "`cluster'", "`strata'", `noreplacement', `sub_size', Sinfo, Cinfo, "`method'")
			}
			set more `actual_more'
			di in gr ")"
			mata: dr_ev_boot("`variance'"=="multiplier", Qreg2_Boot, Qreg2_Results, "`quants'", `level', "`covariance'", "`pointwise'", "`uniform'", "`tests'", "`functional'"=="functional",`obs', `sub_size', `df', "`nodfadjust'"=="nodfadjust")
			mata: mata drop Qreg2_Results Qreg2_Boot
			cap restore
			mata mata drop Cinfo Sinfo
		}
		else{
			restore
			if "`HC1'"=="HC1"{
				mat `covariance'=(`nc'/(`nc'-`ns'+1))*((`obs'-1)/(`obs'-colsof(`coefficients')/`nq'))*(`nc'/(`nc'-1))*`covariance'
			}
		}
		if "`print'"!="noprint"{
			dis
			dis as text _column(0) "Distribution regression"
			dis as text _column(0) "No. of obs." _c
			dis as result _column(20) %-8.0f `obs'
			dis as text _column(0) "Algorithm:" _c
			if "`method'"=="lpm" dis as result _column(20) %-8.0f "linear probability model."
			if "`method'"=="logit" & "`onestep'"=="" dis as result _column(20) %-8.0f "logistic regression."
			if "`method'"=="logit" & "`onestep'"=="onestep" dis as result _column(20) %-8.0f "one-step logistic regression."
			if "`method'"=="probit" & "`onestep'"=="" dis as result _column(20) %-8.0f "probit regression."
			if "`method'"=="probit" & "`onestep'"=="onestep" dis as result _column(20) %-8.0f "one-step probit regression."
			if "`method'"=="cloglog" & "`onestep'"=="" dis as result _column(20) %-8.0f "complementary log-log regression."
			if "`method'"=="cloglog" & "`onestep'"=="onestep" dis as result _column(20) %-8.0f "one-step complementary log-log regression."
			dis as text _column(0) "Variance:" _continue
			if "`variance'"=="novar" dis as result _column(20) %-8.0f "the variance has not been estimated."
			if "`variance'"=="bootstrap" & "`boot_method'"=="empirical" dis as result _column(20) %-8.0f "empirical bootstrap."
			if "`variance'"=="bootstrap" & "`boot_method'"=="weighted" dis as result _column(20) %-8.0f "weighted bootstrap (centered standard exponential weights)."
			if "`variance'"=="bootstrap" & "`boot_method'"=="subsampling" & "`replacement'"=="" dis as result _column(20) %-8.0f "subsampling with replacement (subsamples size is `sub_size')."
			if "`variance'"=="bootstrap" & "`boot_method'"=="subsampling" & "`replacement'"=="noreplacement" dis as result _column(20) %-8.0f "subsampling without replacement (subsamples size is `sub_size')."
			if "`variance'"=="robust" dis as result _column(20) %-8.0f "robust analytical."
			if "`variance'"=="multiplier" dis as result _column(20) %-8.0f "multiplier bootstrap."
			dis
			dis as text "{hline 13}" "{c TT}" "{hline 64}"
			if "`functional'"=="functional"{
				dis as text _column(0) %11s "`dep'" _column(14) "{c |}"  _column(19) %~10s "Coef." _column(29) %~10s "Pointwise" _column(45) %~8s "Pointwise" _column(60) %~21s "Functional"
				dis as text _column(14) "{c |}" _column(29) %~10s "Std. Err." _column(42) %~8s "[`level'% Conf. Int.]" _column(60) %~21s "[`level'% Conf. Int.]"
				dis as text "{hline 13}" "{c BT}" "{hline 64}"
				forvalues i=1/`nq'{
					dis as text _column(0) "y " as result round(`quants'[`i',1],0.001)  
	*				dis as text _column(0) "Ps. R2: " as result round(1-(`mindev'[1,`i']/`rawdev'[1,`i']),0.001) _column(14) "{c |}"
					forvalues j=1/`k'{
							dis as text _column(0) %11s abbrev(word("`conams'",`j'),12) _column(14) "{c |}" as result _column(18) %8.0g (`coefmat'[`j',`i']) _continue
							dis as result _column(29) %8.0g (`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
							dis as result _column(40) %8.0g (`pointwise'[(`i'-1)*`k'+`j',1]) _continue
							dis as result _column(50) %8.0g (`pointwise'[(`i'-1)*`k'+`j',2]) _continue
							dis as result _column(61) %8.0g (`uniform'[(`i'-1)*`k'+`j',1]) _continue
							dis as result _column(71) %8.0g (`uniform'[(`i'-1)*`k'+`j',2])
					}
					if `i'<`nq' dis as text "{hline 13}" "{c +}" "{hline 64}"
					else dis as text "{hline 13}" "{c BT}" "{hline 64}"
				}
				dis
				dis as text _n "Bootstrap inference on the distribution regression process"
				dis as text "{hline 51}" "{c TT}" "{hline 26}"
				dis as text _column(52) "{c |}"  _column(62) %~10s "P-values"
				dis as text _column(0) "Null-hypothesis" _column(40) %11s "Coef." _column(52) "{c |}" _column(55) %~10s "KS-stat." _column(69) %~10s "CMS-stat."
				dis as text "{hline 51}" "{c +}" "{hline 26}"
				local i=1
				foreach null in "No effect: beta(tau)=0 for all taus" "Constant effect: beta(tau)=B for all taus" "Positive effect: beta(tau)>=0 for all taus" "Negative effect: beta(tau)<=0 for all taus"{
					dis as text _column(0) "`null'" as text _column(52) "{c |}" 
					if `i'<5 local max_k=`k'
					else local max_k=`k'-1
					forvalues j=1/`max_k'{
						dis as text _column(30) %20s abbrev(word("`conams'",`j'),20) as text _column(52) "{c |}" as result _column(57) as result %4.3f `tests'[`j',(`i'-1)*2+1] _column(71) as result %4.3f `tests'[`j',(`i'-1)*2+2]
					}
					dis as text _column(30) %20s "all slopes"  as text _column(52) "{c |}" as result _column(57) as result %4.3f `tests'[`k'+1,(`i'-1)*2+1] _column(71) as result %4.3f `tests'[`k'+1,(`i'-1)*2+2]
					local i=`i'+1
				}
				dis as text "{hline 51}" "{c BT}" "{hline 26}"
			}
			else{
				if "`nodfadjust'"==""{
					dis as text _column(0) %11s "`dep'" _column(14) "{c |}"  _column(19) %~10s "Coef." _column(29) %~10s "Std. Err." _column(41) %~8s "t" _column(45) %~8s "P>|t|" _column(59) %~8s "[`level'% Conf. Interval]"
				}
				else{
					dis as text _column(0) %11s "`dep'" _column(14) "{c |}"  _column(19) %~10s "Coef." _column(29) %~10s "Std. Err." _column(41) %~8s "z" _column(45) %~8s "P>|z|" _column(59) %~8s "[`level'% Conf. Interval]"			
				}
				dis as text "{hline 13}" "{c BT}" "{hline 64}"
				forvalues i=1/`nq'{
	*				dis as text _column(0) "Quant. 0" as result round(`quants'[`i',1],0.001)  _column(14) "{c |}" 
					dis as text _column(0) "Threshold: " as result `quants'[`i',1] 
					forvalues j=1/`k'{
							dis as text _column(0) %11s abbrev(word("`conams'",`j'),12) _column(14) "{c |}" as result _column(17) %9.0g (`coefmat'[`j',`i']) _continue
							if "`variance'"!="novar"{
								dis as result _column(28) %9.0g (`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
								dis as result _column(40) %5.2f `coefmat'[`j',`i']/(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
								if "`nodfadjust'"==""{
									dis as result _column(49) %4.3f 2*ttail(`df',abs(`coefmat'[`j',`i']/(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5)) _continue
									dis as result _column(58) %9.0g `coefmat'[`j',`i']-invttail(`df',(100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
									dis as result _column(70) %9.0g `coefmat'[`j',`i']+invttail(`df',(100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5
								}
								else{
									dis as result _column(49) %4.3f 2-2*normal(abs(`coefmat'[`j',`i']/(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5)) _continue
									dis as result _column(58) %9.0g `coefmat'[`j',`i']+invnormal((100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
									dis as result _column(70) %9.0g `coefmat'[`j',`i']-invnormal((100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5							
								}
							}
							else dis
					}
					dis as text "{hline 13}" "{c BT}" "{hline 64}"
				}
			}
		}
		`vv' mat colnames `coefficients' = `conams'
		mat coleq `coefficients' = `eqnams'
		mat rownames `coefficients'=`dep'
		if "`variance'"!="novar"{
			`vv' mat rownames `covariance' = `conams'
			mat roweq `covariance' = `eqnams'
			`vv' mat colnames `covariance' = `conams'
			mat coleq `covariance' = `eqnams'
			ereturn post `coefficients' `covariance', depname(`dep') obs(`obs') dof(`df') esample(`touse')
		}
		else ereturn post `coefficients', depname(`dep') obs(`obs') dof(`df') esample(`touse')
		`vv' mat rownames `coefmat' = `variable_names' _cons
		`vv' mat colnames `coefmat' = `ueqnams'
		ereturn matrix coefmat `coefmat'	
		mat rownames `quants' = `ueqnams'
		mat colnames `quants' = "y"
		ereturn matrix thresholds `quants'
		ereturn local estat_cmd "drprocess_estat"
		ereturn local predict "drprocess_p"
		if "`boot_method'"=="subsampling"{
			ereturn local scalar `"`sub_size'"'
			if `noreplacement'==0 ereturn local replacement "with replacement"
			else ereturn local replacement "without replacement"
		}
		ereturn local bmethod `"`boot_method'"'
		ereturn local vce `"`variance'"'
		if "`onestep'"=="" ereturn local method `"`method'"'
		else ereturn local method `"`method', one-step estimator"'
		ereturn scalar df_m=`k'-1
		if "`cluster'"!=""{
			ereturn local clustvar `"`cluster'"'
			ereturn scalar N_clust=`nc'
		}
		if "`strata'"!=""{
			ereturn local stratvar `"`strata'"'
			ereturn scalar N_strat=`ns'
		}
		if `isweight'==1{
			ereturn local wtype `"`weight'"'
			ereturn local wexp `"`exp'"'
		}
		if "`variance'"=="bootstrap" | "`variance'"=="multiplier"{
			ereturn scalar rep=`reps'
		}
		ereturn local xvar `"`variable_names'"'
		ereturn local varlist `"`varlistt'"'
		ereturn local depvar `"`dep'"'
		if "`functional'"!=""{
			`vv' mat rownames `pointwise' = `conams'
			mat roweq `pointwise' = `eqnams'
			`vv' mat rownames `uniform' = `conams'
			mat roweq `uniform' = `eqnams'
			mat colnames `pointwise'="lower_bound upper_bound"
			mat colnames `uniform'="lower_bound upper_bound"
			tokenize "`varlist'"
			local i=1
			local conam ""
			while "``i''" != "" {
				local conam "`conam' ``i''"
				local i = `i' + 1
			}
			local conam "`conam' _cons"
			mat rownames `tests'=`conam' "all slopes"
			mat colnames `tests'=KS_0 CVM_0 KS_constant CVM_constant KS_pos CVM_pos KS_neg CVM_neg
			ereturn matrix pointwise `pointwise'
			ereturn matrix uniform `uniform'
			ereturn matrix tests `tests'
		}
		ereturn local cmdline `"drprocess `0'"'
		ereturn local title "Distribution regression"
		ereturn local cmd "drprocess"
		ereturn local plotprocess "e(thresholds) Threshold"
	}
end

*Mata function doing the evaluation of the bootstrap results
version 9.2
mata void dr_ev_boot(real scalar score, numeric matrix qte_cov_boot, numeric rowvector qte_cov_def, string scalar quantile, numeric scalar level, string scalar boot_cov, string scalar boot_poit_ci, string scalar boot_unif_ci, string scalar tests, real scalar ptests, real scalar obs, real scalar sub_size, real scalar df, real scalar nodfadjust)
{
	real colvector quants, sel, Kmaxuqf, Kmeanuqf
	real scalar nr, nq, k, i, Kalpha, KSstat, CMSstat, j
	real matrix Vuqf, uniform, test, Kuqf, tempvar, temp, qte_cov_boot1, qte_cov_boot2
	real rowvector seuqf, qte_cov_def1, qte_cov_def2, center, center1, Vuqf1
	quants=st_matrix(quantile)
	nr=rows(qte_cov_boot)
	nq=rows(quants)
	k=cols(qte_cov_def)/nq
	Vuqf=(sub_size/obs):*variance(qte_cov_boot)
	st_matrix(boot_cov,Vuqf)
	Vuqf=diagonal(Vuqf)'
	seuqf=sqrt(Vuqf)
	if(nodfadjust==0) st_matrix(boot_poit_ci,(qte_cov_def'-invttail(df,(100-level)/200):*seuqf',qte_cov_def'+invttail(df,(100-level)/200):*seuqf'))
	if(nodfadjust==1) st_matrix(boot_poit_ci,(qte_cov_def'-invnormal((100-level)/200):*seuqf',qte_cov_def'+invnormal((100-level)/200):*seuqf'))	
	if(ptests==1){
		uniform=J(nq*k,2,0)
		test=J(k+1,8,.)
		//test of no effect, coef by coef, KS
		if(score==0) center=qte_cov_def 
		else center=J(1,cols(qte_cov_def),0)
		Kuqf=((qte_cov_boot-center#J(nr,1,1)):^2:/(Vuqf#J(nr,1,1))):^0.5
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmaxuqf=sqrt(sub_size):*rowmax(Kuqf[.,sel])
			Kalpha=mm_quantile(Kmaxuqf:/sqrt(obs),1,level/100)
			uniform[sel',1]=(qte_cov_def[1,sel]-seuqf[1,sel]*Kalpha)'
			uniform[sel',2]=(qte_cov_def[1,sel]+seuqf[1,sel]*Kalpha)'
			KSstat=sqrt(obs):*max((qte_cov_def[1,sel]:^2:/Vuqf[sel]):^0.5)
			test[i,1]=mean(Kmaxuqf:>=KSstat)
		}
		st_matrix(boot_unif_ci,uniform)
		//test of no effects, coef by coef, CMS
		Kuqf=((qte_cov_boot-center#J(nr,1,1)):^2:/(Vuqf#J(nr,1,1)))
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmeanuqf=sub_size:*mean(Kuqf[.,sel]')
			CMSstat=obs*mean((qte_cov_def[sel]:^2:/Vuqf[sel])')
			test[i,2]=mean(Kmeanuqf':>=CMSstat)
		}
		//Test of no effect for all coefficient (except the constant)
		if(k>1){
			qte_cov_def1=J(1,nq,.)
			qte_cov_boot1=J(nr,nq,.)
			for(i=1;i<=nq;i++){
				sel=(i-1)*k:+(1..(k-1))
				tempvar=variance(qte_cov_boot[.,sel])
				tempvar=invsym(tempvar)
				qte_cov_def1[i]=sqrt(qte_cov_def[sel]*tempvar*qte_cov_def[sel]')
				for(j=1;j<=nr;j++){
					qte_cov_boot1[j,i]=sqrt((qte_cov_boot[j,sel]-center[sel])*tempvar*(qte_cov_boot[j,sel]-center[sel])')
				}
			}
			Kmaxuqf=sqrt(sub_size):*rowmax(qte_cov_boot1)
			KSstat=sqrt(obs)*max(qte_cov_def1)
			test[k+1,1]=mean(Kmaxuqf:>=KSstat)
			Kmeanuqf=sub_size:*mean(qte_cov_boot1':^2)
			CMSstat=obs*mean(qte_cov_def1':^2)
			test[k+1,2]=mean(Kmeanuqf':>=CMSstat)		
		}	
		//test of constant effect, coef by coef, KS
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			qte_cov_def1=qte_cov_def[sel]:-mean(qte_cov_def[sel]')
			qte_cov_boot1=qte_cov_boot[.,sel]-J(1,nq,1)#mean(qte_cov_boot[.,sel]')' 
			if(score==0){
				center1=qte_cov_def1
			} else{
				center1=J(1, cols(qte_cov_def1),0)
			}
			Vuqf1=diagonal(variance(qte_cov_boot1))'
//			Kuqf=((qte_cov_boot1-qte_cov_def1#J(nr,1,1)):^2:/(Vuqf[sel]#J(nr,1,1))):^0.5
			Kuqf=((qte_cov_boot1-center1#J(nr,1,1)):^2:/(Vuqf1#J(nr,1,1))):^0.5
			Kmaxuqf=sqrt(sub_size):*rowmax(Kuqf)
			KSstat=sqrt(obs)*max((qte_cov_def1:^2:/Vuqf1):^0.5)
			test[i,3]=mean(Kmaxuqf:>KSstat)
		//test of constant effects, coef by coef, CMS
//			Kuqf=((qte_cov_boot1-qte_cov_def1#J(nr,1,1)):^2:/(Vuqf[sel]#J(nr,1,1)))
			Kuqf=((qte_cov_boot1-center1#J(nr,1,1)):^2:/(Vuqf1#J(nr,1,1)))
			Kmeanuqf=sub_size:*mean(Kuqf')
			CMSstat=obs*mean((qte_cov_def1:^2:/Vuqf1)')
			test[i,4]=mean(Kmeanuqf':>=CMSstat)
		}
		//test of constant effect, all coef except the constant, KS
		if(k>1){
			qte_cov_def1=colshape(qte_cov_def,k)[.,1..(k-1)]
			qte_cov_def1=qte_cov_def1:-mean(qte_cov_def1)#J(nq,1,1)
			qte_cov_def1=rowshape(qte_cov_def1,1)
			if(score==0){
				center1=qte_cov_def1
			} else{
				center1=J(1,(k-1)*nq,0)
			}
			qte_cov_boot1=J(nr,(k-1)*nq,.)
			for(i=1;i<=nr;i++){
				temp=colshape(qte_cov_boot[i,.],k)[.,1..(k-1)]
				temp=temp:-mean(temp)#J(nq,1,1)
				qte_cov_boot1[i,.]=rowshape(temp,1)
			}
			qte_cov_def2=J(1,nq,.)
			qte_cov_boot2=J(nr,nq,.)
			for(i=1;i<=nq;i++){
				sel=(i-1)*(k-1):+(1..(k-1))
				tempvar=variance(qte_cov_boot1[.,sel])
				tempvar=invsym(tempvar)
				qte_cov_def2[i]=sqrt(qte_cov_def1[sel]*tempvar*qte_cov_def1[sel]')
				for(j=1;j<=nr;j++){
					qte_cov_boot2[j,i]=sqrt((qte_cov_boot1[j,sel]-center1[sel])*tempvar*(qte_cov_boot1[j,sel]-center1[sel])')
				}
			}
			Kmaxuqf=sqrt(sub_size):*rowmax(qte_cov_boot2)
			KSstat=sqrt(obs)*max(qte_cov_def2)
			test[k+1,3]=mean(Kmaxuqf:>=KSstat)
			Kmeanuqf=sub_size:*mean(qte_cov_boot2':^2)
			CMSstat=obs*mean(qte_cov_def2':^2)
			test[k+1,4]=mean(Kmeanuqf':>=CMSstat)		
		}	
		//test of stochastic dominance, coef by coef, KS
		qte_cov_boot1=qte_cov_boot-center#J(nr,1,1)
		qte_cov_boot1=qte_cov_boot1:*(qte_cov_boot1:<=0)
		qte_cov_def1=qte_cov_def:*(qte_cov_def:<=0)
		Kuqf=(qte_cov_boot1:^2:/(Vuqf#J(nr,1,1))):^0.5
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmaxuqf=sqrt(sub_size):*rowmax(Kuqf[.,sel])
			KSstat=sqrt(obs)*max((qte_cov_def1[1,sel]:^2:/Vuqf[sel]):^0.5)
			test[i,5]=mean(Kmaxuqf:>=KSstat)
		}
		//test of stochastic dominance, coef by coef, CMS
		Kuqf=qte_cov_boot1:^2:/(Vuqf#J(rows(qte_cov_boot),1,1))
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmeanuqf=sub_size:*mean(Kuqf[.,sel]')
			CMSstat=obs*mean((qte_cov_def1[1,sel]:^2:/Vuqf[sel])')
			test[i,6]=mean(Kmeanuqf':>=CMSstat)
		}
		//Test of stochastic dominance for all coefficient (except the constant)
		if(k>1){
			qte_cov_def2=J(1,nq,.)
			qte_cov_boot2=J(nr,nq,.)
			for(i=1;i<=nq;i++){
				sel=(i-1)*k:+(1..(k-1))
				tempvar=variance(qte_cov_boot[.,sel])
				tempvar=invsym(tempvar)
				qte_cov_def2[i]=sqrt(qte_cov_def1[sel]*tempvar*qte_cov_def1[sel]')
				for(j=1;j<=nr;j++){
					qte_cov_boot2[j,i]=sqrt(qte_cov_boot1[j,sel]*tempvar*qte_cov_boot1[j,sel]')
				}
			}
			Kmaxuqf=sqrt(sub_size):*rowmax(qte_cov_boot2)
			KSstat=sqrt(obs)*max(qte_cov_def2)
			test[k+1,5]=mean(Kmaxuqf:>=KSstat)
			Kmeanuqf=sub_size:*mean(qte_cov_boot2':^2)
			CMSstat=obs*mean(qte_cov_def2':^2)
			test[k+1,6]=mean(Kmeanuqf':>=CMSstat)		
		}			
		//test of being stochastically dominated, KS
		qte_cov_boot1=qte_cov_boot-center#J(nr,1,1)
		qte_cov_boot1=qte_cov_boot1:*(qte_cov_boot1:>=0)
		qte_cov_def1=qte_cov_def:*(qte_cov_def:>=0)
		Kuqf=(qte_cov_boot1:^2:/(Vuqf#J(nr,1,1))):^0.5
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmaxuqf=sqrt(sub_size):*rowmax(Kuqf[.,sel])
			KSstat=sqrt(obs)*max((qte_cov_def1[1,sel]:^2:/Vuqf[sel]):^0.5)
			test[i,7]=mean(Kmaxuqf:>=KSstat)
		}
		//test of being stochastically dominated, CMS
		Kuqf=qte_cov_boot1:^2:/(Vuqf#J(nr,1,1))
		for(i=1;i<=k;i++){
			sel=(1..nq):*k:+i:-k
			Kmeanuqf=sub_size:*mean(Kuqf[.,sel]')
			CMSstat=obs*mean((qte_cov_def1[1,sel]:^2:/Vuqf[sel])')
			test[i,8]=mean(Kmeanuqf':>=CMSstat)
		}
		//Test of being stochastically dominated for all coefficient (except the constant)
		if(k>1){
			qte_cov_def2=J(1,nq,.)
			qte_cov_boot2=J(nr,nq,.)
			for(i=1;i<=nq;i++){
				sel=(i-1)*k:+(1..(k-1))
				tempvar=variance(qte_cov_boot[.,sel])
				tempvar=invsym(tempvar)
				qte_cov_def2[i]=sqrt(qte_cov_def1[sel]*tempvar*qte_cov_def1[sel]')
				for(j=1;j<=nr;j++){
					qte_cov_boot2[j,i]=sqrt(qte_cov_boot1[j,sel]*tempvar*qte_cov_boot1[j,sel]')
				}
			}
			Kmaxuqf=sqrt(sub_size):*rowmax(qte_cov_boot2)
			KSstat=sqrt(obs)*max(qte_cov_def2)
			test[k+1,7]=mean(Kmaxuqf:>=KSstat)
			Kmeanuqf=sub_size:*mean(qte_cov_boot2':^2)
			CMSstat=obs*mean(qte_cov_def2':^2)
			test[k+1,8]=mean(Kmeanuqf':>=CMSstat)		
		}			
		st_matrix(tests,test)
	}
}

version 9.2
mata real matrix dr_lpm(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quants1, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, real scalar all, string scalar nc, string scalar ns)
{
	real colvector dep, wei, o, quants
	real scalar n, nreg, i
	string rowvector reg2
	real matrix reg, A, coef, cov
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quants1)
	o=order(dep,1)
	reg=reg[o,.]
	dep=dep[o]
	wei=wei[o]
	nreg=rows(quants)
	A=invsym(quadcross(reg,wei,reg))*(reg:*wei)'
	A=mm_colrunsum(A')'
//	for (i=1; i<=rows(A); i++) {
//		A[i,.]=quadrunningsum(A[i,.])
//	}
	if(all==0){
		coef=J(rows(A),nreg,.)
		for(i=1; i<=nreg; i++){
			coef[.,i]=A[.,sum(dep:<=quants[i])]
		}
	}
	else{
		coef=select(A,((dep[2..rows(dep)]\.):!=dep)')
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		o=invorder(o)
		cov=drcov(dep[o], reg[o,.], weight_type, wei[o], coef, quants, variance, cluster, strata, touse, "lpm", nc, ns)
	}
	else cov=0
	st_matrix(results, coef)
	st_matrix(covariance, cov)
}

version 9.2
mata real matrix dr_boot(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar quantile, real matrix start, real scalar reps, string scalar boot_method, string scalar cluster, string scalar strata, real scalar noreplacement, real scalar subsize, real matrix Sinfo, real colvector Cinfo, string scalar method, real scalar all, string scalar tempdep, string scalar tempbweight)
{
	real colvector y, quants, w, o, bweight, wei, b, prob
	real matrix x, boot_res, temp_coef, A
	real scalar index_dep, index_weight, index_touse, n, k, nq, index_temp_bw, i, r, rc
	real rowvector index_reg
	string scalar temp_bw
	index_dep=st_varindex(dep)
	index_reg=st_varindex(tokens(reg))
	index_weight=st_varindex(weight)
	index_touse=st_varindex(touse)
	y=st_data(.,index_dep,index_touse)
	x=st_data(.,(index_reg,index_touse),index_touse)
	n=rows(y)
	k=cols(x)
	if(method=="lpm"){
		o=order(y,1)
		x=x[o,.]
		y=y[o]
	}
	quants=st_matrix(quantile)
	nq=rows(quants)
	w=st_data(.,index_weight,index_touse)
	index_temp_bw=st_addvar("double", temp_bw=st_tempname())
	boot_res=J(reps,nq*k,.)
	if(Cinfo!=. & Sinfo==.){
		subsize=rows(Cinfo)*subsize/n
	}
	else if(Cinfo==. & Sinfo!=.){
		subsize=Sinfo*subsize/n
	}
	else if(Cinfo!=. & Sinfo!=.){
		subsize=Sinfo[.,2]*subsize/n
	}	
	for (r=1; r<=reps; r++){
		if(boot_method=="weighted"){
			bweight=dr_prepare_weights(w, index_touse, cluster, strata, "")
		}
		else{
			if(cluster=="" & strata==""){
				if(noreplacement==0){
					bweight=mm_srswr(subsize, n, 1):*w
				}
				else{
					bweight=mm_srswor(subsize, n, 1):*w
				}
			}
			else{
				bweight=mm_sample(subsize, Sinfo, Cinfo, 1, noreplacement, 1):*w
			}
		}
		st_store((1..n)',index_temp_bw,index_touse,bweight)
		stata("_rmcoll "+reg+" [iweight="+temp_bw+"]",1)
		if(st_global("r(varlist)")!=reg){
			stata(`"dis in red "x" _continue"')
			r=r-1
		}
		else{
			if(method=="lpm"){
				wei=bweight[o]
				A=invsym(quadcross(x,wei,x))*(x:*wei)'
//				for (i=1; i<=rows(A); i++) {
//					A[i,.]=quadrunningsum(A[i,.])
//				}
				A=mm_colrunsum(A')'
				if(all==0){
					temp_coef=J(k,nq,.)
					for(i=1; i<=nq; i++){
						temp_coef[.,i]=A[.,sum(y:<=quants[i])]
					}
				}
				else{
					temp_coef=select(A,((y[2..n]\.):!=y)')
				}
			}
			else if(method=="logit"){
				temp_coef=J(k,nq,.)
				for(i=1;i<=nq;i++){
					b=start[.,i]
					prob=logisticcdf(x*b)
					intlog((y:<=quants[i]), x, bweight, ., b, prob)
					if(colmissing(b)){
						st_store(., tempdep, touse, y:<=quants[i])
						st_store(., tempbweight, touse, bweight)
						rc=_stata("logit "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
						if(rc==0){
							temp_coef[.,i]=st_matrix("e(b)")'
						}
					}
					else temp_coef[.,i]=b
				}
			}
			else if(method=="probit"){
				temp_coef=J(k,nq,.)
				for(i=1;i<=nq;i++){
					b=start[.,i]
					prob=normal(x*b)
					intprob((y:<=quants[i]), x, bweight, ., b, prob)
					if(colmissing(b)){
						st_store(., tempdep, touse, y:<=quants[i])
						st_store(., tempbweight, touse, bweight)
						rc=_stata("probit "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
						if(rc==0){
							temp_coef[.,i]=st_matrix("e(b)")'
						}
					}
					else temp_coef[.,i]=b
				}
			}
			else if(method=="cloglog"){
				temp_coef=J(k,nq,.)
				for(i=1;i<=nq;i++){
					b=start[.,i]
					prob=1:-exp(-exp(x*b))
					intcloglog((y:<=quants[i]), x, bweight, ., b, prob)
					if(colmissing(b)){
						st_store(., tempdep, touse, y:<=quants[i])
						st_store(., tempbweight, touse, bweight)
						rc=_stata("cloglog "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
						if(rc==0){
							temp_coef[.,i]=st_matrix("e(b)")'
						}
					}
					else temp_coef[.,i]=b
				}
			}
			boot_res[r,.]=vec(temp_coef)'
			stata(`"di in gr "." _c"')
		}
	}
	return(boot_res)
}

version 9.2
mata real matrix dr_multiplier(string scalar dep, string scalar reg1, string scalar weight, string scalar touse, string scalar quantile, real matrix start, real scalar reps, string scalar boot_method, string scalar cluster, string scalar strata, real scalar noreplacement, real scalar subsize, real matrix Sinfo, real colvector Cinfo, string scalar method)
{
	real colvector y, quants, w, bweight, prob, resid, ptp, fit, den
	real matrix x, boot_res, score, reg, J
	real scalar index_dep, index_weight, index_touse, n, k, nq, i, r
	real rowvector index_reg
	index_dep=st_varindex(dep)
	index_reg=st_varindex(tokens(reg1))
	index_weight=st_varindex(weight)
	index_touse=st_varindex(touse)
	y=st_data(.,index_dep,index_touse)
	x=st_data(.,(index_reg,index_touse),index_touse)
	n=rows(y)
	k=cols(x)
	quants=st_matrix(quantile)
	nq=rows(quants)
	w=st_data(.,index_weight,index_touse)
	score=J(n,nq*k,.)
	reg=x:*sqrt(w)
	if(method=="lpm") J=invsym(quadcross(reg,reg):/n)
	for (i=1; i<=nq; i++){
		fit=x*start[.,i]
		if(method=="lpm") resid=((y:<=quants[i])-fit)
		else{
			if(method=="logit") prob=logisticcdf(fit)
			else if (method=="probit") prob=normal(fit)
			else if (method=="cloglog") prob=1:-exp(-exp(fit))
			ptp=prob:*(1:-prob)
			if(method=="logit"){
				J=invsym(quadcross(reg,ptp,reg):/n)
				resid=((y:<=quants[i])-prob)
			}
			else{
				if(method=="probit") den=normalden(fit)
				else if (method=="cloglog") den=exp(-exp(fit)):*exp(fit)
				resid=((y:<=quants[i])-prob):*den:/ptp
				J=invsym(quadcross(reg,den:^2:/ptp,reg):/n)
			}
		}
		score[.,(i-1)*k:+(1..k)]=(resid:*x)*J
	}
	boot_res=J(reps,nq*k,.)
	if(Cinfo!=. & Sinfo==.){
		subsize=rows(Cinfo)*subsize/n
	}
	else if(Cinfo==. & Sinfo!=.){
		subsize=Sinfo*subsize/n
	}
	else if(Cinfo!=. & Sinfo!=.){
		subsize=Sinfo[.,2]*subsize/n
	}	
	for (r=1; r<=reps; r++){
		if(boot_method=="exponential" | boot_method=="Gaussian" | boot_method=="wild"){
			bweight=dr_prepare_weights(w, index_touse, cluster, strata, boot_method)
		}
		else{
			if(cluster=="" & strata==""){
				if(noreplacement==0){
					bweight=mm_srswr(subsize, n, 1):*w*n/subsize
				}
				else{
					bweight=mm_srswor(subsize, n, 1):*w*n/subsize
				}
			}
			else{
				bweight=mm_sample(subsize, Sinfo, Cinfo, 1, noreplacement, 1):*w*n/subsize
			}
		}
		boot_res[r,.]=mean(score:*bweight)
		stata(`"di in gr "." _c"')
	}
	return(boot_res)
}

version 9.2
mata real matrix drcov(real colvector dep1, real matrix reg1, string scalar weight_type, real colvector weights, real matrix coef, real colvector tau, string scalar variance, string scalar cluster1, string scalar strata1, string scalar touse, string scalar method, string scalar ncluster, string scalar nstrata)
{
	real scalar k, n, nq, nc, ns, i, q, c, s, j, n_effective
	real matrix resid, reg, J, cov, temp, temp1q, temp1i, temp2q, temp2i
	real colvector cluster, uniq_clust, strata, uniq_strata, residq, residi, fit, prob, ptp, den, index
	real rowvector tempmq, tempmi
	k=cols(reg1)
	n=rows(reg1)
	nq=length(tau)
	variance=J(nq,1,variance)
	resid=J(n,nq,.)
	if(cluster1!=""){
		cluster=st_data(.,cluster1,touse)
		uniq_clust=uniqrows(cluster)
		nc=rows(uniq_clust)
		st_numscalar(ncluster, nc)
	}
	else st_numscalar(ncluster, n)
	if(strata1!=""){
		strata=st_data(.,strata1,touse)
		uniq_strata=uniqrows(strata)
		ns=rows(uniq_strata)
		st_numscalar(nstrata, ns)
	}
	else st_numscalar(nstrata, 1)
	reg=reg1:*sqrt(weights)
	if(method=="lpm"){
		for(i=1;i<=nq;i++) resid[,i]=((dep1:<=tau[i])-reg1*coef[.,i])
		J=invsym(quadcross(reg,reg))
	}
	else{
		J=J(k,nq*k,.)
		for(i=1;i<=nq;i++){
			fit=reg1*coef[.,i]
			if(method=="logit") prob=logisticcdf(fit)
			else if (method=="probit") prob=normal(fit)
			else if (method=="cloglog") prob=1:-exp(-exp(fit))
			ptp=prob:*(1:-prob)
			if(method=="logit"){
				J[.,(i-1)*k:+(1..k)]=invsym(quadcross(reg,ptp,reg))
				resid[,i]=((dep1:<=tau[i])-prob)
			}
			else{
				if(method=="probit") den=normalden(fit)
				else if (method=="cloglog") den=exp(-exp(fit)):*exp(fit)
				resid[,i]=((dep1:<=tau[i])-prob):*den:/ptp
				J[.,(i-1)*k:+(1..k)]=invsym(quadcross(reg,den:^2:/ptp,reg))
			}
		}
	}
	resid=resid:*sqrt(weights)
	cov=J(nq*k,nq*k,.)
	for(q=1; q<=nq; q++){
		for(i=1; i<=q; i++){
			if(cluster1=="" & strata1==""){
				temp=quadcross(reg:*resid[.,i],reg:*resid[.,q])
//				temp=n/(n-k)*temp
			}
			else if(strata1==""){
				temp=J(k,k,0)
				for(c=1; c<=nc; c++){
					index=my_selectindex(cluster:==uniq_clust[c])
					residq=resid[index,q]:*reg[index,.]
					residi=resid[index,i]:*reg[index,.]
					temp=temp:+colsum(residq)'colsum(residi)
				}
//				temp=((n-1)/(n-k))*(nc/(nc-1))*temp				
			}
			else if(cluster1==""){
				temp=J(k,k,0)
				for(s=1; s<=ns; s++){
					index=my_selectindex(strata:==uniq_strata[s])
					temp1q=resid[index,q]:*reg[index,.]
					temp1i=resid[index,i]:*reg[index,.]
					tempmq=mean(temp1q)
					tempmi=mean(temp1i)
					for(j=1;j<=rows(temp1q);j++){
						temp=temp+(temp1q[j,.]-tempmq)'(temp1i[j,.]-tempmi)
					}
				}
//				temp=n/(n-ns*k)*temp
			}
			else{
				temp=J(k,k,0)
				for(s=1; s<=ns; s++){
					uniq_clust=uniqrows(select(cluster,strata:==uniq_strata[s]))
					nc=rows(uniq_clust)
					index=my_selectindex(strata:==uniq_strata[s])
					temp1q=resid[index,q]:*reg[index,.]
					temp1i=resid[index,i]:*reg[index,.]
					tempmq=mean(temp1q)
					tempmi=mean(temp1i)
					for(c=1; c<=nc; c++){
						index=my_selectindex((cluster:==uniq_clust[c]):*(strata:==uniq_strata[s]))
						temp2q=resid[index,q]:*reg[index,.]
						temp2i=resid[index,i]:*reg[index,.]
						temp=temp:+colsum(temp2q:-tempmq#J(rows(temp2q),1,1))'colsum(temp2i:-tempmi#J(rows(temp2i),1,1))
					}
				}
//				temp=(rows(uniqrows((cluster,strata)))-ns)/ns*temp
			}			
			if(method=="lpm") cov[((q-1)*k+1)..(q*k),((i-1)*k+1)..(i*k)]=J*temp*J 
			else cov[((q-1)*k+1)..(q*k),((i-1)*k+1)..(i*k)]=J[.,(q-1)*k:+(1..k)]*temp*J[.,(i-1)*k:+(1..k)]
			if(i!=q) cov[((i-1)*k+1)..(i*k),((q-1)*k+1)..(q*k)]=cov[((q-1)*k+1)..(q*k),((i-1)*k+1)..(i*k)]'
		}
	}
	if(weight_type=="fweight"){
		n_effective=colsum(weights)
		cov=cov*n/n_effective
	}
	return(cov)
}

*Parsing vce
program define dr_VCEParse, rclass
	syntax [anything(name=var)] [, Reps(integer 100) BMethod(string) Cluster(varname numeric) Strata(varname numeric) funcint(string) Subsize(integer -10) noReplacement weight(real 0)  hc1 NODFadjust]
	if "`var'"=="" & "`funcint'"==""{
		local var "robust"
	}
	else if "`var'"==""{
		local var "bootstrap"
	}
	else if "`var'"=="boot" | "`var'"=="boots" | "`var'"=="bootst" | "`var'"=="bootstr" | "`var'"=="bootstra"{
		local var "bootstrap"
	}
	else if "`var'"=="mult" | "`var'"=="multi" | "`var'"=="multip" | "`var'"=="multipl" | "`var'"=="multipli" | "`var'"=="multiplie" | "`var'"=="multiplier" | "`var'"=="score"{
		local var "multiplier"
	}
	else if "`var'"=="r" | "`var'"=="ro" | "`var'"=="rob" | "`var'"=="robu" | "`var'"=="robus"{
		local var "robust"
	}
	else if "`var'"!="novar" & "`var'"!="robust"{
		dis as error "The option vce(`var') is not allowed."
		error 400
	}
	if ("`var'"=="robust" | "`var'"=="novar") & "`funcint'"=="functional"{
		dis as error "Only the (classical or multiplier) bootstrap can be used to estimate the variance when the option functional is activated."
		error 400
	}
	if `reps'<2{
		dis as error "At least two bootstrap replications must be performed. We recommend using at least 100 replications."
		error 400
	}
	if "`var'"=="bootstrap"{
		if "`bmethod'"==""{
			local bmethod "empirical"
		}
		if "`bmethod'"=="subsampling" & `subsize'==-10{
			dis as error "The option subsize must be specified if bmethod(subsampling) is selected."
			error 400
		}
		if `subsize'!=-10 & `subsize'<0{
			dis as error "The option subsize(`subsize') is not allowed."
			error 400
		}
		if "`bmethod'"!="empirical" & "`bmethod'"!="subsampling" & "`bmethod'"!="weighted"{
			dis as error "The option bmethod(`bmethod') is not allowed."
			error 400
		}
		if "`bmethod'"!="subsampling" & (`subsize'!=-10 | "`replacement'"!=""){
			dis as error "The options subsize and noreplacement cannot be used when boot_method(`boot_method')."
			error 400
		}
		local noreplacement=("`replacement'"=="noreplacement")		
	}
	if "`var'"=="multiplier"{
		if "`bmethod'"==""{
			local bmethod "wild"
		}
		if "`bmethod'"=="subsampling" & `subsize'==-10{
			dis as error "The option subsize must be specified if bmethod(subsampling) is selected."
			error 400
		}
		if `subsize'!=-10 & `subsize'<0{
			dis as error "The option subsize(`subsize') is not allowed."
			error 400
		}
		if "`bmethod'"!="empirical" & "`bmethod'"!="subsampling" & "`bmethod'"!="exponential" & "`bmethod'"!="Gaussian" & "`bmethod'"!="wild"{
			dis as error "The option bmethod(`bmethod') is not allowed."
			error 400
		}
		if "`bmethod'"!="subsampling" & (`subsize'!=-10 | "`replacement'"!=""){
			dis as error "The options subsize and noreplacement cannot be used when boot_method(`boot_method')."
			error 400
		}
		local noreplacement=("`replacement'"=="noreplacement")		
	}
	return local var `var'
	return local reps `reps'
	return local boot_method `bmethod'
	return local cluster `cluster'
	return local strata `strata'
	return local sub_size `subsize'
	return local noreplacement `noreplacement'
	if "`hc1'"=="hc1" 	return local HC1 "HC1"
	return local nodfadjust `nodfadjust'
end

*Parsing method
program define dr_methodParse, rclass
	syntax [anything(name=method)] , [ MAX_it(real 100) ONEstep]
	if "`method'"==""{
		local method "logit"
	}
	if "`method'"!="logit" & "`method'"!="probit"  & "`method'"!="lpm" & "`method'"!="cloglog"{
		dis as error "The selected method in the option 'method' has not been implemented."
		error 400
	}
	return local method `method'
	return local max_it `max_it'
	return local onestep `onestep'
end

*bootstrap weights for weighted bootstrap with clustering and/or stratification
version 9.2
mata real colvector dr_prepare_weights(real colvector weight, real scalar touse, string scalar cluster1, string scalar strata1, string scalar method)
{
	real scalar n, nc, c, ns, s
	real colvector bootw, cluster, uniq_clust, clust_w, index, strata, uniq_strata, strata_w
	n=rows(weight)
	if(cluster1=="" & strata1==""){
		if(method=="") bootw=-ln(uniform(n,1)):*weight
		else if(method=="exponential") bootw=(-ln(uniform(n,1)):-1):*weight
		else if(method=="Gaussian") bootw=invnormal(uniform(n,1)):*weight
		else if(method=="wild") bootw=(invnormal(uniform(n,1))/sqrt(2)+(invnormal(uniform(n,1)):^2:-1):/2):*weight
	}
	else if(cluster1!="" & strata1==""){
		bootw=J(n,1,.)
		cluster=st_data(.,cluster1,touse)
		uniq_clust=uniqrows(cluster)
		nc=rows(uniq_clust)
		if(method=="") clust_w=-ln(uniform(nc,1))
		else if(method=="exponential") clust_w=(-ln(uniform(nc,1)):-1)
		else if(method=="Gaussian") clust_w=invnormal(uniform(nc,1))
		else if(method=="wild") clust_w=invnormal(uniform(nc,1))/sqrt(2)+(invnormal(uniform(nc,1)):^2:-1):/2
		for(c=1;c<=nc;c++){
			index=my_selectindex(cluster:==uniq_clust[c])
			bootw[index]=J(rows(index),1,clust_w[c]):*weight[index]
		}
	}
	else if(cluster1=="" & strata1!=""){
		if(method=="") bootw=-ln(uniform(n,1)):*weight
		else if(method=="exponential") bootw=(-ln(uniform(n,1)):-1):*weight
		else if(method=="Gaussian") bootw=invnormal(uniform(n,1)):*weight
		else if(method=="wild") bootw=(invnormal(uniform(n,1))/sqrt(2)+(invnormal(uniform(n,1)):^2:-1):/2):*weight
		strata=st_data(.,strata1,touse)
		uniq_strata=uniqrows(strata)
		ns=rows(uniq_strata)
		for(s=1;s<=ns;s++){
			index=my_selectindex(strata:==uniq_strata[s])
			bootw[index]=bootw[index]:*sum(weight[index])/sum(bootw[index])
		}
	}
	else{
		bootw=J(n,1,.)
		cluster=st_data(.,cluster1,touse)
		strata=st_data(.,strata1,touse)
		uniq_strata=uniqrows(strata)
		ns=rows(uniq_strata)
		if(method=="") strata_w=-ln(uniform(ns,1))
		else if(method=="exponential") strata_w=(-ln(uniform(ns,1)):-1)
		else if(method=="Gaussian") strata_w=invnormal(uniform(ns,1))
		else if(method=="wild") strata_w=invnormal(uniform(ns,1))/sqrt(2)+(invnormal(uniform(ns,1)):^2:-1):/2
		for(s=1;s<=ns;s++){
			uniq_clust=uniqrows(select(cluster,strata:==uniq_strata[s]))
			for(c=1;c<=rows(uniq_clust);c++){
				index=my_selectindex((cluster:==uniq_clust[c]):*(strata:==uniq_strata[s]))
				bootw[index]=J(rows(index),1,strata_w[s]):*weight[index]
			}
		}
		for(s=1;s<=ns;s++){
			index=my_selectindex(strata:==uniq_strata[s])
			bootw[index]=bootw[index]:*sum(weight[index])/sum(bootw[index])
		}
	}
	return(bootw)
}

version 9.2
mata real colvector my_selectindex(real colvector input) return(select((1..rows(input))',input))

version 9.2
mata real matrix dr_boot_1step(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar quantile, real scalar reps, string scalar boot_method, string scalar cluster, string scalar strata, real scalar noreplacement, real scalar subsize, real matrix Sinfo, real colvector Cinfo, string scalar method, real matrix start, string scalar tempdep, string scalar tempbweight)
{
	real colvector y, quants, w, bweight, fit, prob, b, den, pp
	real matrix x, boot_res, temp_coef, temp
	real scalar index_touse, n, k, nq, index_temp_bw, i, r, index1, quants1, rc
	real rowvector index, direction
	string scalar temp_bw
	string rowvector reg1
	index_touse=st_varindex(touse)
	index_temp_bw=st_addvar("double", temp_bw=st_tempname())
	quants=st_matrix(quantile)
	y=st_data(.,dep,touse)
	n=rows(y)
	reg1=tokens(reg)
	if(length(reg1)>0) x=st_data(.,reg1,touse),J(n,1,1)
	else x=J(n,1,1)
	k=cols(x)
	w=st_data(.,weight,touse)
	w=w:/mean(w)
	nq=rows(quants)
	index=(ceil((nq+1)/2)..nq),((ceil((nq+1)/2)-1)..1)
	direction=J(1,nq-ceil((nq+1)/2)+1,1),J(1,ceil((nq+1)/2)-1,-1)
	boot_res=J(reps,nq*k,.)
	if(Cinfo!=. & Sinfo==.){
		subsize=rows(Cinfo)*subsize/n
	}
	else if(Cinfo==. & Sinfo!=.){
		subsize=Sinfo*subsize/n
	}
	else if(Cinfo!=. & Sinfo!=.){
		subsize=Sinfo[.,2]*subsize/n
	}	
	for (r=1; r<=reps; r++){
		if(boot_method=="weighted"){
			bweight=dr_prepare_weights(w, index_touse, cluster, strata, "")
		}
		else{
			if(cluster=="" & strata==""){
				if(noreplacement==0){
					bweight=mm_srswr(subsize, n, 1):*w
				}
				else{
					bweight=mm_srswor(subsize, n, 1):*w
				}
			}
			else{
				bweight=mm_sample(subsize, Sinfo, Cinfo, 1, noreplacement, 1):*w
			}
		}
		st_store(.,index_temp_bw,index_touse,bweight)
		stata("_rmcoll "+reg+" [iweight="+temp_bw+"]",1)
		if(st_global("r(varlist)")!=reg){
			stata(`"dis in red "x" _continue"')
			r=r-1
		}
		else{
			temp_coef=J(k,nq,.)
			for (i=1; i<=rows(quants); i++) {
				index1=index[i]
				quants1=quants[index1]
				if (i>1){
					fit=x*temp_coef[.,index1-direction[i]]				
					if (method=="logit"){
						prob=logisticcdf(fit)					
						temp=cross(x,prob:*(1:-prob):*bweight,x):/n
						temp=luinv(temp)
						temp=temp*mean((prob:-(y:<=quants1)):*x,bweight)'
					}
					else{
						if (method=="probit"){
							prob=normal(fit)
							den=normalden(fit)
						}
						else if (method=="cloglog"){
							prob=1:-exp(-exp(fit))
							den=exp(-exp(fit)):*exp(fit)
						}
						pp=prob:*(1:-prob)
						temp=cross(x,den:^2:*bweight:/pp,x):/n
						temp=luinv(temp)
						temp=temp*mean((prob:-(y:<=quants1)):*den:/pp:*x,bweight)'
					}
					temp_coef[.,index1]=temp_coef[.,index1-direction[i]]-temp
				}
				if (i==1 | colmissing(temp_coef[.,index1])){
					b=start[.,index1]
					if (method=="logit"){
						prob=logisticcdf(x*b)
						intlog((y:<=quants1), x, bweight, ., b, prob)
						if(colmissing(b)){
							st_store(., tempdep, touse, y:<=quants1)
							st_store(., tempbweight, touse, bweight)
							rc=_stata("logit "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
							if(rc==0){
								b=st_matrix("e(b)")'
							}
						}
					}
					else if (method=="probit"){
						prob=normal(x*b)
						intprob((y:<=quants1), x, bweight, ., b, prob)
						if(colmissing(b)){
							st_store(., tempdep, touse, y:<=quants1)
							st_store(., tempbweight, touse, bweight)
							rc=_stata("probit "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
							if(rc==0){
								b=st_matrix("e(b)")'
							}
						}
					}
					else if (method=="cloglog"){
						prob=1:-exp(-exp(x*b))
						intcloglog((y:<=quants1), x, bweight, ., b, prob)
						if(colmissing(b)){
							st_store(., tempdep, touse, y:<=quants1)
							st_store(., tempbweight, touse, bweight)
							rc=_stata("cloglog "+tempdep+" "+reg+" [pweight="+tempbweight+"], asis",1)
							if(rc==0){
								b=st_matrix("e(b)")'
							}
						}
					}
					temp_coef[,index1]=b
				}
			}
			boot_res[r,.]=vec(temp_coef)'
			stata(`"di in gr "." _c"')
		}
	}
	return(boot_res)
}

mata real matrix dr_logit(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quants1, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, o, weio, b, b_start, prob
	real scalar n, k, nreg, i, rc
	string rowvector reg2
	real matrix reg, rego, A, coef, cov
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quants1)
	nreg=rows(quants)
	o=order(dep,1)
	rego=reg[o,.]
	weio=wei[o]
	A=invsym(quadcross(rego,weio,rego))*(rego:*weio)'
	A=mm_colrunsum(A')'
	b_start=J(k,nreg,.)
	for(i=1; i<=nreg; i++){
		b_start[.,i]=A[.,sum(dep:<=quants[i])]
	}
	coef=J(k,nreg,.)
	for(i=1;i<=nreg;i++){
		prob=logisticcdf(reg*b_start[.,i])
		intlog((dep:<=quants[i]), reg, wei, ., b=b_start[.,i], prob)
		if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){
			st_store(.,tempdep,touse,dep:<=quants[i])
			rc=_stata("logit "+tempdep+" "+reg1+" [pweight="+wei1+"], asis",1)
			if(rc==0){
				coef[.,i]=st_matrix("e(b)")'
			}
			else{
				coef[.,i]=J(k,1,0)
				displayas("err")
				printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
			}
		}
		else coef[,i]=b
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "logit", nc, ns)
	}
	else cov=0
	st_matrix(results, coef)
	st_matrix(covariance, cov)
}

version 9.2
mata void intlog(real colvector dep, real matrix reg, real colvector we, real scalar convergence, real colvector b, real colvector prob)
{
//variable declarations
	real scalar objo, objn, it
	real rowvector db
	objo=0
	objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
	db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	it=1
	while(it<100 & sum(abs(db):>1e-8)>0 & abs(objn-objo)>1e-8){	
		objo=objn
		b=b+db'
		it=it+1
		prob=logisticcdf(reg*b)
		objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
		db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	}
	convergence=(it==100)
}

*Mata, logistic distribution
version 9.2
mata real colvector logisticcdf(real colvector x){
	real colvector cdf
	real scalar lm
	cdf=1:/(1:+exp(-x))
	if(colmissing(cdf)){
		lm=log(maxdouble())
		cdf[my_selectindex(x:<(-lm)),1]=J(sum(x:<(-lm)),1,0)
	}
	return(cdf)
}

version 9.2
mata void dr_logit1(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quantsn, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, fit, prob, b
	real scalar n, k, nq, i, index1, quants1, rc
	string rowvector reg2
	real matrix reg, coef, temp, cov
	real rowvector index, direction
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quantsn)
	nq=rows(quants)
	coef=J(k,nq,.)
	index=(ceil((nq+1)/2)..nq),((ceil((nq+1)/2)-1)..1)
	direction=J(1,nq-ceil((nq+1)/2)+1,1),J(1,ceil((nq+1)/2)-1,-1)
	for(i=1; i<=nq; i++) {
		index1=index[i]	
		quants1=quants[index1]
		if(i>1){
			fit=reg*coef[.,index1-direction[i]]
			prob=logisticcdf(fit)
			temp=cross(reg,prob:*(1:-prob):*wei,reg):/n
			temp=luinv(temp)
			temp=temp*mean((prob:-(dep:<=quants1)):*reg,wei)'
			coef[.,index1]=coef[.,index1-direction[i]]-temp
		}
		if(i==1 | colmissing(coef[.,index1])){
			b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants1))	
			prob=logisticcdf(reg*b)
			intlog((dep:<=quants1), reg, wei, ., b, prob)		
			if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){
				st_store(.,tempdep,touse,dep:<=quants1)
				rc=_stata("logit "+tempdep+" "+reg1+" [pweight="+wei1+"], asis",1)
				if(rc==0){
					coef[.,index1]=st_matrix("e(b)")'
				}
				else{
					coef[.,i]=J(k,1,0)
					displayas("err")
					printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
				}
				b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants[i]))
				prob=logisticcdf(reg*b)
			}
			else coef[.,index1]=b
		}
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "logit", nc, ns)
	}
	else cov=0
	st_matrix(results,coef)
	st_matrix(covariance,cov)
}

mata real matrix dr_probit(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quants1, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, o, weio, b_start, b, prob
	real scalar n, k, nreg, i, rc
	string rowvector reg2
	real matrix reg, rego, A, coef, cov
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quants1)
	nreg=rows(quants)
	o=order(dep,1)
	rego=reg[o,.]
	weio=wei[o]
	A=invsym(quadcross(rego,weio,rego))*(rego:*weio)'
	A=mm_colrunsum(A')'
	b_start=J(k,nreg,.)
	for(i=1; i<=nreg; i++){
		b_start[.,i]=A[.,sum(dep:<=quants[i])]
	}
	coef=J(k,nreg,.)
	for(i=1;i<=nreg;i++){
		prob=normal(reg*b_start[.,i])
		intprob((dep:<=quants[i]), reg, wei, ., b=b_start[.,i], prob)
		if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){		
			st_store(.,tempdep,touse,dep:<=quants[i])
			rc=_stata("probit "+tempdep+" "+reg1+" [pweight="+wei1+"], asis",1)
			if(rc==0){
				coef[.,i]=st_matrix("e(b)")'
				b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants[i]))
				prob=normal(reg*b)
			}
			else{
				coef[.,i]=J(k,1,0)
				displayas("err")
				printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
			}
		}
		else coef[,i]=b
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "probit", nc, ns)
	}
	else cov=0
	st_matrix(results, coef)
	st_matrix(covariance, cov)
}

version 9.2
mata void intprob(real colvector dep, real matrix reg, real colvector we, real scalar convergence, real colvector b, real colvector prob)
{
//variable declarations
	real scalar objo, objn, it
	real colvector fit, pp, den
	real rowvector db
	objo=0
	objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
	db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	it=1
	while(it<100 & sum(abs(db):>1e-8)>0 & abs(objn-objo)>1e-8){
		objo=objn
		b=b+db'
		it=it+1
		fit=reg*b
		prob=normal(fit)
		pp=prob:*(1:-prob)
		den=normalden(fit)
		objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
		db=colsum(we:*(reg:*(dep-prob):*den:/pp))*invsym((we:*reg:*den:^2:/pp)'reg)
	}
	convergence=(it==100)
}

version 9.2
mata void dr_probit1(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quantsn, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, fit, prob, b, den, pp
	real scalar n, k, nq, i, index1, quants1, rc
	string rowvector reg2
	real matrix reg, coef, temp, cov
	real rowvector index, direction
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quantsn)
	nq=rows(quants)
	coef=J(k,nq,.)
	index=(ceil((nq+1)/2)..nq),((ceil((nq+1)/2)-1)..1)
	direction=J(1,nq-ceil((nq+1)/2)+1,1),J(1,ceil((nq+1)/2)-1,-1)
	for(i=1; i<=nq; i++) {
		index1=index[i]
		quants1=quants[index1]
		if(i>1){
			fit=reg*coef[.,index1-direction[i]]
			prob=normal(fit)
			den=normalden(fit)
			pp=prob:*(1:-prob)
			temp=cross(reg,den:^2:*wei:/pp,reg):/n
			temp=luinv(temp)
			temp=temp*mean((prob:-(dep:<=quants1)):*den:/pp:*reg,wei)'
			coef[.,index1]=coef[.,index1-direction[i]]-temp
		}
		if(i==1 | colmissing(coef[.,index1])){
			b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants1))
			prob=normal(reg*b)
			intprob((dep:<=quants1), reg, wei, ., b, prob)
			if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){
				st_store(.,tempdep,touse,dep:<=quants1)
				rc=_stata("probit "+tempdep+" "+reg1+" [pweight="+wei1+"], asis",1)
				if(rc==0){
					coef[.,index1]=st_matrix("e(b)")'
				}
				else{
					coef[.,i]=J(k,1,0)
					displayas("err")
					printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
				}
				b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants[i]))
				prob=logisticcdf(reg*b)
			}
			else coef[.,index1]=b
		}
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "probit", nc, ns)
	}
	else cov=0
	st_matrix(results,coef)
	st_matrix(covariance,cov)
}

mata real matrix dr_cloglog(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quants1, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, o, weio, b, prob
	real scalar n, k, nreg, i, rc
	string rowvector reg2
	real matrix reg, rego, A, b_start, coef, cov
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quants1)
	nreg=rows(quants)
	o=order(dep,1)
	rego=reg[o,.]
	weio=wei[o]
	A=invsym(quadcross(rego,weio,rego))*(rego:*weio)'
	A=mm_colrunsum(A')'
	b_start=J(k,nreg,.)
	for(i=1; i<=nreg; i++){
		b_start[.,i]=A[.,sum(dep:<=quants[i])]
	}
	coef=J(k,nreg,.)
	for(i=1;i<=nreg;i++){
		prob=1:-exp(-exp(reg*b_start[.,i]))
		intcloglog((dep:<=quants[i]), reg, wei, ., b=b_start[.,i], prob)	
		if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){		
			st_store(.,tempdep,touse,dep:<=quants[i])
			rc=_stata("cloglog "+tempdep+" "+reg1+" [pweight="+wei1+"], iter(1000) asis",1)
			if(rc==0){
				coef[.,i]=st_matrix("e(b)")'
				b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants[i]))
				prob=1:-exp(-exp(reg*b))
			}
			else{
				coef[.,i]=J(k,1,0)
				displayas("err")
				printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
			}
		}
		else coef[,i]=b
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "cloglog", nc, ns)
	}
	else cov=0
	st_matrix(results, coef)
	st_matrix(covariance, cov)
}

version 9.2
mata void intcloglog(real colvector dep, real matrix reg, real colvector we, real scalar convergence, real colvector b, real colvector prob)
{
//variable declarations
	real scalar objo, objn, it
	real colvector fit, pp, den
	real rowvector db
	objo=0
	objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
	db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	it=1
	while(it<100 & sum(abs(db):>1e-8)>0 & abs(objn-objo)>1e-8){
		objo=objn
		b=b+db'
		it=it+1
		fit=reg*b
		prob=1:-exp(-exp(fit))
		pp=prob:*(1:-prob)
		den=exp(-exp(fit)):*exp(fit)
		objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
		db=colsum(we:*(reg:*(dep-prob):*den:/pp))*invsym((we:*reg:*den:^2:/pp)'reg)
	}
	convergence=(it==100)
}

version 9.2
mata void dr_cloglog1(string scalar dep1, string scalar reg1, string scalar weight_type, string scalar wei1, string scalar touse, string scalar quantsn, string scalar variance, string scalar cluster, string scalar strata, string scalar results, string scalar covariance, string scalar nc, string scalar ns, string scalar tempdep)
{
	real colvector dep, wei, quants, fit, prob, b, den, pp
	real scalar n, k, nq, i, index1, quants1, rc
	string rowvector reg2
	real matrix reg, coef, temp, cov
	real rowvector index, direction
	dep=st_data(.,dep1,touse)
	n=rows(dep)
	reg2=tokens(reg1)
	if(length(reg2)>0) reg=st_data(.,reg2,touse),J(n,1,1)
	else reg=J(n,1,1)
	k=cols(reg)
	wei=st_data(.,wei1,touse)
	quants=st_matrix(quantsn)
	nq=rows(quants)
	coef=J(k,nq,.)
	index=(ceil((nq+1)/2)..nq),((ceil((nq+1)/2)-1)..1)
	direction=J(1,nq-ceil((nq+1)/2)+1,1),J(1,ceil((nq+1)/2)-1,-1)
	for(i=1; i<=nq; i++) {
		index1=index[i]
		quants1=quants[index1]
		if(i>1){
			fit=reg*coef[.,index1-direction[i]]
			prob=1:-exp(-exp(fit))
			den=exp(-exp(fit)):*exp(fit)
			pp=prob:*(1:-prob)
			temp=cross(reg,den:^2:*wei:/pp,reg):/n
			temp=luinv(temp)
			temp=temp*mean((prob:-(dep:<=quants1)):*den:/pp:*reg,wei)'
			coef[.,index1]=coef[.,index1-direction[i]]-temp
		}
		if(i==1 | colmissing(coef[.,index1])){
			b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants1))
			prob=1:-exp(-exp(reg*b))
			intcloglog((dep:<=quants1), reg, wei, ., b, prob)
			if(colmissing(b) | sum(prob:<=epsilon(1))+sum(prob:>=1-epsilon(1))==n){
				st_store(.,tempdep,touse,dep:<=quants1)
				rc=_stata("cloglog "+tempdep+" "+reg1+" [pweight="+wei1+"], asis iter(1000)",1)
				if(rc==0){
					coef[.,index1]=st_matrix("e(b)")'
				}
				else{
					coef[.,i]=J(k,1,0)
					displayas("err")
					printf("{red}WARNING: The distribution regression at the threshold %f could not be estimated. The coefficients at this threshold will take the value of 0 but should be interpreted as missing.\n", quants[i])
				}
				b=invsym(cross(reg,wei,reg))*cross(reg,wei,(dep:<=quants[i]))
				prob=1:-exp(-exp(reg*b))
			}
			else coef[.,index1]=b
		}
	}
	if(variance!="novar" & variance!="boot" & variance!="bootstrap" & variance!="multiplier"){
		cov=drcov(dep, reg, weight_type, wei, coef, quants, variance, cluster, strata, touse, "cloglog", nc, ns)
	}
	else cov=0
	st_matrix(results,coef)
	st_matrix(covariance,cov)
}
