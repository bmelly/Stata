*mdqr: minimum distance quantile regression
*! version 0.0.5  16.11.2022  Blaise Melly

*To do (potentially): (1) parallel processing with the parallel package (or multishell) DONE
*(2) Covariance between coefficients are different quantiles (multi-variate GMM?)
*(4) bootstrap all the quantiles together to obtain valid covariances.

cap prog drop mdqr
program mdqr, eclass byable(recall) sortpreserve
	local stata_version = _caller() 
	version 9.2
*check that moremata is installed
	capt findfile lmoremata.mlib
	if _rc {
      	di as error "-moremata- is required; type {stata ssc install moremata}"
		error 499
	}
*check that qrprocess is installed
	capt findfile qrprocess.ado
	if _rc {
      	di as error `"-qrprocess- is required; type {stata "net install qrprocess, from(https://raw.githubusercontent.com/bmelly/Stata/main/)"}"'
		error 499
	}
*if the command is used without arguments it shows the previous results
	if replay() {
		if "`e(cmd)'"!="mdqr" { 
			error 301 
		} 
		if _by() {
			error 190 
		}
        syntax [, Level(cilevel)]
		ereturn display, level(`level')
 	}
	else {
*syntax
		syntax anything [if] [in] [pweight] , Group(varlist) [Quantiles(numlist >0 <1 sort) Cluster(varname) Est_command(string) est_opts(string) qr_opts(string) BOOTstrap(string) Level(cilevel) noPrint Save_first(name) Load_first(name) Version(string) n_small(integer 1) PARallel]
*separate dependent variable from regressors
		gettoken dep anything : anything
		confirm numeric variable `dep'
*separate exogenous variables
		gettoken exo anything : anything, parse("(")
		local fvops = 0
		if "`exo'" == "(" {
			local exo ""
		}
		else {
			if `stata_version' >= 11{
				version 11.1: fvexpand `exo'
				if "`r(fvops)'" == "true"{
					local fvops = 1
					local exo `r(varlist)'
				}
			}
			else if "`exo'" != ""{ 
				unab exo : `exo'
				confirm numeric variable `exo'
			}
			gettoken away anything : anything, parse("(")
		}
*separate endogenous variables
		gettoken endo anything : anything, parse("=")		
		if "`endo'" == "=" {
			local endo ""
		}
		else {
			if `stata_version' >= 11{
				version 11.1: fvexpand `endo'
				if "`r(fvops)'" == "true"{
					local fvops = 1
					local endo `r(varlist)'
				}
			}
			else if "`endo'" != ""{
				confirm numeric variable `endo' 
			}
			gettoken away anything : anything, parse("=")
		}
*separate instruments		
		gettoken inst anything : anything, parse(")")		
		if "`inst'" == ")" {
			local "`inst'" ""
		}
		else {
			if `stata_version' >= 11{
				version 11.1: fvexpand `inst'
				if "`r(fvops)'" == "true"{
					local fvops = 1
					local inst `r(varlist)'
				}
			}
			else if "`inst'" != ""{
				confirm numeric variable `inst' 
			}
		}	
*to deal with factor variables
		if "`version'" != ""{
			local vv:  di "version `version', missing: "
		}
		else if `fvops' | "`est_command'"=="ivreghdfe"{
			local vv: di "version " string(max(11, `stata_version')) ", missing: " 
		}
*sample definition
		marksample touse
		markout `touse' `dep' `exo' `endo' `inst' `group'
*Generate a simple group variables to simplify iterations	
		tempvar groupvar
		quietly egen `groupvar' = group(`group') if `touse'
		quietly su `groupvar'  if `touse', meanonly
		local ngroup = r(max)
*Identify time-varying exogenous variables
		tempvar n_tv_var test
		quiet gen `n_tv_var' = 1 + `n_small' if `touse'
		quietly gen `test' = .	
		if "`exo'" != ""{
			quietly foreach v of local exo {
				if `fvops' == 1{
					unopvarlist `v'
					local v_unfactored `r(varlist)'
				} 
				else {
					local v_unfactored "`v'"
				}
				bysort `touse' `groupvar' (`v_unfactored') : replace `test' = (`v_unfactored'[1] == `v_unfactored'[_N]) if `touse'
				su `test'  if `touse', meanonly
				if r(min) == 0 {
					local exo_tv "`exo_tv' `v'"
					replace `n_tv_var' = `n_tv_var' + (1 - `test') if `touse'
				}
			}
 		}	
*Identify time-varying endogenous variables.
		if "`endo'" != "" {
			quietly foreach v of local endo { 
				bysort `touse' `groupvar' (`v') : replace `test' = `v'[1] == `v'[_N] if `touse'
				su `test'  if `touse', meanonly 
				if r(min) == 0 {
					local endo_tv "`endo_tv' `v'"
					replace `n_tv_var' = `n_tv_var' + (1 - `test') if `touse'
				} 	
			}
		}
*n_small is at equal to the number of time-varying variables + 1
		tempvar n_by_group		
		quiet egen  `n_by_group' = count(`touse'), by (`groupvar')
		quiet replace `touse' = 0 if `n_by_group' < `n_tv_var'	
*we redo the group variable		
		quiet drop `groupvar'
		quietly egen `groupvar' = group(`group') if `touse'
		quietly su `groupvar'  if `touse', meanonly
		local ngroup = r(max)
*Check that the instruments are linear functions of the 1st stage regressors
*Maybe it is too much to do that.
		if "`inst'" != "" {
			quietly foreach v of local inst { 
				forvalues i = 1/`ngroup' {					
					`vv' regress `v' `exo_tv' `endo_tv' [`weight'`exp'] if `groupvar' == `i' & `touse'
					if e(r2) < 0.999999 {
						display as error "The instrumental variables must be linear functions of the time-varying variables."
						error 498
					} 			
				}
			}
		}		
*clusters
		if "`cluster'" == ""{
			*clustering by group if cluster variable is missing
			local cluster "`groupvar'"
		}
		else {
			*check that the cluster variable is constant within groups.
			quietly bysort `touse' `groupvar' (`cluster') : replace `test' = `cluster'[1] == `cluster'[_N] if `touse'
			quiet su `test'  if `touse', meanonly 
			if r(min) == 0 {
				display as error "The clusters must nest the groups, otherwise the s.e. will be incorrect."
				error 498
			} 			
		}
*The number of instruments must be at least as large as the number of endogenous variables. 
		if wordcount("`endo'") > wordcount("`inst'") {
			dis as error "The model is not identified; there are fewer instrumental variables than endogenous variables."
			error 481
		}
*Estimation command
		if "`est_command'" == "" {
			if "`endo'" == ""{
				local est_command "regress"
			}
			else if wordcount("`endo'") == wordcount("`inst'") {
				if `stata_version' < 10 {
					local est_command "ivreg "
				}
				else {
					local est_command "ivregress 2sls"
				}
			}
			else {
				if `stata_version' < 10 {
					local est_command "ivreg2 "
					local est_opts "`est_opts' gmm "
				}
				else {
					local est_command "ivregress gmm "
				}
			}
		}
*VCE cannot be changed (to make sure that the s.e. are clustered!)
		if strpos("`est_opts'", "vce")>0{
			dis as error "The option vce() cannot be modified with est_opts."
			dis as error "The standard errors must be clustered. You can change the cluser variable with the option cluster."
			dis as error "The bootstrap can be activated with the option bootstrap."
			error 198
		}
*The cluster variable must be specified directly with the option cluster (for simplicity).		
		if (strpos("`bootstrap'", "cluster")>0 | strpos("`bootstrap'", "cl(")>0 | strpos("`bootstrap'", "clu(")>0 | strpos("`bootstrap'", "clus(")>0 | strpos("`bootstrap'", "clust(")>0 | strpos("`bootstrap'", "cluste(")>0){
			dis as error "The cluster variable cannot be specified with the option bootstrap."
			dis as error "You can change the cluser variable with the option cluster."
			error 198
		}		
*weights
		local isweight=("`exp'"!="")
*number of observations
		quietly sum `dep' if `touse'
		local obs=r(N)
*temporary names
		tempname quants covariance coefficients coefmat
		tempvar fitted
*cleaning of the quantiles
*default quantiles: 0.1 0.25 0.5 0.75 0.9
		if "`quantiles'"==""{
			local quantiles "0.1 0.25 0.5 0.75 0.9"
		}
		mata: st_matrix("`quants'", strtoreal(tokens(st_local("quantiles")))')
		local nq=rowsof(`quants')
*1st stage Estimation
		if "`load_first'" == "" {
			forvalues q = 1/`nq'{
				if "`save_first'" == "" {
					tempvar fit`q'
				}
				else {
					local fit`q' "`save_first'`q'"
				}
				quiet gen `fit`q''=.
			}			
			if "`parallel'" == ""{
				forvalues i = 1/`ngroup' {		
					`vv' qrprocess `dep' `exo_tv' `endo_tv' [`weight'`exp'] if `groupvar' == `i', vce(novar) quantile(`quantiles') noprint `qr_opts'
					forvalues q = 1/`nq'{
						quietly predict `fitted' if `groupvar' == `i', equation(#`q')
						quietly replace `fit`q''=`fitted' if `groupvar' == `i'
						quietly drop `fitted'
					}
				}
			}
			else{
				forvalues q = 1/`nq'{
					local fit_list "`fit_list' `fit`q''"
				}
				quietly sort `groupvar'
				quietly parallel, by(`groupvar'): by `groupvar': par_qrprocess, command(`vv' qrprocess `dep' `exo_tv' `endo_tv' [`weight'`exp']) opts(vce(novar) quantile(`quantiles') noprint `qr_opts') nq(`nq') fit(`fit_list')
			}
		}
		else {
			forvalues q = 1/`nq'{
				local fit`q' "`load_first'`q'"
				confirm numeric variable `fit`q''
			}
		}
*2nd stage estimation
*If GMM, by default the cluster weighting matrix is used, if not overruled:
		if "`est_command'"=="ivregress gmm" & strpos("`est_opts'","wmat")==0{
			local est_opts "wmatrix(cluster `cluster') `est_opts'"
		}
*if bootstrap: clustered bootstrap
		if "`bootstrap'"!=""{
			*if a column in bootstrap
			if strpos("`bootstrap'",":")>0 {
				gettoken bootstrap away : bootstrap, parse(":")
			}
			*if no comma in bootstrap
			if strpos("`bootstrap'",",")==0 {
				local bootstrap "`bootstrap', cluster(`cluster'):"
			}
			else {
				local bootstrap "`bootstrap' cluster(`cluster'):"
			}
		} 
		else {
			local est_opts "cluster(`cluster') `est_opts'" 
		}
*This is peculiar but "ivregress gmm" and "ivregress 2sls" check if the instruments
*are colinear with the endogenous regressors. This actually prevents implementing 
*the random effects estimator! There is an option "perfect" to avoid this check.
		if "`est_command'" == "ivregress 2sls" | "`est_command'" == "ivregress gmm" {
			local est_opts "`est_opts' perfect"
		}
*estimation
		forvalues q = 1/`nq'{
			if "`inst'" == "" {
				quietly `vv' `bootstrap' `est_command' `fit`q'' `exo' if `touse' [`weight'`exp'], `est_opts'
			} 
			else {
				quietly `vv' `bootstrap' `est_command' `fit`q'' `exo' (`endo' = `inst') if `touse' [`weight'`exp'], `est_opts'
			}
			if `q' == 1 {
				mata: st_local("k", strofreal(cols(st_matrix("e(b)"))))
				mata: mdqr_coef_mat = J(cols(st_matrix("e(b)")), `nq', .)
				mata: mdqr_var_mat = J(cols(st_matrix("e(b)")) * `nq', cols(st_matrix("e(b)")) * `nq', 0)
			}
			mata: mdqr_coef_mat[.,`q'] = st_matrix("e(b)")'
			mata: mdqr_var_mat[(`k'*(`q'-1)+1)..(`k'*`q'),(`k'*(`q'-1)+1)..(`k'*`q')] = st_matrix("e(V)")'
		}
		mata: st_matrix("`coefmat'", mdqr_coef_mat)
		mata: st_matrix("`coefficients'", vec(mdqr_coef_mat)')
		mata: st_matrix("`covariance'", mdqr_var_mat)
*Degrees of freedom (not important and not well-defined in theory but used by plotprocess)
		if "`e(df_r)'" == ""{
			local df_r = e(N) - e(df_m) -1
		}
		else{
			local df_r = e(df_r)
		} 
*Names of the columns of the matrices.
		local variable_names : colnames e(b)
		capture{
			tokenize "`variable_names'"
			local nq=rowsof(`quants')
			forvalues j = 1/`nq' {
				local i=1
				while "``i''" != "" {
					local conams "`conams' ``i''"
					local eqnams "`eqnams' q`j'"
					local i = `i' + 1
				}
 			}
		}
		if _rc{
			dis in red "Due to the high number of regressors and quantiles the rows and columns of the matrices cannot be named correctly."
			tokenize "`variable_names'"
			local nq=rowsof(`quants')
			local conams ""
			local eqnams ""
			forvalues j = 1/`nq' {
				local i=1
				while "``i''" != "" {
					local conams "`conams' x`i'"
					local eqnams "`eqnams' q`j'"
					local i = `i' + 1
				}
 			}
		}
		`vv' mat colnames `coefficients' = `conams'
		mat coleq `coefficients' = `eqnams'
		mat rownames `coefficients'=`dep'	
		`vv' mat rownames `covariance' = `conams'
		mat roweq `covariance' = `eqnams'
		`vv' mat colnames `covariance' = `conams'
		mat coleq `covariance' = `eqnams'
		`vv' mat rownames `coefmat' = `variable_names'
		`vv' mat colnames `coefmat' = `ueqnams'
		mat rownames `quants' = `ueqnams'
		mat colnames `quants' = "tau"
*Print the results		
		if "`print'"!="noprint"{
			dis
			dis as text _column(0) "Minimum distance quantile regression"
			dis as text _column(0) "No. of obs." _c
			dis as result _column(30) %-8.0f `obs'
			dis as text _column(0) "1st stage estimation : " _c
			dis as result _column(30) %-8.0f "qrprocess"
			dis as text _column(0) "2nd stage estimation : " _c
			dis as result _column(30) %-8.0f "`est_command'"
			dis as text _column(0) "Variance : " _continue
			if "`bootstrap'"=="bootstrap"{
				dis as result _column(30) %-8.0f "clustered bootstrap"
			}
			else {
				dis as result _column(30) %-8.0f "clustered robust"
			}
			dis
			dis as text "{hline 13}" "{c TT}" "{hline 64}"
			dis as text _column(0) %11s "`dep'" _column(14) "{c |}"  _column(19) %~10s "Coef." _column(29) %~10s "Std. Err." _column(41) %~8s "t" _column(45) %~8s "P>|t|" _column(59) %~8s "[`level'% Conf. Interval]"
			dis as text "{hline 13}" "{c +}" "{hline 64}"
			forvalues i=1/`nq'{
				dis as text _column(0) "Quant. " as result "0" as result %-5.0g `quants'[`i',1]  _column(14) "{c |}" 
				forvalues j=1/`k'{
					dis as text _column(0) %11s abbrev(word("`conams'",`j'),12) _column(14) "{c |}" as result _column(17) %9.0g (`coefmat'[`j',`i']) _continue
					dis as result _column(28) %9.0g (`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
					dis as result _column(40) %5.2f `coefmat'[`j',`i']/(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
					dis as result _column(49) %4.3f 2-2*normal(abs(`coefmat'[`j',`i']/(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5)) _continue
					dis as result _column(58) %9.0g `coefmat'[`j',`i']+invnormal((100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5 _continue
					dis as result _column(70) %9.0g `coefmat'[`j',`i']-invnormal((100-`level')/200)*(`covariance'[(`i'-1)*`k'+`j',(`i'-1)*`k'+`j'])^0.5							
				}
				dis as text "{hline 13}" "{c BT}" "{hline 64}"
			}
		}
*Post the results
		ereturn post `coefficients' `covariance', depname(`dep') obs(`obs') esample(`touse')
		ereturn matrix coefmat `coefmat'
		ereturn matrix quantiles `quants'
*		ereturn local estat_cmd "qrprocess_estat"
*		ereturn local predict "qrprocess_p"
		if "`bootstrap'"!=""{
			ereturn local vce "bootstrap"
		}
		else {
			ereturn local vce "cluster"
		}
		ereturn local method `"`est_command'"'
		ereturn local clustvar `"`cluster'"'
		ereturn scalar N_groups=`ngroup'
		ereturn scalar df_r=`df_r'
		if `isweight'==1{
			ereturn local wtype `"`weight'"'
			ereturn local wexp `"`exp'"'
		}
		*drop "_cons" from variable names
		local variable_names = subinword("`variable_names'","_cons","",1)
		ereturn local xvar `"`variable_names'"'
		ereturn local depvar `"`dep'"'
		ereturn local cmdline `"mdqr `0'"'
		ereturn local title "Minimum distance quantile regression"
		ereturn local cmd "mdqr"
		ereturn local plotprocess "e(quantiles) Quantile"
		cap mata: mata drop mdqr_coef_mat
		cap mata: mata drop mdqr_var_mat
	}
end
