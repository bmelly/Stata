*mdqr: minimum distance quantile regression
*! version 0.0.7  19.01.2023  Blaise Melly
*Estimation in Mata

*To do (potentially): *(2) Covariance between coefficients are different quantiles (multi-variate GMM?)
*(4) bootstrap all the quantiles together to obtain valid covariances.

cap prog drop mdqr
program mdqr, eclass byable(recall) sortpreserve
	local stata_version = _caller() 
	version 12
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
*Generate the fitted values for the first stage (objective: get an early error message if variables with the same name already exist).
		forvalues q = 1/`nq'{
				if "`save_first'" == "" {
					tempvar fit`q'
				}
				else {
					local fit`q' "`save_first'`q'"
				}
				quiet gen `fit`q''=.
			}			
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
*			if `stata_version' >= 11{
				version 11.1: fvexpand `exo'
				local exo `r(varlist)'
				if "`r(fvops)'" == "true"{
					local fvops = 1
				}
/*			}
			else if "`exo'" != ""{ 
				unab exo : `exo'
				confirm numeric variable `exo'
			}*/
			gettoken away anything : anything, parse("(")
		}
*separate endogenous variables
		gettoken endo anything : anything, parse("=")		
		if "`endo'" == "=" {
			local endo ""
		}
		else {
*			if `stata_version' >= 11{
				version 11.1: fvexpand `endo'
				if "`r(fvops)'" == "true"{
					local fvops = 1
					local endo `r(varlist)'
				}
/*			}
			else if "`endo'" != ""{
				confirm numeric variable `endo' 
			}*/
			gettoken away anything : anything, parse("=")
		}
*separate instruments		
		gettoken inst anything : anything, parse(")")		
		if "`inst'" == ")" {
			local "`inst'" ""
		}
		else {
*			if `stata_version' >= 11{
				version 11.1: fvexpand `inst'
				if "`r(fvops)'" == "true"{
					local fvops = 1
					local inst `r(varlist)'
				}
/*			}
			else if "`inst'" != ""{
				confirm numeric variable `inst' 
			}*/
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
		tempvar n_tv_var test temp
		quiet gen `n_tv_var' = 1 if `touse'
		quietly gen `test' = .	
		quietly gen `temp' = .
		local potentially_tv "`exo' `endo'"
*if factor variables:  generate temporary variables for the factors
		if `fvops' {
			`vv' _fv_check_depvar `dep'
			`vv'  quiet _rmcoll `potentially_tv ' [`weight'`exp'] if `touse'
			local potentially_tv "`r(varlist)'"
			`vv' fvexpand `potentially_tv ' if `touse'
			local varlistt "`r(varlist)'"
			local potentially_tv  ""
			local i=1
			foreach vn of local varlistt{						
				`vv' _ms_parse_parts `vn'
				if ~`r(omit)'{
					if r(type)=="variable"{
						local potentially_tv  `potentially_tv ' `vn'
					}
					else{
						tempvar factor`i'
						quiet gen `factor`i''=`vn' if `touse'
						local potentially_tv  `potentially_tv ' `factor`i''
					}
					local i=`i'+1
				}
			}
		}
local tv "`potentially_tv'"		
		/*
*time-varying variables		
		if "`potentially_tv'" != ""{
*			local pastv ""
*			local factor = 0
*			foreach v of local potentially_tv {
*				if `fvops' == 1{
*					unopvarlist `v'
*					if "`v'" != "`r(varlist)'"{
*						local factor = `factor' + 1
*						local v `r(varlist)'
*						if "`v'" == "`pastv'"{						
*							quietly replace `n_tv_var' = `n_tv_var' + (1 - `test') if `touse'
*							continue
*						}
*					}
*					else{
*						local factor = 0
*					}
*					local pastv "`v'"
*				}
			foreach v of local potentially_tv {
				quietly bysort `touse' `groupvar' (`v') : replace `test' = (`v'[1] == `v'[_N]) if `touse'
				quietly su `test'  if `touse', meanonly
				if r(min) == 0 {
					local tv "`tv' `v'"
					quietly replace `n_tv_var' = `n_tv_var' + (1 - `test') if `touse'			
				}
			}
		}
		
		tempvar n_by_group new_touse group2
		quiet egen  `n_by_group' = count(`touse'), by (`groupvar')
*		quiet gen `new_touse' = 1 if `n_by_group' >= `n_small' + `n_tv_var' & `touse'
		quiet replace `touse' = 0 if `n_by_group' < `n_small' + min(`n_tv_var', 2)
*we redo the group variable		
*/
		quiet drop `groupvar'
		quietly egen `groupvar' = group(`group') if `touse'
		quietly su `groupvar'  if `touse', meanonly
		local ngroup = r(max)		
*Check that the instruments are linear functions of the 1st stage regressors
*Maybe it is too much to do that.
		if "`inst'" != "" {
			quietly foreach v of local inst { 
				forvalues i = 1/`ngroup' {					
					quietly `vv' regress `v' `tv' [`weight'`exp'] if `groupvar' == `i' & `touse'
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
/*				if `stata_version' < 10 {
					local est_command "ivreg "
				}
				else {*/
					local est_command "ivregress 2sls"
*				}
			}
			else {
/*				if `stata_version' < 10 {
					local est_command "ivreg2 "
					local est_opts "`est_opts' gmm "
				}
				else {*/
					local est_command "ivregress gmm "
*				}
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
*1st stage Estimation
		if "`load_first'" == "" {
			forvalues q = 1/`nq'{
				local fit_list "`fit_list' `fit`q''"
			}
*weights
			if "`weight'"==""{
				tempvar exp1
				quietly gen byte `exp1' = 1 if `touse'
			}
			else{
				gettoken away exp1 : exp
			}	
			if "`parallel'" == ""{
				mata: mdqr_fn("`dep'", "`tv'", "`exp1'", "`touse'", "`groupvar'", "`quants'", 0.9995, 0.00001, 100, "`fit_list'", `n_small')
			}
			else{
				quietly sort `groupvar', stable
				quietly parallel, by(`groupvar'): par_qrprocess,  dep(`dep') reg(`tv') weight(`exp1') touse(`touse') group(`groupvar') quantiles(`quantiles') fitted(`fit_list') n_small(`n_small')
			}
		}
		else {
			forvalues q = 1/`nq'{
				local fit`q' "`load_first'`q'"
				confirm numeric variable `fit`q''
			}
		}
*number of observations
		quietly sum `dep' if `touse'
		local obs=r(N)
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

cap mata mata drop mdqr_fn()
version 12
mata void mdqr_fn(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar group, string scalar quantile, real scalar beta, real scalar small, real scalar max_it, string scalar fitted, real scalar n_small)
{
	real colvector quants, y, w, uy
	real rowvector conv, rdev, mdev
	string rowvector reg1
	real matrix x, coef, cov
	real scalar alpha, i, n
	quants = st_matrix(quantile)
	nq = rows(quants)
	y = st_data(., dep, touse)
	n = rows(y)
	int_touse = J(n, 1, 1)
	reg1 = tokens(reg)
	if(length(reg1) > 0) x = J(n, 1, 1), st_data(., reg1, touse)
	else x = J(n, 1, 1)
	k = cols(x)
	w = st_data(., weight, touse)
	g = st_data(., group, touse)	
	fit = J(n, nq, .)
	ngroup = max(g)
	conv = 0
	coef = 0
	for(i = 1; i <= ngroup; i++){
		temp_g = my_selectindex(g :== i)
		if(rows(temp_g) < n_small+1){
			int_touse[temp_g] = J(rows(temp_g), 1, 0)
		}
		else{
			temp_x_all = x[temp_g, ]
			temp_minmax = colminmax(temp_x_all[, (2..k)])
			select_var = (1, select(2..k, temp_minmax[1,] :< temp_minmax[2,]))	
			hqrdp(temp_x_all[, select_var], temp1=., temp2=., temp3=., p = (1, J(1, cols(select_var)-1, 0)))	
			temp_x = temp_x_all[, select_var][, p[1..sum(abs(diagonal(temp3)):>10^(-9))]]
			if(rows(temp_x) - cols(temp_x) < n_small){
				int_touse[temp_g] = J(rows(temp_x), 1, 0)
			}
			else{
				temp_y = y[temp_g]
				temp_w = w[temp_g]
				for (q = 1; q <= rows(quants); q++) {			
					if(q==1) coef=0			
					coef = int_fnm(temp_y, temp_x, temp_w, quants[q, 1], beta, small, max_it, conv, coef)
					if(colmissing(coef) | conv==0){	
						stata("_qreg "+dep+" "+reg+" if "+group+"=="+strofreal(i)+" [aweight="+weight+"], quantile("+strofreal(quants[q,1])+")", 1)			
						coef=st_matrix("e(b)")'
						coef=coef[(rows(coef), 1..(rows(coef)-1)), 1]						
						fit[temp_g, q] = temp_x_all * coef
						coef=0
					}
					else{
						fit[temp_g, q] = temp_x * coef
					}
				}
			}
		}	
	}
	st_store(., tokens(fitted), touse, fit)
	st_store(., touse, touse, int_touse)
}

cap mata mata drop int_fnm()
version 12
mata real colvector int_fnm(real colvector dep, real matrix X, real colvector wei, real scalar p, real scalar beta, real scalar small, real scalar max_it, real scalar convergence, real colvector start)
{
	real colvector weight, c, b, x, s, y, r, z, w, q, rhs, dy, dx, ds, dz, dw, fx, fs, fw, fz, fp, fd, dxdz, dsdw, xinv, sinv, xi
	real scalar n, gap, it, mu, g
	real matrix A, AQ, invAQ
	weight=wei:/mean(wei)
	n=rows(X)
	A=(X:*weight)
	c=-(dep:*weight)
	b=colsum(A):*(1-p)
	x=J(n,1,1-p)
	s=J(n,1,p)
	if(start==0){
		y = (invsym(cross(A, A))*cross(A,c))
		y[cols(X)]=y[cols(X)]+mm_quantile(c - cross(A' , y),1,p)
	}
	else{
		y=start
	}
	r = c - cross(A' , y)
	r = r + 0.001 * (r :== 0)
	z = r :* (r :> 0)
	w = z - r
	it = 0
	while(it < max_it){
    	it = it + 1
    	q=1:/(z:/x+w:/s)
    	r = z - w
		AQ = A:*sqrt(q)
		rhs=r:*sqrt(q)
		invAQ=invsym(cross(AQ, AQ))
		dy = invAQ*cross(AQ,rhs)
		dx = q :* (A*dy - r)    
		ds = -dx
		dz = -z :* (1 :+ dx:/x)
		dw = -w :* (1 :+ ds:/s)
		fx = bound(x, dx)
		fs = bound(s, ds)
		fw = bound(w, dw)
		fz = bound(z, dz)
		fp = rowmin((fx, fs))
		fd = colmin((fw, fz))
		fp = min((beta * fp\ 1))
		fd = min((beta * fd, 1))
		if(min((fp, fd)) < 1){
			mu = z ' x + w ' s
			g = (z + fd * dz) ' (x + fp * dx) + (w + fd * dw) ' (s + fp * ds)
			mu = mu * (g / mu) ^3 / ( 2* n)
			dxdz = dx :* dz
			dsdw = ds :* dw
			xinv = 1 :/ x
			sinv = 1 :/ s
			xi = mu * (xinv - sinv)
			rhs = rhs + sqrt(q):*(dxdz - dsdw :- xi)
			dy = invAQ*cross(AQ,rhs)
			dx = q :* (A*dy  - r -dxdz + dsdw:+ xi)
			ds = -dx
			dz = mu * xinv - z - xinv :* z :* dx - dxdz
			dw = mu * sinv - w - sinv :* w :* ds - dsdw
 			fx = bound(x, dx)
			fs = bound(s, ds)
			fw = bound(w, dw)
			fz = bound(z, dz)
			fp = rowmin((fx, fs))
			fd = colmin((fw, fz))
			fp = min((beta * fp\ 1))
			fd = min((beta * fd, 1))
		}  
		x = x + fp * dx
		s = s + fp * ds
		gap=fd * dy
		y = y + gap
		if(max(abs(gap)) < small){
			if(c'x-b*y+sum(w)<small){
				break
			}
		}
		w = w + fd * dw
		z = z + fd * dz
	}
	convergence=(it < max_it)
	return(-y)
}

cap mata mata drop bound()
version 12
mata real colvector bound(real colvector x, real colvector dx)
{
	real colvector b, f
	f=(dx:>=0)
	b=f:*1e20-(1:-f):*x:/dx
	return(b)
}

cap mata mata drop my_selectindex()
version 12
mata real colvector my_selectindex(real colvector input) return(select((1..rows(input))', input))
