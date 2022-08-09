program define drprocess_p
	version 6, missing

		/* Step 1:
			place command-unique options in local myopts
			Note that standard options are
			LR:
				Index XB Cooksd Hat 
				REsiduals RSTAndard RSTUdent
				STDF STDP STDR noOFFset
			SE:
				Index XB STDP noOFFset
		*/
	local myopts "Probability Residuals Difference REArranged(numlist min=1 sort) Quantile(numlist >0 <1 min=1 sort)"

		/* Step 2:
			call _propts, exit if done, 
			else collect what was returned.
		*/
	_pred_me "`myopts'" `0'
	if `s(done)' { exit }
	local vtyp  `s(typ)'
	local varn `s(varn)'
	local 0 `"`s(rest)'"'


		/* Step 3:
			Parse your syntax.
		*/
	syntax [if] [in] [, `myopts' noOFFset EQuation(string)]
	if "`rearranged'"!=""{
		local cfg="rearrangement"
	}
	else if "`quantile'"!=""{
		local cfg="quantile"
	}
	else if "`probability'"!=""{
		local cfg "pr"
	}

	if "`equatio'" == "" {	    /* we're version 6 --  7-char locals! */
		tempname b
		mat `b' = get(_b)
		local eqnames : coleq `b'
		gettoken equatio : eqnames
	}

		/* Step 4:
			Concatenate switch options together
		*/
	local type "`residua'`differe'`cfg'"
	

		/* Step 5:
			quickly process default case if you can 
			Do not forget -nooffset- option.
		*/
		
	if "`type'"=="" | "`type'"=="pr" | "`type'"=="xb"{
		tempname thresh	
		mat `thresh'=e(thresholds)
		if "`type'"=="" {
			di in smcl in gr ///
			"(option {bf:pr} assumed; Pr(`e(depvar)'))"
			local type "pr"
		}
		if "`equation'"==""{
			local equation "#1"
		}	
		local nq=wordcount("`equation'")
		if `nq'>1{
			forvalues i=1/`nq'{
				local names "`names' `varn'`i'"
			}
		}
		else{
			local names "`varn'"
		}
		tokenize "`equation'", parse(" ")
		foreach name of local names{
			local eqnum=substr("`1'",2,100)
			local temp=`thresh'[`eqnum',1]
			_predict `vtyp' `name' `if' `in', xb `offset' equation(`1')
			if "`type'"=="xb"{
				label var `name' "linear index at `temp'"
			}
			if "`type'"=="pr"{
				if "`e(method)'"=="logit"{
					quiet replace `name'=1/(1+exp(-`name'))
				}
				else if "`e(method)'"=="probit"{
					quiet replace `name'=normal(`name')
				}
				else if "`e(method)'"=="cloglog"{
					quiet replace `name'=1-exp(-exp(`name'))
				}
				label var `name' "Pr(`e(depvar)'<=`temp')"
			}
			mac shift
		}
		exit
	}

		/* Step 6:
			mark sample (this is not e(sample)).
		*/
	marksample touse
	local method "`e(method)'"
	if "`type'"=="rearrangement" {
		tempname quants
		local nq=wordcount("`rearranged'")
		if `nq'>1{
			forvalues i=1/`nq'{
				quietly gen `vtyp' `varn'`i'=.
				local pn "`pn' `varn'`i'"
			}
			tokenize "`rearranged'", parse(" ")
			local i=1
			while "`1'" != "" {
				matrix `quants'=nullmat(`quants')\(`1')
				mac shift 
				local i=`i'+1
			}
		}
		else{
				quietly gen `vtyp' `varn'=.
				local pn "`varn'"
				matrix `quants'=`rearranged'				
		}
		mata: dr_rearranged("`e(xvar)'", "`touse'", "`quants'", "`pn'", "`method'")
		if `nq'==1 {
			di "Rearranged conditional distribution at the `rearranged' threshold" 
		}
		else{
			dis "Rearranged conditional distribution functions for the following thresholds:"
			forvalues i=1/`nq'{
				local temp=`quants'[`i',1]
				dis "The conditional cdf at `temp' is saved in `varn'`i'"
			}
		}
		exit
	}

		if "`type'"=="quantile" {
			tempname eval
			local nq=wordcount("`quantile'")
			if `nq'>1{
				forvalues i=1/`nq'{
					quietly gen `vtyp' `varn'`i'=.
					local pn "`pn' `varn'`i'"
				}
				tokenize "`quantile'", parse(" ")
				local i=1
				while "`1'" != "" {
					matrix `eval'=nullmat(`eval')\(`1')
					mac shift 
					local i=`i'+1
				}
			}
			else{
				quietly gen `vtyp' `varn'=.
				local pn "`varn'"
				matrix `eval'=`quantile'
			}
			mata: dr_quantile("`e(xvar)'", "`touse'", "`eval'", "`pn'", "`method'")
			if `nq'==1 {
				di "Conditional `quantile'th quantile function." 
			}
			else{
				dis "Conditional quantile functions at the following quantile:"
				forvalues i=1/`nq'{
					local temp=`eval'[`i',1]
					dis "`temp'th quantile is saved in `varn'`i'"
				}
			}
			exit
		}



		/* Step 8:
			handle switch options that can be used in-sample or 
			out-of-sample one at a time.
			Be careful in coding that number of missing values
			created is shown.
			Do all intermediate calculations in double.
		*/

	if "`type'"=="residuals" {
		tempvar pr thresh
		qui _predict double `pr' if `touse', xb `offset' eq(`equatio')
		if "`e(method)'"=="logit"{
			quiet replace `pr'=1/(1+exp(-`pr'))
		}
		else if "`e(method)'"=="probit"{
			quiet replace `pr'=normal(`pr')
		}
		else if "`e(method)'"=="cloglog"{
			quiet replace `pr'=1-exp(-exp(`pr'))
		}
		mat `thresh'=e(thresholds)
		cap mat `thresh'=`thresh'["`equatio'",1]
		if _rc>0{
			mat `thresh'=`thresh'[real(substr("`equatio'",2,100)),1]
		}
		gen `vtyp' `varn' = ((`e(depvar)'<=`thresh'[1,1])-`pr')/sqrt(`pr'*(1-`pr')) if `touse'
		label var `varn' "Pearson residuals:  `equatio'"
		exit		
	}
	tokenize "`equatio'", parse(",") 
	local eq1 `"`1'"'
	local eq2 `"`3'"'
	
	if "`type'"=="difference" {
		tempvar xb1 xb2
		qui _predict double `xb1' if `touse', `offset' eq(`eq1')
		qui _predict double `xb2' if `touse', `offset' eq(`eq2')
		gen `vtyp' `varn' = `xb1' - `xb2' if `touse'
		lab var `varn' "Fitted diff.:  `eq1' - `eq2'"
		exit
	}
	
		/* Step 9:
			handle switch options that can be used in-sample only.
			Same comments as for step 8.
		*/
	*qui replace `touse'=0 if !e(sample)


			/* Step 10.
				Issue r(198), syntax error.
				The user specified more than one option
			*/
	error 198
end

*Mata function doing the rearrangement
version 9.2
mata void dr_rearranged(string scalar reg, string scalar touse, string scalar quant_wished, string scalar out, string scalar method)
{
	real colvector quants, pred_quants, q_index
	real scalar nw, i, n
	real matrix x, coef, fit
	quants=st_matrix("e(thresholds)")
	pred_quants=st_matrix(quant_wished)
	nw=rows(pred_quants)
	q_index=J(nw,1,.)
	for(i=1;i<=nw;i++){
//		if(min(abs(quants:-pred_quants[i]))>0.000001){
//			"The "+strofreal(pred_quants[i])+" threshold was not included in the estimation."
//			"Include this threshold when you call drprocess if you want ot obtaint the predicted values at this threshold."
//			exit(400)
//		}
		q_index[i,1]=sum(quants:<=pred_quants[i])
	}
	n = sum(st_data(., touse, touse))
	if(reg!=""){
		x = get_reg(reg, touse, n), J(n,1,1)
	} else {
		x=J(n,1,1)
	}
	coef=st_matrix("e(coefmat)")
	fit=cross(x',coef)'
	for(i=1;i<=rows(fit);i++){
		fit[i,.]=colmin(colmax((J(1,n,-20)\fit[i,.]))\J(1,n,20))
	}
	if (method=="logit") fit=1:/(1:+exp(-fit))
	else if (method=="probit")  fit=normal(fit)
	else if (method=="cloglog")  fit=1:-exp(-exp(fit))
	for(i=1;i<=n;i++){
		fit[.,i]=sort(fit[.,i],-1)
	}
	st_store(.,tokens(out),touse,fit[q_index,.]')
}

*Mata function calculating the cdf
version 9.2
mata void dr_quantile(string scalar reg, string scalar touse, string scalar eval, string scalar out, string scalar method)
{
	real colvector thresh, quants, evalu
	real scalar nw, n, i, j
	real matrix x, coef, fit
	thresh=st_matrix("e(thresholds)")
	evalu=st_matrix(eval)
	nw=rows(evalu)
	n = sum(st_data(., touse, touse))
	if(reg!=""){
		x = get_reg(reg, touse, n), J(n,1,1)
	} else {
		x=J(n,1,1)
	}
	coef=st_matrix("e(coefmat)")
	fit=cross(x',coef)
	if (method=="logit") fit=1:/(1:+exp(-fit))
	else if (method=="probit")  fit=normal(fit)
	else if (method=="cloglog")  fit=1:-exp(-exp(fit))
	quants=J(n,nw,.)
	for(i=1;i<=n;i++){
		for(j=1;j<=nw;j++){
			quants[i,j]=thresh[max(sum(fit[i,.]:<=evalu[j,1])\1),1]
		}
	}
	st_store(.,tokens(out),touse,quants)
}

mata real matrix get_reg(string scalar input, string scalar touse, real scalar nobs)
{
	string rowvector varlist
	real scalar nreg, i
	real matrix output
	varlist = tokens(input) 
	nreg = cols(varlist)
	output = J(nobs, nreg, .)
	for(i=1; i <= nreg; i++){		
		output[.,i] = st_data(., varlist[1,i], touse)
	}
	return(output)
}
