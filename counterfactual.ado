*Version 1.0.0 01mar2023
*Codes implementing the estimators proposed in Chernozhukov, Fernandez-Val and Melly
*Codes for QTE, pointwise and uniform confidence intervals based on bootstrap and Kolmogorov-Smirnov statistic

program counterfactual, eclass
	version 9.2
	capt findfile lmoremata.mlib
	if _rc {
      	di as error "-moremata- is required; type {stata ssc install moremata} and restart Stata."
		error 499
	}
	if replay() {
		if "`e(cmd)'"!="counterfactual" { 
			error 301 
		} 
		if _by() {
			error 190 
		}
        syntax [, Level(cilevel)]
		ereturn display, level(`level')
 	}
	else {
		syntax varlist [if] [in] [pweight/], [Group(varname) Counterfactual(varlist) Method(string) QLow(real 0.1) QHigh(real 0.9) QStep(real 0.1) Quantiles(numlist >0 <1 sort) NReg(real 100) Reps(integer 100) Level(cilevel) First(real 0.1) Last(real 0.9) noboot noprint SCale(varlist) counterscale(varlist) SAVing(string) CONS_test(string) beta(real 0.9995) small(real 0.00001) max_it(real 50) Censoring(varname) Firstc(real 0.1) Secondc(real 0.05) NSteps(integer 3) RIght est_opts(string)]
		local nreg=round(`nreg')
		if `nreg'<1{
			dis as error "The option nreg must be a strictly positive integer."
			exit
		}
		if "`group'"=="" & "`counterfactual'"==""{
			di as error "One of the options group and counterfactual are required"
			exit
		}
		if "`method'"==""{
			local method "qr"
		}
		if "`method'"!="qr" & "`method'"!="logit" & "`method'"!="probit" & "`method'"!="lpm" & "`method'"!="loc" & "`method'"!="locsca" & "`method'"!="cox" & "`method'"!="cqr"{
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
		if "`method'"=="locsca" & "`scale'"!="" & "`group'"=="" & "`counterscale'"==""{
			dis as error "If the location scale estimator is used with the option scale, then either the group option or the counterscale option must be provided."
			exit
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
		markout `touse' `group' `counterfactual' `scale' `counterscale' `censoring' 
		gettoken dep varlist : varlist
		if "`exp'"==""{
			tempvar exp
			gen `exp'=1
		}
		quiet _rmcoll `varlist' [pw=`exp'] if `touse'
		local varlist `r(varlist)'
		qte_cov_int `dep' `varlist' `if' `in' [pweight=`exp'], group(`group') counterfactual(`counterfactual') method(`method') qlow(`qlow') qhigh(`qhigh') qstep(`qstep') quantiles(`quantiles') nreg(`nreg') scale(`scale') counterscale(`counterscale') beta(`beta') small(`small') max_it(`max_it') censoring(`censoring') firstc(`firstc') secondc(`secondc') nsteps(`nsteps') `right' est_opts(`est_opts')
		tempname coefficients covariance qte distributions obs0 obsc quants constest tests fit ref nreg0
		sca `nreg0'=r(nreg)
		matrix `quants'=r(quants)
		local obs=r(obs)
		sca `obs0'=r(obs0)
		sca `obsc'=r(obsc)
		quietly gen `ref'=`touse'
		if "`group'"!=""{
			quietly replace `ref'=0 if `group'==1
		}
		mata: qte_cov_obs=mm_quantile(st_data(.,"`dep'","`ref'"),st_data(.,"`exp'","`ref'"),qte_cov_quant)'
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
			mata: qte_cov_def=qte_cov_res
			mata: qte_cov_deff=qte_cov_fitted
			mata: qte_cov_defc=qte_cov_counter
			mata: qte_cov_boot=J(0,cols(qte_cov_res)*4,.)
			preserve
			local actual_more=c(more)
			set more off
			di in gr "(bootstrapping " _c
			forvalues i=1/`reps'{
				bsample
				qte_cov_int `dep' `varlist' `if' `in' [pweight=`exp'], group(`group') counterfactual(`counterfactual') method(`method') qlow(`qlow') qhigh(`qhigh') qstep(`qstep') quantiles(`quantiles') nreg(`nreg') scale(`scale') counterscale(`counterscale') beta(`beta') small(`small') max_it(`max_it') censoring(`censoring') firstc(`firstc') secondc(`secondc') nsteps(`nsteps') `right'
				mata: qte_cov_boot=qte_cov_boot\(qte_cov_res,mm_quantile(st_data(.,"`dep'","`touse'"),st_data(.,"`exp'","`touse'"),qte_cov_quant)',qte_cov_fitted,qte_cov_counter)
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
				di in gr "." _c
			}
			set more `actual_more'
			di in gr ")"
			mata: ev_boot(qte_cov_boot, qte_cov_def, qte_cov_obs, qte_cov_deff, qte_cov_defc, qte_cov_quant, `last', `first', `level', "`qte'", "`covariance'", "`tests'", "`constest'","`distributions'","`fit'")
			mata: mata drop qte_cov_boot qte_cov_coef qte_cov_counter qte_cov_def qte_cov_defc qte_cov_deff qte_cov_fitted qte_cov_obs qte_cov_quant qte_cov_res

		}
		else{
			mata: st_matrix("`qte'",(qte_cov_res\J(3,cols(qte_cov_res),.))')
			mata: st_matrix("`distributions'",(qte_cov_obs',J(rows(qte_cov_quant),3,.),qte_cov_fitted',J(rows(qte_cov_quant),3,.),qte_cov_counter',J(rows(qte_cov_quant),3,.)))
			mata: st_matrix("`fit'",(qte_cov_obs'-qte_cov_fitted',J(rows(qte_cov_quant),3,.)))
			mata: st_matrix("`covariance'",J(cols(qte_cov_res),cols(qte_cov_res),0))
			mata: mata drop qte_cov_coef qte_cov_counter qte_cov_fitted qte_cov_obs qte_cov_quant qte_cov_res
		}
		mat `coefficients'=`qte'[....,1]'
		local nq=rowsof(`quants')
		forvalues i=1/`nq'{
			local conam "`conam' q`i'"
		}
		mat colnames `coefficients'=`conam'
		mat rownames `coefficients'=`dep'
		mat colnames `covariance'=`conam'
		mat rownames `covariance'=`conam'
		if "`boot'" == "" {
			ereturn post `coefficients' `covariance', dep(`dep') obs(`obs') esample(`touse')
			ereturn local vce "bootstrap"
		}
		else {
			ereturn post `coefficients', dep(`dep') obs(`obs') esample(`touse')
			ereturn local vce "novar"
		}
		mat rownames `qte'=`conam'
		mat colnames `qte'=qte point_se uniform_lb uniform_ub
		ereturn matrix qte=`qte'
		ereturn local cmd "counterfactual"
		ereturn local plotdeco "e(quantiles) Quantile"
		ereturn scalar obs0=`obs0'
		ereturn scalar obsc=`obsc'
		ereturn scalar nreg=`nreg'
		mat rownames `quants'=`conam'
		mat colnames `quants'=quantile
		ereturn matrix quantiles=`quants'
		mat colnames `distributions'=observed point_se uniform_lb uniform_ub fitted point_se uniform_lb uniform_ub counterfactual poit_se uniform_lb uniform_ub
		mat rownames `distributions'=`conam'
		ereturn matrix distributions=`distributions'
		mat colnames `fit'=error point_se uniform_lb uniform_ub
		mat rownames `fit'=`conam'
		ereturn matrix fit=`fit'
		if "`boot'"==""{
			mat colnames `tests'=KS CMS
			local conam "constant_0"
			if `constest'[1,1]!=0{
				local nct=rowsof(`constest')
				forvalues i=1/`nct'{
					local conam "`conam' constant_`i'"
				}
			}
			local conam "`conam' constant_m stoch_dom_pos stoch_dom_neg fit"
			mat rownames `tests'=`conam'
			ereturn matrix tests=`tests'
		}
	}
*Display the results
*header
	if "`print'"!="noprint"{
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
		dis as result _column(43) %-8.0f `obsc'
		dis
		display as text _n "Quantile Treatment Effects" 
		dis as text "{hline 12}" "{c TT}" "{hline 65}"
		dis as text _column(13) "{c |}"  _column(26) %~10s "Pointwise" _column(40) %~8s "Pointwise" _column(58) %~21s "Functional"
		dis as text _column(0) "Quantile" _column(13) "{c |}" _column(16) %~10s "QTE" _column(26) %~10s "Std. Err." _column(35) %~8s "[`level'% Conf. Interval]" _column(58) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 12}" "{c +}" "{hline 65}"
		local k=1
		tempname results quants
		mat `results'=e(qte)
		mat `quants'=e(quantiles)
		local nq=rowsof(`quants')
*Results for each quantile 
		while `k'<=`nq'{
			dis as text _column(0) %17s `quants'[`k',1] as text _column(13) "{c |}" _column(16) as result %8.0g (`results'[`k',1]) _column (26) as result %8.0g (`results'[`k',2]) _column (37) as result %8.0g (`results'[`k',1]-invnormal(0.5*(1+`level'/100))*`results'[`k',2]) _column(47) as result %8.0g (`results'[`k',1]+invnormal(0.5*(1+`level'/100))*`results'[`k',2]) _column(60) as result %8.0g (`results'[`k',3]) _column(70) as result %8.0g (`results'[`k',4]) 
			if `k'==`nq' {
				dis as text "{hline 12}" "{c BT}" "{hline 65}" 
			}
			local k=`k'+1
		}
*Tests
		if "`boot'"==""{
			mat `tests'=e(tests)
			dis
			dis as text _n "Bootstrap inference on the counterfactual quantile processes"
			dis as text "{hline 47}" "{c TT}" "{hline 30}"
			dis as text _column(48) "{c |}"  _column(60) %~10s "P-values"
			dis as text _column(0) "Null-hypothesis" _column(48) "{c |}" _column(50) %~10s "KS-statistic" _column(65) %~10s "CMS-statistic"
			dis as text "{hline 47}" "{c +}" "{hline 30}"
			dis as text _column(0) %17s "Correct specification of the parametric model" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[rowsof(`tests'),1]) _column(65) as result %8.0g (`tests'[rowsof(`tests'),2]) 
			dis as text _column(0) %17s "No effect: QTE(tau)=0 for all taus" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[1,1]) _column(65) as result %8.0g (`tests'[1,2]) 
			if `constest'[1,1]!=0{
				local nct=rowsof(`constest')
				forvalues i=1/`nct'{
					local temp=`constest'[`i',1]
					dis as text _column(0) %17s "Constant effect: QTE(tau)=`temp' for all taus" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[1+`i',1]) _column(65) as result %8.0g (`tests'[1+`i',2]) 
				}
			}
			dis as text _column(0) %17s "Constant effect: QTE(tau)=QTE(0.5) for all taus" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[rowsof(`tests')-3,1]) _column(65) as result %8.0g (`tests'[rowsof(`tests')-3,2]) 
			dis as text _column(0) %17s "Stochastic dominance: QTE(tau)>0 for all taus" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[rowsof(`tests')-2,1]) _column(65) as result %8.0g (`tests'[rowsof(`tests')-2,2]) 
			dis as text _column(0) %17s "Stochastic dominance: QTE(tau)<0 for all taus" as text _column(48) "{c |}" _column(50) as result %8.0g (`tests'[rowsof(`tests')-1,1]) _column(65) as result %8.0g (`tests'[rowsof(`tests')-1,2]) 
			dis as text "{hline 47}" "{c BT}" "{hline 30}"
		}
	}
end

*Mata function doing the evaluation of the bootstrap results
version 9.2
mata void ev_boot(numeric matrix qte_cov_booti, numeric rowvector qte_cov_def, numeric rowvector qte_cov_obs, numeric rowvector qte_cov_deff, numeric rowvector qte_cov_defc, numeric colvector qte_cov_quant, numeric scalar max, numeric scalar min, numeric scalar level, string scalar results, string scalar covariance, string scalar tests, string scalar contest, string scalar distributions, string scalar fit)
{
	qte_cov_boot=qte_cov_booti[.,1..rows(qte_cov_quant)]
	constest=st_matrix(contest)
	st_matrix(covariance,variance(qte_cov_boot))
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-qte_cov_def#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	sel=(sum(qte_cov_quant:<min)+1)..sum(qte_cov_quant:<=max)
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	lb=qte_cov_def-seuqf*Kalpha
	ub=qte_cov_def+seuqf*Kalpha
	st_matrix(results,(qte_cov_def\seuqf\lb\ub)')
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
	st_matrix(tests,test)
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
//variance of the observed distribution
	qte_cov_boot=qte_cov_booti[.,(rows(qte_cov_quant)+1)..(2*rows(qte_cov_quant))]
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-qte_cov_obs#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	se1=seuqf
	lb1=qte_cov_obs-seuqf*Kalpha
	ub1=qte_cov_obs+seuqf*Kalpha
//variance of the fitted distribution
	qte_cov_boot=qte_cov_booti[.,(2*rows(qte_cov_quant)+1)..(3*rows(qte_cov_quant))]
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-qte_cov_deff#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	se2=seuqf
	lb2=qte_cov_deff-seuqf*Kalpha
	ub2=qte_cov_deff+seuqf*Kalpha
//variance of the counterfactual distribution
	qte_cov_boot=qte_cov_booti[.,(3*rows(qte_cov_quant)+1)..(4*rows(qte_cov_quant))]
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-qte_cov_defc#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	se3=seuqf
	lb3=qte_cov_defc-seuqf*Kalpha
	ub3=qte_cov_defc+seuqf*Kalpha
//variance of the difference between the fitted and observed distributions
	qte_cov_boot=qte_cov_booti[.,(rows(qte_cov_quant)+1)..(2*rows(qte_cov_quant))]-qte_cov_booti[.,(2*rows(qte_cov_quant)+1)..(3*rows(qte_cov_quant))]
	Vuqf=diagonal(variance(qte_cov_boot))'
	Kuqf=((qte_cov_boot-(qte_cov_obs-qte_cov_deff)#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1))):^0.5
	Kmaxuqf=rowmax(Kuqf[.,sel])
	Kalpha=mm_quantile(Kmaxuqf,1,level*0.01)
	seuqf=sqrt(Vuqf)
	se4=seuqf
	lb4=qte_cov_obs-qte_cov_deff-seuqf*Kalpha
	ub4=qte_cov_obs-qte_cov_deff+seuqf*Kalpha
//test of no misspecification, KS
	KSstat=max(((qte_cov_obs-qte_cov_deff)[sel]:^2:/Vuqf[sel]):^0.5)
//test of no misspecification, CMS
	Kuqf=((qte_cov_boot-(qte_cov_obs-qte_cov_deff)#J(rows(qte_cov_boot),1,1)):^2:/(Vuqf#J(rows(qte_cov_boot),1,1)))
	Kmeanuqf=mean(Kuqf[.,sel]')
	CMSstat=mean(((qte_cov_obs-qte_cov_deff)[sel]:^2:/Vuqf[sel])')
	test=test\(mean(Kmaxuqf:>KSstat),mean(Kmeanuqf':>CMSstat))
	st_matrix(tests,test)
	st_matrix(distributions,(qte_cov_obs',se1',lb1',ub1',qte_cov_deff',se2',lb2',ub2',qte_cov_defc',se3',lb3',ub3'))
	st_matrix(fit,(qte_cov_obs'-qte_cov_deff',se4',lb4',ub4'))
}

*Generic function, counterfactual distribution using group==0 to estimate and group==1 to predict
program qte_cov_int, rclas
	version 9.2
	syntax varlist [if] [in] [pweight/] , [Group(varname) Counterfactual(varlist)  Method(string)  QLow(real 0.1) QHigh(real 0.9) QStep(real 0.1) Quantiles(numlist >0 <1 sort) NReg(real 100) scale(varlist) counterscale(varlist) beta(real 0.9995) small(real 0.00001) max_it(real 50) Censoring(varname) Firstc(real 0.1) Secondc(real 0.05) NSteps(integer 3) right est_opts(string)] 
		marksample touse
		markout `touse' `counterfactual' `group'
		tempname quants results obs0 obsc obs
		if "`method'"=="cox"{
			quietly sum `dep' if `touse'
			if r(min)<0{
				dis in red "Only positive dependent variables allowed"
				exit
			}
		}
*if a sequence has been selected: divide by hundred if higher than 1, check that the selected values make sense
		if "`quantiles'"==""{
			if `qlow'>=1{
				local qlow=`qlow'/100
				local qhigh=`qhigh'/100
				local qstep=`qstep'/100
			}
			if `qlow'>=1 | `qhigh'>=1 | `qstep'>=1 | `qlow'<=0 | `qhigh'<=0 | `qstep'<=0 | `qstep'>(`qhigh'-`qlow') | (((`qigh'-`qlow')/`qstep')-round((`qigh'-`qlow')/`qstep'))~=0{
				di in red `"Options qlow, qhigh and qstep are invalid"'
				exit
			}
			local i=0
			while `i'<=((`qhigh'-`qlow')/`qstep'){
				local temp=(`qlow'+`i'*`qstep')
				matrix `quants'=nullmat(`quants')\(`temp')
				local i=`i'+1
			}
		}
*if the quantiles have been given directly: check that they are numbers and between 0 and 1 (or between 1 and 100)
		else{
			tokenize "`quantiles'", parse(" ")
			local i=1
			while "`1'" != "" {
				matrix `quants'= nullmat(`quants') \ (`1')
				mac shift 
				local i=`i'+1
			}
		}
		gettoken dep varlist : varlist
		tempvar ref counter nreg0
		sca `nreg0'=`nreg'
		if "`group'"!=""{
			quietly gen `ref'=`touse'
			quietly replace `ref'=0 if `group'==1
			quietly gen `counter'=`touse'
			quietly replace `counter'=0 if `group'==0
			quietly sum `ref'
			sca `obs0'=r(sum)
			quietly sum `counter'
			sca `obsc'=r(sum)
			sca `obs'=`obs0'+`obsc'
			if "`method'"=="qr"{
				local qlow = 0.5/`nreg0'
				local qhigh = 1-0.5/`nreg0'
				local qstep = 1/`nreg0'
				qrprocess `dep' `varlist' if `ref' == 1 [pw=`exp'], vce(novar) qlow(`qlow') qhigh(`qhigh') qstep(`qstep') noprint `est_opts'
				mata: qte_cov_coef = st_matrix("e(quantiles)")'\st_matrix("e(coefmat)")
				mata: qte_cov_coef = est_qr("`dep'", "`varlist'", "`exp'","`ref'","`nreg0'",`beta',`small',`max_it')
				mata: qte_cov_fitted = rqpred("`quants'", "`varlist'", "`exp'","`ref'",qte_cov_coef)
				mata: qte_cov_counter = rqpred("`quants'", "`varlist'", "`exp'","`counter'",qte_cov_coef)
			}
			if "`method'"=="cqr"{
				if "`right'"==""{	
					local right=0
				}
				else{
					local right=1
				}
				mata: qte_cov_coef=est_cqr("`dep'", "`censoring'", "`varlist'", "`exp'", "`ref'", "`nreg0'", `firstc', `secondc', `nsteps', `right', `beta',`small',`max_it')
				mata: qte_cov_fitted=rqpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef)
				mata: qte_cov_counter=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef)
			}
			if "`method'"=="logit"{
				mata: qte_cov_coef=est_logit("`dep'","`varlist'","`exp'","`ref'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef,1)
				mata: qte_cov_counter=distpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef,1)
			}
			if "`method'"=="probit"{
				mata: qte_cov_coef=est_probit("`dep'","`varlist'","`exp'","`ref'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef,2)
				mata: qte_cov_counter=distpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef,2)
			}
			if "`method'"=="lpm"{
				mata: qte_cov_coef=est_lpm("`dep'","`varlist'","`exp'","`ref'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef,3)
				mata: qte_cov_counter=distpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef,3)
			}
			if "`method'"=="loc"{
				mata: qte_cov_coef=est_loc("`dep'","`varlist'","`exp'","`ref'","`nreg0'")
				mata: qte_cov_fitted=rqpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef)
				mata: qte_cov_counter=rqpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef)
			}
			if "`method'"=="locsca"{
				if "`scale'"==""{
					local scale "`varlist'"
				}
				mata: qte_cov_coef=est_locsca("`dep'","`varlist'","`scale'","`exp'","`ref'","`nreg0'")
				mata: qte_cov_fitted=lspred("`quants'","`varlist'","`scale'","`exp'","`ref'",qte_cov_coef)
				mata: qte_cov_counter=lspred("`quants'","`varlist'","`scale'","`exp'","`counter'",qte_cov_coef)
			}
			if "`method'"=="cox"{
				mata: qte_cov_coef=est_cox("`dep'","`varlist'","`exp'","`ref'")
				mata: qte_cov_fitted=coxpred("`quants'","`varlist'","`exp'","`ref'",qte_cov_coef)
				mata: qte_cov_counter=coxpred("`quants'","`varlist'","`exp'","`counter'",qte_cov_coef)
			}
		} 
		else{
			quietly sum `reg' if `touse'
			sca `obs0'=r(N)
			quietly sum `counterfactual' if `touse'
			sca `obsc'=r(N)
			quietly sum `dep' if `touse'
			sca `obs'=r(N)
			if "`method'"=="qr"{
				local qlow = 0.5/`nreg0'
				local qhigh = 1-0.5/`nreg0'
				local qstep = 1/`nreg0'
				qrprocess `dep' `varlist' if `touse' == 1 [pw=`exp'], vce(novar) qlow(`qlow') qhigh(`qhigh') qstep(`qstep') noprint `est_opts'
				mata: qte_cov_coef = st_matrix("e(quantiles)")'\st_matrix("e(coefmat)")
				mata: qte_cov_fitted=rqpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef)
				mata: qte_cov_counter=rqpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef)
			}
			if "`method'"=="cqr"{
				if "`right'"==""{	
					local right=0
				}
				else{
					local right=1
				}
				mata: qte_cov_coef=est_cqr("`dep'", "`censoring'", "`varlist'", "`exp'", "`touse'", `nreg0', `firstc', `secondc', `nsteps', `right', `beta',`small',`max_it')
				mata: qte_cov_fitted=rqpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef)
				mata: qte_cov_counter=rqpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef)
			}
			if "`method'"=="logit"{
				mata: qte_cov_coef=est_logit("`dep'","`varlist'","`exp'","`touse'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef,1)
				mata: qte_cov_counter=distpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef,1)
			}
			if "`method'"=="probit"{
				mata: qte_cov_coef=est_probit("`dep'","`varlist'","`exp'","`touse'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef,2)
				mata: qte_cov_counter=distpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef,2)
			}
			if "`method'"=="lpm"{
				mata: qte_cov_coef=est_lpm("`dep'","`varlist'","`exp'","`touse'","`nreg0'")
				mata: qte_cov_fitted=distpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef,3)
				mata: qte_cov_counter=distpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef,3)
			}
			if "`method'"=="loc"{
				mata: qte_cov_coef=est_loc("`dep'","`varlist'","`exp'","`touse'","`nreg0'")
				mata: qte_cov_fitted=rqpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef)
				mata: qte_cov_counter=rqpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef)
			}
			if "`method'"=="locsca"{
				if "`scale'"==""{
					local scale "`varlist'"
				}
				if "`counterscale'"==""{
					local counterscale "`counterfactual'"
				}
				mata: qte_cov_coef=est_locsca("`dep'","`varlist'","`scale'","`exp'","`touse'","`nreg0'")
				mata: qte_cov_fitted=lspred("`quants'","`varlist'","`scale'","`exp'","`touse'",qte_cov_coef)
				mata: qte_cov_counter=lspred("`quants'","`counterfactual'","`counterscale'","`exp'","`touse'",qte_cov_coef)
			}
			if "`method'"=="cox"{
				mata: qte_cov_coef=est_cox("`dep'","`varlist'","`exp'","`touse'")
				mata: qte_cov_fitted=coxpred("`quants'","`varlist'","`exp'","`touse'",qte_cov_coef)
				mata: qte_cov_counter=coxpred("`quants'","`counterfactual'","`exp'","`touse'",qte_cov_coef)
			}
		}
		mata: qte_cov_res=qte_cov_counter-qte_cov_fitted
		mata: qte_cov_quant=st_matrix("`quants'")
		return scalar obs=`obs'
		return scalar obs0=`obs0'
		return scalar obsc=`obsc'
		return matrix quants `quants'
		return scalar nreg=`nreg0'
end

*Mata function doing the unconditional estimation using cox
version 9.2
mata numeric matrix coxpred(string scalar quants, string scalar varlist, string scalar weights, string scalar touse, numeric matrix Coef1)
{
	Quants=st_matrix(quants)
	Reg=st_data(.,tokens(varlist),touse)
	Wei=st_data(.,weights,touse)
	beta=Coef1[1..cols(Reg),1]
	t=Coef1[(cols(Reg)+1)..rows(Coef1),1]
	S0=Coef1[(cols(Reg)+1)..rows(Coef1),2]'
	index=exp(Reg*beta)
	Sc=(S0#J(rows(Reg),1,1)):^(index#J(1,cols(S0),1))
	Sc=1:-Sc
	F=mean(Sc,Wei)
	RQ_deco_ReT=mm_quantile(t,F'-(0\F[1..(length(F)-1)]'),Quants)'
	return(RQ_deco_ReT)
}
*Mata function doing the unconditional estimation using locsca
version 9.2
mata numeric matrix lspred(string scalar quants, string scalar varlist, string scalar scale1, string scalar weights, string scalar touse, numeric matrix Coef1)
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

* Clean Distribution Function                                                       
mata real cleandist(real matrix v){
	n = rows(v)
	m=cols(v)
	for( i=1 ; i<=n ; i++ ){
		for( j=1 ; j<=m ; j++ ){
			if( v[i,j] < 0 ){
				v[i,j] = 0
			}
			if( v[i,j] > 1 ){
				v[i,j] = 1
			}
		}
	}
	return(v)
}

mata:
numeric getquantile(numeric colvector y, numeric colvector F, numeric colvector TAU_)
{
	real colvector Q
	NumRows = rows(y)
	NumTau = rows(TAU_)
	for( i=1 ; i<=NumTau ; i++ ){
		Q = Q \ y[min((colsum( F :<= TAU_[i] )+1,NumRows))]
	}
	return(Q)
}
end

*Mata function doing the unconditional estimation using logit, probit or lpm
version 9.2
mata numeric matrix distpred(string scalar quants,string scalar varlist,string scalar weights,string scalar touse,numeric matrix Coef1,method)
{
	Quants=st_matrix(quants)
	Reg=st_data(.,tokens(varlist),touse)
	Wei=st_data(.,weights,touse)
	ys=Coef1[1,.]
	Coef=Coef1[2..rows(Coef1),1..(cols(Coef1)-1)]
	Pred=cross(Reg'\J(1,rows(Reg),1),Coef)
	if(method==1) Pred=invlogit(Pred)
	if(method==2) Pred=normal(Pred)
//	if(method==3) Pred=cleandist(Pred)
	Pred=mean(Pred,Wei)'\1
	Pred=sort(Pred, 1)
	RQ_deco_ReT=getquantile(ys',Pred,Quants)
	return(RQ_deco_ReT')
}

*Mata function doing the unconditional estimation using qr
version 9.2
mata numeric matrix rqpred(string scalar quants, string scalar varlist, string scalar weights, string scalar touse, numeric matrix Coef1)
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
mata numeric matrix est_qr(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar nquant1, numeric scalar beta, numeric scalar small, numeric scalar max_it)
{
	y=st_data(.,dep,touse)
	x=st_data(.,tokens(reg),touse)
	n=rows(x)
	x=x,J(n,1,1)
	k=cols(x)
	w=st_data(.,weight,touse)
	nquant=st_numscalar(nquant1)
	coef=J(cols(x)+1,nquant,.)
	coef[1,.]=(0.5/nquant:+(0..(nquant-1)):/nquant)
	for (i=1; i<=nquant; i++) {
		coef[2..(k+1),i]=rq_fnm(x, y, w, coef[1,i], beta, small, max_it, dep, reg, weight, touse)
	}
	return(coef)
}

*interior QR
mata real vector rq_fnm(numeric matrix X, numeric colvector dep, numeric colvector weight, numeric scalar p, numeric scalar beta, numeric scalar small, numeric scalar max_it, string scalar depo, string scalar rego, string scalar weighto, string scalar touse)
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
mata real vector bound(numeric vector x, numeric vector dx)
{
	b = J(1,length(x),1e20)
	f=select(1..length(x),dx:<0)
	b[f] = -x[f] :/ dx[f]
	return(b)
}

*Mata function doing the conditional estimation using probit
version 9.2
mata numeric matrix est_logit(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	if(nreg<.){
		depeval=uniqrows(mm_quantile(dep,wei,(1..nreg):/nreg:-0.5/nreg)')
	}
	else{
		depeval=uniqrows(dep)
	}
	nreg=rows(depeval)
	st_numscalar(nreg1,nreg)
	coef=J(cols(reg)+1,0,.)
	idx = st_addvar("double", st_tempname())
	for (i=1; i<=nreg; i++) {
		level=depeval[i]
		st_store(.,idx,touse,dep:<=level)
		stata("logit "+st_varname(idx)+" "+reg1+" [iweight="+wei1+"] if "+touse+"==1, asis",1)
		coef=coef,(level\st_matrix("e(b)")')
	}
	coef=coef[.,select(1..nreg,colmissing(coef):==0)]
	coef=coef,(max(dep)\J(cols(reg),1,.))
	return(coef)
}

*Mata function doing the conditional estimation using probit
version 9.2
mata numeric matrix est_probit(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	if(nreg<.){
		depeval=uniqrows(mm_quantile(dep,wei,(1..nreg):/nreg:-0.5/nreg)')
	}
	else{
		depeval=uniqrows(dep)
	}
	if(max(dep)==max(depeval)){		
		nreg=rows(depeval)-1	
	} 
	else nreg=rows(depeval)
	st_numscalar(nreg1,nreg)
	coef=J(cols(reg)+1,0,.)
	idx = st_addvar("double", st_tempname())
	for (i=1; i<=nreg; i++) {
		level=depeval[i]
		st_store(.,idx,touse,dep:<=level)
		stata("probit "+st_varname(idx)+" "+reg1+" [iweight="+wei1+"] if "+touse+"==1, asis",1)
		coef=coef,(level\st_matrix("e(b)")')
	}
	coef=coef,(max(dep)\J(cols(reg),1,.))
	return(coef)
}

*Mata function doing the conditional estimation using linear probability model
version 9.2
mata numeric matrix est_lpm(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	if(nreg<.){
		depeval=uniqrows(mm_quantile(dep,wei,(1..nreg):/nreg:-0.5/nreg)')
	}
	else{
		depeval=uniqrows(dep)
	}
	if(max(dep)==max(depeval)){		
		nreg=rows(depeval)-1	
	} 
	else nreg=rows(depeval)
	st_numscalar(nreg1,nreg)
	coef=J(cols(reg)+1,0,.)
	for (i=1; i<=nreg; i++) {
		level=depeval[i]
		coef=coef,(level\invsym(cross(reg,wei,reg))*cross(reg,wei,dep:<=level))
	}
	coef=coef,(max(dep)\J(cols(reg),1,.))
	return(coef)
}

*Mata function doing the conditional estimation using location model
version 9.2
mata numeric matrix est_loc(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse, string scalar nreg1)
{
	dep=st_data(.,dep1,touse)
	reg=st_data(.,tokens(reg1),touse)
	reg=reg,J(rows(reg),1,1)
	wei=st_data(.,wei1,touse)
	nreg=st_numscalar(nreg1)
	coef=J(cols(reg)+1,nreg,.)
	coef[1,.]=(0.5/nreg:+(0..(nreg-1)):/nreg)
	beta=invsym(cross(reg,wei,reg))*cross(reg,wei,dep)
	resid=dep:-reg*beta
	coef[2..rows(coef),.]=beta#J(1,nreg,1)
	coef[cols(reg)+1,.]=coef[cols(reg)+1,.]+mm_quantile(resid,wei,coef[1,.])
	return(coef)
}

*Mata function doing the conditional estimation using location model
version 9.2
mata numeric matrix est_locsca(string scalar dep1, string scalar reg1, string scalar scale1, string scalar wei1, string scalar touse, string scalar nreg1)
{
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
	return(coef)
}

*Mata function doing the conditional estimation using location model
version 9.2
mata numeric matrix est_cox(string scalar dep1, string scalar reg1, string scalar wei1, string scalar touse)
{
	stata("preserve",1)
	stata("sort "+dep1,1)
	stata("stset "+dep1+" [pweight="+wei1+"]",1)
	S0 = st_tempname()
	stata("stcox "+reg1+" if "+touse+" , basesurv("+S0+")",1)
	coef=st_matrix("e(b)")'
	t=st_data(.,"_t",touse)
	S0=st_data(.,S0,touse)
	stata("stset, clear",1)
	coef=(coef,J(rows(coef),1,.))\(t,S0)
	return(coef)
}

*Mata function doing the conditional estimation using cqr
mata numeric matrix est_cqr(string scalar dep, string scalar censoring, string scalar reg, string scalar weight, string scalar touse, string scalar nquant1, numeric scalar firstc, numeric scalar secondc, numeric scalar nsteps, numeric scalar right, numeric scalar beta, numeric scalar small, numeric scalar max_it)
{
	c=st_data(.,censoring,touse)
	y=st_data(.,dep,touse)
	x=st_data(.,tokens(reg),touse)
	x=x,J(rows(x),1,1)
	w=st_data(.,weight,touse)
	nquant=st_numscalar(nquant1)
	quants=(0.5/nquant:+(0..(nquant-1)):/nquant)'
	if(right==0) coef=est_cqrl(y, c, x, w, quants, firstc, secondc, nsteps, beta, small, max_it, dep, reg, weight, touse)
	if(right==1){
		y=-y
		c=-c	
		coef=-est_cqrl(y, c, x, w, 1:-quants, firstc, secondc, nsteps, beta, small, max_it, dep, reg, weight, touse)
		coef[1,.]=quants'
	}
	return(coef)
}

mata numeric matrix est_cqrl(numeric colvector y, numeric colvector c, numeric matrix x, numeric colvector w, numeric colvector quants, numeric scalar c1, numeric scalar c2, numeric scalar nsteps, numeric scalar beta, numeric scalar small, numeric scalar max_it, string scalar dep, string scalar reg, string scalar weight, string scalar touse)
{
	coef=J(cols(x),rows(quants),.)
	ncensored=(y:>c)
	idx = st_addvar("double", st_tempname())
	st_store(.,idx,touse,ncensored)
	stata("logit "+st_varname(idx)+" "+reg+" [iweight="+weight+"] if "+touse+"==1",1)
	temp=st_matrix("e(b)")'
	pred=invlogit(x*temp)
	nq=rows(quants)
	coef=J(cols(x),nq,.)
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
		coef[.,i]=temp
	}
	coef=quants'\coef
	return(coef)
}
