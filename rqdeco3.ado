*! version 1.1 27jan2010 Blaise Melly
*  version 1.0: 26jan2010 
*  version 1.1: 27jan2010, names of the columns of the matrices

program rqdeco3, rclass byable(recall)
	version 9.2
*1 replay
	if replay(){
		if "`r(cmd)'" != "rqdeco3" {
			error 301
		}
		if _by() {
			error 190 
		}
		syntax [, Level(cilevel) noPRint]
		tempname results se quants
		local nquantreg=r(nquantreg)
		local obs=r(obs)
		local obs0=r(obs0)
		local obs1=r(obs1)
		matrix `results'=r(results)
		matrix `se'=r(se)
		matrix `quants'=r(quants)
		local nq=rowsof(`quants')
	}
*2 estimation
	if !replay(){
		syntax varlist(min=2) [if] [in] [pw] , by(varname) [ NQuantreg(integer 100) QLow(real 0.1) QHigh(real 0.9) QStep(real 0.1) Quantiles(string) vce(string) Level(cilevel) reps(integer 50) noPRint SAving(string)]
*2.1 Consistency checks
		local vce=substr("`vce'",1,1)
*2.1.a if no method for variance has been selected: no veriance estimation
		if "`vce'"==""{
			local vce="n"
		}
*2.1.b if bootstrap: check that at least 2 replications have been chosen
		if "`vce'"=="b"{
			if (`reps')<2 {
				di in red `"The number of replication must be higher or equal to 2"'
				error
			}
		}
*2.1.c Check that either bootstrap or no variance estimation has been selected
		else if "`vce'"~="n"{
			di in red `"The option vce can take only the values "bootstrap" or "none""'
			exit
		}
*2.1.d Prepare the quantiles at which we estimate the decomposition
		tempname quants nquantreg1 results se
*if a sequence has been selected: divide by hundred if higher than 1, check that the selected values make sense
		if "`quantiles'"==""{
			if `qlow'>=1{
				local qlow=`qlow'/100
				local qhigh=`qhigh'/100
				local qstep=`qstep'/100
			}
			if `qlow'>=1 | `qhigh'>=1 | `qstep'>=1 | `qlow'<=0 | `qhigh'<=0 | `qstep'<=0 | `qstep'>(`qhigh'-`qlow') | (((`qhigh'-`qlow')/`qstep')-round((`qhigh'-`qlow')/`qstep'))~=0{
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
			SetQ `quantiles'
			local quantiles="`r(quants)'"
			tokenize "`quantiles'", parse(" ")
			local i=1
			while "`1'" != "" {
				matrix `quants'=nullmat(`quants')\(`1')
				mac shift 
				local i=`i'+1
			}
		}
		sca `nquantreg1'=`nquantreg'
*2.1.e error message if one or more than 2 groups
		quietly summarize `by' `if' `in'
		if r(min)==r(max) {
			di as error "`by' is a constant"
			exit 198
		}
		capture assert `by'==0 | `by'==1 | `by'==.
		if _rc != 0 {
			di as error "`by' must take on two values, 0 or 1"
			exit 198
		}
*2.1.f error message if level not between 0 and 100
		if `level' <= 0 | `level' >= 100 {
			di as err "level() must be between 0 and 100"
			exit 198
		}
*2.2 Prepare the estimation
		marksample touse
		quietly: sum `by' if `touse'
		local obs=r(N)
		quietly: sum if `touse' & `by'==0
		local obs0=r(N)
		quietly: sum if `touse' & `by'==1
		local obs1=r(N)
		if "`vce'"=="b"{
			di in gr "Fitting base model"
		}
*2.3 call the estimation command
		rqdeco_c `dep' `varlist' if `touse' [`weight'`exp'], by(`by') nq(`nquantreg') ql(`qlow') qh(`qhigh') qs(`qstep') q(`quantiles') obs(`obs')
		mat `results'=r(results)
		mat colnames `results'=quantile fitted_1 fitted_x1m1r0 fitted_x1m0r0 fitted_0 total_differential residuals median characteristics
		local nq=rowsof(`quants')
*2.4 Variance estimation
*2.4.a: no variance estimation: input missing
		if "`vce'"=="n"{
			local k 1
			while `k'<=`nq'{
				mat `se'=nullmat(`se')\(.,.,.)
				local k=`k'+1
			}
		}
*2.4.b: bootstrap
		if "`vce'"=="b"{
	*preserve and keep only useful observations and variables
			preserve
			quietly keep if `touse'
			tokenize `exp'
			local exp1 "`2'" /* exp1 is the name of the weighting variable */
			quietly keep `dep' `varlist' `by' `exp1'
	*prepare bootstraping
			if `"`saving'"'=="" {
				tempfile BOOTRES
				local filetmp "yes"
			}
			else {
				_prefix_saving `saving'
				local BOOTRES	`"`s(filename)'"'
				if "`double'" == "" {
					local double	`"`s(double)'"'
				}
				local every	`"`s(every)'"'
				local replace	`"`s(replace)'"'
			}
			tempvar bwt
			qui gen double `bwt' = .
			tempfile data
			tempname res res1 res2 res3 res4
			di in gr "(bootstrapping " _c
			local j 1
			quietly save "`data'"
			local actual_more=c(more)
			set more off
	*loop for each replication
			while `j'<=`reps'{
				quietly use "`data'", clear
				bsampl_w `bwt'
				if "`weight'"~=""{
					replace `bwt'=`bwt'*`exp1'
				}
				capture{
					rqdeco_c `dep' `varlist' [pw=`bwt'], by(`by') nq(`nquantreg') ql(`qlow') qh(`qhigh') qs(`qstep') q(`quantiles')  obs(`obs')
				}
				local rc = _rc
				if (`rc'==0) {
					quietly drop _all
					mat `res'=r(results)
					mat `res1'=vec(`res'[1...,6])'
					quietly svmat `res1', names("X1_M1_R1_")
					mat `res2'=vec(`res'[1...,7])'
					quietly svmat `res2', names("X1_M1_R0_")
					mat `res3'=vec(`res'[1...,8])'
					quietly svmat `res3', names("X1_M0_R0_")
					mat `res4'=vec(`res'[1...,9])'
					quietly svmat `res4', names("X0_M0_R0_")
					if( `j'==1){
						quietly save "`BOOTRES'", `replace'
					}
					else{
						quietly append using "`BOOTRES'"
						quietly save "`BOOTRES'", replace
					}
					di in gr "." _c
					local j=`j'+1
				}
				else {
					if _rc == 1 {
						 exit 1 
					}
					di in red "*" _c
				}
			}
	* evaluate the boostrap results	
			set more `actual_more'
			di in gr ")"
			di
			qui use "`BOOTRES'", clear
			tempname se
			mata: rqboot("`quants'","`se'")
			mat colnames `se'=se_total_differential se_residuals se_median se_characteristics
			use "`data'"
		}
	}
*3. Display the results
*3.1 header
	display as text _n "Decomposition of differences in distribution using quantile regression" 
	dis as text _column(10) "Total number of observations" _c
	dis as result _column(54) %8.0f `obs'
	dis as text _column(12) "Number of observations in group 0" _c
	dis as result _column(54) %8.0f `obs0'
	dis as text _column(12) "Number of observations in group 1" _c
	dis as result _column(54) %8.0f `obs1'
	dis
	dis as text _column(10) "Number of quantile regressions estimated" _c
	di as result _column(54) %8.0f `nquantreg'
	di
	if "`vce'"=="b"{
		di as text "The variance has been estimated by bootstraping the results " as result `reps'  as text " times"
	}
	else if "`vce'"=="n"{
		di as text "The variance has not been computed." _new "Use the option vce if you want to compute it."
	}
	dis
	if "`print'"!="noprint"{
		dis as text "{hline 17}" "{c TT}" "{hline 59}"
		dis as text _column(0) %17s "Component" _column(18) "{c |}" _column(21) %~10s "Effects" _column(31) %~10s "Std. Err." _column(42) %~8s "t" _column(50) %~8s "P>|t|" _column(58) %~21s "[`level'% Conf. Interval]"
		dis as text "{hline 17}" "{c +}" "{hline 59}"
		local k 1
*3.2 Results for each quantile 
		while `k'<=`nq'{
			dis as text "Quantile " as result %-8.7g `quants'[`k',1] as text _column(18) "{c |}"
			dis as text _column(0) %17s "Raw difference" as text _column(18) "{c |}" _column(21) as result %8.0g (`results'[`k',6]) _column (31) as result %8.0g (`se'[`k',1]) _column (42) as result %3.2f (`results'[`k',6]/`se'[`k',1]) _column(51) as result %4.3f chi2tail(1,(`results'[`k',6]/`se'[`k',1])^2) _column(59) as result %8.0g (`results'[`k',6]-invnormal(0.5*(1+`level'/100))*`se'[`k',1]) _column(70) as result %8.0g (`results'[`k',6]+invnormal(0.5*(1+`level'/100))*`se'[`k',1])
			dis as text _column(0) %17s "Residuals" as text _column(18) "{c |}" _column(21) as result %8.0g (`results'[`k',7]) _column (31) as result %8.0g (`se'[`k',2]) _column (42) as result %3.2f (`results'[`k',7]/`se'[`k',2]) _column(51) as result %4.3f chi2tail(1,(`results'[`k',7]/`se'[`k',3])^2) _column(59) as result %8.0g (`results'[`k',7]-invnormal(0.5*(1+`level'/100))*`se'[`k',2]) _column(70) as result %8.0g (`results'[`k',7]+invnormal(0.5*(1+`level'/100))*`se'[`k',2])
			dis as text _column(0) %17s "Median" as text _column(18) "{c |}" _column(21) as result %8.0g (`results'[`k',8]) _column (31) as result %8.0g (`se'[`k',3]) _column (42) as result %3.2f (`results'[`k',8]/`se'[`k',3]) _column(51) as result %4.3f chi2tail(1,(`results'[`k',8]/`se'[`k',3])^2)  _column(59) as result %8.0g (`results'[`k',8]-invnormal(0.5*(1+`level'/100))*`se'[`k',3]) _column(70) as result %8.0g (`results'[`k',8]+invnormal(0.5*(1+`level'/100))*`se'[`k',3])
			dis as text _column(0) %17s "Characteristics" as text _column(18) "{c |}" _column(21) as result %8.0g (`results'[`k',9]) _column (31) as result %8.0g (`se'[`k',4]) _column (42) as result %3.2f (`results'[`k',9]/`se'[`k',4]) _column(51) as result %4.3f chi2tail(1,(`results'[`k',9]/`se'[`k',4])^2)  _column(59) as result %8.0g (`results'[`k',9]-invnormal(0.5*(1+`level'/100))*`se'[`k',4]) _column(70) as result %8.0g (`results'[`k',9]+invnormal(0.5*(1+`level'/100))*`se'[`k',4])
			if `k'<`nq' {
				dis as text "{hline 17}" "{c +}" "{hline 59}" 
			} 
			else {
				dis as text "{hline 17}" "{c BT}" "{hline 59}" 
			}
			local k=`k'+1
		}
	}
*3.3 Return the results
	return matrix results=`results'
	return matrix se=`se'
	return local cmd "rqdeco3"
	return scalar obs=`obs'
	return scalar obs0=`obs0'
	return scalar obs1=`obs1'
	return scalar nquantreg=`nquantreg'
	return matrix quants=`quants'
end

*Programs taken from sqreg to check that the inputed quantiles make sense
program define SetQ, rclass
	local orig "`*'"
	tokenize "`*'", parse(" ,")
	while "`1'" != "" {
		FixNumb "`orig'" `1'
		ret local quants "`return(quants)' `r(q)'"
		mac shift 
		if "`1'"=="," {
			mac shift
		}
	}
end
program define FixNumb , rclass
	local orig "`1'"
	mac shift
	capture confirm number `1'
	if _rc {
		Invalid "`orig'" "`1' not a number"
	}
	if `1' >= 1 {
		ret local q = `1'/100
	}
	else 	ret local q `1'
	if `return(q)'<=0 | `return(q)'>=1 {
		Invalid "`orig'" "`return(q)' out of range"
	}
end
program define Invalid
	di in red "quantiles(`1') invalid"
	if "`2'" != "" {
		di in red "`2'"
	}
	exit 198
end

*Estimation of the decomposition
program rqdeco_c, rclass
	syntax varlist [if] [pw], by(varname) nq(integer) obs(integer) [QLow(real 0.1) QHigh(real 0.9) QStep(real 0.1) Quantiles(string)]
	marksample touse
	tempname coefr coefo obs1 nq1 q1 quants isw results
	quietly sum `by' if `touse'
	sca `obs1'=r(N)
	sca `nq1'=`nq'
*put the quantiles in a matrix
	if "`quantiles'"==""{
		if `qlow'>=1{
			local qlow=`qlow'/100
			local qhigh=`qhigh'/100
			local qstep=`qstep'/100
		}
		local i=0
		while `i'<=((`qhigh'-`qlow')/`qstep'){
			local temp=(`qlow'+`i'*`qstep')
			matrix `quants'=nullmat(`quants')\(`temp')
			local i=`i'+1
		}
	}
	else{
		SetQ `quantiles'
		local quantiles="`r(quants)'"
		tokenize "`quantiles'", parse(" ")
		local i=1
		while "`1'" != "" {
			matrix `quants'=nullmat(`quants')\(`1')
			mac shift 
			local i=`i'+1
		}
	}
	sca `isw'=("`weight'"!="")
	if `isw'==1{
		local aw "aw"
*if weights: generate interaction between by and touse (needed by mata function rqdest)
		tempvar s0 s1
		gen `s0'=`touse'*(1-`by')
		gen `s1'=`touse'*`by'
*if weights: exp1 contains the name of the weight
		tokenize `exp'
		local exp1 "`2'"	
	}
*estimate the quantile regressions
	forvalues i=1/`nq'{
		local estq=1/(2*`nq')+(`i'-1)/`nq'
		quietly: _qreg `varlist' if `by'==0 & `touse' [`aw'`exp'], quantile(`estq') nolog
		matrix `coefr'=nullmat(`coefr'),e(b)'
		quietly: _qreg `varlist' if `by'==1 & `touse' [`aw'`exp'], quantile(`estq') nolog
		matrix `coefo'=nullmat(`coefo'),e(b)'
	}
	quietly: _qreg `varlist' if `by'==0 & `touse' [`aw'`exp'], quantile(0.5) nolog
	matrix `coefr'=nullmat(`coefr'),e(b)'
	quietly: _qreg `varlist' if `by'==1 & `touse' [`aw'`exp'], quantile(0.5) nolog
	matrix `coefo'=nullmat(`coefo'),e(b)'
*separate dependent variables from regressors
	gettoken dep varlist : varlist
*call mata
	mata: rqdest("`obs1'","`quants'","`nq1'","`varlist'","`by'","`isw'","`exp1'","`touse'","`coefr'","`coefo'","`results'","`s0'","`s1'")
*return the matrix of results
	return matrix results=`results'
end

*Mata function doing the estimation
version 9.2
mata void rqdest(string scalar obs,string scalar quants,string scalar nquantreg1,string scalar varlist, string scalar by,string scalar isw, string scalar exp1,string scalar touse,string scalar coefr, string scalar coefo,string scalar results,string scalar s0,string scalar s1)
{
//read the data into Mata
	RQ_deco_ObS=st_numscalar(obs)
	RQ_deco_QaN=st_matrix(quants)
	RQ_deco_NqA=st_numscalar(nquantreg1)
	RQ_deco_ReG=st_data(.,tokens(varlist),touse)
	RQ_deco_By=st_data(.,by,touse)
	if(st_numscalar(isw)==0){
		RQ_deco_WeR=RQ_deco_WeO=1
	}
	else{
		RQ_deco_We=st_data(.,exp1,s0)
		RQ_deco_WeR=RQ_deco_We[.,J(1,RQ_deco_NqA,1)]
		RQ_deco_We=st_data(.,exp1,s1)
		RQ_deco_WeO=RQ_deco_We[.,J(1,RQ_deco_NqA,1)]
	}
	RQ_deco_Cr=st_matrix(coefr)
	mediancr=RQ_deco_Cr[.,cols(RQ_deco_Cr)]
	RQ_deco_Cr=RQ_deco_Cr[.,1..(cols(RQ_deco_Cr)-1)]
	RQ_deco_Co=st_matrix(coefo)
	medianco=RQ_deco_Co[.,cols(RQ_deco_Co)]
	RQ_deco_Co=RQ_deco_Co[.,1..(cols(RQ_deco_Co)-1)]
	RQ_deco_Cmorr=medianco:+RQ_deco_Cr:-mediancr
//generate the fitted values
	RQ_deco_Po=cross(RQ_deco_ReG'\J(1,RQ_deco_ObS,1),RQ_deco_Co)
	RQ_deco_Pr=cross(RQ_deco_ReG'\J(1,RQ_deco_ObS,1),RQ_deco_Cr)
	RQ_deco_Pmorr=cross(RQ_deco_ReG'\J(1,RQ_deco_ObS,1),RQ_deco_Cmorr)
//calculate the quantiles of interest
	RQ_deco_ReT=RQ_deco_QaN
//first using coefficients from group with by==1 and X from group with by==1
	RQ_deco_ReT=RQ_deco_ReT,bm_quantile(vec(select(RQ_deco_Po,RQ_deco_By)),vec(RQ_deco_WeO),RQ_deco_QaN)
//second using median coefficients from group with by==1, residual coefficients from by==0, and X from by==1
	RQ_deco_ReT=RQ_deco_ReT,bm_quantile(vec(select(RQ_deco_Pmorr,RQ_deco_By)),vec(RQ_deco_WeO),RQ_deco_QaN)
//third using coefficients from by==0, and X from by==1
	RQ_deco_ReT=RQ_deco_ReT,bm_quantile(vec(select(RQ_deco_Pr,RQ_deco_By)),vec(RQ_deco_WeO),RQ_deco_QaN)
//fourth using coefficients from by==0, and X from by==0
	RQ_deco_ReT=RQ_deco_ReT,bm_quantile(vec(select(RQ_deco_Pr,RQ_deco_By:==0)),vec(RQ_deco_WeR),RQ_deco_QaN)
//calculate the effects of interest
//total difference
	RQ_deco_ReT=RQ_deco_ReT,RQ_deco_ReT[.,2]-RQ_deco_ReT[.,5]
//effects of changes in residuals
	RQ_deco_ReT=RQ_deco_ReT,RQ_deco_ReT[.,2]-RQ_deco_ReT[.,3]
//effects of changes in median coefficients
	RQ_deco_ReT=RQ_deco_ReT,RQ_deco_ReT[.,3]-RQ_deco_ReT[.,4]
//effects of changes in charcteristics
	RQ_deco_ReT=RQ_deco_ReT,RQ_deco_ReT[.,4]-RQ_deco_ReT[.,5]
//return the results
	st_matrix(results,RQ_deco_ReT)
}

*Mata function evaluating the bootstrap
version 9.2
mata void rqboot(string scalar quants,string scalar se)
{
//read the data into Mata
	RQ_deco_QaN=st_matrix(quants)
	RQ_deco_BoOt=st_data(.,.)
//calculate the standard errors
	RQ_deco_Se1=sqrt(diagonal(variance(RQ_deco_BoOt,1)))
	RQ_deco_Se=J(0,4,.)
	nq=rows(RQ_deco_QaN)
	for(i=1;i<=nq;i++) RQ_deco_Se=RQ_deco_Se\(RQ_deco_Se1[i],RQ_deco_Se1[nq+i],RQ_deco_Se1[2*nq+i],RQ_deco_Se1[3*nq+i])
//return the results
	st_matrix(se,RQ_deco_Se)
}

*mata function calculating weighted quantiles
version 9.2
mata real colvector bm_quantile(real colvector x, real colvector w, real colvector p)
{
//read the data into Mata
	if(w==1) w=J(rows(x),1,1)
	o=order(x,1)
	w=w[o]
	x=x[o]
	w=w/sum(w)
	i=1
	j=1
	a=0
	q=J(0,1,.)
	while (j<=rows(p)){
		while(a<p[j]){
			a=a+w[i]
			i++
		}
		q=q\x[i]
		j++
	}
	return(q)
}
