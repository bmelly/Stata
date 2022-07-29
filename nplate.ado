*! version 2.1.2  3mar2010  Markus Froelich and Blaise Melly
*1.0.1: add proportion of compliers
*2.0.1: correct estimation of the variance
*2.1.0: more information provided by the command
*2.1.1: higher-order kernels, check that the treatment is binary, check that the dependent variable is between 0 and 1 for the logit
*2.1.2: check that the treatment variable is between 0 and 1 if treat_lin is not activated.

program nplate, eclass
	version 9.2
	capt findfile lmoremata.mlib
	if _rc {
      	di as error "-moremata- is required; type {stata ssc install moremata}"
		exit
	}
	capt findfile kdens.mata
	if _rc {
      	di as error "-kdens- is required; type {stata ssc install kdens}"
		exit
	}
	syntax anything(name=0) [if] [in] [, Continuous(varlist) Dummy(varlist) Unordered(varlist) Kernel(string) Bandwidth(real -1) Lambda(real 1) treat_lin dep_log Variance mata_opt LEvel(cilevel) trim(real 0.01) inst_lin] 
	if "`continuous'"!=""{
		unab continuous : `continuous'
	}
	if "`unordered'"!=""{
		unab unordered : `unordered'
	}
	if "`dummy'"!=""{
		unab dummy : `dummy'
	}
	gettoken dependent 0 : 0, parse("(")
	gettoken weg 0 : 0, parse("(")
	gettoken 0 : 0, parse(")")
	gettoken treatment 0 : 0, parse("=")
	gettoken weg 0 : 0, parse("=")
	local instrument="`0'"
	marksample touse
	markout `touse' `dependent' `treatment' `instrument' `continuous' `dummy' `unordered'
	if "`treatment'"==""{
		dis in red "A treatment variable must be provided."
		exit
	}
	if "`instrument'"==""{
		dis in red "An instrumental variable must be provided."
		exit
	}
	if "`mata_opt'"=="mata_opt" {
		capt findfile optimize.mata
     		if _rc{
			di as error "The option mata_opt can only be used with Stata 10 or newer."
			exit
		}
	}
*check that there are one enodogen and one instrument
	tokenize `treatment'
	local treatment="`1'"
	if "`2'"!=""{
		di in red "Only one treatment variable may be specified"
		exit
	}
	tokenize `instrument'
	local instrument="`1'"
	if "`2'"!=""{
		di in red "Only one instrumental variable may be specified"
		exit
	}
*check that the instrument is binary 0 or 1 
	quietly tab `instrument'  if `touse'
	if r(r)!=2{
		di in red "The instrument, `instrument', is not binary."
		di in red "You can apply Theorem 8 of Froelich (2007) by defining a new instrument that is 0 for the lowest value, 1 for the highest value and missing otherwise."
		di in red "Use nplate again with this new instrument."
		exit
	}
	quietly sum `instrument'  if `touse'
	if r(min)!=0 | r(max)!=1{
		di in red "The instrument, `instrument', is not a 0/1 variable."
		exit
	}
	if "`treat_lin'"==""{
		quietly sum `treatment'
		if r(max)>1 | r(min)<0{
			dis as error "The range of the treatment variable must be between 0 and 1 when the option treat_lin is not activated."
			exit
		}
	}
	if "`dep_log'"!=""{
		quietly sum `dependent'
		if r(max)>1 | r(min)<0{
			dis as error "The range of the dependent variable must be between 0 and 1 when the option dep_logit is activated."
			exit
		}
	}
	quietly summarize `touse'
	local obs=r(sum)
* check that there is no singularity within the continuous variables
	if "`continuous'"!=""{
		_rmcollright `continuous'
		local continuous "`r(varlist)'"
	}
* check that there is no singularity within the dummy variables
	if "`dummy'"!=""{
		_rmcollright `dummy'
		local dummy "`r(varlist)'"
	}
	if "`unordered'"!=""{	
		foreach x in `unordered'{
			tempvar `x'
			quietly tabulate `x' if `touse', generate(``x'')
			local nc=r(r)
			forvalues i=2/`nc'{	
				local listu1 "`listu1' `x'`i'"
			}
			drop ``x'`index''1
			unab temp:``x'`index''*
			local listu "`listu' `temp'"
		}
	}
	if ("`dummy'"!="" & "`continuous'"!="") | ("`dummy'"!="" & "`unordered'"!="") | ("`continuous'"!="" & "`unordered'"!=""){
		quietly _rmcollright `dummy' `continuous' `listu'
		if "`r(dropped)'"!=""{
			di in red "The covariates are multicollinear."
			exit
		}
	}
	if `lambda'<0 | `lambda'>1{
		di as error "The options lambda must be between 0 and 1"
		exit
	}
*if kernel is missing, set kernel to epan2
	if "`kernel'"==""{
		local kernel "epan2"
	} 
*Calculate the regression functions
	tempname m1 m0 mu1 mu0
	local continuous1 "`continuous'"
	if "`continuous'"==""{
		local continuous "empty"
	}
	local dummy1 "`dummy'"
	if "`dummy'"==""{
		local dummy "empty"
	}
	local unordered1 "`unordered'"
	local listu1 "`listu'"
	if "`unordered'"==""{
		local unordered "empty"
		local listu "empty"
	}
	if (("`dummy'"=="empty" & "`unordered'"=="empty") | `lambda'==1) & ("`continuous'"=="empty" | `bandwidth'<=0){
		if "`dep_log'"=="dep_log"{
			quietly logit `dependent' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==0 
			quietly predict `m0' if `touse'
			quietly logit `dependent' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==1
			quietly predict `m1' if `touse'
		}
		else{
			quietly regress `dependent' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==0 
			quietly predict `m0' if `touse'
			quietly regress `dependent' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==1
			quietly predict `m1' if `touse'
		}
		if "`treat_lin'"=="treat_lin"{
			quietly regress `treatment' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==0 
			quietly predict `mu0' if `touse'
			quietly regress `treatment' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==1
			quietly predict `mu1' if `touse'
		}
		else{
			quietly logit `treatment' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==0 
			quietly predict `mu0' if `touse'
			quietly logit `treatment' `dummy1' `continuous1' `listu1' if `touse' & `instrument'==1
			quietly predict `mu1' if `touse'
		}
		if "`variance'"=="variance"{
			tempname temp pi vy1 cvy1d1 vd1 vy0 cvy0d0 vd0
			if "`inst_lin'"=="inst_lin"{
				quietly regress `instrument' `dummy1' `continuous1' `listu1' if `touse'
				quietly predict `pi' if `touse'
			}
			else{
				quietly logit `instrument' `dummy1' `continuous1' `listu1' if `touse'
				quietly predict `pi' if `touse'
			}
			if "`dep_log'"=="dep_log"{
				quietly generate `vy1'=`m1'*(1-`m1') if `touse'
				quietly generate `vy0'=`m0'*(1-`m0') if `touse'
			}
			else{
				quietly generate `temp'=(`dependent'-`m1')^2 if `touse' & `instrument'==1
				quietly regress `temp' `dummy1' `continuous1' `listu1'
				quietly predict `vy1' if `touse'
				quietly drop `temp'	
				quietly generate `temp'=(`dependent'-`m0')^2 if `touse' & `instrument'==0
				quietly regress `temp' `dummy1' `continuous1' `listu1'
				quietly predict `vy0' if `touse'
				quietly drop `temp'	
			}
			quietly generate `temp'=(`dependent'-`m1')*(`treatment'-`mu1') if `touse' & `instrument'==1
			quietly regress `temp' `dummy1' `continuous1' `listu1'
			quietly predict `cvy1d1' if `touse'
			quietly drop `temp'
			quietly generate `temp'=(`dependent'-`m0')*(`treatment'-`mu0') if `touse' & `instrument'==0
			quietly regress `temp' `dummy1' `continuous1' `listu1'
			quietly predict `cvy0d0' if `touse'
			quietly generate `vd1'=`mu1'*(1-`mu1') if `touse'
			quietly generate `vd0'=`mu0'*(1-`mu0') if `touse'
		}
	}
	else{
		tempname temp
		quietly generate `m0'=.
		quietly generate `m1'=.
		quietly generate `mu0'=.
		quietly generate `mu1'=.
		if "`dep_log'"=="dep_log"{
			quietly gen `temp'=`touse'*(1-`instrument')
			if "`mata_opt'"=="mata_opt"{
				mata: loclog1("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m0'")
			}
			else {
				mata: loclog("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m0'")
			}
			quietly replace `temp'=`touse'*`instrument'
			if "`mata_opt'"=="mata_opt"{
				mata: loclog1("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m1'")
			}
			else {
				mata: loclog("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m1'")
			}
		}
		else{
			quietly gen `temp'=`touse'*(1-`instrument')
			mata: loclin("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m0'")
			quietly replace `temp'=`touse'*`instrument'
			mata: loclin("`dependent'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`m1'")
		}
		if "`treat_lin'"==""{
			quietly replace `temp'=`touse'*(1-`instrument')
			if "`mata_opt'"=="mata_opt"{
				mata: loclog1("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu0'")
			}
			else {
				mata: loclog("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu0'")
			}
			quietly replace `temp'=`touse'*`instrument'
			if "`mata_opt'"=="mata_opt"{
				mata: loclog1("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu1'")
			}
			else {
				mata: loclog("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu1'")
			}
		}
		else{
			quietly replace `temp'=`touse'*(1-`instrument')
			mata: loclin("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu0'")
			quietly replace `temp'=`touse'*`instrument'
			mata: loclin("`treatment'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`mu1'")
		}
		if "`variance'"=="variance"{
			tempname temp1 pi vy1 cvy1d1 vd1 vy0 cvy0d0 vd0
			quiet gen `pi'=.
			quiet gen `vy1'=.
			quiet gen `vy0'=.
			quiet gen `vd1'=.
			quiet gen `vd0'=.
			quiet gen `cvy1d1'=.
			quiet gen `cvy0d0'=.
			if "`inst_lin'"=="inst_lin"{
				mata: loclin("`instrument'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`touse'","`touse'","`pi'")
			}
			else{
				if "`mata_opt'"=="mata_opt"{
					mata: loclog1("`instrument'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`touse'","`touse'","`pi'")
				}
				else{
					mata: loclog("`instrument'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`touse'","`touse'","`pi'")
				}
			}
			if "`dep_log'"=="dep_log"{
				quietly replace `vy1'=`m1'*(1-`m1') if `touse'
				quietly replace `vy0'=`m0'*(1-`m0') if `touse'
			}
			else{
				quietly generate `temp1'=(`dependent'-`m1')^2 if `touse' & `instrument'==1
				quietly replace `temp'=`touse'*`instrument'
				mata: loclin("`temp1'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`vy1'")
				quietly drop `temp1'	
				quietly generate `temp1'=(`dependent'-`m0')^2 if `touse' & `instrument'==0
				quietly replace `temp'=`touse'*(1-`instrument')
				mata: loclin("`temp1'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`vy0'")
				quietly drop `temp1'	
			}
			quietly generate `temp1'=(`dependent'-`m1')*(`treatment'-`mu1') if `touse' & `instrument'==1
			quietly replace `temp'=`touse'*`instrument'
			mata: loclin("`temp1'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`cvy1d1'")
			quietly drop `temp1'	
			quietly generate `temp1'=(`dependent'-`m0')*(`treatment'-`mu0') if `touse' & `instrument'==0
			quietly replace `temp'=`touse'*(1-`instrument')
			mata: loclin("`temp1'","`continuous'","`dummy'","`unordered'","`listu'","`kernel'",`bandwidth',`lambda',"`temp'","`touse'","`cvy0d0'")
			quietly drop `temp1'	
			quietly replace `vd1'=`mu1'*(1-`mu1') if `touse'
			quietly replace `vd0'=`mu0'*(1-`mu0') if `touse'
		}
	}
	tempname a b c d e f g h i j k l res var late
	quietly sum `m1' if `touse'
	sca `a'=r(mean)
	quietly sum `m0' if `touse'
	sca `b'=r(mean)
	quietly sum `mu1' if `touse'
	sca `c'=r(mean)
	quietly sum `mu0' if `touse'
	sca `d'=r(mean)
	sca `res'=(`a'-`b')/(`c'-`d')
*estimate the variance
	if "`variance'"=="variance"{
		drop `temp'
		quietly gen `temp'=`vy1'/`pi' if `touse' & `pi'>`trim'
		quietly sum `temp' if `touse'
		sca `e'=r(mean)
		drop `temp'
		quietly gen `temp'=`cvy1d1'/`pi' if `touse' & `pi'>`trim'
		quietly sum `temp' if `touse'
		sca `f'=r(mean)
		drop `temp'
		quietly gen `temp'=`vd1'/`pi' if `touse' & `pi'>`trim'
		quietly sum `temp' if `touse'
		sca `g'=r(mean)
		drop `temp'
		quietly gen `temp'=`vy0'/(1-`pi') if `touse' & `pi'<(1-`trim')
		quietly sum `temp' if `touse'
		sca `h'=r(mean)
		drop `temp'
		quietly gen `temp'=`cvy0d0'/(1-`pi') if `touse' & `pi'<(1-`trim')
		quietly sum `temp' if `touse'
		sca `i'=r(mean)
		drop `temp'
		quietly gen `temp'=`vd0'/(1-`pi') if `touse' & `pi'<(1-`trim')
		quietly sum `temp' if `touse'
		sca `j'=r(mean)
		quietly replace `temp'=(`m1'-`m0'-`res'*`mu1'-`res'*`mu0')^2 if `touse'
		quietly sum `temp' if `touse'
		sca `k'=r(mean)
		sca `l'=(`c'-`d')^2
		mat `var'=(1/`l')*(`e'-2*`res'*`f'+`g'+`h'-2*`res'*`i'+`j'+`k')/`obs'
		mat rownames `var'=late
		mat colnames `var'=late
		mat `late'=`res'
		mat rownames `late'=late
		mat colnames `late'=late
		ereturn post `late' `var', esample(`touse')
	} 
	else{
		mat `late'=`res'
		mat rownames `late'=late
		mat colnames `late'=late
		ereturn post `late', esample(`touse')
	}
*display results
	dis
	dis in green "Nonparametric IV estimation of local average treatment effects"
	dis "Estimator suggested in Froelich (2007)"
	dis
	dis "Dependent variable:" _column(30) "`dependent'"
	dis "Treatment variable:" _column(30) "`treatment'"
	dis "Control variable(s):" _column(30) "`continuous1' `dummy1' `unordered1'"
	dis "Number of observations:" _column(30)  "`obs'"
	quietly sum `mu1'
	local a=r(mean)
	quietly sum `mu0'
	local complier=round(`a'-r(mean),0.001)
	dis "Proportion of compliers:" _column(30) "`complier'"
	dis
	if `bandwidth'<=0{
		local bandwidth "infinity"
	}
	if "`dep_log'"=="dep_log"{
		dis "E[Y|X,Z] is estimated by local logit with h=`bandwidth' and lambda=`lambda'."
	}
	else{
		dis "E[Y|X,Z] is estimated by local linear regression with h=`bandwidth' and lambda=`lambda'."
	}
	if "`treat_lin'"=="treat_lin"{
		dis "E[D|X,Z] is estimated by local linear regression with h=`bandwidth' and lambda=`lambda'."
	}
	else{
		dis "E[D|X,Z] is estimated by local logit regression with h=`bandwidth' and lambda=`lambda'."
	}
	if "`variance'"=="variance"{
		dis 
		dis "Estimation of the standard errors:"
		if "`inst_lin'"=="inst_lin"{
			dis "E[Z|X] is estimated by local linear regression with h=`bandwidth' and lambda=`lambda'."
		}
		else{
			dis "E[Z|X] is estimated by local logit regression with h=`bandwidth' and lambda=`lambda'."
		}
	}
	dis
	ereturn display, level(`level')
end

*Higher order kernels
mata real colvector fm_kern(string scalar name, real colvector u)
{
	if(name=="epanechnikov_o3" | name=="epanechnikov_o4"){
		w=(3/4):*(15/8:-7/8*5:*u:^2):*(1:-u:^2):*(u:^2:<1)
	} else if(name=="epanechnikov_o5" | name=="epanechnikov_o6"){ 
		w=(3/4):*(175/64:-105/32*5:*u:^2+231/320*25:*u:^4):*(1:-u:^2):*(u:^2:<1)
	} else if(name=="gaussian_o3" | name=="gaussian_o4"){
		w=(1/2):*(3:-u:^2):*normalden(u)
	} else if(name=="gaussian_o5" | name=="gaussian_o6"){
		w=(1/8):*(15:-10:*u:^2+u:^4):*normalden(u)
	} else if(name=="gaussian_o7" | name=="gaussian_o8"){
		w=(1/48):*(105:-105:*u:^2:+21:*u:^4-u:^6):*normalden(u)
	} else{
		w=mm_kern(name,u)
	}
	return(w)
}

*Mata function calculating mixed kernel and returning the regressors with a constant in the first column
version 9.2
mata void mkernel(regd, regc, regu, ev, evu, string scalar kernel, real scalar band, real scalar lambda,real scalar n,real scalar nd, real scalar nc,real scalar nu, w, reg)
{
	if(nd+nu>0) regdt=regd:-ev[1..(nd+nu)] 
	else regdt=regd
	if(nc>0) regct=(regc:-ev[(nd+nu+1)..(nd+nu+nc)])
	else regct=regc
	if(nu>0) regut=regu:-evu
	else regut=regu
	w=J(n,1,1)
	if(nc>0 & band>0) for(i=1;i<=nc;i++) w=w:*fm_kern(kernel,regct[.,i]:/band)
	if((nd+nu)>0 & lambda<1) w=w:*(lambda:^((nd+nu):-rowsum(regdt:==0)))
	if(nd>0) reg=select((J(n,1,1),regdt[.,1..nd],regut,regct),w)
	else reg=select((J(n,1,1),regut,regct),w)
} 

*Mata function estimating the propensity score by local linear regression
version 9.2
mata void loclin(string scalar dep, string scalar continuous, string scalar dummy, string scalar unordered, string scalar unord_list, string scalar kernel, real scalar bandwidth, real scalar lambda, string scalar touse, string scalar touse1, string scalar out)
{
	y=st_data(.,dep,touse)
	n_use=rows(y)
	if(continuous~="empty") {
		xc=st_data(.,tokens(continuous),touse)
		evc=st_data(.,tokens(continuous),touse1)
		xct=xc*luinv(cholesky(variance(xc)))'
		evct=evc*luinv(cholesky(variance(xc)))'
		n=rows(evc)
	}  
	if(dummy~="empty"){
		xd=st_data(.,tokens(dummy),touse) 
		evd=st_data(.,tokens(dummy),touse1) 
		n=rows(evd)
	}
	if(unordered~="empty"){
		xur=st_data(.,tokens(unord_list),touse)
		xuk=st_data(.,tokens(unordered),touse)	
		evur=st_data(.,tokens(unord_list),touse1)
		evuk=st_data(.,tokens(unordered),touse1)	
		n=rows(evuk)
	}
	if(continuous=="empty"){
		xct=J(n_use,0,0)
		evct=J(n,0,0)
	}
	if(dummy=="empty"){
		xd=J(n_use,0,0)
		evd=J(n,0,0)
	}
	if(unordered=="empty"){
		xur=xuk=J(n_use,0,0)
		evur=evuk=J(n,0,0)
	}
	nc=cols(xc)
	nd=cols(xd)
	nu=cols(xuk)
	pred=J(n,1,.)
	for(i=1; i<=n; i++){
		h1=bandwidth
		l1=lambda
		mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		dettemp=det(reg'reg)
		while(dettemp<1e-7 & h1<100){
			h1=h1*1.05
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
			dettemp=det(reg'reg)
		}
		if(dettemp<1e-7){
			h1=1e10
			l1=1
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		}
		yt=select(y,w)
		wt=select(w,w)
		pred[i,1]=(invsym((reg:*wt)'reg)*(reg:*wt)'yt)[1,1]
	}
	st_store(.,out,touse1,pred)
}

*Mata, logistic distribution
version 9.2
mata real colvector logisticcdf(x) return(1:/(1:+exp(-x)))

*Mata, logistic density
version 9.2
mata real colvector logisticpdf(x)
{
temp=exp(-x)
return(temp:/(1:+temp)^2)
}

*Mata, objective function of the weighted logit estimator
version 9.2
mata void lnwlogit(todo, p,  x, y, w,  lnf, S, H)
{
	prob=logisticcdf(x*p')
	lnf=w:*log(y:*prob:+(1:-y):*(1:-prob))
	if (todo >= 1) {
		S=w:*(x:*(y-prob))
		if (todo==2) {
			H=-(w:*x:*prob:*(1:-prob))'x
		}
	}
}

*Mata function estimating the propensity score by local logit, optimization using Stata 10 optimizer
version 9.2
mata void loclog(string scalar dep, string scalar continuous, string scalar dummy, string scalar unordered, string scalar unord_list, string scalar kernel, real scalar bandwidth, real scalar lambda, string scalar touse, string scalar touse1, string scalar out)
{
	y=st_data(.,dep,touse)
	n_use=rows(y)
	if(continuous~="empty") {
		xc=st_data(.,tokens(continuous),touse)
		evc=st_data(.,tokens(continuous),touse1)
		xct=xc*luinv(cholesky(variance(xc)))'
		evct=evc*luinv(cholesky(variance(xc)))'
		n=rows(evc)
	}  
	if(dummy~="empty"){
		xd=st_data(.,tokens(dummy),touse) 
		evd=st_data(.,tokens(dummy),touse1) 
		n=rows(evd)
	}
	if(unordered~="empty"){
		xur=st_data(.,tokens(unord_list),touse)
		xuk=st_data(.,tokens(unordered),touse)	
		evur=st_data(.,tokens(unord_list),touse1)
		evuk=st_data(.,tokens(unordered),touse1)	
		n=rows(evuk)
	}
	if(continuous=="empty"){
		xct=J(n_use,0,0)
		evct=J(n,0,0)
	}
	if(dummy=="empty"){
		xd=J(n_use,0,0)
		evd=J(n,0,0)
	}
	if(unordered=="empty"){
		xur=xuk=J(n_use,0,0)
		evur=evuk=J(n,0,0)
	}
	nc=cols(xc)
	nd=cols(xd)
	nu=cols(xuk)
	pred=J(n,1,0)
	S = optimize_init()
	optimize_init_evaluator(S, &lnwlogit())
	optimize_init_evaluatortype(S, "v2")
	optimize_init_conv_maxiter(S, 300)
	optimize_init_verbose(S, 0)
	optimize_init_tracelevel(S, "none")
	for(i=1; i<=n; i++){
		h1=bandwidth
		l1=lambda
		mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		dettemp=det(reg'reg)
		while(dettemp<1e-7 & h1<100){
			h1=h1*1.05
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
			dettemp=det(reg'reg)
		}
		if(dettemp<1e-7){
			h1=1e10
			l1=1
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		}
		ret=1
		while(ret!=0){
			yt=select(y,w)
			wt=select(w,w)
			optimize_init_params(S,((invsym((reg:*wt)'reg)*(reg:*wt)'yt))')
			optimize_init_argument(S, 1, reg)
			optimize_init_argument(S, 2, yt)
			optimize_init_argument(S, 3, wt)
			ret = _optimize(S)
			if(ret!=0){
				if(h1<100){
					h1=h1*1.05
				}
				else{
					h1=1e10
					l1=1
				}
				mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
			}
		}
		pred[i,1]=logisticcdf(optimize_result_params(S)[1])
	}
	st_store(.,out,touse1,pred)
}

version 9.2
mata real rowvector intlog(numeric colvector dep, numeric matrix reg ,numeric colvector we, convergence)
{
	objo=0
	b=((invsym((reg:*we)'reg)*(reg:*we)'dep))'
	prob=logisticcdf(reg*b')
	objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
	db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	it=1
	while(it<100 & sum(abs(db):>1e-8)>0 & abs(objn-objo)>1e-8){
		objo=objn
		b=b+db
		it=it+1
		prob=logisticcdf(reg*b')
		objn=colsum(we:*log(dep:*prob:+(1:-dep):*(1:-prob)))
		db=colsum(we:*(reg:*(dep-prob)))*invsym((we:*reg:*prob:*(1:-prob))'reg)
	}
	convergence=(it==100)
	return(b)
}

*Mata function estimating the propensity score by local logit, optimization using self written codes
version 9.2
mata void loclog1(string scalar dep, string scalar continuous, string scalar dummy, string scalar unordered, string scalar unord_list, string scalar kernel, real scalar bandwidth, real scalar lambda, string scalar touse, string scalar touse1, string scalar out)
{
	y=st_data(.,dep,touse)
	n_use=rows(y)
	if(continuous~="empty") {
		xc=st_data(.,tokens(continuous),touse)
		evc=st_data(.,tokens(continuous),touse1)
		xct=xc*luinv(cholesky(variance(xc)))'
		evct=evc*luinv(cholesky(variance(xc)))'
		n=rows(evc)
	}  
	if(dummy~="empty"){
		xd=st_data(.,tokens(dummy),touse) 
		evd=st_data(.,tokens(dummy),touse1) 
		n=rows(evd)
	}
	if(unordered~="empty"){
		xur=st_data(.,tokens(unord_list),touse)
		xuk=st_data(.,tokens(unordered),touse)	
		evur=st_data(.,tokens(unord_list),touse1)
		evuk=st_data(.,tokens(unordered),touse1)	
		n=rows(evuk)
	}
	if(continuous=="empty"){
		xct=J(n_use,0,0)
		evct=J(n,0,0)
	}
	if(dummy=="empty"){
		xd=J(n_use,0,0)
		evd=J(n,0,0)
	}
	if(unordered=="empty"){
		xur=xuk=J(n_use,0,0)
		evur=evuk=J(n,0,0)
	}
	nc=cols(xc)
	nd=cols(xd)
	nu=cols(xuk)
	pred=J(n,1,0)
	for(i=1; i<=n; i++){
		h1=bandwidth
		l1=lambda
		mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		dettemp=det(reg'reg)
		while(dettemp<1e-7 & h1<100){
			h1=h1*1.05
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
			dettemp=det(reg'reg)
		}
		if(dettemp<1e-7){
			h1=1e10
			l1=1
			mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
		}
		convergence=.
		while(convergence!=0){
			convergence=.
			yt=select(y,w)
			wt=select(w,w)
			bt = intlog(yt,reg,wt,convergence)
			if(convergence!=0){
				if(h1<100){
					h1=h1*1.05
				}
				else{
					h1=1e10
					l1=1
				}
				mkernel((xd,xuk),xct,xur,(evd[i,.],evuk[i,.],evct[i,.]),evur[i,.],kernel,h1,l1,n_use,nd,nc,nu,w=.,reg=.)
			}
		}
		pred[i,1]=logisticcdf(bt[1])
	}
	st_store(.,out,touse1,pred)
}
