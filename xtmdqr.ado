*xtmdqr: panel data minimum distance quantile regression
*! version 0.0.1  03.08.2022  Blaise Melly

*Ideas: (1) an option to add a linear trend? This is trivial but I am
*sure that many users will try to add indicator variables for periods...

program xtmdqr, eclass byable(recall) sortpreserve
	version 9.2
*if the command is used without arguments it shows the previous results
	if replay() {
		if "`e(cmd)'"!="xtmdqr" { 
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
		syntax varlist [if] [in] [pweight] [, BE FE RE Quantiles(numlist >0 <1 sort) Cluster(varname) qr_opts(string) BOOTstrap(string) Level(cilevel) noPrint Save_first(name) Load_first(name) n_small(integer 1)]
		if ("`be'"!="") + ("`fe'"!="") + ("`re'"!="") > 1 { 
			di in red "choose only one of be, fe, or re"
			error 198 
		}
		if ("`fe'" == "") & ("`be'" == "") {
			local re "re"
		}
*separate dependent variable from regressors
		gettoken dep reg : varlist
		if _caller() >= 11{
			version 11.1: fvexpand `reg'
			if "`r(fvops)'" == "true"{
				local fvops = 1
				local reg `r(varlist)'
			}
		}
		else{ 
			confirm numeric variable `reg'
		}
		quietly xtset
		local group `r(panelvar)'
		local i=1
		quietly foreach v of local reg {
			tempname mean_`i'
			bysort `group': egen `mean_`i'' = mean(`v')
			if "`be'"==""{
				tempname demean_`i'
				generate `demean_`i'' = `v' - `mean_`i''
			}
			if "`be'" == "be" | "`re'" == "re" {
				local inst "`inst' `mean_`i''" 
			}
			if "`fe'" == "fe" | "`re'" == "re" {
				local inst "`inst' `demean_`i''"
			}
			local i = `i' + 1 
		}
		mdqr `dep' (`reg' = `inst') `if' `in' [`weight'`exp'], group(`group') quantiles(`quantiles') cluster(`cluster') qr_opts(`qr_opts') bootstrap(`bootstrap') level(`level') `print' save_first(`save_first') load_first(`load_first') n_small(`n_small')
		ereturn local cmd "xtmdqr"
	}
end
