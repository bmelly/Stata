cap prog drop par_qrprocess
program define par_qrprocess, byable(recall)
	syntax anything [if] [in] [pweight], command(string) opts(string) fit(namelist) nq(integer) group_var(varname) ext_touse(varname)
	marksample touse
	count if !missing(`group_var') & `touse'
	if r(N)>0{
		capture `command' `anything' [`weight'`exp'] if `touse', `opts'
		if _rc==0{
			tokenize `fit'
			forvalues q = 1/`nq'{
				tempvar temp
				predict `temp' if `touse', equation(#`q')
				replace ``q'' = `temp' if `touse'
			}
		}
		else{
			replace `ext_touse' = 0  if `touse'
		}
	}
 end
 