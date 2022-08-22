cap prog drop par_qrprocess
program define par_qrprocess, byable(recall)
	syntax , command(string) opts(string) fit(namelist) nq(numlist)
	marksample touse
	`command' if `touse', `opts'
	tokenize `fit'
	forvalues q = 1/`nq'{
		tempvar temp
		quietly predict `temp' if `touse', equation(#`q')
		quiet replace ``q'' = `temp' if `touse'
	}
 end
 