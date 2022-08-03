*plotdeco: plot quantile decompositions
*! version 0.0.1  03.08.2022  Blaise Melly

program plotdeco
*version control
	version 9.2
	*the last estimates must exist en allow for plotdeco
	if "`e(plotdeco)'" == ""{
		if "`e(cmd)'" == "" {
			error 301
		}
		else{ 
			dis in red "plotdeco cannot be used after `e(cmd)'."
			error 332
		}
	}
	syntax [anything(name=namelist)] [, Pointwise Uniform Both None lcolor(string) ucolor(string) pcolor(string) legend(string) ytitle(string) xtitle(string) title(string) OTHER_graph_options(string) COMBINE_options(string) Level(string) compare(string) compare_ci compare_color(string)]
	*put the quantiles in a temp variable
	tempname quantile quantiles
	mat `quantile' = e(quantiles)
	local nq=rowsof(`quantile')
	quiet svmat `quantile', name(`quantiles')
	*default plot
	if "`namelist'"=="" {
		if wordcount("`e(plotdeco)'") > 1{
			local namelist "`e(plotdeco)' all"
		}
		else local namelist "`e(plotdeco)'"
	}
	*clean the names of the plot to match the name of the matrices
	tokenize "`namelist'"
	*number of temporary variables for the fitted values
	local n_estimates = 1
	*number of plots
	local kplot=wordcount("`namelist'")
	forvalues i = 1/`kplot'{
		if "``i''" != "all"{
			if "``i''" == "tot" | "``i''" == "total" {
				local i "total_difference"
			}
			if "``i''" == "char" {
				local i "characteristics" 
			}
			if "``i''" == "coef" {
				local i "coefficients" 
			}
			if "``i''" == "resid" {
				local i "residuals" 
			}
			confirm matrix e(``i'')
			local plotlist "`plotlist' ``i''"
		}
		else {
			local n_estimates = wordcount("`e(plotdeco)'")
			local plotlist "`plotlist' all"
		}
	}
	local namelist "`plotlist'"
	*generate the temporary variables for the points estimates
	forvalues i = 1/`n_estimates'{
		tempvar estimates`i'
		quiet gen `estimates`i''=.
	}
	*by defaut, the confidence intervals are plotted if they have been estimated
	if "`pointwise'"=="" & "`uniform'"=="" & "`both'"=="" & "`none'"==""{
		if "`e(vce)'"!="novar"{
			local pointwise "pointwise"
			local uniform "uniform"
		}
	}
	if	"`both'"=="both"{
		local pointwise "pointwise"
		local uniform "uniform"
	}
	*the level can be changed only for the pointwise CI
	if "`level'"!="" & "`uniform'"=="uniform"{
		dis in red "The option level cannot be specified for uniform confidence bands."
		dis in red "The option level must be specified when `e(cmd)' is called."
		error 400
	}
	*Are the CI available?
	if "`pointwise'"=="pointwise" | "`uniform'" == "uniform"{
		if "`e(vce)'"=="novar"{
			dis in red "The confidence intervals cannot be plotted because the s.e. have not been estimated by `e(cmd)'."
			dis in red "Call `e(cmd)' to obtain the confidence intervals (do not deactivate the bootstrap)."
			local pointwise ""
			local uniform ""
		}
	}
	*defaut colors
	if "`lcolor'"==""{
		local lcolor "black"
	}
	if "`pcolor'"==""{
		local pcolor "gs8"
	}	
	if "`ucolor'"==""{
		local ucolor "gs13"
	}
	*generate temporary variables for the pointwise CI
	if "`pointwise'" == "pointwise"{
		tempvar pointl pointu
		quie gen `pointl'=.
		quie gen `pointu'=.
		if "`level'"==""{
			local level = c(level)
		}
	}
	*generate the temporary variables for the uniform CI
	if "`uniform'"=="uniform"{
		tempvar unifl unifu
		quie gen `unifl'=.
		quie gen `unifu'=.
	}
	*defaut values
	if "`xtitle'" == ""{
		local xtitle "Quantile"
	}
	if "`xtitle'"=="off"{
		local xtitle ""
	}
	if "`ytitle'"=="off"{
		local ytitle "Quantile effect"
	}
	*if legend is off: no legend
	local legend_all `"`legend'"'
	*defaut for legend (when not off):
	*no legend when several panels for the panel with a single component (title is enough)
	if `kplot'>1 & `"`legend'"'==""{
		local legend "off"
	}
	*extract the titles
	if `kplot'>1 & wordcount("`title'")>=`kplot'{
		local different_titles="diff"
		forvalues i=1/`kplot'{
			gettoken title`i' title : title , parse("||")
			gettoken trash title:title , parse("||")
			gettoken trash title:title , parse("||")
		}
	}
	*default legend for the single components
	if `"`legend'"'=="" & ("`pointwise'"=="pointwise" | "`uniform'"=="uniform"){
		if "`pointwise'"=="pointwise" & "`uniform'"=="uniform"{
			local legend "order(3 "Estimates" 2 "Pointwise CI" 1 "Uniform CB") rows(1)"
		}
		else if "`pointwise'"=="pointwise"{
			local legend "order(2 "Estimates" 1 "Pointwise CI") rows(1)"
		}
		else{
			local legend "order(2 "Estimates" 1 "Uniform CB") rows(1)"
		}
	}
/*	if "`compare'"!=""{
		tempname keep_est
		quietly estimates store `keep_est'
		gettoken compare_fct compare_opt : compare, parse(",")
		if "`e(wtype)'"!=""{
			local compare_weight "`e(wtype)'=`e(wexp)'"
		}
		if _caller() >= 11{
			version 11.1: quietly `compare_fct' `e(depvar)' `e(xvar)' [`compare_weight'] `compare_opt'
		}
		else{
			version 9.2: quietly `compare_fct' `compare_fct' `e(xvar)' [`compare_weight'] `compare_opt'
		}
		if "`compare_color'"==""{
			local compare_color "red*0.6"
		}
	}*/
	*loop over the panels
	forvalues i=1/`kplot'{
		local temp_var = word("`namelist'",`i')
		*titles
		if "`different_titles'"=="diff"{
			local temptitle `title`i''
		}
		else if "`title'"==""{
			if "`temp_var'" == "characteristics" local temptitle "Characteristics"
			if "`temp_var'" == "coefficients" local temptitle "Coefficients"
			if "`temp_var'" == "residuals" local temptitle "Residuals"
			if "`temp_var'" == "total_difference" local temptitle "Observed differences"
			if "`temp_var'" == "all" local temptitle "Quantile decomposition"
		}
		else if "`title'"=="off"{
			local temptitle ""
		}
		else local temptitle "`title'"
		*plot for a single component
		if "`temp_var'" != "all"{
			*point estimate
			mata: plotdeco_Temp_Mat = st_matrix("e(`temp_var')")
			mata: st_store((1..`nq')', "`estimates1'", plotdeco_Temp_Mat[.,1])
			local estgraph "(line `estimates1' `quantiles', lcolor(`lcolor'))"
			*pointwise CI
			if "`pointwise'"=="pointwise"{
				mata: st_store((1..`nq')', ("`pointl'", "`pointu'"), (plotdeco_Temp_Mat[., 1] + invnormal((100-`level')/200) * plotdeco_Temp_Mat[., 2], plotdeco_Temp_Mat[., 1] - invnormal((100-`level')/200) * plotdeco_Temp_Mat[., 2]))
				local pointgraph "(rarea `pointl' `pointu' `quantiles', fcolor(`pcolor') lcolor(`pcolor'))"
			}
			*uniform CI
			if "`uniform'"=="uniform"{
				mata: st_store((1..`nq')', ("`unifl'", "`unifu'"), plotdeco_Temp_Mat[., (3, 4)])
				local unifgraph "(rarea `unifl' `unifu' `quantiles', fcolor(`ucolor') lcolor(`ucolor'))"
			}
			*legend
			local templegend `"`legend'"'
		}
		*plot for all components
		else {
			tokenize "`e(plotdeco)'"
			local estgraph ""
			*loop over the different lines
			forvalues j = 1/`n_estimates'{
				mata: st_store((1..`nq')', "`estimates`j''", st_matrix("e(``j'')")[.,1])
				local estgraph "`estgraph' (line `estimates`j'' `quantiles')"
				if "``j''" == "total_difference" local templegend_all `"`templegend_all' `j' "tot." "' 
				else if "``j''" == "coefficients" local templegend_all `"`templegend_all' `j' "coef." "'
				else if "``j''" == "characteristics" local templegend_all `"`templegend_all' `j' "char." "'
				else if "``j''" == "residuals" local templegend_all `"`templegend_all' `j' "resid." "'
				else local templegend_all `"`templegend_all' `j' "``j''" "'
			}
			if `"`legend_all'"' == ""{
				if `kplot' == 1 {
					local templegend `"order(`templegend_all') rows(1)"'
				}
				else {
					local templegend `"order(`templegend_all') rows(1) symxsize(*0.3) colgap(*0.5) keygap(*0.5)"'
				}
			}
			else {
				local templegend `"`legend_all'"'
			}
			local pointgraph ""
			local unifgraph ""
		}
/*		if "`compare'"!=""{
			local ols_coef = _b[`temp_var']
			local olsgraph "(scatteri `ols_coef' 0 `ols_coef' 1, recast(line) lcolor(`compare_color'))"
		}	
		if "`compare_ci'"!=""{
			if "`level'"==""{
				local level=c(level)
			}
			local ols_lb = _b[`temp_var']+invnormal((100-`level')/200)*_se[`temp_var']
			local ols_ub = _b[`temp_var']+invnormal(1-(100-`level')/200)*_se[`temp_var']
			local olscigraph "(scatteri `ols_lb' 0 `ols_lb' 1, recast(line) lcolor(`compare_color') lpattern(dash)) (scatteri `ols_ub' 0 `ols_ub' 1, recast(line) lcolor(`compare_color') lpattern(dash))"
		}
*/		
		if `kplot'==1{
			twoway `unifgraph' `pointgraph' `estgraph' `olsgraph' `olscigraph', ytitle("`ytitle'") xtitle("`xtitle'") title("`temptitle'") legend(`templegend') `other_graph_options'
		}
		else{
			tempname graph`i'
			local graphcombine "`graphcombine' `graph`i''"
			twoway `unifgraph' `pointgraph' `estgraph' `olsgraph' `olscigraph', ytitle("`ytitle'") xtitle("`xtitle'") title("`temptitle'") legend(`templegend') `other_graph_options' nodraw name(`graph`i'')
		}
	}
	if `kplot'>1{
		graph combine `graphcombine', `combine_options'
	}
/*	if "`compare'"!=""{
		quietly estimates restore `keep_est'
	}*/
	cap mata: mata drop plotdeco_Temp_Mat
end	
