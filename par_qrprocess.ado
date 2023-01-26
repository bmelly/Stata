cap prog drop par_qrprocess
program define par_qrprocess, byable(recall)
	syntax, dep(varname) reg(varlist) weight(varname) touse(varname) group(varname) quantiles(numlist) fitted(namelist) n_small(integer)
	tempname quants
	mata: st_matrix("`quants'", strtoreal(tokens(st_local("quantiles")))')
	mata: parallel_fn("`dep'", "`reg'", "`weight'", "`touse'", "`group'", "`quants'", 0.9995, 0.00001, 100, "`fitted'", `n_small')
end
 
cap mata mata drop parallel_fn()
version 12
mata void parallel_fn(string scalar dep, string scalar reg, string scalar weight, string scalar touse, string scalar group, string scalar quantile, real scalar beta, real scalar small, real scalar max_it, string scalar fitted, real scalar n_small)
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
	unique_g = uniqrows(g)
	fit = J(n, nq, .)
	ngroup = max(g)
	conv = 0
	coef = 0
	for(j = 1; j <= rows(unique_g); j++){
		i = unique_g[j]
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
tokens(fitted)
rows(fit)
cols(fit)	
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
