function rleg, x, c
	;return legendre polynomials (m=0)
	;x is indep variable
	;c is coeffs
	nx = n_elements(x)
	nc = n_elements(c)
	res = dblarr(nx)
	for i=0, nc-1 do begin
		res += c[i] * legendre(x,i)
	endfor
	return,res
end
