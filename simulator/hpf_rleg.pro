;+
; NAME:
;  hpf_rleg
;
; PURPOSE:
;
;	A convenience function to implement a variable order legendre polynomial
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	out = hpf_rleg(x,c)
;
; INPUTS:
;
;	x: The independent variable (must be -1 < x < 1)
;
;	c: The coefficients
;
; OUTPUTS:
;	
;	The legendre function result
;
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_rleg, x, c
	;return legendre polynomials (m=0)
	;x is indep variable
	;c is coeffs
	nx = n_elements(x)
	nc = n_elements(c)
	res = dblarr(nx)
	for i=0, nc-1 do begin
		res += c[i] * legendre(x,i,/double)
	endfor
	return,res
end
