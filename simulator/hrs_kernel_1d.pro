;+
; NAME:
;  hrs_kernel_1d
;
; PURPOSE:
;
;	Make a kernel based on a mean extracted image of the orders on an HRS exposure
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	out = hrs_kernel_1d(npix,width)
;
; INPUTS:
;
;	npix: Number of pixels of the output array
;
;	width: Width of the kernel in pixels (corresponds roughly to the width of the tophat)
;
; OUTPUTS:
;
;	The 1-d kernel
;	
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-12-2014
;-


function hrs_kernel_1d, npix, width
	t1 = [35454.459d,35367.913,34272.269,28004.848,8753.4973,743.31283,146.71539,25.814310,7.9165832,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000]
	t1 = t1/max(t1)
	
	x1 = dindgen(n_elements(t1))
	
	p0 = 6.5
	
	scale = width / 6.5
	
	x2 = x1
	
	x2[1] = scale
	x2[2:*] = (scale - 1d) + x2[2:*]
	
	;try scaling whole thing
	x2 = x1 * scale
	

	if npix mod 2 eq 0 then begin
		print,'NEED ODD # PIXELS'
		return,!values.f_nan
	endif
	
	xs = abs(dindgen(npix) - npix/2)
	
	ys = interpol(t1,x2,xs)
	ngtm = where(xs gt max(xs),nn)
	if nn ge 1 then ys[ngtm] = 0d
	
	tot = total(ys,/double)
	
	ys /= tot
	
	return, ys
	
end
