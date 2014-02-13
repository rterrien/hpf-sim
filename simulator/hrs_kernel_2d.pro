;+
; NAME:
;  hrs_kernel_2d
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
;	out = hrs_kernel_2d(npix_x,npix_y,width)
;
; INPUTS:
;
;	npix_x: Number of x pixels of the output array
;
;	npix_y: Number of y pixels of the output array (must be integer mult of x)
;
;	width: Width of the kernel in pixels (corresponds roughly to the width of the tophat)
;
; OUTPUTS:
;
;	The 2-d kernel
;	
; KEYWORD PARAMETERS: 
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-12-2014
;-


function hrs_kernel_2d, npix_x, npix_y, width

	;implementing as a circle, then rebin to stretch for now.

	npix = npix_x < npix_y
	
	if (npix_x > npix_y) mod npix ne 0 then begin
		print,'dimensions must be integer multiples'
		return,!values.f_nan
	endif

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
	
	
	
	xs = abs(dindgen(npix) - npix/2)
	;ys = abs(dindgen(npix_y) - npix_y/2)
	xs2 = rebin(xs,npix,npix)
	ys2 = rebin(1#xs,npix,npix)
	
	dd = sqrt(xs2^2. + ys2^2.)
	
	;xx = dblarr(npix_x,npix_x)
	
	zs = interpol(t1,x2,dd)
	ngtm = where(dd gt max(x2),nn)
	if nn ge 1 then zs[ngtm] = 0d
	
	zs_final = rebin(zs,npix_x,npix_y)
	
	tot = total(zs_final,/double)
	
	zs_final = zs_final/tot
	
	return, zs_final
	
end
