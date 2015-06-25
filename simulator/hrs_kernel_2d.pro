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
;		03-13-14: RCT modified to stretch inner 6.5 pixels, maintain edge shapes
;-


function hrs_kernel_2d, npix_x, npix_y, width

	;implementing as a circle, then rebin to stretch if needed

	npix = npix_x < npix_y
	
	if (npix_x > npix_y) mod npix ne 0 then begin
		print,'dimensions must be integer multiples'
		return,!values.f_nan
	endif

	; Input kernel as read off from an average of beams on an HRS image
	;t1a = [35454.459d,35367.913,34272.269,28004.848,8753.4973,743.31283,146.71539,25.814310,7.9165832,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000]
	
	t1a = [35454.459d,35367.913,34272.269,28004.848,8753.4973,743.31283,146.71539,25.814310,7.9165832,1.0000000,0.1000000,0.0100000,0.0010000,0.0001000]
	
	t1a = t1a/max(t1a)
	
	x1a = dindgen(n_elements(t1a))
	
	; Assume fiber width is 6.5 pixels (info from chad)
	p0 = 6.5
	scale = width / p0
	
	; Upsample to get to 6.5/2 = 3.25, quarter pixels
	x1 = dindgen(n_elements(t1a)*4)/4.
	t1 = interpol(t1a,x1a,x1)
	
	x2 = x1
	
	; Stretch inner 3.25 pixels to make top hat wider
	x2[0:13] = x2[0:13]*scale
	xoffset = x2[14]*scale - x2[14]
	x2[14:*] = x2[14:*] + xoffset
	
	; Generate distance (from center) array
	xs = abs(dindgen(npix) - npix/2)

	xs2 = rebin(xs,npix,npix)
	ys2 = rebin(1#xs,npix,npix)
	
	dd = sqrt(xs2^2. + ys2^2.)
	
	; Calculate value for each distance
	
	zs = interpol(t1,x2,dd)
	ngtm = where(dd gt max(x2),nn)
	if nn ge 1 then zs[ngtm] = 0d
	
	; stretch if needed
	zs_final = rebin(zs,npix_x,npix_y)
	
	; normalize
	tot = total(zs_final,/double)
	
	zs_final = zs_final/tot

	return, zs_final
	
end
