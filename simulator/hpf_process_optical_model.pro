;+
; NAME:
;  hpf_process_optical_model
;
; PURPOSE:
;
;  Process the general parameters from the optical model into an array of x/y/wavelength for the pixels and pixel edges
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_process_optical_model(optical_params)
;
; INPUTS:
;
;	optical_params: The optical_params structure defined in initialize_optical_params
;
;	
; OUTPUTS:
;	
;	The structure containing x/y/wvl according to the echellogram.
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-

function hpf_process_optical_model, optical_params

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	pixel_size = 18d-3 ;mm
	
	restore, optical_params.model_file
	
	xs = model.xs
	ys = model.ys
	ws = model.ws

	norders = n_elements(model.orders)
	
	nfibers = optical_params.nfibers
	buffer = 100.
	upfactor = optical_params.projection_upsample
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;MODIFY ORDER SHAPES IF NECESSARY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if keyword_set(optical_params.orders_shape ne 0) then begin
		hpf_mod_echel, xs, ys, ws, xs1, ys1, ws1, ver=optical_params.orders_shape
		xs_old = xs
		ys_old = ys
		ws_old = ws
		xs = xs1
		ys = ys1
		ws = ws1
	endif	

	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PROCESS ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	fs_base = dblarr((2048 + buffer)*upfactor,norders,nfibers)
	
	xs_base = ((dindgen((2048+buffer)*upfactor) - (1024.+buffer/2d)*upfactor)/upfactor)
	ws_base = dblarr((2048 + buffer)*upfactor,norders)
	
	fs_base_noup = dblarr((2048 + buffer),norders,nfibers)
	xs_base_noup = ((dindgen((2048+buffer)) - (1024.+buffer/2d)))
	ws_base_noup = dblarr((2048 + buffer),norders)

	
	;find left and right limits of each pixel for binning the spectrum
	xs_base_left = xs_base - 0.5
	xs_base_right = xs_base + 0.5
	ws_base_left = dindgen((2048 + buffer) * upfactor,norders)
	ws_base_right = dindgen((2048 + buffer) * upfactor,norders)
	
	ys_base = dblarr((2048 + buffer) * upfactor,norders)
	xs_base_2d = dblarr((2048 + buffer) * upfactor,norders)
	
	xs_base_left_noup = xs_base_noup - 0.5
	xs_base_right_noup = xs_base_noup + 0.5
	ws_base_left_noup = dindgen((2048 + buffer),norders)
	ws_base_right_noup = dindgen((2048 + buffer),norders)
	ys_base_noup = dblarr((2048 + buffer),norders)
	xs_base_2d_noup = dblarr((2048 + buffer),norders)


	;fill in the wavelengths that correspond to each bin and the edge of each bin
	;note that wavelength depends only on x-position
	for i=0, norders-1 do begin
		;for each order, there are 2048 x pixels
		;derive the wl for each pixel based on the x/wl map alone
		w1a = mpfitfun('poly',xs_pix[*,i],ws[*,i],1d,[0d,0d,0d],yfit=w1b,/quiet)
		w1c = poly(xs_base,w1a)
		ws_base[*,i] = w1c
		dw1 = (ws_base[1:*,i] - ws_base[*,i])/2d
		dex = interpol(dw1,xs_base[0:-2],xs_base[-1])
		dw2 = [dw1,dex]
		ws_base_right[*,i] = ws_base[*,i] + dw2
		ws_base_left[*,i] = ws_base[*,i] - dw2
		
		y1a = mpfitfun('poly',xs_pix[*,i],ys_pix[*,i],1d,[0d,0d,0d],yfit=y1b,/quiet)
		y1c = poly(xs_base,y1a)
		ys_base[*,i] = y1c

		;replicate the x base array for convolving and plotting
		xs_base_2d[*,i] = xs_base
		
		w1c_noup = poly(xs_base_noup,w1a)
		ws_base_noup[*,i] = w1c_noup
		dw1_noup = (ws_base_noup[1:*,i] - ws_base_noup[*,i])/2d
		dex_noup = interpol(dw1_noup,xs_base_noup[0:-2],xs_base_noup[-1])
		dw2_noup = [dw1_noup,dex_noup]
		ws_base_right_noup[*,i] = ws_base_noup[*,i] + dw2
		ws_base_left_noup[*,i] = ws_base_noup[*,i] - dw2
		
		y1c_noup = poly(xs_base_noup,y1a)
		ys_base_noup[*,i] = y1c_noup

		;replicate the x base array for convolving and plotting
		xs_base_2d_noup[*,i] = xs_base_noup

	endfor
	
	out = {xs:xs, ys:ys, ws:ws, xs_pix:xs_pix, ys_pix:ys_pix, fs_base:fs_base, xs_base:xs_base, ws_base:ws_base, xs_base_left:xs_base_left, xs_base_right:xs_base_right, ws_base_left:ws_base_left, ws_base_right:ws_base_right, ys_base:ys_base, xs_base_2d:xs_base_2d, fs_base_noup:fs_base_noup, xs_base_noup:xs_base_noup, ws_base_noup:ws_base_noup, xs_base_left_noup:xs_base_left_noup, xs_base_right_noup:xs_base_right_noup, ws_base_left_noup:ws_base_left_noup, ws_base_right_noup:ws_base_right_noup, ys_base_noup:ys_base_noup, xs_base_2d_noup:xs_base_2d_noup}
	
	
	return,out
	
end
