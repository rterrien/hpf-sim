;; fiber_projection
;; HPF DETECTOR SIMULATOR
;; This routine does the projection of fibers onto the detector
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; HZPFSIM_IMG
;;
;; NOTES
;; 
;;
;; parameters:
;; w - spectrum wavelengths in microns
;; f - spectrum fluxes
;; res - resolution
;; pixel_sampling - number of pixels per resolution element
;; wlimg - image of wavelengths
;; specimg - image of fluxes
;; warray - wavelength array (for extraction)
;; fiber_size - fiber size in pixels
;; fiber_gap - gap between fibers in pixels
;; fiber_fractions - dblarr(7) with fiber power fractions

pro fiber_projection_mod, w, f, res, pixel_sampling, wlimg, specimg, calw=calw, calf=calf, diag_out = diag_out, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, fiber_fractions = fiber_fractions, fiber_scale = fiber_scale, orders_gaps = orders_gaps, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, warray = warray


	pixel_size = 18d-3 ;mm
	
	if n_elements(optical_model) eq 0 then optical_model = 'ramsey_513'
	case optical_model of
		'barnes_1212': optical_model_file = 'model_barnes_1212.sav'
		'ramsey_513': optical_model_file = 'model_ramsey_513.sav'
		else: optical_model_file = 'model_ramsey_513.sav'
	endcase
	
	restore,optical_model_file
	
	
	
	xs = model.xs
	ys = model.ys
	ws = model.ws
	

	mod_echel, xs, ys, ws, xs1, ys1, ws1
	xs_old = xs
	ys_old = ys
	ws_old = ws
	xs = xs1
	ys = ys1
	ws = ws1


	
	norders = n_elements(model.orders)
	
	;make the warray, although this is not *required* for extraction now
	warray=MAKE_ARRAY(2048, norders, /Double, Value=!Values.F_NAN)
	
	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PROCESS ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;_base arrays are interpolated versions of the Barnes input to match what is required
	;for the array
	buffer=100
	upfactor = 1
	
	fs_base = dblarr((2048 + buffer)*upfactor,norders)
	fs_base2 = fs_base
	
	xs_base = ((dindgen((2048+buffer)*upfactor) - (1024.+buffer/2)*upfactor)/upfactor)
	ws_base = dblarr((2048 + buffer)*upfactor,norders)
	
	;find left and right limits of each pixel for binning the spectrum
	xs_base_left = xs_base - 0.5
	xs_base_right = xs_base + 0.5
	ws_base_left = dindgen((2048 + buffer) * upfactor,norders)
	ws_base_right = dindgen((2048 + buffer) * upfactor,norders)
	
	ys_base = dblarr((2048 + buffer) * upfactor,norders)
	xs_base_2d = dblarr((2048 + buffer) * upfactor,norders)

	;fill in the wavelengths that correspond to each bin and the edge of each bin
	;note that wavelength depends only on x-position
	for i=0, norders-1 do begin
		;for each order, there are 2048 x pixels
		;derive the wl for each pixel based on the x/wl map alone
		ws_base[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base,/spline)
		;same for left and right limits of each pixel
		ws_base_left[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_left,/spline)
		ws_base_right[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_right,/spline)
		;find the y for each x based on the x/y map alone
		ys_base[*,i] = interpol(ys_pix[*,i],xs_pix[*,i],xs_base,/spline)
		;replicate the x base array for convolving and plotting
		xs_base_2d[*,i] = xs_base
	endfor

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PROJECTION PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	calflag = calw ne !null
	
	diag_out = ''
	
	;make sure we have all the fiber dimensions
	
	
	if fiber_scale eq !null or fiber_core_um eq !null or fiber_cladding_um eq !null or fiber_buffer_um eq !null then stop
	
	;translate fiber dimensions into pixels
	fiber_core_pix = fiber_core_um * fiber_scale
	fiber_cladding_pix = fiber_cladding_um * fiber_scale
	fiber_buffer_pix = fiber_buffer_um * fiber_scale
	
	fiber_outside_thickness = (fiber_buffer_pix - fiber_core_pix)/2d
	
	if fiber_fractions eq !null then fiber_fractions = replicate(1d/7d,7)
	
	diag_out = diag_out + string(13B) + 'Fiber Scale [pix / um]: ' + string(fiber_scale,format='(G6.3)')
	diag_out = diag_out + string(13B) + 'Fiber Core um: ' + string(fiber_core_um,format='(D5.2)') + ' pixels: ' + string(fiber_core_pix,format='(D5.2)')
	diag_out = diag_out + string(13B) + 'Fiber Cladding um: ' + string(fiber_cladding_um,format='(D5.2)') + ' pixels: ' + string(fiber_cladding_pix,format='(D5.2)')
	diag_out = diag_out + string(13B) + 'Fiber Buffer um: ' + string(fiber_buffer_um,format='(D5.2)') + ' pixels: ' + string(fiber_buffer_pix,format='(D5.2)')
	
	no = n_elements(f[*,0])
	
	subimg_height = round(9.*fiber_buffer_pix + 10.*fiber_buffer_pix + 20)
	
	;subimg_height = round(9.*fiber_buffer_pix + 10.*fiber_buffer_pix + 20) ; remove extra cal sep fiber
	
	if keyword_set(orders_gaps) then begin
		gap = orders_gaps
		lambdalow = orders_lambdalow
		lambdahigh = orders_lambdahigh
		n1 = n_elements(gap)
		n2 = n_elements(lambdalow)
		n3 = n_elements(lambdahigh)
		if n1 ne n2 or n1 ne n3 then stop
	endif else begin
		print,'WARNING: MAKE SURE THE OPTICAL MODEL INPUTS ARE SET'
		stop
	endelse

	nw = n_elements(w)
	
	array_pix = 2048L

	if wlimg eq !null then wlimg = dblarr(array_pix,array_pix)
	if specimg eq !null then specimg = dblarr(array_pix,array_pix)
	
	gap = ys_pix[0,1:*] - ys_pix[0,*]
	
	;y0pos = ROUND((array_pix - ROUND(TOTAL(gap)))/2.) ;;origin for star fiber
	y0pos = round(subimg_height) + 20
	cgap=ROUND(TOTAL(gap, /Cumulative)) ;;cumulative gap for positioning other fibers relative to origin
	
	norders=N_ELEMENTS(lambdalow)
	warray=MAKE_ARRAY(array_pix, norders, /Double, Value=!Values.F_NAN)
	

	gap = [0,reform(gap)]
	if calflag then begin
		calflsr = dblarr(no,array_pix) ;resampled cal fiber
		calwlsr = dblarr(no,array_pix) ;cal wavelengths
		cal_li0 = value_locate(calw,lambdalow) ;same as above but for cal fiber
		cal_li1 = value_locate(calw,lambdahigh)
	endif
	li0=value_locate(w,lambdalow) ;what array index does lambdalow correspond to?
	li1=value_locate(w, lambdahigh) ;what array index does lambdahigh correspond to?

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;FILL IN ORDERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for order=0, no-1 do begin
		print,'order',order
		;select out wavelength and flux arrays with appropriate wavelength limits
		subfls = (reform(f[order,*]))[li0[order]:li1[order]] ;fluxes for the order
		subwls = reform(w[li0[order]:li1[order]]) ;wls for the order
		
	
		;rebin corresponding to the resolution and pixel sampling defined above
		rebin_ryan, pixel_sampling, res, subwls, subfls, wr, fr, min(subwls,/nan), max(subwls,/nan)
		;if diag then printf,diaglun,'n pixels for order ',n_elements(wr)
		diag_out = diag_out + string(13B) + 'n pixels for order ' + string(order,format='(I2)') + ' '+ string(n_elements(wr),format='(I6)')
		
		;take 5 pixels off edges for safety
		wr = wr[4:-6]
		fr = fr[4:-6]
		
		;then just take first 2048 pixels to put into detector
		wsub = wr[0:2047]
		fsub = fr[0:2047]
		
		;figure out total height of sub image
		subimg_height = round(9.*fiber_buffer_pix + 10.*fiber_buffer_pix + 20)
		subimg = dblarr(2048,subimg_height)
		subwlimg = dblarr(2048,subimg_height)
		
		
		
		;upsample along xd direction for blurring
		highsample_factor = 18d
		subimg_highsample = dblarr(2048,subimg_height*highsample_factor)
		
		;define origins for each fiber, including the calibration fibers on top and bottom
		fiber_0s = round(indgen(7)*(2*fiber_buffer_pix) + 3*fiber_buffer_pix + fiber_outside_thickness + 10)
		cal_0s = round([fiber_outside_thickness+10+fiber_buffer_pix,max(fiber_0s) + 2*fiber_buffer_pix])

		;also translate these into the corresponding positions for the highsample array
		fiber_0s_highsample = round(highsample_factor*(indgen(7)*(2*fiber_buffer_pix) + 3*fiber_buffer_pix + fiber_outside_thickness + 10d))
		cal_0s_highsample = round(highsample_factor*([fiber_outside_thickness+10+fiber_buffer_pix,max(fiber_0s) + 2*fiber_buffer_pix]))
		
		
		
		;loop through and fill in flux into fibers (and wavelength values
		for i=0, 6 do begin
			fline = fsub * fiber_fractions[i] / round(fiber_core_pix)
			flines = dblarr(2048,round(fiber_core_pix))
			wlines = dblarr(2048,round(fiber_core_pix))
			for j=0, round(fiber_core_pix)-1 do begin 
				flines[0,j] = fline
				wlines[0,j] = wsub
			endfor
			subimg[0,fiber_0s[i]] = flines
			subwlimg[0,fiber_0s[i]] = wsub
		endfor
		
		;repeat the flux filling for the highsample array 
		for i=0, 6 do begin 
			fline = fsub * fiber_fractions[i] / round(fiber_core_pix*highsample_factor) 
			flines = dblarr(2048,round(fiber_core_pix*highsample_factor)) 
			for j=0, round(fiber_core_pix*highsample_factor)-1 do begin 
				flines[0,j] = fline 
			endfor
			subimg_highsample[0,fiber_0s_highsample[i]] = flines
		endfor

		
		;repeat everything for the calibration array if cal flux is provided
		if calflag then begin
			subcalwls = reform(calw[cal_li0[order]:cal_li1[order]])
			subcalfls = (reform(calf[order,*]))[cal_li0[order]:cal_li1[order]]
			rebin_ryan, pixel_sampling, res, subcalwls, subcalfls, calwr, calfr, min(subcalwls,/nan), max(subcalwls,/nan)
			calfr = calfr[4:-6]
			calwr = calwr[4:-6]
			calfsub = calfr[0:2047]
			calwsub = calwr[0:2047]
			cline = calfsub / 2d / round(fiber_core_pix)
			clines = dblarr(2048,round(fiber_core_pix))
			cwlines = dblarr(2048,round(fiber_core_pix))
			for j=0, round(fiber_core_pix)-1 do begin 
				clines[0,j] = cline
				cwlines[0,j] = calwsub
			endfor
			for j=0, 1 do begin
				subimg[0,cal_0s[j]] = clines
				subwlimg[0,cal_0s[j]] = cwlines
			endfor
			
			cline = calfsub / 2d / round(fiber_core_pix*highsample_factor)
			clines = dblarr(2048,round(fiber_core_pix*highsample_factor))
			for j=0, round(fiber_core_pix*highsample_factor)-1 do begin 
				clines[0,j] = cline
			endfor
			for j=0, 1 do begin
				subimg_highsample[0,cal_0s_highsample[j]] = clines
			endfor
			
		endif

		;blur on the x-dispersion direction
		
	
		
		xs = indgen(subimg_height)
		xs_highsample = indgen(subimg_height*highsample_factor)
		for i=0, 2047 do begin
			if i mod 100 eq 0 then print,'column ',i
			farr = reform(subimg_highsample[i,*])
			;form kernel with 20X FWHM
			kernel_width = 19 * highsample_factor
			kernel_center = floor(kernel_width / 2)
			gx = dindgen(kernel_width)
			kernel = gaussian(gx,[1.,kernel_center,(1d*highsample_factor)/2.3548d])
			cfarr = convol(farr,kernel,/normalize,/edge_truncate)
			ys = rebin(cfarr,subimg_height) * 1d * highsample_factor
			subimg[i,*] = ys
		endfor

		;project subimg into img
		;order = i
		
		p1 = round(-(y0pos + cgap[order] - subimg_height/2))
		p2 = p1 - subimg_height + 1 
		p1_positive = 2047 + p1
		p2_positive = 2047 + p2
		;;print,'p1 ',p1
		;;print,'p2 ',p2
		;specimg[0,round(-(y0pos + cgap[order] - subimg_height/2))] += subimg
		specimg[*,p2_positive:p1_positive] += subimg
		;wlimg[0,round(-(y0pos + cgap[order] - subimg_height/2))] = subwlimg
		wlimg[*,p2_positive:p1_positive] = subwlimg

		warray[*,order] = wsub
		
	endfor
	

end


