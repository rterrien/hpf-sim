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

pro fiber_projection, w, f, res, pixel_sampling, wlimg, specimg, calw=calw, calf=calf, diag_out = diag_out, fiber_size = fiber_size, fiber_gap = fiber_gap, fiber_fractions = fiber_fractions


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PROJECTION PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	calflag = calw ne !null
	
	diag_out = ''
	
	if fiber_size eq !null then fiber_size = pixel_sampling
	if fiber_gap eq !null then fiber_gap = pixel_sampling
	if fiber_fractions eq !null then fiber_fractions = replicate(1d/7d,7)
	diag_out = diag_out + string(13B) + 'Fiber Projection Size: ' + string(fiber_size,format='(D5.2)')
	diag_out = diag_out + string(13B) + 'Fiber Gap Size: ' + string(fiber_gap,format='(D5.2)')
	diag_out = diag_out + string(13B) + 'Fiber power fractions: ' + string(fiber_fractions,format='(D5.2)')

	
	no = n_elements(f[*,0])
	
	gap=[0.0, 2.79, 2.67, 2.55, 2.45, 2.35, 2.25, 2.16, 2.08, 2.00, 1.93, 1.86, 1.79, 1.73, 1.67, 1.62, 1.56] * 1000 / 18. ;;pixels, center to center
	ordernum = LINDGEN(17)+46
	lambdalow=[13173, 12893, 12624, 12366, 12119, 11881, 11653, 11433, 11221, 11017, 10821, 10631, 10448, 10270, 10099, 9934, 9773]/1d4
	lambdahigh=[13390, 13105, 12832, 12570, 12319, 12078, 11845, 11622, 11407, 11199, 10999, 10806, 10620, 10440, 10266, 10098, 9935]/1d4


	nw      = n_elements(w)
	nypix = 77
	nxarr = 2048L
	
	if wlimg eq !null then wlimg = dblarr(nxarr,nxarr)
	if specimg eq !null then specimg = dblarr(nxarr,nxarr)
	
	y0pos= ROUND((nxarr - ROUND(TOTAL(gap)))/2. + 50) ;what does this do?
	y0poscal=y0pos + nypix + 10 ;;origin of the cal fiber
	cgap=ROUND(TOTAL(gap, /Cumulative))
	norders=N_ELEMENTS(lambdalow)
	warray=MAKE_ARRAY(nxarr, norders, /Double, Value=!Values.F_NAN)
	
	
	nxarr = 2048
	calflsr = dblarr(no,nxarr) ;resampled cal fiber
	calwlsr = dblarr(no,nxarr) ;cal wavelengths
	li0=value_locate(w,lambdalow) ;what array index does lambdalow correspond to?
	li1=value_locate(w, lambdahigh) ;what array index does lambdahigh correspond to?

	cal_li0 = value_locate(calw,lambdalow) ;same as above but for cal fiber
	cal_li1 = value_locate(calw,lambdahigh)
	

	for order=0, no-1 do begin
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
		subimg_height = 9.*fiber_size + 10.*fiber_gap + 20
		subimg = dblarr(2048,subimg_height)
		subwlimg = dblarr(2048,subimg_height)
		
		;define fiber locations on sub image
		fiber_0s = indgen(7)*(fiber_size + fiber_gap) + fiber_size + 2*fiber_gap + 10
		cal_0s = [10,max(fiber_0s)+fiber_size+2*fiber_gap]
		
		
		;loop through and fill in flux into fibers
		for i=0, 6 do begin
			fline = fsub * fiber_fractions[i] / pixel_sampling
			flines = dblarr(2048,pixel_sampling)
			wlines = dblarr(2048,pixel_sampling)
			for j=0, pixel_sampling-1 do begin 
				flines[0,j] = fline
				wlines[0,j] = wsub
			endfor
			subimg[0,fiber_0s[i]] = flines
			subwlimg[0,fiber_0s[i]] = wsub
			
		endfor
		
		if calflag then begin
			subcalwls = reform(calw[cal_li0[order]:cal_li1[order]])
			subcalfls = (reform(calf[order,*]))[cal_li0[order]:cal_li1[order]]
			rebin_ryan, pixel_sampling, res, subcalwls, subcalfls, calwr, calfr, min(subcalwls,/nan), max(subcalwls,/nan)
			calfr = calfr[4:-6]
			calwr = calwr[4:-6]
			calfsub = calfr[0:2047]
			calwsub = calwr[0:2047]
			cline = calfsub / 2d
			clines = dblarr(2048,pixel_sampling)
			cwlines = dblarr(2048,pixel_sampling)
			for j=0, pixel_sampling-1 do begin 
				clines[0,j] = cline
				cwlines[0,j] = calwsub
			endfor
			for j=0, 1 do begin
				subimg[0,cal_0s[j]] = clines
				subwlimg[0,cal_0s[j]] = cwlines
			endfor
			
		endif
		
		;blur subimg along x-disp axis
		for i=0, 2047 do begin
			farr = reform(subimg[i,*])
			;form kernel with 20X FWHM
			kernel_width = 19
			kernel_center = floor(kernel_width / 2)
			gx = dindgen(kernel_width)
			kernel = gaussian(gx,[1.,kernel_center,pixel_sampling/2.3548d])
			cfarr = convol(farr,kernel,/normalize,/edge_truncate)
			subimg[i,*] = cfarr
		endfor

		
		
		;project subimg into img
		specimg[0,round(-(y0pos + cgap[order] - subimg_height/2))] += subimg
		wlimg[0,round(-(y0pos + cgap[order] - subimg_height/2))] = subwlimg

		
	endfor
	
end


