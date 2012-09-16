;; simple_slit_projection
;; HPF DETECTOR SIMULATOR
;; This routine does the baseline projection of a dispersed spectrum into the cross-dispersion
;; direction
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

pro simple_slit_projection, w, f, res, pixel_sampling, wlimg, specimg, calw=calw, calf=calf, diag_out = diag_out, warray = warray


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PROJECTION PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	diag_out = ''
	
	calflag = calw ne !null
	
	no = n_elements(f[*,0])
	
	gap=[0.0, 2.79, 2.67, 2.55, 2.45, 2.35, 2.25, 2.16, 2.08, 2.00, 1.93, 1.86, 1.79, 1.73, 1.67, 1.62, 1.56] * 1000 / 18. ;;pixels, center to center
	ordernum = LINDGEN(17)+46
	lambdalow=[13173, 12893, 12624, 12366, 12119, 11881, 11653, 11433, 11221, 11017, 10821, 10631, 10448, 10270, 10099, 9934, 9773]/1d4
	lambdahigh=[13390, 13105, 12832, 12570, 12319, 12078, 11845, 11622, 11407, 11199, 10999, 10806, 10620, 10440, 10266, 10098, 9935]/1d4


	nw      = n_elements(w)
	nypix = 22.
	nxarr = 2048L

	if wlimg eq !null then wlimg = dblarr(nxarr,nxarr)
	if specimg eq !null then specimg = dblarr(nxarr,nxarr)
	
	y0pos= ROUND((nxarr - ROUND(TOTAL(gap)))/2.) ;what does this do?
	y0poscal=y0pos + nypix + 10 ;;origin of the cal fiber
	cgap=ROUND(TOTAL(gap, /Cumulative))
	norders=N_ELEMENTS(lambdalow)
	warray=MAKE_ARRAY(nxarr, norders, /Double, Value=!Values.F_NAN)

	
	nxarr = 2048
	calflsr = dblarr(no,nxarr) ;resampled cal fiber
	calwlsr = dblarr(no,nxarr) ;cal wavelengths
	li0=value_locate(w,lambdalow) ;what array index does lambdalow correspond to?
	li1=value_locate(w, lambdahigh) ;what array index does lambdahigh correspond to?

	if calflag then begin
		cal_li0 = value_locate(calw,lambdalow) ;same as above but for cal fiber
		cal_li1 = value_locate(calw,lambdahigh)
	endif

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

		;create order images
		farr = dblarr(2048,nypix)
		warr = dblarr(2048,nypix)
		
		calfarr = dblarr(2048,nypix)
		
		
		
		;fill in order images
		for i=0, nypix-1 do begin
			farr[*,i] = fsub / double(nypix)
			warr[*,i] = wsub
		endfor
		
		;project order images into actual array
		specimg[0,-(y0pos + cgap[order] - nypix/2)] = farr
		wlimg[0,-(y0pos + cgap[order] - nypix/2)] = warr
		
		;store wavelengths
		warray[0,order] = wsub
		
		
		
		;same as above but for cal fiber
		if calflag then begin
			subcalwls = reform(calw[cal_li0[order]:cal_li1[order]])
			subcalfls = (reform(calf[order,*]))[cal_li0[order]:cal_li1[order]]
			rebin_ryan, pixel_sampling, res, subcalwls, subcalfls, calwr, calfr, min(subcalwls,/nan), max(subcalwls,/nan)
			calfr = calfr[4:-6]
			calwr = calwr[4:-6]
			calfsub = calfr[0:2047]
			calwsub = calwr[0:2047]
			calfarr = dblarr(2048,nypix)
			for i=0, nypix-1 do calfarr[*,i] = calfsub / double(nypix)
			specimg[0,-(y0poscal + cgap[order] - nypix/2)] = calfarr
			
		endif
		
		
	endfor

end
    
	
