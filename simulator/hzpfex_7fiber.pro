;; hzpfex_7fiber.pro (and nonan.pro)
;; HPF DETECTOR SIMULATOR
;; This is the extraction code
;; It takes data from the image files, adds up with pre-defined orders
;; (profile fitting will be implemented in the future?)
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; extract
;;
;; NOTES
;; 
;;
;; parameters:
;; infile - the input image file
;; inwlfile - the input wavelength image
;; outfile - output file
;; varfile - do not use
;; tellfile - telluric spectrum .sav file for a simple telluric correction

PRO nonan, array

;;Removes nan's in a 2-d array by averaging over their
;;neighbors
;;from chad's nonaninf.pro

  s=SIZE(array)
  
  junk = where(finite(array), complement=indx)

  ;;there is nothing to do
  IF indx[0] EQ -1 THEN RETURN

  FOR i=0,N_ELEMENTS(indx)-1 DO BEGIN
     xy=WHEREUNWRAP(array,indx[i])
     x1=(xy[0]-1) > 0
     x2=(xy[0]+1) < (s[1]-1)
     y1=(xy[1]-1) > 0
     y2=(xy[1]+1) < (s[2]-1)
     array[xy[0],xy[1]]=MEAN(array[x1:x2,y1:y2],/NaN)
  ENDFOR

END



PRO HZPFEX_7FIBER, infile, inwlfile, outfile, varfile = varfile, tellfile = tellfile, crimg = crimg, diag_output = diag_output, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = orders_gaps, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um


	diag_output = ''
	
	diag_output = diag_output + '***********************************************'
	diag_output = diag_output + string(13B) + 'BEGIN HZPFEX_7FIBER RESULTS'
	diag_output = diag_output + string(13B) + '**********************************************'
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;OPTICAL MODEL DEFINITIONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	array_pix    = 2048           ; array size
	slitl_pix   = 22              ; slit length		   
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
		
	if fiber_scale eq !null or fiber_core_um eq !null or fiber_cladding_um eq !null or fiber_buffer_um eq !null then stop
	
	fiber_core_pix = fiber_core_um * fiber_scale
	fiber_cladding_pix = fiber_cladding_um * fiber_scale
	fiber_buffer_pix = fiber_buffer_um * fiber_scale
	
	fiber_outside_thickness = (fiber_buffer_pix - fiber_core_pix)/2d
	

	
	
	diag_output = diag_output + string(13B) + 'Fiber Scale [pix / um]: ' + string(fiber_scale,format='(G6.3)')
	diag_output = diag_output + string(13B) + 'Fiber Core um: ' + string(fiber_core_um,format='(D5.2)') + ' pixels: ' + string(fiber_core_pix,format='(D5.2)')
	diag_output = diag_output + string(13B) + 'Fiber Cladding um: ' + string(fiber_cladding_um,format='(D5.2)') + ' pixels: ' + string(fiber_cladding_pix,format='(D5.2)')
	diag_output = diag_output + string(13B) + 'Fiber Buffer um: ' + string(fiber_buffer_um,format='(D5.2)') + ' pixels: ' + string(fiber_buffer_pix,format='(D5.2)')

	
	subimg_height = round(9.*fiber_buffer_pix + 10.*fiber_buffer_pix + 20)

	array_pix = 2048L
	norders=n_elements(lambdalow)
	y0pos = round(subimg_height) + 20
	cgap=ROUND(TOTAL(gap, /Cumulative))
	
	
	

  
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ IN IMAGE AND WAVELENGTH FILES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	image=MRDFITS(infile, 0, hdr)
	
	imsize=SIZE(image, /Dimensions)
	IF (imsize[0] NE array_pix) || (imsize[1] NE array_pix) THEN BEGIN
		PRINT, 'Error, image '+infile+' is wrong size'
		RETURN
	ENDIF

	;;Read in the wavelength files
    restore,inwlfile
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REPLACE BAD PIXELS IF THEY EXIST
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;only replace pixels on orders - prevent edge pixels from being averaged with
	;interorder pixels
	a = where(~finite(image),na)
	
	diag_output = diag_output + string(13B) + 'Number Bad Pixels: '+string(na,format='(I7)')
	
	if na gt 0 then begin
		for order=0, norders-1 do begin
			p1 = round(-(y0pos + cgap[order] - subimg_height/2))
			p2 = p1 - subimg_height + 1 
			p1_positive = 2047 + p1
			p2_positive = 2047 + p2
			subimg = image[*,p2_positive:p1_positive]
			nonan,subimg
			image[*,p2_positive:p1_positive] = subimg
		endfor
	endif
  
  

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONSTRUCT SUB-SPEC ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	subspec=MAKE_ARRAY(array_pix, 2, norders, /Double, Value=0d)
	calspec=MAKE_ARRAY(array_pix, 2, norders, /Double, Value=0d) ;cfb added for cal fiber
	
	
	;if varfile ne !null then subvar = make_array(array_pix,norders,/Double,value=!values.f_nan)
	FOR order=0, norders-1 DO BEGIN
		;select out order image defined by projection routine
		p1 = round(-(y0pos + cgap[order] - subimg_height/2))
		p2 = p1 - subimg_height + 1 
		p1_positive = 2047 + p1
		p2_positive = 2047 + p2
		
		subimg = image[*,p2_positive:p1_positive]
		
		fiber_0s = round(indgen(7)*(2*fiber_buffer_pix) + 3*fiber_buffer_pix + fiber_outside_thickness + 10)
		cal_0s = round([fiber_outside_thickness+10+fiber_buffer_pix,max(fiber_0s) + 2*fiber_buffer_pix])

		wls = reform(warray[*,order])
		;binsizes
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]]

		subspec[*,0,order] = wls
		calspec[*,0,order] = wls
		;plot,subimg[200,*],ps=10
		;vline,round(fiber_0s-2),co=fsc_color('red')
		;vline,round(fiber_0s+fiber_core_pix+1),co=fsc_color('red')
		
		
		for i=0, 6 do begin
			;0s are fiber/cal_0s and + fiber_core_pix
			;pad this with 1 pixel each way for extraction? as a first guess
			;fiber_img = subimg[*,fiber_0s[i]-2 : round(fiber_0s[i] + fiber_core_pix + 1)]
			fiber_img = subimg[*,fiber_0s[i]-1 : round(fiber_0s[i] + fiber_core_pix)]
			fiber_collapsed = total(fiber_img,2,/double)
			subspec[*,1,order] += fiber_collapsed
			;plot,fiber_img[200,*],ps=10
			;stop
		endfor
		
		for i=0, 1 do begin
			cal_img = subimg[*,cal_0s[i]-2 : round(cal_0s[i] + fiber_core_pix + 1)]
			cal_collapsed = total(cal_img,2,/double)
			calspec[*,1,order] += cal_collapsed
			;stop
		endfor
		
		;divide by binsize to get flux/wavelength
		subspec[*,1,order] /= double(xsi)
		calspec[*,1,order] /= double(xsi)
		
		
	ENDFOR

	

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;LOAD TELLURIC SPECTRUM IF NECESSARY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if keyword_set(tellfile) then begin
		tellspec = mrdfits(tellfile)
		tellwl = reform(tellspec[0,*])
		tellfl = reform(tellspec[1,*])
		nans = where(~finite(tellfl),nnans)
		if nnans eq 0 then stop
		tellfl[nans] = 1.d
		diag_output = diag_output + string(13B) + 'Telluric Correction File: '+tellfile
	endif else begin
		diag_output = diag_output + string(13B) + 'Telluric Correction File: None'
	endelse
	
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REASSEMBLE THE ORDERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;Loop backwards since the last order is the shortest wavelength  
  	FOR order=norders-1, 0, -1 DO BEGIN
		subfinite=WHERE(FINITE(subspec[*, 0, order]))
		subwave=reform(subspec[subfinite, 0, order])
		subflux=reform(subspec[subfinite, 1, order])
		subcalwave=reform(calspec[subfinite,0,order])
    	subcalflux=reform(calspec[subfinite,1,order])

		
		;telluric correct if necessary
		;testing - remove below
		subflux2 = subflux
		;end testing
		if keyword_set(tellfile) then begin
			tellfl_1 = interpol(tellfl,tellwl,subwave)
			subflux = subflux / tellfl_1
			bad = where(tellfl_1 lt .5,nbad)
			if nbad ne 0 then subflux[bad] = !values.f_nan
			tellcor = 1
			sxaddpar,hdr,'tellcor',tellcor  
		endif
		;continuum normalize
		nel = n_elements(subwave)
		boxcarmax,subwave,subflux,nel/50.,wlsm,ms
		wlfit = ((subwave - min(subwave))/ (max(subwave) - min(subwave)))*2.d - 1.d
		wlfitin = wlfit[where(subflux eq subflux)]
		fit = svdfit(wlfitin,ms,5,/legendre,/double)
		fl_nvals = double(rleg(wlfit,fit))
		subflux = subflux / fl_nvals

		;;Three possibilities for order overlap
		if order eq norders-1 then begin
			wave = subwave
			flux = subflux
			calwave=calspec[subfinite,0,order]
			calflux=calspec[subfinite,1,order]
		endif else begin
			CASE 1 OF
			wave[-1] GT subwave[0]: BEGIN
				inew=VALUE_LOCATE(subwave, wave[-1])
				wave=[wave, subwave[inew:*]]
				flux=[flux, subflux[inew:*]]
				calwave=[calwave,subcalwave[inew:*]]
        		calflux=[calflux,subcalflux[inew:*]]
				if dovar then var = [var, subvariance[inew:*]]
			END
			wave[-1] EQ subwave[0]: BEGIN
				wave=[wave, subwave[1:*]]
				flux=[flux, subflux[1:*]]
				calwave=[calwave,subcalwave[1:*]]
				calflux=[calflux,subcalflux[1:*]]
				if dovar then var = [var, subvariance[1:*]]
			END
			wave[-1] LT subwave[0]: BEGIN
				delta0=wave[-1]-wave[-2]
				delta1=subwave[1]-subwave[0]
				delta=MEAN([delta0, delta1])
				nfill=FLOOR((subwave[0]-wave[-1])/delta)
				wfill=(dindgen(nfill)+1)*delta + wave[-1]
				wave=[wave, wfill, subwave]
				flux=[flux, MAKE_ARRAY(nfill, /Double, Value=!Values.F_NaN), subflux]
				calwave=[calwave,wfill,subcalwave]
        		calflux=[calflux,wfill,subcalflux]
			END
			ENDCASE
		endelse
	ENDFOR

	mwrfits, [1d#wave, 1d#flux], outfile, hdr,/create
	mwrfits,[1d#calwave,1d#calflux],outfile ;output cal spec to ext 1 of same file
	
	diag_output = diag_output + string(13B) + 'Mean Flux: '+string(mean(flux,/nan),format='(G10.5)')
	diag_output = diag_output + string(13B) + 'StDev Flux: '+string(stddev(flux,/nan),format='(G10.5)')
	diag_output = diag_output + string(13B) + 'Mean Flux in Cal: '+string(mean(calflux,/nan),format='(G10.5)')
	
	diag_output += string(13B) + '***********************************************'
	diag_output += string(13B) + 'END HZPFEX_7FIBER RESULTS'
	diag_output += string(13B) + '**********************************************'

  
END                                                        