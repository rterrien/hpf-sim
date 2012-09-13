;; hzpfex_basic2.pro (and nonan.pro)
;; HPF DETECTOR SIMULATOR
;; This is the extraction code
;; It takes data from the image files, collaposes along hard-coded columns,
;; and outputs a sum-extraction
;; 
;; CALLS
;; hzpfex_basic2
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



PRO HZPFEX_BASIC2, infile, inwlfile, outfile, varfile = varfile, tellfile = tellfile, crimg = crimg, diag_output = diag_output


	diag_output = ''
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;OPTICAL MODEL DEFINITIONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	array_pix    = 2048           ; array size
	slitl_pix   = 22              ; slit length		   
	gap=[0.0, 2.79, 2.67, 2.55, 2.45, 2.35, 2.25, $
	   2.16, 2.08, 2.00, 1.93, 1.86, 1.79, 1.73, $
	   1.67, 1.62, 1.56] * 1000 / 18. ;;pixels, center to center
	ordernum = LINDGEN(17)+46
	lambdalow=[13173, 12893, 12624, 12366, 12119, 11881, 11653, 11433, $
			 11221, 11017, 10821, 10631, 10448, 10270, 10099, 9934, 9773]/1d4 
	lambdahigh=[13390, 13105, 12832, 12570, 12319, 12078, 11845, 11622, $
			  11407, 11199, 10999, 10806, 10620, 10440, 10266, 10098, 9935]/1d4
			  
	norders=n_elements(lambdalow)
	y0pos = ROUND((array_pix - ROUND(TOTAL(gap)))/2.)
	y0poscal=y0pos + slitl_pix + 10 ;;origin of the cal fiber
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
			subpos=(-(y0pos+cgap[order]-slitl_pix/2))
			subimg=image[0:*, subpos:subpos+slitl_pix-1]
			nonan,subimg
			image[0:*, subpos:subpos+slitl_pix-1] = subimg
		endfor
	endif
  
;;  if keyword_set(crimg) then begin
;;  	crname = 'cr46/'+file_basename(outfile,'.fits') + '_crimg.fits'
;;  	mwrfits,image,crname,/create
;;  endif
  
  

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONSTRUCT SUB-SPEC ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	subspec=MAKE_ARRAY(array_pix, 2, norders, /Double, Value=!Values.F_NAN)
	calspec=MAKE_ARRAY(array_pix, 2, norders, /Double, Value=!Values.F_NAN) ;cfb added for cal fiber
	;if varfile ne !null then subvar = make_array(array_pix,norders,/Double,value=!values.f_nan)
	FOR order=0, norders-1 DO BEGIN
		subpos=(-(y0pos+cgap[order]-slitl_pix/2))
		subimg=image[0:*, subpos:subpos+slitl_pix-1]
		;;handle cal fiber, 8/3/2012, cfb
		subcalpos=(-(y0poscal+cgap[order]-slitl_pix/2))
		subcal=image[0:*, subcalpos:subcalpos+slitl_pix-1]
		;figure out binsizes to get flux/wavelength
		wls = reform(warray[*,order])
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]]
		wl1 = warray[*,order]
		fl1 = total(subimg,2,/double)/double(xsi) 
		cl1 = total(subcal,2,/double) ;; cfb added
		subspec[*, 0, order]=wl1
		subspec[*, 1, order]=fl1
		;;added to handle cal fiber, 8/3/2012, cfb
		calspec[*,0,order]=wl1
		calspec[*,1,order]=cl1
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

  ;if dovar then $
  	;WRITESPEC, [1#wave, 1#flux, 1#var], outfile, '0.0 0.0 0.0 '+outfile, lformat='D15.13' else $
 	;WRITESPEC, [1#wave, 1#flux], outfile, '0.0 0.0 0.0 '+outfile, lformat='D15.13'
 	;WRITESPEC, [1#wave, 1#flux, 1#var], outfile, '0.0 0.0 0.0 '+outfile, lformat='D20.18' else $
 	;WRITESPEC, [1#wave, 1#flux], outfile, '0.0 0.0 0.0 '+outfile, lformat='D20.18'
 	;mwrfits, [1#wave, 1#flux, 1#var], outfile, hdr,/create else $
	mwrfits, [1d#wave, 1d#flux], outfile, hdr,/create
	mwrfits,[1d#calwave,1d#calflux],outfile ;output cal spec to ext 1 of same file
	
	diag_output = diag_output + string(13B) + 'Mean Flux: '+string(mean(flux,/nan),format='(G10.5)')
	diag_output = diag_output + string(13B) + 'StDev Flux: '+string(stddev(flux,/nan),format='(G10.5)')
	diag_output = diag_output + string(13B) + 'Mean Flux in Cal: '+string(mean(calflux,/nan),format='(G10.5)')
	

  
END                                                        