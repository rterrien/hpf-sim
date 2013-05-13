;; extract.pro
;; HPF DETECTOR SIMULATOR
;; This is the second-level script which runs the basic spectra extraction
;; 
;; CALLS
;; hzpfex_basic2
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; script_all_XX
;;
;; NOTES
;; 
;;
;; parameters:
;; prefix_fl - fluence directory, so it can find the wavelength images
;; prefix_im - the image directory
;; prefix_sp - the spectra (output) directory
;; var - do not use
;; tellfile - a .sav file with the telluric spectrum, which the extraction uses to 
;; 	correct the contaminated spectrum
;; crimg - do not use

pro extract, prefix_fl, prefix_im, prefix_sp, var = var, tellfile=tellfile, crimg = crimg, diagfile = diagfile, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = orders_gaps, projection_type = projection_type, orders_norders = orders_norders, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, straight_orders = straight_orders

	diag = diagfile ne !null
	if diag then begin
		openu,diaglun,diagfile,/get_lun,/append
		printf,diaglun,string(13B)+string(13B)+'############# BEGIN EXTRACT.PRO ##############'
		printf,diaglun,'Run at: ',systime()
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;FIGURE OUT ORDER LOCATIONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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


	if ~keyword_set(projection_type) then projection_type = 'simple' 
	if diag then printf,diaglun,string(13B)+'Projection Type: '+projection_type

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;FIGURE OUT LOCATIONS OF INPUTS / OUTPUTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	image_files = file_search(prefix_im+'*_image.fits',/fully_qualify_path)
	;if keyword_set(var) then var_files = file_search(prefix_im+'*_variance.fits')
	names = file_basename(image_files,'_image.fits')
	wl_files = prefix_fl+strmid(names,0,10)+'.fits_wl.sav'
	;32 for btsettl
	;9 for others
	if diag then begin
		printf,diaglun,string(13B)+'IMAGE FILES ------------------'
		for di=0, n_elements(image_files)-1 do printf,diaglun,image_files[di]
		
		printf,diaglun,string(13B)+'NAMES ------------------------'
		for di=0, n_elements(names)-1 do printf,diaglun,names[di]
		
		printf,diaglun,string(13B)+'WL FILES ---------------------'
		for di=0, n_elements(wl_files)-1 do printf,diaglun,wl_files[di]
	endif


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;LOOP THROUGH AND EXTRACT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	nf = n_elements(image_files)
	
	for i=0, nf-1 do begin
		wname = wl_files[i]
		fname = image_files[i]
		oname = prefix_sp+names[i]+'_spec.fits'
		;if keyword_set(var) then vname = var_files[i]
		if diag then begin
			printf,diaglun,'image: ',fname
			printf,diaglun,'wave:  ',wname
			printf,diaglun,'out: ',oname
		endif
		;if keyword_set(var) then hzpfex_basic,fname,wname,oname,varfile=vname else $
		;	hzpfex_basic,fname,wname,oname
		if keyword_set(crimg) then crflag=1 else crflag=0
		
		case projection_type of
			'simple': hzpfex_basic2,fname,wname,oname,tellfile=tellfile,crimg=crflag,diag_output = diag_output, orders_lambdalow = lambdalow, orders_lambdahigh = lambdahigh, orders_gaps = gap
			'7_fibers': hzpfex_7fiber,fname,wname,oname,tellfile=tellfile,crimg=crflag,diag_output = diag_output, orders_lambdalow = lambdalow, orders_lambdahigh = lambdahigh, orders_gaps = gap, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um
			'tram1': hzpfex_tram2, fname, wname, oname, tellfile = tellfile, crimg = crflag, diag_output = diag_output, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = gap, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um,nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, straight_orders = straight_orders
			else: stop
		endcase
		
		if diag then begin
			printf,diaglun,string(13B)+'## BEGIN HZPFEX_BASIC2 REPORT ##'
			printf,diaglun,diag_output
			printf,diaglun,string(13B)+'## END HZPFEX_BASIC2 REPORT ##'
		endif
		
	endfor
	
	if diag then begin
		printf,diaglun,string(13B)+string(13B)+'############# END EXTRACT.PRO ##############'
		free_lun,diaglun
	endif


end
	