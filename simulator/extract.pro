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

pro extract, prefix_fl, prefix_im, prefix_sp, var = var, tellfile=tellfile, crimg = crimg, diagfile = diagfile

	diag = diagfile ne !null
	if diag then begin
		openu,diaglun,diagfile,/get_lun,/append
		printf,diaglun,string(13B)+string(13B)+'############# BEGIN EXTRACT.PRO ##############'
		printf,diaglun,'Run at: ',systime()
	endif
		

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
		hzpfex_basic2,fname,wname,oname,tellfile=tellfile,crimg=crflag,diag_output = diag_output
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
	