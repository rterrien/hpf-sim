;; test_mask_sb2
;; HPF DETECTOR SIMULATOR
;; This routine prepares the masks and runs the mask CCFs and organizes the results
;; 
;; CALLS
;; mask_cc3
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; script_all_xx
;;
;; NOTES
;; 
;;
;; parameters:
;; prefix_sp - location for the spectra files
;; vel_file - location of the velocity .sav file (not used - replaced by header)
;; sav_file - filename to save rvs
;; var - do not use
;; tellfile - telluric .fits file, if included it passes to mask_cc3, which excludes
;;   mask points near telluric features
;; maskfile - file which has a mask (output from my mask creation routine which stores
;;   all 4 parameters of fit for each line, so take elem 1 and every 4th after to get centers
;; weights - file which has weights for the mask points above. 


pro test_mask_sb2, prefix_sp, vel_file, sav_file, var=var, tellfile = tellfile, maskfile = maskfile, weightsfile = weightsfile, diagfile = diagfile
	
	diag = diagfile ne !null
	if diag then begin
		openu,diaglun,diagfile,/get_lun,/append
		printf,diaglun,string(13B)+string(13B)+'############# BEGIN TEST_MASK_SB2.PRO ##############'
		printf,diaglun,'Run at: ',systime()
	endif


	c=299792458.d ;m/s
	mask_width = 3d3 ;4 km/s mask widths
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SET MASKS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;;;;;;;;;MASK FOR HAUSCHILDT ******************************************
	;load the mask
	;;	restore,'res.sav'
	;;	restore,'res_weights.sav'
	;;	;res and weights
	;;	;get first element and every 4th thereafter (line centers)
	;;	inds = indgen(n_elements(res)/4)*4 + 1
	;;	mask_midpoints = res[inds]
	;;	;use 4km/s mask widths
	;;	widths = 4.d3 / c * mask_midpoints
	;;	lp = mask_midpoints - widths ;left sides of mask regions
	;;	rp = mask_midpoints + widths ;right sides of mask regions
	;;;;;;;;;;MASK FOR HAUSCHILDT ******************************************
		
	;;;;;;;;;;MASK FOR SINGLE GAUSSIAN ******************************************
	;;mask_midpoints = [1.1314000d]
	;;	weights = [1.d]
	;;	;use 4km/s mask widths
	;;	widths = 2.d3 / c * mask_midpoints
	;;	lp = mask_midpoints - widths ;left sides of mask regions
	;;	rp = mask_midpoints + widths ;right sides of mask regions
	;;;;;;;;;;MASK FOR SINGLE GAUSSIAN ******************************************
	
	;;;;;;;;;;MASK FOR BT34 ******************************************
	if keyword_set(maskfile) then begin
		restore,maskfile
		restore,weightsfile
	endif else begin
		;load the mask
		restore,'btmask2.sav'
		restore,'btmask_weights2.sav'
		maskfile = 'btmask2.sav'
	endelse
	
	if diag then printf,diaglun,'Mask File: '+maskfile
	
	;res and weights
	;get first element and every 4th thereafter (line centers)
	inds = indgen(n_elements(res)/5)*5 + 1
	mask_midpoints = res[inds]
	;use 4km/s mask widths
	widths = mask_width / c * mask_midpoints
	lp = mask_midpoints - widths ;left sides of mask regions
	rp = mask_midpoints + widths ;right sides of mask regions
	
	if diag then printf,diaglun,'Mask Width(m/s): ',mask_width
	;;;;;;;;;;MASK FOR BT34 ******************************************
	
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;LOOP THROUGH SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	files = file_search(prefix_sp+'*.fits')
	vels = dblarr(n_elements(files))
	sns = dblarr(n_elements(files))
	rvs = dblarr(n_elements(files))

	for i=0, n_elements(files)-1 do begin
		spec = mrdfits(files[i],0,hdr)
		wl = reform(spec[0,*])
		fl = reform(spec[1,*])
		vel1 = sxpar(hdr,'vrad')
		sns1 = sqrt(sxpar(hdr,'medcts'))
		vels[i] = vel1
		sns[i] = sns1
		
		if diag then begin
			printf,diaglun,'CCF for #: ',string(i,format='(I5)')
			printf,diaglun,'file: '+files[i]
			printf,diaglun,'RV(m/s): '+string(vel1,format='(D+10.3)')
		endif

		;there used to be a normalization here, but since there is an order-by-order
		;smoothing in the extraction code, it is not necessary
		fln = fl 
		
		;ccf recover velocities
		mask_cc3,wl,fln,rv,lp,rp,weights,rverr,/nan,tellfile=tellfile, diag_output = diag_output
		
		if diag then begin
			printf,diaglun,string(13B)+'******* mask_cc3.pro OUTPUT BEGIN ********'
			printf,diaglun,diag_output
			printf,diaglun,string(13B)+'******* mask_cc3.pro OUTPUT END **********'
		endif
				
		rvs[i] = rv
	
		;save every 10th result incase something goes wrong	
		if i mod 10 eq 0 then save,rvs,filename=sav_file+'_tmp' ;this is slow, so save progress
		
	endfor
	
	save,rvs,filename=sav_file
	
	dv = vels - rvs
	
	save,sns,filename='sns_'+sav_file
	save,dv,filename='dv_'+sav_file
	
	if diag then begin
		printf,diaglun,string(13B)+string(13B)+'############# END TEST_MASK_SB2.PRO ##############'
		free_lun,diaglun
	endif

	
	
end