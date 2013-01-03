;; expose.pro
;; HPF DETECTOR SIMULATOR
;; This is the second-level script that controls the detector
;; Exposes the fluence files and stores images
;; 
;; CALLS
;; initialize_detector
;; set_wlim_detector
;; set_qewl_detector
;; reset_detector
;; read_detector
;; set_fluence_detector
;; 
;; 
;; DIRECTLY MODIFIES
;; det.op_mask
;; 
;; CALLED BY
;; script_all_XX
;;
;; NOTES
;; 
;;
;; parameters:
;; prefix_fl - flux directory - string
;; prefix_out - output directory - string
;; var - don't use
;; exptime - exposure time in seconds - scalar
;; iparam - a detector structure, include parameters to replace defaults
;; badpixelmask - the filename of a bytarr .sav file
;; outnum - a number to put in the image filename (usually exposure time)
;; filerange - string expression to feed file_search to use a subset of fluence files



pro expose, prefix_fl, prefix_out, var=var, exptime = exptime, iparam = iparam, badpixelmask = badpixelmask, outnum=outnum, filerange = filerange, diagfile = diagfile

	common seeds, seed
	
	diag = diagfile ne !null
	if diag then begin
		openu,diaglun,diagfile,/get_lun,/append
		printf,diaglun,string(13B)+string(13B)+'############# BEGIN EXPOSE.PRO ##############'
		printf,diaglun,'Run at: ',systime()
	endif
	
	;include necessary code for detector
	;!path = '/Volumes/RAID/HZPF/simulator/rvtest_20120403/pro/:'+!path
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;FIGURE OUT LOCATIONS OF FLUENCE FILES AND OUTPUT FILES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	fluence_dir = prefix_fl
	if keyword_set(filerange) then fluence_files = file_search(filerange,/fully_qualify_path) else $
		fluence_files = file_search(fluence_dir+'*[0-999].fits',/fully_qualify_path)
	wlim_files = file_search(fluence_dir+'*wlimg.fits',/fully_qualify_path)
	wl_files = fluence_files + '_wl.sav'
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DIAGNOSTIC PRINTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if diag then begin
		printf,diaglun,string(13B)+'FLUENCE FILES ------------------'
		for di=0, n_elements(fluence_files)-1 do printf,diaglun,fluence_files[di]
		
		printf,diaglun,string(13B)+'WLIM FILES ---------------------'
		for di=0, n_elements(wlim_files)-1 do printf,diaglun,wlim_files[di]
		
		printf,diaglun,string(13B)+'WL FILES -----------------------'
		for di=0, n_elements(wl_files)-1 do printf,diaglun,wl_files[di]
		
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CREATE THE DETECTOR
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	initialize_detector, 2048., det, parameters=iparam, diag_output = diag_output
	if diag then begin
		printf,diaglun,string(13B)+'## BEGIN INITIALIZE_DETECTOR REPORT ##'
		printf,diaglun,diag_output
		printf,diaglun,string(13B)+'## END INITIALIZE_DETECTOR REPORT ##'
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;LOAD THE QE FUNCTIONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.qe_flag then begin
		qe_function_lookup,det,qearr
		qe_wl = reform(qearr[0,*])
		qe_qe = reform(qearr[1,*])
	endif
	
	if diag and det.qe_flag then begin
		printf,diaglun,string(13B)+'QE CHOICE: ',det.qe_loc
		printf,diaglun,string(13B)+'QE WL: ',qe_wl
		printf,diaglun,string(13B)+'QE   : ',qe_qe
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SET BADPIXEL MASK
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.op_mask_flag then begin
		bpfile = 'badpixelmask1.sav'
		if keyword_set(badpixelmask) then bpfile = badpixelmask
		restore,bpfile
		det.op_mask = mask
	endif

	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;LOOP THROUGH FLUENCES, EXPOSE AND STORE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	nf = n_elements(fluence_files)
	
	if diag then printf,diaglun,string(13B)+'# Fluence Files :',nf
		

	for i=0, nf-1 do begin
		npix = det.n
		
		;get fluence
		flu = mrdfits(fluence_files[i],0,hdr)
		
		;get wavelength image
		wlim = mrdfits(wlim_files[i],0)
		set_wlim_detector,det,wlim
		
		if det.qe_flag then set_qewl_detector,det,qearr
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;SET GLASS FILTER ARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		if det.filter_flag then begin
			glass_filter_detector,det
			if diag then begin
				printf,diaglun,string(13B)+'GLASS FILTER INCLUDED'
				printf,diaglun,'GLASS TYPE: ',det.filter_glass
				printf,diaglun,'GLASS THICKNESS (mm): ',det.filter_thick
				printf,diaglun,'AR COATING: ',det.filter_ar_coating_flag
			endif
		endif else begin
			if diag then printf,diaglun,'Glass Filter Not Included'
		endelse

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;SET THERMAL BACKGROUND ARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		bgdiag = ''
		if det.bg_flag then begin
			thermal_bg_detector,det,diag_output = bgdiag
			printf,diaglun,bgdiag
			if diag then begin
				printf,diaglun,string(13B)+'THERMAL BACKGROUND INCLUDED'
				printf,diaglun,'Thermal Background Mean: ',mean(det.bgarr)
				printf,diaglun,'Thermal Background StDev: ',stddev(det.bgarr)
			endif
		endif else begin
			if diag then printf,diaglun,'Thermal Background Not Included'
		endelse

		
		;reset detector
		reset_detector,det
		
		;read f1
		read_detector,det,f1

		;set fluence
		set_fluence_detector,det,flu
		
		;let time pass
		timestep_detector,det,double(exptime)
		;timestep_detector,det,.01

		;read f2
		read_detector,det,f2

		cf = f2 - f1
		
		;figure out mean sn
		meansn = det.signal_cts / sqrt( det.signal_cts + det.thermal_cts + double(det.read_noise)^2. )

		;store result
		h = []
		sxaddpar,h,'exptime',exptime
		get_juldate,jd
		vrad = sxpar(hdr,'vrad')
		sxaddpar,h,'jd',jd
		sxaddpar,h,'vrad',vrad
		sxaddpar,h,'dark_current',det.dark_current
		sxaddpar,h,'qe_flag',uint(det.qe_flag)
		sxaddpar,h,'persist_flag',uint(det.persist_flag)
		sxaddpar,h,'ipc_flag',uint(det.ipc_flag)
		sxaddpar,h,'read_noise',det.read_noise
		sxaddpar,h,'reset_noise',det.reset_noise
		sxaddpar,h,'photon_noise_flag',uint(det.photon_noise_flag)
		sxaddpar,h,'op_mask_flag',uint(det.op_mask_flag)
		sxaddpar,h,'well_depth',det.well_depth
		sxaddpar,h,'flatten',uint(det.flat_flag)
		sxaddpar,h,'meancts',det.signal_cts
		sxaddpar,h,'meansn',meansn
		sxaddpar,h,'signalcts',det.signal_cts
		sxaddpar,h,'thermcts',det.thermal_cts
		sxaddpar,h,'ad_flag',det.ad_flag
		

		outfile = prefix_out + file_basename(fluence_files[i],'.fits')+'_'+outnum+'_image.fits'
		
		if diag then printf,diaglun,string(13B)+'Expose Output: '+outfile
		
		mwrfits,cf,outfile,h,/create
		
		
	endfor
	
	if diag then begin
		printf,diaglun,string(13B)+'Last Header -------'
		printf,diaglun,h
		printf,diaglun,string(13B)+string(13B)+'############# END EXPOSE.PRO ##############'
		free_lun,diaglun
	endif
	
	
end