;+
; NAME:
;  hpf_expose_cds
;
; PURPOSE:
;
;  Step the detector through a CDS exposure and associated timestep.
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_expose_cds, det_params, fluence, exptime, infile, outfile
;
; INPUTS:
;
;	det_params: structure defined by initialize_det_params
;
;	exptime: the exposure time in seconds
;
;	infile: the fits file with the fluence
;
; 	outfile: the location for the fits output image
;
; OUTPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_expose_cds, det_params, exptime, infile, outfile


	common seeds, seed

	fluence = mrdfits(infile,0,hdr)

	
	;if det.qe_flag then set_qewl_detector,det,qearr
	
	;reset detector
	hpf_reset_detector,det_params
	
	;read f1
	hpf_read_detector,det_params,f1

	;set fluence
	hpf_set_fluence_detector,det_params,fluence
	
	;let time pass
	hpf_timestep_detector,det_params,double(exptime)
	;timestep_detector,det,.01

	;read f2
	hpf_read_detector,det_params,f2

	cf = f2 - f1
	
	;figure out mean sn
	meansn = det_params.signal_cts / sqrt( det_params.signal_cts + det_params.thermal_cts + double(det_params.read_noise)^2. )

	;store result
	h = hdr
	sxaddpar,h,'exptime',exptime
	get_juldate,jd
	sxaddpar,h,'jd_exp',jd
	sxaddpar,h,'dark_current',det_params.dark_current
	sxaddpar,h,'qe_flag',uint(det_params.qe_flag)
	sxaddpar,h,'persist_flag',uint(det_params.persist_flag)
	sxaddpar,h,'ipc_flag',uint(det_params.ipc_flag)
	sxaddpar,h,'read_noise',det_params.read_noise
	sxaddpar,h,'reset_noise',det_params.reset_noise
	sxaddpar,h,'photon_noise_flag',uint(det_params.photon_noise_flag)
	sxaddpar,h,'op_mask_flag',uint(det_params.op_mask_flag)
	sxaddpar,h,'well_depth',det_params.well_depth
	sxaddpar,h,'flatten',uint(det_params.flat_flag)
	sxaddpar,h,'meancts',det_params.signal_cts
	sxaddpar,h,'meansn',meansn
	sxaddpar,h,'signalcts',det_params.signal_cts
	sxaddpar,h,'thermcts',det_params.thermal_cts
	sxaddpar,h,'ad_flag',det_params.ad_flag
	

	;outfile = prefix_out + file_basename(fluence_files[i],'.fits')+'_'+outnum+'_image.fits'
	
	;if diag then printf,diaglun,string(13B)+'Expose Output: '+outfile
	
	mwrfits,cf,outfile,h,/create


	

	
	
end