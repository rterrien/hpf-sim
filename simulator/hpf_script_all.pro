;+
; NAME:
;  hpf_script_all
;
; PURPOSE:
;
;	The top-level script that calls up all the necessary functions to do a HPF simulation
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	hpf_script_all, n
;
; INPUTS:
;
;	n: The index of the run (for filenames)
;
; OUTPUTS:
;	
; KEYWORD PARAMETERS:
;
;	datdir/outdir: Manually set the directories for data and output
;
;	nvels: The number of images to simulate (each with its own velocity) (default=2)
;
;	skip_fluence/expose/extract/mask: Skip the indicated step
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_script_all, index, datdir = datdir, outdir = outdir, nvels = nvels, skip_fluence = skip_fluence, skip_expose = skip_expose, skip_extract = skip_extract, skip_mask = skip_mask, invel = invel, start = start, finish = finish

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; BASIC OVERALL PATHS and PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	print,memory()/1d6

	ind = string(index,format='(I02)')

	remove_outputs = 0

	exposure_time = 300d

	if ~keyword_set(nvels) then nvels = 2
	vels=(randomu(seed,nvels,/uniform)-0.5) * 60d * 1000d
	if n_elements(invel) ne 0 then begin
		vels = invel
		nvels = n_elements(invel)
	endif

	num = indgen(nvels)
	
	if n_elements(start) ne 0 then begin
		nvels = finish - start + 1
		num = indgen(nvels) + start
		vels=(randomu(seed,nvels,/uniform)-0.5) * 60d * 1000d
	endif

	if ~(keyword_set(datdir) and keyword_set(outdir)) then begin
		computer = where_am_i()
		case computer of
		-1: begin
			print,'need datdir and outdir paths'
			stop
		end
		0: begin
			datdir = '/Volumes/RAID/HZPF/simulator/current_20130101/data/'
			outdir = '/Volumes/RAID/HZPF/simulator/current_20130101/result/'
		end
		1: begin
			datdir = '~/research/hzpf_sims/prelim/sbtest/test5/data/'
			outdir = '~/research/hzpf_sims/prelim/sbtest/test5/out/'
		end
		2: begin
			datdir = '/gpfs/scratch/rct151/simulator/data/'
			outdir = '/gpfs/scratch/rct151/simulator/out/'
		end
		endcase
	endif

	diag_output = outdir+'diag_output_'+ind+'.txt'

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; FORM NECESSARY DIRECTORIES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	fldir = datdir+'fl'+ind+'/'
	imdir = datdir+'im'+ind+'/'
	spdir = datdir+'sp'+ind+'/'
	crdir = datdir+'cr'+ind+'/'

	spawn,'mkdir ' + fldir
	spawn,'mkdir ' + imdir
	spawn,'mkdir ' + spdir
	spawn,'mkdir ' + crdir

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; DETECTOR PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	init_det_param={ $
		dark_current:0.00d, $ ;0.05
		qe_flag:0b, $
		qe_loc:'support/h2rg_2.5qe_24micron', $
		ad_flag:1b, $
		persist_flag:0b, $
		ipc_flag:0b, $
		ipc_mean:0.01, $
		ipc_sd:0.002, $
		read_noise:18L, $ ;18
		reset_noise:0L, $
		photon_noise_Flag:1b, $
		op_mask_flag:0b, $
		well_depth:1d30, $
		flat_flag:0b, $
		bg_flag:0b, $
		bg_temp:170.d, $
		filter_flag:0b, $
		filter_glass:'kzfsn5short',$
		filter_thick:6d,$
		filter_ar_coating_flag:1b,$
		optical_model_version:2,$
		optical_model_rangechoice:'red'}
		
	det_params = hpf_initialize_det_params(init_params = init_det_param)
	
	print,'MEMORY AFTER DET INIT'
	print,memory()/1d6
		
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; SPECTRA PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	 init_optical_params = { wl:ptrarr(3,/allocate_heap), $
    fl:ptrarr(3,/allocate_heap), $
    shift_wl:ptrarr(3,/allocate_heap), $
    spec_file:['support/bt_34_extended.fits','~/work/spectra/ffp_35_300_10nwpermode_5000.0_001.00.fits',''], $ ;'~/work/spectra/ffp_30_80.fits'   support/bt_34_extended.fits
    type:['STAR','CAL','FLAT'],$
    rv:0d, $
    rv_type:'RELATIVISTIC', $
    tellcontam_flag:[0,0,1], $
    tellcontam_file:['support/tellspec.fits','','support/tell_vzyj_unsmoothed_nocont.fits'], $
    upsample_factor:[12,1,1], $ ;12,1,1
    filter:[1d,1d-9,1d], $ ;1d-2 for ffp_30_80
    jmag:[9d,!values.f_nan,!values.f_nan], $
    flat:[0,0,1], $
    normalize_output:[1,0,1], $ ; 1,0,1
    output_per_wavelength:[1,1,0], $ ; 1,1,0
    lfc_lims:[[.813d,.935d],[.964d,1.105d],[1.159d,1.284d]],$
    lfc_lims_flag:1 }
	
	spec_params = hpf_initialize_spec_params(init_params = init_optical_params)
	print,'MEMORY AFTER SPEC INIT'
	print,memory()/1d6
	hpf_load_specs, spec_params, det_params
	print,'MEMORY AFTER LOAD SPEC'
	print,memory()/1d6
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; OPTICAL PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	 init_optical_params = {fiber_core_um:300, $ ;300 100/300/300
    fiber_cladding_um:489.51d, $ ;400
    fiber_buffer_um:489.51d, $ ;500
    fiber_extra_sep_um:0d , $;200d, $
    fiber_scale:0.0286d, $ ;3 pix/ 100um
    nfibers:3, $
    slitwidth_um:100d, $;100d, $90
    projection_upsample:2d, $ ;4x4
    convol_upsample: 2d, $
    orders_shape: 0, $
    orders_inty:0,$
    wlimg_file:'', $
    model_file:'support/model_v1_nov30_2014.sav', $ ;support/model_020714_shift.sav or model_072014_shift.sav or support/model_072014_shift.sav
    kernel_type:'step3_tilt', $ ;or gaussian or hrs or step3 or circlehat or step3_tilt
    warp_style:'polyclip', $ ;polyclip or tri_surf
    extract_width_mod:2, $
    linear_warp:0, $
    bypass_warp:0, $
    extra_orders_blue:2, $
    extra_orders_red:2, $
    reimager_kernel:0 ,$
    slit_angle:0d} ;degrees

	optical_params = hpf_initialize_optical_params(init_params = init_optical_params)
	print,'MEMORY AFTER LOAD SPEC'
	print,memory()/1d6
	optical_params.wlimg_file = fldir + 'wlimg.fits'
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; PROJ PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	proj_params = hpf_process_optical_model(optical_params)
	print,'MEMORY AFTER PROC OPT'
	print,memory()/1d6
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; MASK PARAMETERS
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	init_mask_params = {mask_file:'support/btmask2.sav',$
    weights_file:'support/btmask_weights2.sav',$
    width_ms:3.d3,$
    excl_nan:1,$
    ccf1_velrange:80000d,$
    ccf1_nvel:501d,$
    ccf2_velrange:20000d,$ ;20000
    ccf2_nvel:501d,$
    ccf1_nterms:5,$
    ccf2_nterms:4,$
    excl_tell:0,$
    telluric_spec:'support/tellspec3_detected.fits',$
    tell_excl_level:0.98d }
	
    mask_params = hpf_initialize_mask(spec_params, init_params = init_mask_params)

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   SINGLE ROUND
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	outs = fldir+'test1fl'+string(num,format='(I03)')+'.fits'
	
	; Create Fluence array
	
	if ~keyword_set(skip_fluence) then begin
	
		for i=0, nvels-1 do begin
		
			hpf_apply_velshift,spec_params,vels[i]

			print,'MEMORY BEFORE TRAM PROJ #',i
			print,memory()/1d6
			if i eq 0 then hpf_tram_project, spec_params, optical_params, proj_params, det_params, outs[i], /generate_wlimg else $
				hpf_tram_project, spec_params, optical_params, proj_params, det_params, outs[i]
		
		endfor
	endif
	
		
	; Expose the fluence arrays
	
	if ~keyword_set(skip_expose) then begin
	
		fluence_files = file_search(fldir+'*[0-999].fits',/fully_qualify_path,count = nfiles)
		for i=0, nfiles-1 do begin
			outfile = imdir + file_basename(fluence_files[i],'.fits')+'_image.fits'
			hpf_expose_cds, det_params, exposure_time, fluence_files[i], outfile
		endfor
	endif
	
	; Extract the exposures
	
	if ~keyword_set(skip_extract) then begin
	
		images = file_search(imdir+'*[0-999]_image.fits',/fully_qualify_path,count=nfiles)
		for i=0, nfiles-1 do begin	
			outfile = spdir + file_basename(images[i],'_image.fits')+'_spec.fits'
			hpf_extract, spec_params, optical_params, proj_params, det_params, images[i], outfile
		endfor

	endif
		
	; Do the mask analysis
	
	if ~keyword_set(skip_mask) then begin
		
		spec_files = file_search(spdir+'*_spec.fits',/fully_qualify_path,count=nfiles)
	
		measured_rvs = dblarr(nfiles)
		
		known_rvs = dblarr(nfiles)
	
		for i=0, nfiles-1 do begin
			spec = mrdfits(spec_files[i],0,h)
			spec_struct = {wl:double(reform(spec[0,*])), fl:double(reform(spec[1:*,*]))}
			known_rvs[i] = sxpar(h,'RV')
			
			mask_res = hpf_mask_rv(spec_params, mask_params, spec_struct)
			measured_rvs[i] = mask_res.rv_out
		endfor
	
		overall_scatter = stddev(measured_rvs - known_rvs)
		outfile = outdir+'results_'+string(index,format='(I03)')+'_'+string(systime(/julian),format='(D13.5)')+'.txt'
		openw,1,outfile
		printf,1,index
		printf,1,'overall scatter: ',overall_scatter
		close,1
		
		out_rvs = outdir+'result_rvs_'+string(index,format='(I03)')+'_'+string(systime(/julian),format='(D13.5)')+'.fits'
	
		res = {measured_rvs:measured_rvs, known_rvs:known_rvs}
		mwrfits,res,out_rvs,/create
		print,'OVERALL SCATTER: ',overall_scatter
	endif
	
		
	stop

end
