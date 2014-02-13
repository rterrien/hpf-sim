;+
; NAME:
;  hpf_initialize_det_params
;
; PURPOSE:
;
;  Initialize the structure that contains parameters related to the detector array
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_initialize_det_params()
;
; INPUTS:
;
; OUTPUTS:
;	
;	A structure containing the detector parameters
;
; KEYWORD PARAMETERS:
;
;	init_params: A structure containing a subset of the parameters listed below, which will supersede the default parameters in this file
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_initialize_det_params, init_params = init_params

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;INITIALIZE STRUCTURE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	m=4
	n=2048
	det = {n:n, $
		   m: m, $ 						;reference pixels
		   dark_current: 0.05d, $		;e-/s (double)
		   read_noise: 30L, $			;e- (long)
		   reset_noise: 18L, $			;e- (long)
		   photon_noise_flag: 1b, $		;flag to include photon noise in the timestep routine
		   outputs:32, $				;number of outputs (int)
		   ad_flag:1b, $				;flag to include A/D conversion on array readout
		   op_mask:bytarr(n,n),$		;operability mask (byte)
		   op_mask_flag:1b,$			;flag to include pixel operability
		   qe_flag:1b, $				;flag to include QE
		   qe_loc:'', $					;either '17', '25', or a filename with QE(wl)
		   qewl:ptr_new(/allocate_heap), $
		   gain: 1.0, $					;gain photons -> e- (float)
		   arr:dblarr(n,n), $			;pixel array (doubles)
		   var:dblarr(n,n), $			;variance array (doubles)
		   time:0., $					;current detector time (float)
		   fluence:dblarr(n,n), $		;current fluence onto detector (doubles)
		   wlim:dblarr(n,n), $			;wavelengths for each pixel
		   persist_time: 100., $		;persistence decay constant (s) (float)
		   persist_percent: 1., $		;persistence percentage (float)
		   persist_flag: 1b, $			;flag to include persistence
		   persist:dblarr(n,n), $		;"persisted" charge (doubles)
		   well_depth:80000D, $			;well depth (e-) (long)
		   ipc_flag: 1b, $				;flag to include IPC
		   ipc_mean: 0.02, $			;mean of IPC coeffs (float)
		   ipc_sd: 0.002, $				;sd of IPC coeffs (float)
		   ipc_coeff:fltarr(n,n,4), $	;IPC coefficients for [top,right,bottom,left] neighbors (float)
		   flat_flag:0b, $				;flag to flatten the array in the expose routine
		   bg_flag:0b, $				;flag if BG is used
		   bg_temp:0d, $				;BG temp for 1-component model
		   bgarr:dblarr(n,n), $			;array of BG fluence
		   signal_cts:0d, $				;an estimate of the mean # of photons from the signal
		   thermal_cts:0d, $			;an estimate of the mean BG counts
		   filter_flag:0b, $			;flag to include a glass filter
		   filter_thick:0d, $			;thickness of the glass filter in mm
		   filter_glass:'', $			;either 'kzfsn5' or 'pk50'
		   filter_ar_coating_flag:0b, $	;flag to use antireflective coating on glass
		   filter_arr:dblarr(n,n),$ 	;the wavelength, transmission array for the filter
		   filter_spec:dblarr(2,1000)} ;filter spectrum, for the non-dispersed bg flux

	det_names = tag_names(det)
	;Set init_params if they are input
	if (size(init_params))[2] eq 8 then begin
		par_names = tag_names(init_params)
		for i=0, n_elements(par_names) - 1 do begin
			ind = where(det_names eq par_names[i],nind)
			if nind eq 0 then begin
				;if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter not in det struct: '+par_names[i]
			endif else begin
				s_det = size(det.(ind))
				s_par = size(init_params.(i))
				if s_det[1] ne s_par[1] then begin
					;if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter dimension mismatch: '+ par_names[i]
				endif else begin
					det.(ind) = init_params.(i)
				endelse
			endelse
		endfor
	endif else begin
		if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'No init init_params'
	endelse
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SET BADPIXEL MASK
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.op_mask_flag then begin
		bpfile = 'support/badpixelmask1.sav'
		if keyword_set(badpixelmask) then bpfile = badpixelmask
		restore,bpfile
		det.op_mask = mask
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;INITIALIZE QE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.qe_flag then begin
		qe_function_lookup,det,qearr
		*det.qewl = qearr
	endif
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;INITIALIZE IPC
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	det.ipc_coeff = randomu(seed,n,n,4,/normal) * det.ipc_sd + det.ipc_mean


	return, det
	
end
	
