;; initialize_detector
;; HPF DETECTOR SIMULATOR
;; This code initializes the detector structure
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;; all parameters in the detector structure (see below)
;;
;; CALLED BY
;; expose
;;
;; NOTES
;; 
;;
;; parameters:
;; n_pixels - number of pixels to include in detector
;; output_struct - name for the output structure
;; parameters - parameters to use for detector (defaults are used if not replaced)
;; verbose - keyword to set verbose output


pro initialize_detector, n_pixels, output_struct, parameters = parameters, verbose=verbose, diag_output = diag_output

	diag_output = ''
	;common block for random draws
	common seeds, seed
	
	verbose=1
	
	if output_struct ne !NULL then begin
		print, "warning: output struct already defined"
		if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'warning: output structure already defined'
	endif
	
	;Init structure
	m=4
	n=n_pixels
	det = {n:n_pixels, $
		   m: m, $ 						;reference pixels
		   dark_current: 0.05d, $		;e-/s (double)
		   read_noise: 30L, $			;e- (long)
		   reset_noise: 18L, $			;e- (long)
		   photon_noise_flag: 1b, $		;flag to include photon noise in the timestep routine
		   outputs:32, $				;number of outputs (int)
		   ad_flag:1b, $				;flag to include A/D conversion on array readout
		   op_mask:bytarr(n,n),$		;operability mask (byte)
		   op_mask_flag:1b,$			;flag to include pixel operability
		   qe:dblarr(n,n), $			;quantum efficiency (doubles)
		   qe_flag:1b, $				;flag to include QE
		   qe_loc:'', $					;either '17', '25', or a filename with QE(wl)
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
		   well_depth:80000L, $			;well depth (e-) (long)
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
		   filter_spec:dblarr(2,1000)}	;filter spectrum, for the non-dispersed bg flux
				
		   
	;initialize QE
	det.qe[*] = 1.
	
	;initialze ipc coeffecients
	det.ipc_coeff = randomu(seed,n,n,4,/normal) * det.ipc_sd + det.ipc_mean
	
	det_names = tag_names(det)
	
	;Set parameters if they are input
	if (size(parameters))[2] eq 8 then begin
		par_names = tag_names(parameters)
		for i=0, n_elements(par_names) - 1 do begin
			ind = where(det_names eq par_names[i],nind)
			if nind eq 0 then begin
				if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter not in det struct: '+par_names[i]
			endif else begin
				s_det = size(det.(ind))
				s_par = size(parameters.(i))
				if s_det[1] ne s_par[1] then begin
					if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter dimension mismatch: '+ par_names[i]
				endif else begin
					det.(ind) = parameters.(i)
				endelse
			endelse
		endfor
	endif else begin
		if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'No init parameters'
	endelse

	;store output
	output_struct = det
	
end
	
		
				
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
