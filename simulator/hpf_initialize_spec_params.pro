;+
; NAME:
;  hpf_initialize_spec_params
;
; PURPOSE:
;
;  Initialize the structure that contains parameters related to the spectra
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_initialize_spec_params()
;
; INPUTS:
;	
; OUTPUTS:
;	
;	A structure containing parameters for the spectra
;
; KEYWORD PARAMETERS:
;
;	init_params: A structure containing a subset of the parameters listed below, which will supersede the default parameters in this file
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_initialize_spec_params, init_params = init_params

	if n_elements(init_params) ne 0 then nspec = n_elements(init_params.spec_file) else nspec = 3
	
	out = { wl:ptrarr(3,/allocate_heap), $
		fl:ptrarr(3,/allocate_heap), $
		shift_wl:ptrarr(3,/allocate_heap), $
		spec_file:['support/bt_34_extended.fits','~/work/spectra/ffp_35_300_10nwpermode_5000.0_001.00.fits',''], $ ;'~/work/spectra/ffp_30_80.fits'
		type:['STAR','CAL','FLAT'],$
		rv:0d, $
		rv_type:'RELATIVISTIC', $
		tellcontam_flag:[0,0,0], $
		tellcontam_file:['support/tellspec2_detected.fits','',''], $
		upsample_factor:[12,1,1], $
		filter:[1d,1d-9,1], $ ;1d-2 for ffp_30_80
		jmag:[9d,!values.f_nan,!values.f_nan], $
		flat:[0,0,1], $
		normalize_output:[0,0,0], $
		output_per_wavelength:[1,1,0], $ 
		lfc_lims:make_array(2,3,/double,value=!values.f_nan), $
		lfc_lims_flag:0 ,$
		minl:.7 $ ;limits of input spectrum in microns
		maxl:1.4 }
	
	out_names = tag_names(out)
	;Set init_params if they are input
	if (size(init_params))[2] eq 8 then begin
		par_names = tag_names(init_params)
		for i=0, n_elements(par_names) - 1 do begin
			ind = where(out_names eq par_names[i],nind)
			if nind eq 0 then begin
				;if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter not in det struct: '+par_names[i]
			endif else begin
				s_out = size(out.(ind))
				s_par = size(init_params.(i))
				if s_out[1] ne s_par[1] then begin
					;if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'Parameter dimension mismatch: '+ par_names[i]
				endif else begin
					out.(ind) = init_params.(i)
				endelse
			endelse
		endfor
	endif else begin
		;if keyword_set(verbose) then diag_output = diag_output + string(13B) + 'No init init_params'
	endelse
	
	return,out
end
		
	