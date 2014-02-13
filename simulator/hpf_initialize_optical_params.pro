;+
; NAME:
;  hpf_initialize_optical_params
;
; PURPOSE:
;
;  Initialize the structure that contains parameters related to the optical configuration
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_initialize_optical_params()
;
; INPUTS:
;	
; OUTPUTS:
;	
;	A structure containing parameters for the optical configuration
;
; KEYWORD PARAMETERS:
;
;	init_params: A structure containing a subset of the parameters listed below, which will supersede the default parameters in this file
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_initialize_optical_params, init_params = init_params

	if n_elements(init_params) ne 0 then nspec = n_elements(init_params.spec_file) else nspec = 3
	
	out = {fiber_core_um:300d, $
		fiber_cladding_um:400d, $
		fiber_buffer_um:500d, $
		fiber_extra_sep_um:0d , $;200d, $
		fiber_scale:0.03d, $ ;3 pix/ 100um
		nfibers:3, $
		slitwidth_um:100d, $
		projection_upsample:1d, $
		orders_shape: 0, $
		wlimg_file:'', $
		model_file:'support/model_020714_shift.sav', $
		kernel_type:'hrs'} ;or gaussian
		
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

		