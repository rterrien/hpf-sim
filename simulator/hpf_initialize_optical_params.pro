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

	;if n_elements(init_params) ne 0 then nspec = n_elements(init_params.spec_file) else nspec = 3
	
	out = {fiber_core_um:300, $ ;300 100/300/300
		fiber_cladding_um:489.51d, $ ;400
		fiber_buffer_um:489.51d, $ ;500
		fiber_extra_sep_um:0d , $;200d, $
		fiber_scale:0.0286d, $ ;3 pix/ 100um
		nfibers:3, $
		slitwidth_um:100d, $;100d, $90
		projection_upsample:4d, $
		convol_upsample: 4d, $
		orders_shape: 0, $
		orders_inty:0,$
		wlimg_file:'', $
		model_file:'support/model_072014_shift.sav', $ ;support/model_020714_shift.sav
		kernel_type:'step3_tilt', $ ;or gaussian or hrs or step3 or circlehat or step3_tilt
		warp_style:'polyclip', $ ;polyclip or tri_surf
		extract_width_mod:2, $
		linear_warp:0, $
		bypass_warp:0, $
		extra_orders_blue:2, $
		extra_orders_red:2, $
		reimager_kernel:0 ,$
		slit_angle:0d, $
		pixel_size:18d-3, $
		n_pixels:2048d } ;degrees
		
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

		