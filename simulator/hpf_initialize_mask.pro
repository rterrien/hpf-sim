;+
; NAME:
;  hpf_initialize_mask
;
; PURPOSE:
;
;	Initialize the CCF mask and associated parameters for the mask analysis
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_initialize_mask(spec_params)
;
; INPUTS:
;
;	spec_params: A spec_params structure as described in hpf_initialize_spec_params
;
; OUTPUTS:
;	
;	A structure containing the mask/CCF parameters
;
; KEYWORD PARAMETERS:
;
;	init_params: A structure containing a subset of the parameters listed below, which will supersede the default parameters in this file
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_initialize_mask, spec_params, init_params = init_params

	if n_elements(spec_params) ne 0 then nspec = n_elements(spec_params.spec_file) else nspec = 3
	
	c=299792458.d ;m/s

	out = {mask_file:'support/btmask2.sav',$
		weights_file:'support/btmask_weights2.sav',$
		width_ms:3d3,$
		excl_nan:1,$
		excl_tell:0,$
		ccf1_velrange:80000d,$
		ccf1_nvel:501d,$
		ccf2_velrange:20000d,$ ;20000
		ccf2_nvel:501d,$
		ccf1_nterms:5,$
		ccf2_nterms:5}
		
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

	
	restore,out.mask_file
	restore,out.weights_file
	
	;res and weights
	
	;get first element and every 4th thereafter (line centers)
	inds = indgen(n_elements(res)/5)*5 + 1
	mask_midpoints = res[inds]
	;use 4km/s mask widths
	widths = out.width_ms / c * mask_midpoints
	lp = mask_midpoints - widths ;left sides of mask regions
	rp = mask_midpoints + widths ;right sides of mask regions
	out = jjadd_tag(out,'lp',lp,/array_tag)
	out = jjadd_tag(out,'rp',rp,/array_tag)
	out = jjadd_tag(out,'widths',widths,/array_tag)
	out = jjadd_tag(out,'weights',weights,/array_tag)	

	;mask points are in microns
	
	return,out
	
end
