;+
; NAME:
;  hpf_echel_scale_disp
;
; PURPOSE:
;
;  Scale the dispersion of the echellogram (ie stretch the orders in wavelength)
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_echel_add_orders(model,scale,orig_sampling=orig_sampling)
;
; INPUTS:
;
;	model = optical model structure, with xs,ys,ws (n x norders) and orders (norders)
;
;	sampling = number of pixels per slitwidth (using 3pix/100um)
;	
;	
; OUTPUTS:
;	
;	out = modified model structure
;
; KEYWORD PARAMETERS:
;
;	orig_sampling
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 04-28-14
;-


function hpf_echel_scale_disp, model, sampling, orig_sampling = orig_sampling

	if n_elements(orig_sampling) eq 0 then orig_sampling = 2.79
	
	scale = (orig_sampling - sampling) / orig_sampling
	
	xs_new = model.xs * (1d + scale)
	
	return,model_out
	
	
	
end
		

	