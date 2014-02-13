;+
; NAME:
;  hpf_apply_velshift
;
; PURPOSE:
;
;  Apply a velocity shift to the wavelength array of a spectrum. Only apply RV shift to "STAR" spectra.
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_apply_velshift, spec_params, vrad
;
; INPUTS:
;
;	spec_params: A spec_params structure as defined in initialize_spec_params 	
;
;	vrad: Radial velocity in m/s
;	
; OUTPUTS:
;	
;	The spec_params.shift_wl tag is modified to implement the RV shift
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_apply_velshift, spec_params, vrad
	
	c = 2.99792458d8	; m/s
	
	for i=0, n_elements(spec_params.spec_file)-1 do begin
		if spec_params.type[i] ne 'STAR' then continue
		case strupcase(strcompress(spec_params.rv_type,/remove_all)) of
		'RELATIVISTIC': begin
			beta=double(vrad/c)
			*(spec_params.shift_wl[i])=double(*(spec_params.wl[i])*SQRT( (1d + beta) / (1d - beta) ))
		end
		'NEWTONIAN': begin
			*(spec_params.shift_wl[i]) = double(*(spec_params.wl[i])*(1d + vrad/c))
		end
		else: message,'UNKNOWN RV STYLE'
		endcase
		spec_params.rv[i] = vrad
	endfor
	
end