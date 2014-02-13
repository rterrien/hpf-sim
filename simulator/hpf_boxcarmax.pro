;+
; NAME:
;  hpf_boxcarmax
;
; PURPOSE:
;
;  Return the boxcar sliding maximum for normalizing the spectrum.
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_boxcarmax(wl,fl,pix,wl_out,ms_out)
;
; INPUTS:
;
;	wl: wavelength array (units arbitrary)
;
;	fl: flux array (units arbitrary)
;
;	pix: half-width of boxcar in pixels
;
; OUTPUTS:
;	
;	ms_out: the array of maxima
;
; KEYWORD PARAMETERS:
;
;	keep_nan: if an element of the input flux is NaN, then keep NaN in the output array
;
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_boxcarmax, wl, fl, pix, wl_out, ms_out, keep_nan = keep_nan
	nel = n_elements(fl)
	ms = dblarr(nel)
	for i=0, nel-1 do begin
		if ~finite(fl[i]) then ms[i] = !values.f_nan else begin
			temp = max(fl[0 > (i - pix) : (nel-1) < (i + pix)],/nan)
			ms[i] = temp
		endelse
	endfor
	if keyword_set(keep_nan) then begin
		wl_out = wl
		fl_out = ms
		return
	endif
	g = where(finite(ms))
	wl_out = wl[g]
	ms_out = ms[g]
end
